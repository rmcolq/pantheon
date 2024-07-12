#!/usr/bin/env python

import json
import csv
from collections import defaultdict
from pathlib import Path
from mako.lookup import TemplateLookup
import datetime as dt
from datetime import date
import argparse
from mako.template import Template
from mako.runtime import Context
from mako.exceptions import RichTraceback
from io import StringIO
import os
from Bio import SeqIO
import sys

from kraken_utils import KrakenReport, KrakenCombined

def infer_index(columns, index=None):
    # preference for specified index if it is a column
    # then in order "barcode","sample","id" or "sample_id" if they are a column
    # finally falls back on 1st column
    if index and index not in columns:
        index = None

    if not index:
        index = columns[0]
        for key in ["barcode","sample","id", "sample_id"]:
            if key in columns:
                index = key
                break

    return index

def get_summary_info(list_csvs, relative_directory, group="group"):

    summary_info = defaultdict(lambda:defaultdict(str))

    index = None
    keys = set()
    for csv_file in list_csvs:

        with open(csv_file, 'r') as in_csv:
            try:
                reader = csv.DictReader(in_csv, delimiter=",", quotechar='\"', dialect = "unix")
            except:
                reader = csv.DictReader(in_csv, delimiter="\t", quotechar='\"', dialect = "unix")

            index = infer_index(reader.fieldnames, index=None)

            for row in reader:
                sample = row[index]
                summary_info[sample].update(row)
                summary_info[sample]["index"] = sample
                keys.update(row.keys())
                #if "sample_report" not in row and "sample_report" not in summary_info[sample]:
                #    summary_info[sample]["sample_report"] = None

    file_keys = [key for key in keys if key.endswith("path") or key.endswith("report") or key.endswith("file")]

    filepath_key = None
    for key in ["kraken_file", "kraken_report", "kreport"]:
        if key in file_keys:
            filepath_key = key
            break
    if not filepath_key:
         sys.exit("One of the input files needs a 'kraken_report' column giving the path to the kraken report. Near match columns found are %s" %file_keys)

    incomplete = []
    for sample in summary_info:
        for key in file_keys:
            if key in summary_info[sample]:
                if not summary_info[sample][key].startswith("/") and not summary_info[sample][key].startswith("s3:"):
                    summary_info[sample][key] = os.path.join(relative_directory, summary_info[sample][key])
                summary_info[sample][key] = Path(summary_info[sample][key])
        summary_info[sample]["filepath"] = summary_info[sample][filepath_key]

        #if group in keys and group not in summary_info[sample]:
        #    print("No group found for sample %s: removing from output" %sample)
        #    incomplete.append(sample)

    #for sample in incomplete:
    #    del summary_info[sample]

    return summary_info

def update_summary_info_column(summary_info, column, column_dict):
    for sample in summary_info:
        if sample in column_dict:
            summary_info[sample][column] = column_dict[sample]
        else:
            summary_info[sample][column] = None

def update_summary_info_columns(summary_info, new_dict):
    for sample in summary_info:
        summary_info[sample].update(new_dict[sample])


def groups_from_info(summary_info, group_column="group", group_default="case"):
    groups = defaultdict(list)
    group_map = {}

    for sample in summary_info:

        group = group_default
        if group_column in summary_info[sample].keys():
            group = summary_info[sample]["group"]

        if group.lower() in ["negative control", "negative controls", "negative_control", "negative_controls", "nc", "negative"]:
            group = "negative_control"
        elif group.lower() in ["control", "controls"]:
            group = "control"

        groups[group].append(sample)
        group_map[sample] = group


    return groups, group_map

def get_sample_counts(summary_info, list_taxa=[]):
    dict_kraken_reports = {sample:str(summary_info[sample]["filepath"]) for sample in summary_info}
    kc = KrakenCombined(dict_kraken_reports=dict_kraken_reports, list_taxa=list_taxa)
    for sample in summary_info:
        sample_info = kc.get_sample_info(sample)
        new_info = {k:sample_info[k] for k in ('classified','unclassified','host','total') if k in sample_info}
        summary_info[sample].update(new_info)
    kc.filter_taxa(min_percent_of_samples=5, min_number_of_samples=2, min_percent_of_reads=0, by_group=True)
    kc.taxa = kc.taxa[kc.taxa["pass"]]
    return kc

def scores_from_df(grouping, df):
    result = defaultdict(lambda: defaultdict(float))
    has_counts = df>0
    for category in ["cases", "controls", "negative_controls"]:
        freqs = has_counts[grouping[category]].sum(axis=1)/max(len(grouping[category]), 1)
        for taxon,freq in freqs.items():
            result[taxon][f"{category[:-1]}_frequency"] = freq
        max_count = df[grouping[category]].max(axis=1).fillna(0).astype(int)
        for taxon,count in max_count.items():
            result[taxon][f"{category[:-1]}_max_read_count"] = int(count)
    for taxon in result:
        result[taxon]["score"] = result[taxon]["case_frequency"] - result[taxon]["control_frequency"] - result[taxon]["negative_control_frequency"]
    return result

def get_group_scores(cases, controls, negative_controls, kraken_combined):
    grouping = {"cases": cases, "controls": controls, "negative_controls": negative_controls}
    group_scores = {"direct": scores_from_df(grouping, kraken_combined.ucounts), "downstream":scores_from_df(grouping, kraken_combined.counts)}
    return group_scores


def get_scores(kraken_combined, groups):
    group_scores = {}

    if "negative_control" in groups:
        negative_controls = groups["negative_control"]
    else:
        negative_controls = []

    for group in groups:
        cases = groups[group]

        if "control" in groups:
            controls = groups["control"]
        else:
            controls = []
            for other in groups:
                if other not in [group, "negative_control"]:
                    controls.extend(groups[other])

        #if group in ["all", "control", "negative_control"]:
        #    continue
        group_scores[group] = get_group_scores(cases, controls, negative_controls, kraken_combined)

    # Add a summary group for all
    group = "all"
    cases = [n for n in groups["all"] if n not in negative_controls]
    controls = []
    group_scores[group] = get_group_scores(cases,controls,negative_controls,kraken_combined)

    return group_scores

def make_data_dict(kraken_combined, group_scores, group_map, min_read_count=10):
    data_list = []

    num_taxa_added = 0
    for taxon_id, taxon_row in  kraken_combined.taxa.iterrows():
        if not taxon_row["domain"]:
            continue

        # first pass filter
        all_score_zero = True
        for sample, sample_row in kraken_combined.metadata.iterrows():
            if group_scores[group_map[sample]]["downstream"][taxon_id]["case_max_read_count"] < min_read_count:
                continue
            elif group_scores[group_map[sample]]["downstream"][taxon_id]["score"] != 0:
                all_score_zero = False
        if all_score_zero:
            continue

        num_taxa_added += 1
        for sample, sample_row in kraken_combined.metadata.iterrows():
            new_dict = {k:taxon_row[k] for k in ('rank','domain','name') if k in taxon_row}
            new_dict["taxon"] = new_dict["name"]
            new_dict["simple_rank"] = new_dict["rank"][0]
            new_dict["taxon_id"] = str(taxon_id)

            new_dict["sample"] = sample
            new_dict["group"] = group_map[sample]
            new_dict["direct_count"] = kraken_combined.ucounts.loc[taxon_id,sample]
            new_dict["downstream_count"] = kraken_combined.counts.loc[taxon_id,sample]
            new_dict["sample_total"] = sample_row["total"]

            new_dict["direct"] = dict(group_scores[group_map[sample]]["direct"][taxon_id])
            new_dict["direct"]["num_reads"] = kraken_combined.ucounts.loc[taxon_id,sample]
            new_dict["downstream"] = dict(group_scores[group_map[sample]]["downstream"][taxon_id])
            new_dict["downstream"]["num_reads"] = kraken_combined.counts.loc[taxon_id,sample]
            data_list.append(new_dict)

    return data_list

def make_output_report(report_to_generate, template, version, run, data_for_report= {"groups_null":[], "groups_all":[] ,"heatmap_data":"", "summary_table":""}):
    template_dir = os.path.abspath(os.path.dirname(__file__))
    mylookup = TemplateLookup(directories=[template_dir]) #absolute or relative works
    mytemplate = mylookup.get_template(template)

    buf = StringIO()

    ctx = Context(buf,
                    date = date.today(),
                    version = version,
                    run=run,
                    data_for_report = data_for_report)

    try:
        mytemplate.render_context(ctx)

    except:
        traceback = RichTraceback()
        for (filename, lineno, function, line) in traceback.traceback:
            print("File %s, line %s, in %s" % (filename, lineno, function))
            print(line, "\n")
        print("%s: %s" % (str(traceback.error.__class__.__name__), traceback.error))

    with open(report_to_generate, 'w') as fw:
        print("Generating: " + f"{report_to_generate}")
        fw.write(buf.getvalue())


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-i","--input", help="One or more TSV/CSV files containing a 'barcode'/'sample' column. Must include a 'filepath' column specifying the kraken report and ideally also 'group' and 'sample_report' columns", required=True, nargs='+')

    parser.add_argument("-p","--prefix", help="HTML output prefix ", default="pantheon")
    parser.add_argument("-r","--run", help="Run name", default="Pantheon")

    parser.add_argument("-t","--template", help="HTML template for report", default="heatmap_report.mako.html")
    parser.add_argument("-m","--min_reads", type=int, help="Threshold for min number reads", default=5)
    parser.add_argument("--version", help="Pantheon version", default="__version__")
    parser.add_argument("--relative_directory", help="The directory where the pipeline was launched", default=".")

    args = parser.parse_args()

    summary_info = get_summary_info(args.input, args.relative_directory)
    kraken_combined = get_sample_counts(summary_info)

    groups, group_map = groups_from_info(summary_info)
    group_scores = get_scores(kraken_combined, groups)

    data_list = make_data_dict(kraken_combined, group_scores, group_map, args.min_reads)

    outfile = args.prefix + "_heatmap_report.html"


    data_for_report = {}
    data_for_report["groups_null"] = [null] + [i for i in list(groups.keys()) if i != "" and i != "all"]
    data_for_report["groups_all"] = ["all"] + [i for i in list(groups.keys()) if i != "" and i != "all"]
    data_for_report["heatmap_data"] = data_list
    data_for_report["summary_table"] = [summary_info[sample] for sample in summary_info]

    out_columns = ["sample","barcode","group"]
    final_columns = ["status","warnings"]
    additional_columns = [c for c in data_for_report["summary_table"][0].keys() if c not in out_columns and c not in final_columns and not c.endswith("path") and not c.endswith("report")]
    out_columns.extend(additional_columns)
    out_columns.extend(final_columns)
    data_for_report["summary_table_header"] = out_columns

    make_output_report(outfile, args.template, args.version, args.run, data_for_report)


if __name__ == "__main__":
    main()