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

def get_summary_info(list_csvs):

    summary_info = defaultdict(lambda:defaultdict(str))

    index = None
    keys = set()
    for csv_file in list_csvs:

        with open(csv_file, 'r') as in_csv:
            try:
                reader = csv.DictReader(in_csv, delimiter=",", quotechar='\"', dialect = "unix")
            except:
                reader = csv.DictReader(in_csv, delimiter="\t", quotechar='\"', dialect = "unix")

            if not index:
                for key in ["barcode","sample","id"]:
                    if key in reader.fieldnames:
                        index = key
                        break

            if not index or index not in reader.fieldnames:
                sys.exit("At a minimum, input file needs a 'barcode', 'sample' or 'id' column. Existing columns in %s are %s" %(csv_file, reader.fieldnames))

            for row in reader:
                sample = row[index]
                summary_info[sample].update(row)
                keys.update(row.keys())
                if "sample_report" not in row and "sample_report" not in summary_info[sample]:
                    summary_info[sample]["sample_report"] = None

    file_keys = [key for key in keys if key.endswith("path") or key.endswith("report") or key.endswith("file")]

    filepath_key = None
    for key in ["filepath", "kraken_file"]:
        if key in file_keys:
            filepath_key = key
            break

    if not filepath_key:
         sys.exit("One of the input files needs a 'filepath' column giving the path to the kraken report. Near match columns found are %s" %file_keys)

    incomplete = []
    for sample in summary_info:
        for key in file_keys:
            if key in summary_info[sample]:
                summary_info[sample][key] = Path(summary_info[sample][key])
        summary_info[sample]["filepath"] = summary_info[sample][filepath_key]
        if "group" in keys and "group" not in summary_info[sample]:
            print("No group found for sample %s: removing from output" %sample)
            incomplete.append(sample)
    for sample in incomplete:
        del summary_info[sample]

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


def groups_from_info(summary_info):
    groups = defaultdict(lambda: defaultdict(str))
    group_map = defaultdict()


    for sample in summary_info:
        group = "case"
        if "group" in summary_info[sample].keys():
            group = summary_info[sample]["group"]

        if group.lower() in ["negative control", "negative controls", "negative_control", "negative_controls", "nc", "negative"]:
            group = "negative_control"
        elif group.lower() in ["control", "controls"]:
            group = "control"

        filepath = summary_info[sample]["filepath"]

        groups[group][sample] = filepath
        groups["all"][sample] = filepath
        group_map[sample] = group

    return groups, group_map

def get_sample_counts(samples,sample_to_filepath,list_taxons=None,sample_counts=None, taxon_info=None, totals=None, count_info=None):
    if not sample_counts:
        sample_counts = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    if not taxon_info:
        taxon_info = defaultdict(lambda: defaultdict(int))
    if not totals:
        totals = {}
    if not count_info:
        count_info = defaultdict(lambda:{"human":0, "unclassified":0, "classified":0, "num_taxa":0})

    for sample in samples:
        file = sample_to_filepath[sample]
        with file.open('r') as csv_in:
            domain = None
            for line in csv_in:
                if line.startswith("%"):
                    continue
                try:
                    percentage, num_clade_root, num_direct, rank, ncbi, name = line.strip().split('\t')
                except:
                    percentage, num_clade_root, num_direct, a, b, rank, ncbi, name = line.strip().split('\t')

                name = name.lstrip()

                if name.startswith("Homo"):
                    count_info[sample]["human"] += int(num_direct)
                    continue
                elif name == "unclassified":
                    count_info[sample]["unclassified"] += int(num_direct)
                    continue
                elif name in ['root']:
                    continue
                elif list_taxons and name not in list_taxons:
                    continue

                count_info[sample]["classified"] += int(num_direct)
                count_info[sample]["num_taxa"] += 1

                if sample not in totals:
                    totals[sample] = 0
                totals[sample] += int(num_direct)

                sample_counts[name]["direct"][sample] = int(num_direct)
                sample_counts[name]["downstream"][sample] = int(num_clade_root)
                taxon_info[name]["taxon_ncbi"] = ncbi
                taxon_info[name]["simple_taxon_rank"] = rank[0]
                taxon_info[name]["taxon_rank"] = rank
                if rank == "D":
                    domain = name
                taxon_info[name]["domain"] = domain

    return sample_counts, taxon_info, totals, count_info

def get_scores(sample_counts, taxon_info, groups):
    group_scores = {}

    if "negative_control" in groups:
        negative_controls = groups["negative_control"]
    else:
        negative_controls = []

    for group in groups:
        #if group in ["all", "control", "negative_control"]:
        #    continue
        group_scores[group] = {"direct": {}, "downstream":{}}
        cases = groups[group]
        if "control" in groups:
            controls = groups["control"]
        else:
            controls = [n for n in groups["all"] if n not in groups[group]]

        print("group", group)
        print(cases)
        print(controls)
        print(negative_controls)

        num_cases = len(cases)
        num_controls = max(len(controls),1)
        num_negative_controls = max(len(negative_controls),1)

        print(num_cases, num_controls, num_negative_controls)

        for taxon in taxon_info:
            group_scores[group]["direct"][taxon] = {}
            group_scores[group]["downstream"][taxon] = {}

            for read_type in ["direct", "downstream"]:

                case_freq = sum([1 for case in cases if sample_counts[taxon][read_type][case] > 0]) / float(num_cases)
                control_count = sum([1 for control in controls if sample_counts[taxon][read_type][control] > 0])
                control_freq = control_count/float(num_controls)
                negative_control_count = sum([1 for control in negative_controls if sample_counts[taxon][read_type][control] > 0])
                negative_control_freq = negative_control_count/float(num_negative_controls)


                group_scores[group][read_type][taxon]['case_frequency'] = case_freq
                group_scores[group][read_type][taxon]['control_frequency'] = control_freq
                group_scores[group][read_type][taxon]['negative_control_frequency'] = negative_control_freq
                group_scores[group][read_type][taxon]['score'] = case_freq - control_freq - negative_control_freq
                #print(case_freq, control_freq, negative_control_freq, case_freq - control_freq - negative_control_freq)

                group_scores[group][read_type][taxon]['case_max_read_count'] = max([sample_counts[taxon][read_type][case] for case in cases])
                control_counts = [sample_counts[taxon][read_type][control] for control in controls]
                control_counts.append(0)
                group_scores[group][read_type][taxon]['control_max_read_count'] = max(control_counts)
    return group_scores

def make_data_dict(taxon_info, group_scores, sample_counts, group_map, totals, min_read_count=10):
    data_list = []

    for taxon in taxon_info:
        for sample in sample_counts[taxon]["direct"]:
            if group_scores[group_map[sample]]["direct"][taxon]["case_max_read_count"] < min_read_count and group_scores[group_map[sample]]["direct"][taxon]["case_max_read_count"] < min_read_count:
                continue

            new_dict = dict(taxon_info[taxon])

            new_dict["taxon"] = taxon
            new_dict["sample"] = sample
            new_dict["group"] = group_map[sample]
            new_dict["direct_count"] = sample_counts[taxon]["direct"][sample]
            new_dict["downstream_count"] = sample_counts[taxon]["downstream"][sample]
            new_dict["sample_total"] = totals[sample]

            new_dict["direct"] = group_scores[group_map[sample]]["direct"][taxon]
            new_dict["direct"]["num_reads"] = sample_counts[taxon]["direct"][sample]
            new_dict["downstream"] = group_scores[group_map[sample]]["downstream"][taxon]
            new_dict["downstream"]["num_reads"] = sample_counts[taxon]["downstream"][sample]
            data_list.append(new_dict)

    return data_list

def make_output_report(report_to_generate, template, run, data_for_report= {"heatmap_data":"", "summary_table":""}):
    template_dir = os.path.abspath(os.path.dirname(__file__))
    mylookup = TemplateLookup(directories=[template_dir]) #absolute or relative works
    mytemplate = mylookup.get_template(template)

    buf = StringIO()

    ctx = Context(buf,
                    date = date.today(),
                    version = "__version__",
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

    parser.add_argument("-t","--template", help="HTML template for report", default="pantheon_report.mako.html")
    parser.add_argument("-m","--min_reads", type=int, help="Threshold for min number reads", default=10)

    args = parser.parse_args()



    data_for_report = {}
    #
    #data_for_report["summary_table"] = []

    summary_info = get_summary_info(args.input)

    groups, group_map = groups_from_info(summary_info)
    sample_counts,taxon_info,totals,count_info = get_sample_counts(list(groups["all"].keys()), groups["all"])

    update_summary_info_column(summary_info, "group", group_map)
    update_summary_info_columns(summary_info, count_info)

    group_scores = get_scores(sample_counts, taxon_info, groups)
    data_list = make_data_dict(taxon_info, group_scores, sample_counts, group_map, totals, args.min_reads)

    outfile = args.prefix + "_report.html"

    data_for_report = {}
    data_for_report["heatmap_data"] = data_list
    data_for_report["summary_table"] = [summary_info[sample] for sample in summary_info]

    out_columns = ["sample","barcode","group"]
    final_columns = ["status","warnings"]
    additional_columns = [c for c in data_for_report["summary_table"][0].keys() if c not in out_columns and c not in final_columns and not c.endswith("path") and not c.endswith("report")]
    out_columns.extend(additional_columns)
    out_columns.extend(final_columns)
    data_for_report["summary_table_header"] = out_columns

    make_output_report(outfile, args.template, args.run, data_for_report)


if __name__ == "__main__":
    main()