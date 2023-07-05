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

def groups_from_csv(csv_file, summary_table=[], summary_table_headers=[]):
    groups = defaultdict(lambda: defaultdict(str))
    group_map = defaultdict()

    with open(csv_file, 'r') as in_csv:
        try:
            reader = csv.DictReader(in_csv, delimiter=",", quotechar='\"', dialect = "unix")
        except:
            reader = csv.DictReader(in_csv, delimiter="\t", quotechar='\"', dialect = "unix")

        for col in ["sample_id", "filepath"]:
            if col not in reader.fieldnames:
                sys.exit("At a minimum, input file needs a sample_id and filepath column. Existing columns are %s" %reader.fieldnames)

        group = "case"
        barcode = None
        for row in reader:
            id = row["sample_id"]
            filepath = Path(row["filepath"])

            if "barcode" in row:
                barcode = row["barcode"]
            if "group" in row:
                group = row["group"]

            if group.lower() in ["negative control", "negative controls", "negative_control", "negative_controls", "nc", "negative"]:
                group = "negative_control"
            elif group.lower() in ["control", "controls"]:
                group = "control"

            groups[group][id] = filepath
            groups["all"][id] = filepath
            group_map[id] = group
            table_row = {"sample":id, "barcode":barcode, "group":group}
            for header in summary_table_headers:
                if header not in table_row:
                    table_row[header] = None
            summary_table.append(table_row)
    return groups, group_map

def get_sample_counts(samples,sample_id_to_filepath,list_taxons=None,sample_counts=None, taxon_info=None, totals=None):
    if not sample_counts:
        sample_counts = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    if not taxon_info:
        taxon_info = defaultdict(lambda: defaultdict(int))
    if not totals:
        totals = {}

    for sample_id in samples:
        file = sample_id_to_filepath[sample_id]
        with file.open('r') as csv_in:
            for line in csv_in:
                try:
                    percentage, num_clade_root, num_direct, rank, ncbi, name = line.strip().split('\t')
                except:
                    percentage, num_clade_root, num_direct, a, b, rank, ncbi, name = line.strip().split('\t')

                name = name.lstrip()

                if name in ['Homo sapiens', 'unclassified', 'root']:
                    continue
                elif list_taxons and name not in list_taxons:
                    continue

                if sample_id not in totals:
                    totals[sample_id] = 0
                totals[sample_id] += int(num_direct)

                sample_counts[name]["direct"][sample_id] = int(num_direct)
                sample_counts[name]["downstream"][sample_id] = int(num_clade_root)
                taxon_info[name]["taxon_ncbi"] = ncbi
                taxon_info[name]["simple_taxon_rank"] = rank[0]
                taxon_info[name]["taxon_rank"] = rank

    return sample_counts, taxon_info, totals

def get_scores(sample_counts, taxon_info, groups):
    group_scores = {}

    if "negative_control" in groups:
        negative_controls = groups["negative_control"]
    else:
        negative_controls = []

    for group in groups:
        if group in ["all", "control", "negative_control"]:
            continue
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
            if group_scores["direct"][taxon]["case_max_read_count"] < min_read_count and group_scores["direct"][taxon]["case_max_read_count"] < min_read_count:
                continue

            new_dict = dict(taxon_info[taxon])

            new_dict["taxon"] = taxon
            new_dict["sample"] = sample
            new_dict["group"] = group_map[sample]
            new_dict["direct_count"] = sample_counts[taxon]["direct"][sample]
            new_dict["downstream_count"] = sample_counts[taxon]["downstream"][sample]
            new_dict["sample_total"] = totals[sample]

            new_dict["direct"] = group_scores["direct"][taxon]
            new_dict["direct"]["num_reads"] = sample_counts[taxon]["direct"][sample]
            new_dict["downstream"] = group_scores["downstream"][taxon]
            new_dict["downstream"]["num_reads"] = sample_counts[taxon]["downstream"][sample]
            data_list.append(new_dict)

    return data_list

def make_output_report(report_to_generate, template, group, data_for_report= {"heatmap_data":"", "summary_table":""}):
    template_dir = os.path.abspath(os.path.dirname(__file__))
    mylookup = TemplateLookup(directories=[template_dir]) #absolute or relative works
    mytemplate = mylookup.get_template(template)

    buf = StringIO()

    ctx = Context(buf,
                    date = date.today(),
                    version = "__version__",
                    group=group,
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

    parser.add_argument("--input", help="CSV file with sample_id,group,filepath", required=True)

    parser.add_argument("--prefix", help="HTML output prefix ", default="pantheon")

    parser.add_argument("--template", help="HTML template for report", default="pantheon_report.mako.html")
    parser.add_argument("--min_reads", type=int, help="Threshold for min number reads", default=10)

    args = parser.parse_args()

    data_for_report = {"heatmap_data":""}
    data_for_report["summary_table_header"] = ["sample","barcode","group","status","warnings"]
    data_for_report["summary_table"] = []

    groups, group_map = groups_from_csv(args.input, data_for_report["summary_table"], data_for_report["summary_table_header"])
    sample_counts,taxon_info,totals = get_sample_counts(list(groups["all"].keys()), groups["all"])
    group_scores = get_scores(sample_counts, taxon_info, groups)

    relevant_groups = [g for g in groups if g not in ["all", "control", "negative_control"]]
    for group in relevant_groups:
        outfile = args.prefix
        if len(relevant_groups) > 1:
            outfile += "_" + group
        outfile += "_report.html"

        data_list = make_data_dict(taxon_info, group_scores[group], sample_counts, group_map, totals, args.min_reads)
        data_for_report["heatmap_data"] = data_list

        make_output_report(outfile, args.template, group, data_for_report)


if __name__ == "__main__":
    main()