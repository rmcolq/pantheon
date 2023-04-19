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

def groups_from_csv(csv_file):
    groups = defaultdict(lambda: defaultdict(str))
    group_map = defaultdict()

    with open(csv_file, 'r') as in_csv:
        header = True
        for line in in_csv:
            if header:
                header = False
                continue
            id,group,filepath = line.strip().split(",")
            filepath = Path(filepath)
            groups[group][id] = filepath
            groups["all"][id] = filepath
            group_map[id] = group
    return groups, group_map

def get_sample_counts(samples,sample_id_to_filepath,list_taxons=None,sample_counts=None, sample_clade_counts=None, totals=None):
    if not sample_counts:
        sample_counts = defaultdict(lambda: defaultdict(int))
    if not sample_clade_counts:
        sample_clade_counts = defaultdict(lambda: defaultdict(int))
    if not totals:
        totals = {}

    for sample_id in samples:
        file = sample_id_to_filepath[sample_id]
        with file.open() as csv_in:
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

                sample_counts[name][sample_id] = int(num_direct)
                sample_clade_counts[name][sample_id] = int(num_clade_root)
                sample_counts[name]["ncbi_taxon"] = ncbi
                sample_counts[name]["simple_rank"] = rank[0]
                sample_counts[name]["rank"] = rank

    return sample_counts, sample_clade_counts, totals

def add_scores(sample_counts, cases, controls):
    num_cases = len(cases)
    num_controls = max(len(controls),1)
    for taxon in sample_counts:
        case_freq = sum([1 for case in cases if sample_counts[taxon][case] > 0]) / float(num_cases)
        control_count = sum([1 for control in controls if sample_counts[taxon][control] > 0])
        control_freq = control_count/float(num_controls)
        sample_counts[taxon]['score'] = case_freq - control_freq
        sample_counts[taxon]['case_frequency'] = case_freq
        sample_counts[taxon]['control_frequency'] = control_freq
        sample_counts[taxon]['case_max_read_count'] = max([sample_counts[taxon][case] for case in cases])
        control_counts = [sample_counts[taxon][control] for control in controls]
        control_counts.append(0)
        sample_counts[taxon]['control_max_read_count'] = max(control_counts)
    return sample_counts

def make_csv_string(samples, group_map, sample_counts, sample_clade_counts, totals, min_read_count=10):
    csv_string = ""

    header = "sample,group,taxon,taxon_ncbi,taxon_rank,simple_taxon_rank,num_reads,num_clade_reads,sample_total,score,case_frequency,control_frequency,case_max_read_count,control_max_read_count,clade_score,clade_case_frequency,clade_control_frequency,clade_case_max_read_count,clade_control_max_read_count\\n"
    csv_string += header
    for taxon in sample_counts:
        if sample_counts[taxon]["case_max_read_count"] < min_read_count and sample_clade_counts[taxon]["case_max_read_count"] < min_read_count:
            continue
        common_data_raw = "%f,%f,%f,%i,%i"%(sample_counts[taxon]["score"],sample_counts[taxon]["case_frequency"],sample_counts[taxon]["control_frequency"],sample_counts[taxon]["case_max_read_count"],sample_counts[taxon]["control_max_read_count"])
        common_data_clade = "%f,%f,%f,%i,%i"%(sample_clade_counts[taxon]["score"],sample_clade_counts[taxon]["case_frequency"],sample_clade_counts[taxon]["control_frequency"],sample_clade_counts[taxon]["case_max_read_count"],sample_clade_counts[taxon]["control_max_read_count"])

        for sample in samples:
                csv_string += "%s,%s,%s,%s,%s,%s,%i,%i,%i,%s,%s\\n" %(sample, group_map[sample], taxon, sample_counts[taxon]["ncbi_taxon"], sample_counts[taxon]["rank"], sample_counts[taxon]["simple_rank"], sample_counts[taxon][sample], sample_clade_counts[taxon][sample], totals[sample], common_data_raw, common_data_clade)
    return csv_string

def make_output_report(report_to_generate, template, data_for_report= {"HEATMAP_DATA":""}):

    # collate data for tables in the report
    print("have data")

    template_dir = os.path.abspath(os.path.dirname(__file__))

    mylookup = TemplateLookup(directories=[template_dir]) #absolute or relative works
    print(mylookup)
    print("found template")
    mytemplate = mylookup.get_template(template)
    print(mytemplate)
    buf = StringIO()

    ctx = Context(buf,
                    date = date.today(),
                    run_name="test_runname",
                    version = "__version__",
                    sample="sample",
                    heatmap_data = data_for_report["HEATMAP_DATA"])
    print("context")

    try:
        mytemplate.render_context(ctx)
        print("render")

    except:
        print("failed")
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

    parser.add_argument("--input", help="CSV file with ID,group,filepath", required=True)

    parser.add_argument("--output", help="HTML output file summarizing counts", default="report.html")
    parser.add_argument("--template", help="HTML template for report", default="pantheon_report.mako")
    parser.add_argument("--min_reads", type=int, help="Threshold for min number reads", default=10)

    args = parser.parse_args()

    data_for_report = {}

    groups, group_map = groups_from_csv(args.input)
    sample_counts,sample_clade_counts,totals = get_sample_counts(list(groups["all"].keys()), groups["all"])
    for group in groups:
        if group == "all":
            continue
        controls = [n for n in groups["all"] if n not in groups[group]]
        sample_counts = add_scores(sample_counts, groups[group], controls)
        sample_clade_counts = add_scores(sample_clade_counts, groups[group], controls)
        csv_string = make_csv_string(list(groups["all"].keys()), group_map, sample_counts, sample_clade_counts, totals, args.min_reads)

    data_for_report["HEATMAP_DATA"] = csv_string
    make_output_report(args.output, args.template, data_for_report)


if __name__ == "__main__":
    main()