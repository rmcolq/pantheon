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
import json


def parse_counts(kraken_file):
    counts = defaultdict(int)
    with open(kraken_file,'r') as csv_in:
        domain = "None"
        for line in csv_in:
            if line.startswith("%"):
                continue
            try:
                percentage, num_clade_root, num_direct, rank, ncbi, name = line.strip().split('\t')
            except:
                percentage, num_clade_root, num_direct, a, b, rank, ncbi, name = line.strip().split('\t')

            name = name.lstrip()
            num_direct = int(num_direct)
            num_clade_root = int(num_clade_root)

            if name.startswith("Homo"):
                counts["human"] += int(num_direct)
            elif name == "unclassified":
                counts["unclassified"] += int(num_direct)
            else:
                counts["classified"] += int(num_direct)
                if num_direct > 0:
                    counts["taxa"] += 1
    return counts

def get_input_info(csv_file, relative_directory):
    summary_info = defaultdict(lambda:defaultdict(str))
    with open(csv_file, 'r') as in_csv:
        first = in_csv.readline()
        in_csv.seek(0)
        if "sample_id" in first:
            reader = csv.DictReader(in_csv)
        else:
            reader = csv.DictReader(in_csv, fieldnames=["sample_id", "reads", "kraken_report"])

        if "reads" in reader.fieldnames:
            reads_key = "reads"
        else:
            reads_key = "sample_report"

        for row in reader:
            sample_id, reads, kraken_report = row["sample_id"], row[reads_key], row["kraken_report"]
            if sample_id == "example":
                continue
            reads = reads.replace("/mnt/c","")
            if not reads.startswith("/") and not reads.startswith("s3:"):
                reads = os.path.join(relative_directory, reads)
            kraken_report = kraken_report.replace("/mnt/c","")
            if not kraken_report.startswith("/") and not kraken_report.startswith("s3:"):
                kraken_report = os.path.join(relative_directory, kraken_report)
            summary_info[sample_id] = {
                "sample_id":sample_id,
                "reads": reads,
                "kraken_report": kraken_report
            }
            summary_info[sample_id].update(parse_counts(kraken_report))
    return summary_info

def get_extract_info(extract_summary):
    extract_info = {}
    with open(extract_summary, 'r') as f:
        for line in f:
            data = json.loads(line)
            sample_id = list(data.keys())[0]
            if sample_id == "example":
                continue
            extract_info[sample_id] = defaultdict(lambda:defaultdict(str))
            for taxon_dict in data[sample_id]:
                taxid = taxon_dict["taxon"]
                extract_info[sample_id][taxid] = taxon_dict

    return extract_info

def update_summary_with_extract_info(summary_info, extract_info, keys=[]):
    for sample_id in summary_info:
        summary_info[sample_id]["extracted"] = 0
        if sample_id in extract_info:
            for taxon in extract_info[sample_id]:
                if extract_info[sample_id][taxon]["filenames"] != [] and extract_info[sample_id][taxon]["report_count"] > 0:
                    summary_info[sample_id]["extracted"] += 1
        for key in keys:
            if key not in summary_info[sample_id]:
                summary_info[sample_id][key] = 0

def get_assembly_info(assembly_summary):
    assembly_info = []
    with open(assembly_summary, 'r') as in_csv:
        for line in in_csv:
            taxid, sample_id, reference, read_count = line.strip().split(',')
            if taxid == "TaxonID":
                continue
            assembly_info.append({
            'taxid': taxid,
            'sample_id': sample_id,
            'reference': reference,
            'read_count': int(read_count)
            })
    return assembly_info

def update_summary_with_assembly_info(summary_info, assembly_info, min_read_count):
    for sample_id in summary_info:
        summary_info[sample_id]["assembled"] = 0
        summary_info[sample_id]["assemblies"] = {}
    for entry in assembly_info:
        if entry["read_count"] > min_read_count:
            summary_info[entry["sample_id"]]["assembled"] += 1
            if entry["taxid"] in summary_info[entry["sample_id"]]["assemblies"]:
                summary_info[entry["sample_id"]]["assemblies"][entry["taxid"]] = ",".join([summary_info[entry["sample_id"]]["assemblies"][entry["taxid"]], entry["reference"]])
            else:
                summary_info[entry["sample_id"]]["assemblies"][entry["taxid"]] = entry["reference"]

def get_reference_info(reference_summary):
    reference_info = defaultdict(lambda:defaultdict(str))
    with open(reference_summary, 'r') as in_csv:
        for line in in_csv:
            name, taxid, num_references = line.strip().split(',')
            if name == "name":
                continue
            reference_info[taxid] = {"name":name, "taxid":taxid, "num_refs": num_references}

    return reference_info

def add_taxon_name_to_assembly_info(assembly_info, reference_info):
    for entry in assembly_info:
        entry.update({"taxon_name": reference_info[entry["taxid"]]["name"]})


def make_output_report(report_to_generate, template, run, version, data_for_report):
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

    parser.add_argument("-i","--input", help="CSV of input files")
    parser.add_argument("-r","--reference_summary", help="CSV of references and taxon id")
    parser.add_argument("-e","--extract_summary", help="JSON with count/qualities of read sets per sample_id")
    parser.add_argument("-a","--assembly_summary", help="CSV of counts")

    parser.add_argument("-t","--template", help="HTML template for report", default="assembly_report.mako.html")
    parser.add_argument("--run", help="Run name", default = "pantheon")
    parser.add_argument("--min_read_count", help="Min read count for assembly", type=int)
    parser.add_argument("--version", help="Pantheon version", default="__version__")
    parser.add_argument("--relative_directory", help="The directory where the pipeline was launched", default=".")

    args = parser.parse_args()


    data_for_report = {}

    summary_info = get_input_info(args.input, args.relative_directory)
    extract_info = get_extract_info(args.extract_summary)
    update_summary_with_extract_info(summary_info, extract_info, keys=["classified", "unclassified", "human", "taxa", "extracted"])

    ref_info = get_reference_info(args.reference_summary)
    assembly_info = get_assembly_info(args.assembly_summary)
    update_summary_with_assembly_info(summary_info, assembly_info, args.min_read_count)
    add_taxon_name_to_assembly_info(assembly_info, ref_info)

    data_for_report["summary_table_header"] = ["sample_id", "reads", "kraken_report", "classified", "unclassified", "human", "taxa", "extracted", "assembled", "assemblies"]
    data_for_report["summary_table"] = [summary_info[sample_id] for sample_id in summary_info]
    for key in data_for_report["summary_table_header"]:
        for sample_id in summary_info:
            if key not in summary_info[sample_id]:
                print(sample_id, key, summary_info[sample_id])

    data_for_report["reference_table_header"] = ["name", "taxid", "num_refs"]
    data_for_report["reference_table"] = [ref_info[ref] for ref in ref_info]

    data_for_report["assembly_table_header"] = ["taxon_name", "taxid", "sample_id", "reference", "read_count"]
    data_for_report["assembly_table"] = assembly_info

    outfile = "%s_assembly_report.html" %args.run
    make_output_report(outfile, args.template, args.run, args.version, data_for_report)

if __name__ == "__main__":
    main()