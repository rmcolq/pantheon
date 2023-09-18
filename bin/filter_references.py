#!/usr/bin/env python

import mappy as mp
import argparse
from Bio import SeqIO
from collections import defaultdict
import re

def parse_accessions(summary_file):
    accessions = []
    with open(summary_file, 'r') as f:
        for line in f:
            accessions.append(line.split()[0])
    return accessions

def accession_in_list(acc, accessions):
    short_acc = acc.split(".")[0]
    return (short_acc in accessions)


def parse_name(s, strain):
    s = " ".join(s.split()[1:])
    for d in ["polymerase", "nucleocapsid", "neuraminidase", "hemagglutinin"]:
        s = re.split(d, s, flags=re.IGNORECASE)[0]
    for d in ["ha", "pa", "pb1", "pb2", "na"]:
        s = s.replace(" %s gene for" %d, "")
    s = s.strip()
    if strain:
        s = " ".join([s,strain])
    return s

def combine_references(ref_dict, refs, sep, output):
    with open(output, "w") as f:
        for key in ref_dict:
            ref_seqs = [str(refs[acc].seq) for acc in ref_dict[key] if acc]
            f.write(">%s\n%s\n" %(key.replace(" ","_"), sep.join(ref_seqs)))

def collect_refs(ref_file, summary_file, sep, output):
    ref_dict = defaultdict(list)
    ref_index = {}
    reserve = defaultdict(list)
    refs = SeqIO.parse(ref_file, 'fasta')
    accessions = parse_accessions(summary_file)

    for record in refs:
        strain = None
        if "strain" in record.description:
            strain = record.description.split("strain ")[1]
            strain = "(%s)" %strain
        d = record.description.split(",")[0]
        prefix = d.split(" segment ")[0]
        acc = prefix.split()[0]
        name = parse_name(prefix, strain)

        if not accession_in_list(acc, accessions):
            continue

        ref_index[acc] = record
        ref_dict[name].append(acc)

    for name in ref_dict:
        ref_dict[name].sort(key=lambda s: ref_index[s].description, reverse=True)
        ref_dict[name].sort(key=lambda s: len(ref_index[s]), reverse=True)

    combine_references(ref_dict, ref_index, sep, output)


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-r","--reference", help="Reference FASTA", required=True)
    parser.add_argument("-s", "--summary", help="Summary TSV of references downloaded", required=True)
    parser.add_argument("-o", "--output", help="Outfile name", default="filtered_references.fa")

    parser.add_argument("--segment_sep", help="String to use as seperator when concatenating segments", default="NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN")
    args = parser.parse_args()

    collect_refs(args.reference, args.summary, args.segment_sep, args.output)

if __name__ == "__main__":
    main()
