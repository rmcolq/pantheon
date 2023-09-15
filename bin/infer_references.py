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

def combine_references(ref_dict, refs, sep):
    with open("intermediate_references.fa", "w") as f:
        for key in ref_dict:
            ref_seqs = [str(refs[acc].seq) for acc in ref_dict[key] if acc]
            f.write(">%s\n%s\n" %(key.replace(" ","_"), sep.join(ref_seqs)))

def collect_refs(ref_file, summary_file, sep):
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

    combine_references(ref_dict, ref_index, sep)


def map_to_refs(query):
    counts = defaultdict(int)
    a = mp.Aligner("intermediate_references.fa")  # load or build index
    if not a:
        raise Exception("ERROR: failed to load/build index")

    for name, seq, qual in mp.fastx_read(query): # read a fasta/q sequence
        for hit in a.map(seq): # traverse alignments
            counts[hit.ctg] += 1
            #print("{}\t{}\t{}\t{}\t{}".format(name, hit.ctg, hit.r_st, hit.r_en, hit.cigar_str))
            break
    return counts

def output_references(counts, min_count):
    refs = SeqIO.index("intermediate_references.fa", 'fasta')
    with open("filtered_references.fa", "w") as f:
        for ref in counts:
            print("Reference %s has %i mapped reads" %(ref, counts[ref]))
            if counts[ref] > min_count:
                SeqIO.write(refs[ref], f, 'fasta')
            else:
                print("Excluded: this is smaller than min_count %i" %min_count)

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-r","--reference", help="Reference FASTA", required=True)
    parser.add_argument("-s", "--summary", help="Summary TSV of references downloaded", required=True)

    parser.add_argument("-q","--query", help="Read FASTQ.GZ", required=True)

    parser.add_argument("--segment_sep", help="String to use as seperator when concatenating segments", default="NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN")
    parser.add_argument("--min_count", help="Minimum number of mapped reads to continue with reference", default=50, type=int)
    args = parser.parse_args()

    collect_refs(args.reference, args.summary, args.segment_sep)
    counts = map_to_refs(args.query)
    output_references(counts, args.min_count)

if __name__ == "__main__":
    main()
