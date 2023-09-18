#!/usr/bin/env python

import mappy as mp
import argparse
from Bio import SeqIO
from collections import defaultdict
import re

def map_to_refs(query, reference):
    counts = defaultdict(int)
    a = mp.Aligner(reference)  # load or build index
    if not a:
        raise Exception("ERROR: failed to load/build index")

    for name, seq, qual in mp.fastx_read(query): # read a fasta/q sequence
        for hit in a.map(seq): # traverse alignments
            counts[hit.ctg] += 1
            #print("{}\t{}\t{}\t{}\t{}".format(name, hit.ctg, hit.r_st, hit.r_en, hit.cigar_str))
            break
    return counts

def output_references(counts, min_count, references, output, taxid, barcode):
    refs = SeqIO.index(references, 'fasta')
    summary = open("summary.csv","w")
    summary.write("TaxonID,BarcodeID,Reference,ReadCount\n")
    with open(output, "w") as f:
        for ref in counts:
            summary.write("%s,%s,%s,%i\n" %(taxid, barcode, ref, counts[ref]))
            print("Reference %s has %i mapped reads" %(ref, counts[ref]))
            if counts[ref] > min_count:
                SeqIO.write(refs[ref], f, 'fasta')
            else:
                print("Excluded: this is smaller than min_count %i" %min_count)
    summary.close()

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-r","--reference", help="Reference FASTA", required=True)
    parser.add_argument("-q","--query", help="Read FASTQ.GZ", required=True)
    parser.add_argument("-o","--output", help="FASTA outfile", default="subset_references.fa")

    parser.add_argument("--taxid", help="Reference TaxonID", required=True)
    parser.add_argument("--barcode", help="Barcode", required=True)

    parser.add_argument("--min_count", help="Minimum number of mapped reads to continue with reference", default=50, type=int)
    args = parser.parse_args()

    counts = map_to_refs(args.query, args.reference)
    output_references(counts, args.min_count, args.reference, args.output, args.taxid, args.barcode)

if __name__ == "__main__":
    main()
