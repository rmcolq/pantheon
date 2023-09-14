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
    with open("combined_references.fa", "w") as f:
        for key in ref_dict:
            ref_seqs = [str(refs[acc].seq) for acc in ref_dict[key] if acc]
            f.write(">%s\n%s\n" %(key.replace(" ","_"), "".join(ref_seqs)))

def collect_refs(ref_file, summary_file, sep):
    ref_dict = {}
    reserve = defaultdict(list)
    refs = SeqIO.index(ref_file, 'fasta')
    accessions = parse_accessions(summary_file)

    for acc in refs:
        strain = None
        if "strain" in refs[acc].description:
            strain = refs[acc].description.split("strain ")[1]
            strain = "(%s)" %strain
        d = refs[acc].description.split(",")[0]
        prefix = d.split(" segment ")[0]
        acc = prefix.split()[0]
        name = parse_name(prefix, strain)

        if not accession_in_list(acc, accessions):
            continue

        if "segment" in d:
            try:
                seg_num = d.split(" segment ")[1]
                seg_num = int(seg_num.split()[0].split(',')[0])
            except:
                print("Can't parse the segment number from %s so ignoring" %refs[acc].description)
                continue

            if name not in ref_dict:
                ref_dict[name] = []
            if type(ref_dict[name]) == str:
                reserve[name].append(ref_dict[name])
                ref_dict[name] = []
            if len(ref_dict[name]) < seg_num:
                ref_dict[name] = ref_dict[name] + [None]*(seg_num - len(ref_dict[name]))

            ref_dict[name][seg_num-1] = acc

        elif name in ref_dict:
            reserve[name].append(acc)
        else:
            ref_dict[name] = acc

    for name in reserve:
        if type(ref_dict[name]) == str:
            ref_dict[name] = [ref_dict[name]]
        for i in range(len(ref_dict[name])):
            if not ref_dict[name][i] and len(reserve[name]) > 0:
                ref_dict[name][i] = reserve[name].pop()
        try:
            ref_dict[name].extend(reserve[name])
        except:
            print("ERROR", name, ref_dict[name], reserve[name])

    for name in ref_dict:
        print(name, ref_dict[name])
    combine_references(ref_dict, refs, sep)


def map_to_refs(query):
    a = mp.Aligner("combined_references.fa")  # load or build index
    if not a:
        raise Exception("ERROR: failed to load/build index")

    for name, seq, qual in mp.fastx_read(query): # read a fasta/q sequence
        for hit in a.map(seq): # traverse alignments
            print("{}\t{}\t{}\t{}\t{}".format(name, hit.ctg, hit.r_st, hit.r_en, hit.cigar_str))
            break


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("-r","--reference", help="Reference FASTA", required=True)
    parser.add_argument("-s", "--summary", help="Summary TSV of references downloaded", required=True)

    parser.add_argument("-q","--query", help="Read FASTQ.GZ", required=True)

    parser.add_argument("--segment_sep", help="String to use as seperator when concatenating segments", default="NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN")

    args = parser.parse_args()

    collect_refs(args.reference, args.summary, args.segment_sep)
    map_to_refs(args.query)

if __name__ == "__main__":
    main()
