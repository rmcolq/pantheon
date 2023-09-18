#!/usr/bin/env python

#################
#
#   NOTE that the read counts/percentages are taken from the bracken reestimated file, but the reads themselves
#    are extracted based on kraken classifications because bracken does not provide classifications at the read level
#
###############

import sys
import os
import gzip
from Bio import SeqIO
import argparse
import json
from datetime import datetime
from collections import defaultdict


def mean(l):
    if len(l) == 0:
        return 0
    return sum(l) / len(l)


def median(l):
    if len(l) % 2 == 0:
        i = (len(l)) / 2
    else:
        i = (len(l) + 1) / 2
    i = int(i)
    l = sorted(l)
    return l[i]

def load_from_taxonomy(taxonomy_dir):
    taxonomy = os.path.join(taxonomy_dir, "nodes.dmp")
    parents = {}
    children = defaultdict(list)
    try:
        with open(taxonomy, "r") as f:
            for line in f:
                fields = line.split("\t|\t")
                tax_id, parent_tax_id = fields[0], fields[1]
                parents[tax_id] = parent_tax_id
                children[parent_tax_id].append(tax_id)
    except:
        sys.stderr.write(
            "ERROR: Could not find taxonomy nodes.dmp file in %s" %taxonomy_dir)
        sys.exit(2)
    return parents, children

def translate_names(taxonomy_dir, taxon_names):
    taxon_ids = {}
    taxonomy = os.path.join(taxonomy_dir, "names.dmp")
    try:
        with open(taxonomy, "r") as f:
            for line in f:
                fields = line.strip().split("\t|")
                taxon_id, name, name_type = fields[0].strip(), fields[1].strip(), fields[3].strip()
                if name not in taxon_names and taxon_id not in taxon_names:
                    continue
                if name_type != "scientific name":
                    continue
                taxon_ids[taxon_id] = name

    except:
        sys.stderr.write(
            "ERROR: Could not find taxonomy names.dmp file in %s" %taxonomy_dir)
        sys.exit(2)
    print("Translations", taxon_ids)
    return taxon_ids


def parse_depth(name):
    parse_name = name.split(" ")
    depth = 0
    for i in parse_name:
        if i != "":
            break
        depth += 1
    depth = int(depth / 2)
    return depth


def load_report_file(report_file, max_human=None):
    entries = {}
    # parses a kraken or bracken file
    with open(report_file, "r") as f:
        for line in f:
            if line.startswith("% of Seqs"):
                continue
            try:
                (
                    percentage,
                    num_clade_root,
                    num_direct,
                    raw_rank,
                    ncbi,
                    name,
                ) = line.strip().split("\t")
            except:
                (
                    percentage,
                    num_clade_root,
                    num_direct,
                    ignore1,
                    ignore2,
                    raw_rank,
                    ncbi,
                    name,
                ) = line.strip().split("\t")
            percentage = float(percentage)
            num_clade_root = int(num_clade_root)
            num_direct = int(num_direct)
            if num_direct > num_clade_root:
                num_direct, num_clade_root = num_clade_root, num_direct
            name = name.strip()
            rank = raw_rank[0]

            entries[ncbi] = {
                "percentage": percentage,
                "count": num_direct,
                "count_descendants": num_clade_root,
                "raw_rank": raw_rank,
                "rank": rank,
                "name": name
            }

    sys.stdout.write("FOUND %i TAXA IN BRACKEN REPORT\n" % len(entries))
    return entries

def get_taxon_list(taxon, include_parents, parents, include_children, children):
    taxon_list = [taxon]
    if include_parents:
        lookup = taxon
        while (lookup in parents and lookup != "1") :
            lookup = parents[lookup]
            if lookup != "1":
                taxon_list.append(lookup)

    if include_children:
        lookup = [taxon]
        while len(lookup) > 0:
            child = lookup.pop()
            if child != taxon:
                taxon_list.append(child)
            lookup.extend(children[child])
    return taxon_list

def infer_entry(taxon, name, children, report_entries, parent_rank):
    lookup = children[taxon]
    report_entries[taxon] = {
                            "percentage": 0,
                            "count": 0,
                            "count_descendants": 0,
                            "raw_rank": "",
                            "rank": "",
                            "name": name
                            }
    child_ranks = []
    for child in children[taxon]:
        if child in report_entries:
            #print(report_entries[child])
            report_entries[taxon]["count_descendants"] += report_entries[child]["count_descendants"]
            report_entries[taxon]["percentage"] += report_entries[child]["percentage"]
            child_ranks.append(report_entries[child]["raw_rank"])
    #print(child_ranks)
    if len(child_ranks) > 0:
        child_rank = min(child_ranks)
        #print(child_rank)
        if len(child_rank) == 2:
            rank = child_rank[0]
        else:
            rank = parent_rank[child_rank]
        #print(rank)
        report_entries[taxon]["rank"] = rank


def get_taxon_id_lists(
    report_entries,
    parents,
    children,
    names=[],
    target_ranks=[],
    min_count=None,
    min_count_descendants=None,
    min_percent=None,
    top_n=None,
    include_parents=False,
    include_children=False,
):
    lists_to_extract = {}
    for taxon in report_entries:
        entry = report_entries[taxon]
        if len(target_ranks) > 0 and entry["rank"] not in target_ranks:
            continue
        if min_count and entry["count"] < min_count:
            continue
        if min_count_descendants and entry["count_descendants"] < min_count_descendants:
            continue
        if min_percent and entry["percentage"] < min_percent:
            continue
        if len(names) > 0 and entry["name"] not in names and taxon not in names:
            continue

        lists_to_extract[taxon] = get_taxon_list(taxon, include_parents, parents, include_children, children)

    for id in names:
        if id not in lists_to_extract:
            lists_to_extract[id] = [id]

    sys.stdout.write("SELECTED %i TAXA TO EXTRACT\n" % len(lists_to_extract))

    if top_n and len(lists_to_extract) > top_n:
        X = list(lists_to_extract.keys())
        Y = [report_entries[x]["percentage"] for x in X]
        ordered = [x for _, x in sorted(zip(Y, X))]
        to_delete = ordered[top_n:]
        for taxon in to_delete:
            del lists_to_extract[taxon]
        sys.stdout.write("REDUCED TO %i TAXA TO EXTRACT\n" % len(lists_to_extract))

    return lists_to_extract


def check_read_files(reads):
    if reads[-3:] == ".gz":
        read_file = gzip.open(reads, "rt")
        zipped = True
    else:
        read_file = open(reads, "rt")
        zipped = False
    first = read_file.readline()
    if len(first) == 0:
        sys.stderr.write("ERROR: sequence file's first line is blank\n")
        sys.exit(1)
    if first[0] == ">":
        filetype = "fasta"
    elif first[0] == "@":
        filetype = "fastq"
    else:
        sys.stderr.write("ERROR: sequence file must be FASTA or FASTQ\n")
        sys.exit(1)
    return filetype, zipped


def parse_kraken_assignment_line(line):
    line_vals = line.strip().split("\t")
    if len(line_vals) < 5:
        return -1, ""
    if "taxid" in line_vals[2]:
        temp = line_vals[2].split("taxid ")[-1]
        tax_id = temp[:-1]
    else:
        tax_id = line_vals[2]

    read_id = line_vals[1]
    if tax_id == "A":
        tax_id = 81077
    else:
        tax_id = tax_id
    return tax_id, read_id


def extract_taxa(
    report_entries,
    lists_to_extract,
    kraken_assignment_file,
    reads1,
    reads2,
    prefix
):
    # open read files
    filetype, zipped = check_read_files(reads1)
    s_file1 = SeqIO.index(reads1, filetype)
    if reads2:
        s_file2 = SeqIO.index(reads2, filetype)

    # open output files
    outfile_handles = {}
    out_counts = {}
    quals = {}
    lens = {}

    keys = {}
    for taxon in lists_to_extract:
        out_counts[taxon] = 0
        quals[taxon] = []
        lens[taxon] = []
        sys.stdout.write(
            "INCLUDING PARENTS/CHILDREN, HAVE %i TAXA TO INCLUDE IN READ FILES for %s\n"
            % (len(lists_to_extract[taxon]), taxon)
        )
        if len(lists_to_extract[taxon]) == 0:
            continue

        out_prefix = "_".join([prefix, taxon])
        for key in lists_to_extract[taxon]:
            if key not in keys:
                keys[key] = []
            keys[key].append(taxon)
        if reads2:
            outfile_handles[taxon] = {
                1: open("%s_1.%s" % (out_prefix, filetype), "w"),
                2: open("%s_2.%s" % (out_prefix, filetype), "w"),
            }
            print(
                "opening %s_1.%s and %s_2.%s"
                % (out_prefix, filetype, out_prefix, filetype)
            )
        else:
            outfile_handles[taxon] = open("%s.%s" % (out_prefix, filetype), "w")
            print("opening %s.%s" % (out_prefix, filetype))


    with open(kraken_assignment_file, "r") as kfile:
        for line in kfile:
            tax_id, read_id = parse_kraken_assignment_line(line)
            if tax_id in keys:
                if reads2:
                    if read_id in s_file1 and read_id in s_file2:
                        read1 = s_file1[read_id]
                        read2 = s_file2[read_id]
                    else:
                        sys.stderr.write(
                            "ERROR: read id %s not found in read files\n" % read_id
                        )
                        sys.exit(1)

                    for taxon in keys[tax_id]:
                        SeqIO.write(read1, outfile_handles[taxon][1], filetype)
                        SeqIO.write(read2, outfile_handles[taxon][2], filetype)
                        out_counts[taxon] += 2
                        quals[taxon].append(
                            median(read1.letter_annotations["phred_quality"])
                        )
                        quals[taxon].append(
                            median(read2.letter_annotations["phred_quality"])
                        )
                        lens[taxon].append(len(read1))
                        lens[taxon].append(len(read2))

                else:
                    if read_id in s_file1:
                        read = s_file1[read_id]
                    else:
                        sys.stderr.write(
                            "ERROR: read id %s not found in read file\n" % read_id
                        )
                        sys.exit(1)

                    for taxon in keys[tax_id]:
                        SeqIO.write(read, outfile_handles[taxon], filetype)
                        out_counts[taxon] += 1
                        quals[taxon].append(
                            median(read.letter_annotations["phred_quality"])
                        )
                        lens[taxon].append(len(read))
    if reads2:
        for handle_dict in outfile_handles:
            outfile_handles[handle_dict][1].close()
            outfile_handles[handle_dict][2].close()
    else:
        for handle in outfile_handles:
            if outfile_handles[handle]:
                outfile_handles[handle].close()

    summary = []
    for taxon in lists_to_extract:
        out_prefix = "_".join([prefix, taxon])
        if reads2:
            summary.append(
                {
                    "human_readable": report_entries[taxon]["name"],
                    "taxon": taxon,
                    "tax_level": report_entries[taxon]["rank"],
                    "report_count": report_entries[taxon]["count_descendants"],
                    "filenames": [
                        "%s_1.%s" % (out_prefix, filetype),
                        "%s_2.%s" % (out_prefix, filetype),
                    ],
                    "qc_metrics": {
                        "num_reads": out_counts[taxon],
                        "avg_qual": mean(quals[taxon]),
                        "mean_len": mean(lens[taxon]),
                    },
                }
            )
        else:
            summary.append(
                {
                    "human_readable": report_entries[taxon]["name"],
                    "taxon": taxon,
                    "tax_level": report_entries[taxon]["rank"],
                    "report_count": report_entries[taxon]["count_descendants"],
                    "filenames": [
                        "%s.%s" % (out_prefix, filetype),
                    ],
                    "qc_metrics": {
                        "num_reads": out_counts[taxon],
                        "avg_qual": mean(quals[taxon]),
                        "mean_len": mean(lens[taxon]),
                    },
                }
            )

    return out_counts, summary

def extend_summary(summary, report_entries):
    for taxon in ['9606', '0']:
        if taxon in report_entries:
            summary.append(
                {
                    "human_readable": report_entries[taxon]["name"],
                    "taxon": taxon,
                    "tax_level": report_entries[taxon]["rank"],
                    "report_count": report_entries[taxon]["count_descendants"],
                    "report_percentage": report_entries[taxon]["percentage"],
                    "filenames": [],
                    "qc_metrics": {},
                }
            )

    return summary

# Main method
def main():
    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-k",
        dest="kraken_assignment_file",
        required=True,
        help="Kraken assignment file to parse",
    )
    parser.add_argument(
        "-r",
        dest="report_file",
        required=False,
        help="Kraken or Bracken file of taxon relationships and quantities",
    )
    parser.add_argument(
        "-t",
        dest="taxonomy",
        required=False,
        help="Taxonomy directory containing the nodes.dmp file. If not provided will infer from report file but this may lead to fewer reads extracted",
    )
    parser.add_argument(
        "-s",
        "-s1",
        "-1",
        dest="reads1",
        required=True,
        help="FASTA/FASTQ File containing the raw reads.",
    )
    parser.add_argument(
        "-s2",
        "-2",
        dest="reads2",
        default="",
        help="2nd FASTA/FASTQ File containing the raw reads (paired).",
    )

    parser.add_argument(
        "-p",
        "--prefix",
        dest="prefix",
        default="",
        help="Prefix for output files",
    )

    parser.add_argument(
        "--taxid",
        dest="taxid",
        required=False,
        nargs="*",
        default=[],
        help="List of taxonomy ID[s] or names to extract (space-delimited) - each to their own file",
    )
    parser.add_argument(
        "--rank", dest="rank", required=False, nargs="*", help="Rank(s) to extract"
    )
    parser.add_argument(
        "--max_human",
        dest="max_human",
        required=False,
        type=int,
        help="Maximum human reads to allow",
    )
    parser.add_argument(
        "--min_count",
        dest="min_count",
        required=False,
        type=int,
        help="Minimum direct read count",
    )
    parser.add_argument(
        "--min_count_descendants",
        dest="min_count_descendants",
        required=False,
        type=int,
        help="Minimum read count at taxon level or descendants",
    )
    parser.add_argument(
        "--min_percent",
        dest="min_percent",
        required=False,
        type=float,
        help="Minimum percentage of reads e.g 4",
    )
    parser.add_argument(
        "--n",
        dest="top_n",
        required=False,
        type=int,
        help="Maximum number of taxa to extract (top n)",
    )
    parser.add_argument(
        "--include_parents",
        dest="include_parents",
        required=False,
        action="store_true",
        default=False,
        help="Include reads classified at parent levels of the specified taxids",
    )
    parser.add_argument(
        "--include_children",
        dest="include_children",
        required=False,
        action="store_true",
        default=False,
        help="Include reads classified more specifically than the specified taxids",
    )
    parser.set_defaults(append=False)

    args = parser.parse_args()

    if not args.report_file:
        sys.stderr.write(
            "ERROR: require at least one report file from bracken or kraken\n"
        )
        sys.exit(1)

    # Start Program
    now = datetime.now()
    time = now.strftime("%m/%d/%Y, %H:%M:%S")
    sys.stdout.write("PROGRAM START TIME: " + time + "\n")

    rank_dict = {
        "kingdom": "K",
        "domain": "D",
        "phylum": "P",
        "class": "C",
        "order": "O",
        "family": "F",
        "genus": "G",
        "species": "S",
        "K": "K",
        "D": "D",
        "P": "P",
        "C": "C",
        "O": "O",
        "F": "F",
        "G": "G",
        "S": "S",
    }
    parent_rank = {
            "K": "R",
            "D": "K",
            "P": "D",
            "C": "P",
            "O": "C",
            "F": "O",
            "G": "F",
            "S": "G",
        }
    if args.rank:
        target_ranks = [rank_dict[r] for r in args.rank]
    else:
        target_ranks = []
    print(target_ranks)

    # get taxids to extract
    all_ids = []
    for id in args.taxid:
        all_ids.extend(id.split(","))
    id_dict = translate_names(args.taxonomy, all_ids)
    all_ids = id_dict.keys()
    sys.stdout.write("Restricting to taxa [%s]\n" %(",".join(all_ids)))

    parent, children = None, None
    parent, children = load_from_taxonomy(args.taxonomy)

    report_entries = load_report_file(args.report_file, args.max_human)
    for taxon in all_ids:
        if taxon not in report_entries:
            infer_entry(taxon, id_dict[taxon], children, report_entries, parent_rank)
    if '9606' not in report_entries:
        infer_entry('9606', "Homo sapiens", children, report_entries, parent_rank)

    lists_to_extract = get_taxon_id_lists(
        report_entries,
        parent,
        children,
        names=all_ids,
        target_ranks=target_ranks,
        min_count=args.min_count,
        min_count_descendants=args.min_count_descendants,
        min_percent=args.min_percent,
        top_n=args.top_n,
        include_parents=args.include_parents,
        include_children=args.include_children,
        )

    out_counts, summary = extract_taxa(
        report_entries,
        lists_to_extract,
        args.kraken_assignment_file,
        args.reads1,
        args.reads2,
        args.prefix
    )

    now = datetime.now()
    time = now.strftime("%m/%d/%Y, %H:%M:%S")
    sys.stdout.write("PROGRAM END TIME: " + time + "\n")

    sys.stdout.write("READ COUNTS: \n")

    for taxon in out_counts:
        sys.stdout.write("%s: %i\n" % (taxon, out_counts[taxon]))

    summary = extend_summary(summary, report_entries)
    out_summary = "_".join([args.prefix, "summary.json"])
    with open(out_summary, "w") as f:
        json.dump(summary, f)

    out_human = "_".join([args.prefix, "human.txt"])
    with open(out_human, "w") as f:
        f.write("%i" %report_entries['9606']["count_descendants"])

    sys.exit(0)


if __name__ == "__main__":
    main()
