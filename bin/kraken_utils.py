from collections import defaultdict
import pandas as pd
import glob
import os
import sys

class KrakenEntry:
    def __init__(self, row=None, domain=None, hierarchy=[]):
        self.taxon_id = 0
        self.name = "unclassified"
        self.rank = "U"
        self.depth = 0
        self.count = 0 # inclusive count
        self.ucount = 0 # unique_count
        self.domain = domain
        self.parent = None
        self.children = set()
        self.sibling_rank = 0
        self.hierarchy = hierarchy
        if row is not None:
            self.add_row(row)

    def print(self):
        print(f"{self.taxon_id},{self.name},{self.rank},{self.depth},{self.count},{self.ucount},{self.domain},{self.parent},{self.children},{self.sibling_rank},{self.hierarchy}")

    def parse_depth(self, name):
        parse_name = name.split(" ")
        depth = 0
        for i in parse_name:
            if i != "":
                break
            depth += 1
        depth = int(depth / 2)
        return depth

    def add_row(self, row):
        self.taxon_id = row["Taxonomy ID"]
        self.name = row["Scientific Name"].strip()
        self.depth = self.parse_depth(row["Scientific Name"])
        self.rank = row["Rank"]
        self.count = row["Clades"] # inclusive count
        self.ucount = row["Taxonomies"] # unique_count
        self.hierarchy = self.hierarchy[:self.depth]

    def add_parent(self, parent):
        self.parent = parent

    def add_child(self, child):
        self.children.add(child)

    def set_sibling_rank(self, rank):
        self.sibling_rank = rank

    def to_taxon_dict(self):
        return {"taxon_id": self.taxon_id,
                "name": self.name,
                "rank": self.rank,
                "depth": self.depth,
                "domain": self.domain,
                "parent": self.parent,
                "children": self.children,
                "sibling_rank": self.sibling_rank,
                "hierarchy": self.hierarchy}

class KrakenReport:
    def __init__(self, file_name=None, host_id=9606):
        self.entries = defaultdict(KrakenEntry)
        self.total = 0
        self.unclassified = 0
        self.classified = 0
        self.host = 0
        self.domains = defaultdict(int)
        self.num_taxa = 0
        if file_name:
            self.load_df(file_name)
            self.unclassified = self.entries[0].count
            self.classified = self.entries[1].count
            self.total = self.classified + self.unclassified
            self.host = self.entries[host_id].count
            self.num_taxa = len(self.entries)
            self.domains = self.get_domains()

    def add_parent_child(self, parent_id, child_id):
        self.entries[child_id].add_parent(parent_id)
        self.entries[parent_id].add_child(child_id)

    def set_sibling_ranks(self):
        for entry_id, entry in self.entries.items():
            if entry.sibling_rank > 0 or entry.parent in [None, 1, 131567]:
                continue
            if entry.rank in ["D", "R", "R1"]:
                entry.set_sibling_rank(1)
            elif len(self.entries[entry.parent].children) == 1:
                entry.set_sibling_rank(1)
            else:
                sibling_dict = {i:self.entries[i].count for i in self.entries[entry.parent].children}
                sorted_counts = sorted(sibling_dict.values(), reverse=True)
                for i,c in sibling_dict.items():
                    rank = sorted_counts.index(c) + 1
                    self.entries[i].set_sibling_rank(rank)

    def check_sibling_ranks(self):
        for entry_id in self.entries:
            if self.entries[entry_id].sibling_rank == 0:
                print(entry_id)

    def load_df(self, file_name):
        if not os.path.isfile(file_name):
            return
        with open(file_name, 'r') as f:
            names = f.readline().strip().split('\t')
            header = 0
            if not names[0].startswith("%"):
                names = ["% of Seqs", "Clades", "Taxonomies", "Rank", "Taxonomy ID", "Scientific Name"]
                header = None
        df = pd.read_csv(file_name, sep="\t", skipinitialspace = False, header=header, names=names)
        hierarchy = []
        domain = None
        for i,row in df.iterrows():
            try:
                if row["Rank"] == "D":
                    domain = row["Scientific Name"].strip()
                    self.domains[domain] = row["Taxonomy ID"]
                entry = KrakenEntry(row=row, domain=domain, hierarchy=hierarchy)

            except:
                print(f"Found badly formatted row:\n{row}")
                print(f"Quitting load of {file_name}")
                break

            self.entries[entry.taxon_id] = entry
            hierarchy = entry.hierarchy.copy()
            if len(hierarchy) > 0:
                self.add_parent_child(hierarchy[-1], entry.taxon_id)
            if entry.taxon_id > 0:
                hierarchy.append(entry.taxon_id)
        self.set_sibling_ranks()
        #self.check_sibling_ranks()

    def get_host_count(self, host_id=9606):
        self.host = self.entries[host_id].count
        return self.host

    def get_classified(self):
        self.classified = sum([self.entries[taxon_id].ucount for taxon_id in self.entries])
        return self.classified

    def get_domains(self):
        self.domains = set([self.entries[taxon_id].domain for taxon_id in self.entries])
        return self.domains

    def restrict(self, list_taxa=[], include_children=False):
        print("Restricting taxa")
        if len(list_taxa) == 0:
            return

        all_taxa = set(list_taxa)
        if include_children:
            to_expand = set(list_taxa)
            while len(to_expand) > 0:
                taxon_id = to_expand.pop()
                all_taxa.add(taxon_id)
                to_expand.add(self.entries[taxon_id].children)

        reject_taxa = [taxon_id for taxon_id in self.entries if taxon_id not in all_taxa]
        for taxon_id in reject_taxa:
            del self.entries[taxon_id]

        self.classified = get_classified()
        self.total = self.classified + self.unclassified
        self.domains = self.get_domains()
        self.num_taxa =len(self.entries)

    def get_percentage(self, taxon_id, domain=None):

        if domain and self.entries[taxon_id].domain!=domain:
            return 0.0

        denominator = self.classified
        if domain:
            denominator = self.entries[self.domains[domain]].count

        count=self.entries[taxon_id].count
        return float(count)/denominator

    def to_source_target_df(self, max_rank=None, domain=None, trace_ids=[]):
        records = []
        ignore = set()
        skip = set()
        for entry_id, entry in self.entries.items():

            # we don't want the connections to cellular organisms, root etc
            if entry.sibling_rank == 0:
                continue

            # filter by domain where required
            if domain and entry.domain!=domain:
                continue

            # filter by rank when specified
            if max_rank:
                if entry.sibling_rank > max_rank:
                    ignore.add(entry_id)
                    if entry_id in trace_ids:
                        print(entry_id, 1)
                    continue
                elif entry.parent in ignore:
                    ignore.add(entry_id)
                    if entry_id in trace_ids:
                        print(entry_id, 2)
                    continue

            # filter if an intermediate rank
            if entry.rank not in ["K","D","D1","D2","P","C","O","F","G","S", "S1", "S2"]:
                skip.add(entry_id)
                if entry_id in trace_ids:
                    print(entry_id, 3)
                continue

            index = 1
            while index < len(entry.hierarchy) and entry.hierarchy[-index] in skip:
                index += 1
            source_id = entry.hierarchy[-index]
            if entry_id in trace_ids:
                print(entry_id, 4, source_id)
            records.append({"source":self.entries[source_id].name, "target": entry.name, "value":entry.count, "percentage":self.get_percentage(entry_id, domain=domain)})

        all_links = pd.DataFrame.from_records(records)
        all_links.to_csv("source_target.csv")

        print(len(records), len(ignore), len(skip))

        return records

    def to_df(self, sample_id="sample_id", ranks=[], type="count"):
        if not ranks or len(ranks) == 0:
            taxon_ids =  [e for e in self.entries.keys()]
        else:
            taxon_ids =  [e for e in self.entries.keys() if self.entries[e].rank in ranks]
        if type=="count":
            return pd.DataFrame({sample_id:[self.entries[e].count for e in taxon_ids]}, index=taxon_ids)
        elif type=="ucount":
            return pd.DataFrame({sample_id:[self.entries[e].ucount for e in taxon_ids]}, index=taxon_ids)
        else:
            sys.exit(f"Required type must be in ['count', 'ucount'], user specified {type}")

class KrakenCombined:
    def __init__(self, list_kraken_reports=[], dict_kraken_reports={}, list_taxa=[], metadata=None, counts=None, taxa=None, ranks=[]):
        self.metadata = pd.DataFrame(columns=["group","pass","classified","unclassified","host","total"])
        if metadata:
            self.metadata = pd.read_csv(metadata, index_col=0)

        self.counts = pd.DataFrame()
        if counts:
            self.counts = pd.read_csv(counts, index_col=0)

        self.ucounts = pd.DataFrame()

        self.taxa = pd.DataFrame(columns=["name","rank","depth","domain","parent","children","sibling_rank","hierarchy"])
        if taxa:
            self.taxa = pd.read_csv(taxa, index_col=0)

        self.load_from_list_kraken_reports(list_kraken_reports, list_taxa, ranks)
        self.load_from_dict_kraken_reports(dict_kraken_reports, list_taxa, ranks)

    def load_from_list_kraken_reports(self, list_kraken_reports, list_taxa=[], ranks=[]):
        for file_name in list_kraken_reports:
            k = KrakenReport(file_name)
            if len(list_taxa) > 0:
                k.restrict(list_taxa, include_children=True)
            sample_id = file_name.split("/")[-1].split(".")[0]
            self.add_kraken_report_to_metadata(k, sample_id)
            self.add_kraken_report_to_counts(k, sample_id, ranks=ranks)
            self.add_kraken_report_to_taxa(k)

    def load_from_dict_kraken_reports(self, dict_kraken_reports, list_taxa=[], ranks=[]):
        for sample_id,file_name in dict_kraken_reports.items():
            k = KrakenReport(file_name)
            if len(list_taxa) > 0:
                k.restrict(list_taxa, include_children=True)
            self.add_kraken_report_to_metadata(k, sample_id)
            self.add_kraken_report_to_counts(k, sample_id, ranks=ranks)
            self.add_kraken_report_to_taxa(k)

    def add_kraken_report_to_metadata(self, kraken_report, sample_id):
        update = {"group": "default", "pass":True, "classified":kraken_report.classified, "unclassified":kraken_report.unclassified, "host": kraken_report.host, "total":kraken_report.total}
        self.metadata.loc[sample_id] = update

    def add_kraken_report_to_counts(self, kraken_report, sample_id, ranks=[]):
        self.counts = self.counts.join(kraken_report.to_df(sample_id, ranks=ranks, type="count"), how='outer')
        self.counts.fillna(0, inplace=True)
        self.counts = self.counts.astype('int64')

        self.ucounts = self.ucounts.join(kraken_report.to_df(sample_id, ranks=ranks, type="ucount"), how='outer')
        self.ucounts.fillna(0, inplace=True)
        self.ucounts = self.ucounts.astype('int64')

    def add_kraken_report_to_taxa(self, kraken_report):
        for taxon_id in kraken_report.entries:
            taxon_dict = kraken_report.entries[taxon_id].to_taxon_dict()
            if taxon_dict["rank"] == "U":
                continue
            self.taxa.loc[taxon_id] = taxon_dict

    def load_metadata(self, metadata_csv):
        with open(metadata_csv, "r") as f:
            for line in f:
                sample_id, group = line.strip().split(",")
                self.metadata.loc[sample_id, "group"] = group

    def normalize_counts(self, by="total"):
        self.normed_counts = self.counts
        for sample_id in self.normed_counts.columns:
            self.normed_counts[sample_id] = self.normed_counts[sample_id]/self.metadata.loc[sample_id, "total"]
        if by == "total":
            return
        x = self.normed_counts.T
        y = clr(x)

    def filter_samples(self, min_count=10):
        self.metadata["pass"] = self.metadata["pass"].fillna(False)
        num_taxa_in_samples = self.counts.astype(bool).sum(axis=0)
        for i,c in num_taxa_in_samples.items():
            self.metadata.loc[i,"pass"] = (c >= min_count)

    def get_taxa_info(self, samples=[], inplace=False):
        if len(samples) > 0:
            counts = self.counts[samples]
        else:
            counts = self.counts
        taxa = counts.astype(bool).sum(axis=1).to_frame("count_samples")
        taxa["percent_samples"] = taxa["count_samples"]/counts.shape[1]*100
        taxa["min"] = counts.min(axis=1)
        taxa["max"] = counts.max(axis=1)
        taxa["mean"] = counts.mean(axis=1)
        taxa["median"] = counts.median(axis=1)
        taxa["count_reads"] = counts.sum(axis=1)
        overall_total = counts.sum().sum()
        taxa["percent_reads"] = taxa["count_reads"]/overall_total*100
        if inplace:
            self.taxa = self.taxa.join(taxa)
        return taxa

    def get_filter(self, min_percent_of_samples=5, min_number_of_samples=2, min_percent_of_reads = 0.01, samples=[]):
        taxa = self.get_taxa_info(samples=samples)
        pass_percent = taxa["percent_samples"] >= min_percent_of_samples
        pass_count = taxa["count_samples"] >= min_number_of_samples
        pass_percent_taxa = taxa["percent_reads"] >= min_percent_of_reads
        return pass_percent & pass_count & pass_percent_taxa

    def filter_taxa(self, min_percent_of_samples=5, min_number_of_samples=2, min_percent_of_reads = 0.01, by_group=False):
        self.get_taxa_info(inplace=True)
        self.taxa["pass"] = False
        if by_group:
            groups = self.metadata['group'].unique()
            print(groups)
            for group in groups:
                samples = self.metadata[self.metadata["group"] == group].index.to_list()
                filter = self.get_filter(min_percent_of_samples=min_percent_of_samples, min_number_of_samples=min_number_of_samples, min_percent_of_reads=min_percent_of_reads, samples=samples)
                self.taxa["pass"] = self.taxa["pass"] | filter
        else:
            self.taxa["pass"] = self.get_filter(min_percent_of_samples=min_percent_of_samples, min_number_of_samples=min_number_of_samples, min_percent_of_reads=min_percent_of_reads)

    def get_taxa(self):
        return self.taxa[self.taxa["pass"]]

    def get_samples(self):
        return self.metadata[self.metadata["pass"]]

    def get_sample_info(self, sample_id):
        return self.metadata.loc[sample_id].to_dict()

    def save(self, prefix):
        self.metadata.to_csv(f"{prefix}.metadata.csv")
        self.counts.to_csv(f"{prefix}.counts.csv")
        self.taxa.to_csv(f"{prefix}.taxa.csv")