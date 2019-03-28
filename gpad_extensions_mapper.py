# from ontobio.io.gafparser import GafParser
from ontobio.io.gpadparser import GpadParser
from ontobio.ontol_factory import OntologyFactory
# from prefixcommons import curie_util
from ontobio.ecomap import EcoMap
from ontobio.util.go_utils import GoAspector
from ontobio.rdfgen.assoc_rdfgen import prefix_context
import json
import csv
import os
import argparse
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--filename")
parser.add_argument("-d", "--dir")
parser.add_argument("-o", "--out_file")
parser.add_argument("-l", "--leftovers_out_file")
parser.add_argument("-p", "--pattern")
parser.add_argument("-q", "--pattern_outfile")
parser.add_argument("-r", "--pattern_sourcefile")

ontology_prefixes = []
# with open("go_context.jsonld") as gcf:
#     context_prefixes = json.load(gcf)
# for k, v in context_prefixes['@context'].items():
for k, v in prefix_context.items():
    if v.startswith("http://purl.obolibrary.org/obo/"):
        ontology_prefixes.append(k)
ontology_prefixes.append("Pfam")

mod_prefixes = [
    "PomBase",
    "UniProtKB",
]

acceptable_evidence_codes = [
    "EXP",
    "IDA",
    "IPI",
    "IMP",
    "IGI",
    "IEP"
]
ecomap = EcoMap()

gaf_indices = {}
gpad_indices = {

}

class CachedGoAspector(GoAspector):

    def __init__(self, cache_filepath):
        GoAspector.__init__(self, OntologyFactory().create("go"))
        self.cache_filepath = cache_filepath
        self.aspect_lookup = {
            "F": [],
            "P": [],
            "C": []
        }
        aspect_file = Path(self.cache_filepath)
        if not aspect_file.is_file():
            print("Creating aspect_lookup cache file: {}".format(self.cache_filepath))
            self.write_cache()
        else:
            with open(cache_filepath) as af:
                try:
                    self.aspect_lookup = json.load(af)
                except json.decoder.JSONDecodeError:
                    print("Corrupt aspect_lookup cache file: {} - Recreating...".format(self.cache_filepath))

    def write_cache(self):
        with open(self.cache_filepath, "w+") as af:
            json.dump(self.aspect_lookup, af)

    def _is_biological_process(self, go_term):
        bp_root = "GO:0008150"
        return go_term in self.aspect_lookup["P"] or go_term == bp_root

    def _is_molecular_function(self, go_term):
        mf_root = "GO:0003674"
        return go_term in self.aspect_lookup["F"] or go_term == mf_root

    def _is_cellular_component(self, go_term):
        cc_root = "GO:0005575"
        return go_term in self.aspect_lookup["C"] or go_term == cc_root

    def go_aspect(self, term):
        if self._is_molecular_function(term):
            return 'F'
        elif self._is_biological_process(term):
            return 'P'
        elif self._is_cellular_component(term):
            return 'C'
        else:
            aspect = super(CachedGoAspector, self).go_aspect(term)
            if aspect:
                self.aspect_lookup[aspect].append(term)
                return aspect
            return None



def format_extensions(ext_json):
    # print(ext_json)
    ej = json.loads(ext_json)
    ret = []
    # print(ej["union_of"][0]["intersection_of"])
    for i in ej["union_of"][0]["intersection_of"]:
        ret.append("{}({})".format(i["property"], i["filler"]))
    # print(ret)
    return ",".join(ret)

def get_relation_and_term(ext):
    ext_parts = ext.split("(")
    relation = ext_parts[0]  # i["property"]
    try:
        ext_term = ext_parts[1].split(")")[0]  # i["filler"]
    except:
        print(ext_parts)
        # [' dppa4']
        # $ grep dppa4 resources/mgi.gpa.test.gpa
        # MGI	MGI:2141165	acts_upstream_of_or_within	GO:0010628	MGI:MGI:5913752|PMID:20699224	ECO:0000315	MGI:MGI:6108355		20190124	MGI	has_regulation_target(Dppa2, dppa4, piwil2),occurs_in(EMAPA:16036)
        ext_term = ext_parts[1].split(")")[0]  # i["filler"]
    return relation, ext_term

def filter_evi_codes(annots):
    filtered = []
    for a in annots:
        if ecomap.ecoclass_to_coderef(a[5])[0] in acceptable_evidence_codes:
            filtered.append(a)
    return filtered

def filter_has_extension(annots):
    filtered = []
    for a in annots:
        if a[10] != "":
            filtered.append(a)
    return filtered

def sum_combos(extension_counts, combo_list):
    cur_sum = 0
    for c in combo_list:
        if c in extension_counts:
            cur_sum += extension_counts[c]
    return cur_sum

def violates_combo_rule(extension_counts, list_of_combo_lists, max_allowed):
    for combo_list in list_of_combo_lists:
        if sum_combos(extension_counts, combo_list) > max_allowed:
            return True
    return False

def following_rules(extension_list, aspect):
    ext_counts = {}
    for e in extension_list:
        if e in ext_counts:
            ext_counts[e] += 1
        else:
            ext_counts[e] = 1

    function_singles_only = [
        "occurs_in(GO:C)",
        "occurs_in(CL)",
        "occurs_in(UBERON)",
        "occurs_in(EMAPA)",
        "has_input(geneID)",
        "has_input(CHEBI)",
        "has_direct_input(geneID)",
        "has_direct_input(CHEBI)",
        "happens_during(GO:P)",
        "part_of(GO:P)",
        "has_regulation_target(geneID)",
        "activated_by(CHEBI)",
        "inhibited_by(CHEBI)"
    ]
    # component_singles_only = [
    #     "occurs_in(GO:C)",
    #     "occurs_in(CL)",
    #     "occurs_in(UBERON)",
    #     "occurs_in(EMAPA)"
    # ]
    component_singles_only = [
        "part_of(GO:C)",
        "part_of(CL)",
        "part_of(UBERON)",
        "part_of(EMAPA)"
    ]
    process_singles_only = [
        "occurs_in(GO:C)",
        "occurs_in(CL)",
        "occurs_in(UBERON)",
        "occurs_in(EMAPA)",
        "has_input(geneID)",
        "has_input(CHEBI)",
        "has_direct_input(geneID)",
        "has_direct_input(CHEBI)",
        "part_of(GO:P)"
    ]
    combos_to_check_for = [
        ["occurs_in(UBERON)", "occurs_in(EMAPA)"]
    ]
    if aspect == "F":
        for s in function_singles_only:
            if s in ext_counts and ext_counts[s] > 1:
                return False
        for ek in ext_counts.keys():
            if ek not in function_singles_only:
                return False    # unrecognised relation-term combo
        combos_to_check_for.append(["has_input(geneID)", "has_input(CHEBI)",
                                    "has_direct_input(geneID)", "has_direct_input(CHEBI)"])
    if aspect == "C":
        for s in component_singles_only:
            if s in ext_counts and ext_counts[s] > 1:
                return False
        for ek in ext_counts.keys():
            if ek not in component_singles_only:
                return False    # unrecognised relation-term combo
    if aspect == "P":
        for s in process_singles_only:
            if s in ext_counts and ext_counts[s] > 1:
                return False
        for ek in ext_counts.keys():
            if ek not in process_singles_only:
                return False    # unrecognised relation-term combo
        combos_to_check_for.append(["has_input(geneID)", "has_input(CHEBI)",
                                    "has_direct_input(geneID)", "has_direct_input(CHEBI)"])
    if violates_combo_rule(ext_counts, combos_to_check_for, 1):
        return False
    return True

def filter_no_rules_broken(annots):
    filtered = []
    for a in annots:
        if not following_rules(a):
            filtered.append(a)
    return filtered


class ExtensionsMapper():
    def __init__(self):
        self.go_aspector = CachedGoAspector("resources/aspect_lookup.json")

    def extensions_list(self, intersection_extensions):
        ext_list = []
        # Assuming these extensions are already separated by comma
        intersection_extensions = self.dedupe_extensions(intersection_extensions)
        for i in intersection_extensions:
            relation, ext_term = i['property'], i['filler']
            term_prefix = ext_term.split(":")[0]
            # is prefix for ontology or mod?
            if term_prefix in ontology_prefixes:
                go_aspect = None
                ont = ""
                if term_prefix == "GO":
                    # need to find aspect of GO term
                    # go_aspect = get_go_aspect(ext_term)
                    go_aspect = self.go_aspector.go_aspect(ext_term)
                    if go_aspect is not None:
                        ont = "GO:" + go_aspect
                    else:
                        print("No aspect found for term: {}".format(ext_term))
                else:
                    ont = term_prefix
                ext_list.append("{}({})".format(relation, ont))
            else:
                ext_list.append("{}(geneID)".format(relation))
        # ordering of extension relations may be inconsistent so sort
        return sorted(ext_list, key=str.lower)

    # Do we need this? Does this work? WE TOTALLY NEED THIS!!! At least maybe as an entry point to find what bucket an
    # annotation falls in.
    # HANDLE ASSOCIATION OBJECT_EXTENSIONS STRUCTURE
    def annot_following_rules(self, annot, aspect):
        ext_list = []
        #TODO for a in annot.split("|"):
        # extensions = annot.split(",")
        # annot["object_extensions"]
        # ext_list = self.extensions_list(extensions)
        ext_list = self.extensions_list(annot)
        # Standardize key - ex:
        #   ["part_of(GO:aspect),part_of(UBERON:)"]
        return following_rules(ext_list, aspect)

    def dedupe_extensions(self, extensions):
        new_extensions = []
        for i in extensions:
            if i not in new_extensions:
                new_extensions.append(i)
        return new_extensions

d = [
    "GeneDB_tsetse",
    "gonuts",
    "reactome",
    "isoform",
    "goa_pdb",
    "goa_uniprot_all.gaf",

]


if __name__ == "__main__":
    args = parser.parse_args()

    filenames = []
    # data = []
    # fname = "/Users/ebertdu/go/go-pombase/gene_association.pombase-06-2018"
    # data = GafParser().parse(fname, skipheader=True)
    if args.filename is not None:
        filenames.append(args.filename)
        # data = GafParser().parse(args.filename, skipheader=True)
    elif args.dir is not None:
        for fname in os.listdir(args.dir):
            # print("Loading file:", fname)
            nono_in_fname = False
            for nono in d:
                if nono in fname:
                    nono_in_fname = True
            if fname.endswith(".tsv") or nono_in_fname:
                continue
            filenames.append(args.dir + fname)
            # data = data + GafParser().parse(fname, skipheader=True)

    # all_dict = {}
    extensions_mapper = ExtensionsMapper()
    gpad_parser = GpadParser()
    print("Creating extension dictionary...")
    ext_dict = {}
    ext_dict['F'] = {}
    ext_dict['P'] = {}
    ext_dict['C'] = {}
    for fname in filenames:
        with open(fname) as f:
            data = []
            print("Loading file:", fname)
            for l in f.readlines():
                if not l.startswith("!"):
                    parts = l.split("\t")
                    # if parts[15] != "" and parts[6] in acceptable_evidence_codes:
                    data.append(parts)
            print(len(data))
            data = filter_has_extension(data)
            print(len(data))
            data = filter_evi_codes(data)
            print("Total GPAD count:", len(data))

            for g in data:
                go_term = g[3]
                aspect = extensions_mapper.go_aspector.go_aspect(go_term)
                ontobio_extensions = gpad_parser._parse_full_extension_expression(g[10])
                ontobio_extensions = extensions_mapper.dedupe_extensions(ontobio_extensions)
                # ontobio_pattern = {
                #     'union_of': [
                #         {
                #             'intersection_of': [
                #                 {'property': 'part_of', 'filler': 'CL:0000678'},
                #                 {'property': 'part_of', 'filler': 'EMAPA:16525'}
                #             ]
                #         },
                #         {
                #               'intersection_of': [
                #                   {'property': 'part_of', 'filler': 'CL:0000678'},
                #                   {'property': 'part_of', 'filler': 'EMAPA:16525'}
                #               ]
                #         }
                #     ]
                # }
                for onto_ext in ontobio_extensions:
                    ext_list = extensions_mapper.extensions_list(onto_ext['intersection_of'])
                    if not following_rules(ext_list, aspect):
                        ext_key = ",".join(ext_list)
                        # if ext_key == "part_of(CL),part_of(EMAPA),part_of(EMAPA)":
                        #     print(ext_list)
                        if ext_key not in ext_dict[aspect]:
                            ext_dict[aspect][ext_key] = [g]
                        else:
                            ext_dict[aspect][ext_key].append(g)
                # Standardize key - ex:
                #   ["part_of(GO:aspect),part_of(UBERON:)"]

            for aspect in ['F','P','C']:
                max_count = 0
                top_k = None
                example_v = None
                for k, v in ext_dict[aspect].items():
                    if len(v) > max_count:
                        max_count = len(v)
                        # print(max_count)
                        top_k = k
                        example_v = v[0]
                # print("Aspect: ", aspect, " Count: ", max_count, " - ", format_extensions(top_k))
                # print(example_v)

            # with open(fname + ".tsv", 'w') as f:
            #     writer = csv.writer(f, delimiter='\t')
            #     for aspect in ['F', 'P', 'C']:
            #         max_count = 0
            #         top_k = None
            #         example_v = None
            #         for k, v in ext_dict[aspect].items():
            #             writer.writerow([aspect, len(v), k])
                    # print("Aspect: ", aspect, " Count: ", max_count, " - ", format_extensions(top_k))
                    # print(example_v)
            # all_dict = {**all_dict, **ext_dict}
            print("ext_dict length:", len(ext_dict))

    def assigner_count(annots, assign):
        count = 0
        for a in annots:
            if a[9] == assign:
                count += 1
        return count

    cols = ['Aspect', 'Total count', 'Extension']
    all_assigners = []
    for aspect in ['F', 'P', 'C']:
        for k, v in ext_dict[aspect].items():
            for a in v:
                if a[9] not in all_assigners:
                    all_assigners.append(a[9])
                    cols.append(a[9])

    def parse_pattern_sourcefile(sourcefile):
        parsed_patterns = []
        with open(sourcefile) as sf:
            for pl in sf.readlines():
                pl = pl.rstrip()
                parsed_patterns.append(pl)
        return parsed_patterns

    class GpadWriter:
        def __init__(self, gpad_file):
            self.gpad_file = gpad_file
            self.writer = csv.writer(pattern_gpad, delimiter="\t")

        def writerow(self, row):
            self.writer.writerow(row)

    class WriterCollection:
        def __init__(self):
            self.writers = {}

        def set_writer(self, writer, writer_name):
            self.writers[writer_name] = writer


    patterns = []
    writers = WriterCollection()
    if args.pattern or args.pattern_sourcefile:
        if args.pattern_sourcefile:
            patterns = parse_pattern_sourcefile(args.pattern_sourcefile)
            # writers = {}
            for patt in patterns:
                pattern_outfile = "{}.gpad".format(patt)
                pattern_gpad = open(pattern_outfile, "w+")
                pattern_gpad_writer = GpadWriter(pattern_gpad)
                writers.set_writer(pattern_gpad_writer, patt)
        else:
            patterns.append(args.pattern)
            pattern_outfile = "{}.gpad".format(args.pattern)
            if args.pattern_outfile:
                pattern_outfile = args.pattern_outfile
            pattern_gpad = open(pattern_outfile, "w+")
            pattern_gpad_writer = GpadWriter(pattern_gpad)
            WriterCollection.set_writer(pattern_gpad_writer, args.pattern)
    out_file = "all.tsv"
    if args.out_file:
        out_file = args.out_file
    with open(out_file, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(cols)
        for aspect in ['F', 'P', 'C']:
            max_count = 0
            top_k = None
            example_v = None
            for k, v in ext_dict[aspect].items():
                row_to_write = [aspect, len(v), k]
                if len(patterns) > 0:
                    for patt in patterns:
                        if k == patt:
                            for a in v:
                                writers.writers[patt].writerow(a[0:len(a)-1])
                for assigner in all_assigners:
                    a_count = assigner_count(v, assigner)
                    row_to_write.append(a_count)
                writer.writerow(row_to_write)
    for writer_name in writers.writers:
        writers.writers[writer_name].gpad_file.close()

    wanted_gafs = []
    leftovers_out_file = "leftovers.gpad"
    if args.leftovers_out_file:
        leftovers_out_file = args.leftovers_out_file
    with open(leftovers_out_file, 'w') as wf:
        for k in ext_dict['P'].keys():
            if 'has_regulation_target' in k:
                wanted_gafs = wanted_gafs + ext_dict['P'][k]
                for g in ext_dict['P'][k]:
                    wf.write("\t".join(g))
