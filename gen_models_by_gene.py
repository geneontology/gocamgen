from gocamgen.gocamgen import AssocGoCamModel
from ontobio.io.gpadparser import GpadParser
from ontobio.ontol_factory import OntologyFactory
from ontobio.ecomap import EcoMap
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('-g', '--gpad_file', help="Filepath of GPAD source with annotations to model", required=True)
parser.add_argument('-s', '--specific_gene', help="If specified, will only translate model for annotations "
                                                  "to this specific gene")
parser.add_argument('-n', '--max_model_limit', help="Only translate specified number of models. Mainly for testing.")
parser.add_argument('-m', '--mod', help="MOD rules to follow for filtering and translating.")


class FilterRule():
    def __init__(self, unwanted_evidence_codes=['IEA', 'IBA'],
                 unwanted_evi_code_ref_combos=[('IKR', 'PMID:21873635')],
                 required_attributes=[]):
        self.unwanted_evidence_codes = unwanted_evidence_codes
        self.unwanted_evi_code_ref_combos = unwanted_evi_code_ref_combos
        self.required_attributes = required_attributes


class WBFilterRule(FilterRule):
    def __init__(self):
        FilterRule.__init__(self)


class MGIFilterRule(FilterRule):
    def __init__(self):
        FilterRule.__init__(self, required_attributes=[{"provided_by": ["MGI"]}])
        self.unwanted_evidence_codes.append('ISO')

mod_filter_map = {
    'WB': WBFilterRule,
    'MGI': MGIFilterRule
}

args = parser.parse_args()

if args.mod in mod_filter_map:
    filter_rule = mod_filter_map[args.mod]()
else:
    filter_rule = FilterRule()
ecomap = EcoMap()

def translate_to_model(gene, assocs, ont):
    model = AssocGoCamModel(gene, assocs, ont)

    model.translate()

    out_filename = "{}.ttl".format(gene.replace(":", "_"))
    model.write(out_filename)

    print("Model for {} written to {}".format(gene, out_filename))
    return model

def validate_line(assoc):
    evi_code = ecomap.ecoclass_to_coderef(assoc["evidence"]["type"])[0]
    if evi_code in filter_rule.unwanted_evidence_codes:
        return False
    if len(filter_rule.unwanted_evi_code_ref_combos) > 0:
        references = assoc["evidence"]["has_supporting_reference"]
        for evi_ref_combo in filter_rule.unwanted_evi_code_ref_combos:
            if evi_ref_combo[0] == evi_code and evi_ref_combo[1] in references:
                return False
    if len(filter_rule.required_attributes) > 0:
        meets_requirement = False
        for attr in filter_rule.required_attributes:
            # a[attr] is dict
            for k in attr.keys():
                if a[k] in attr[k]:
                    meets_requirement = True
        if not meets_requirement:
            return False
    return True

gpad_parser = GpadParser()
assocs = gpad_parser.parse(args.gpad_file, skipheader=True)

assocs_by_gene = {}
for a in assocs:
    # validation function
    if not validate_line(a):
        continue
    subject_id = a["subject"]["id"]
    if subject_id in assocs_by_gene:
        assocs_by_gene[subject_id].append(a)
    else:
        assocs_by_gene[subject_id] = [a]
print("{} distinct genes".format(len(assocs_by_gene)))

ont = OntologyFactory().create("go")

if args.specific_gene:
    # python3 gen_models_by_gene.py --gpad_file wb.gpad --specific_gene WB:WBGene00003609
    if args.specific_gene not in assocs_by_gene:
        print("ERROR: specific gene {} not found in filtered annotation list".format(args.specific_gene))
    else:
        model = translate_to_model(args.specific_gene, assocs_by_gene[args.specific_gene], ont)
else:
    count = 0
    for gene in assocs_by_gene:
        model = translate_to_model(gene, assocs_by_gene[gene], ont)
        count += 1
        if args.max_model_limit and count == args.max_model_limit:
            break

