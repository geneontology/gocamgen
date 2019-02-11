from gocamgen.gocamgen import AssocGoCamModel
from gpad_extensions_mapper import ExtensionsMapper
from ontobio.io.gpadparser import GpadParser
from ontobio.ontol_factory import OntologyFactory
from ontobio.ecomap import EcoMap
import argparse
from os import path

parser = argparse.ArgumentParser()
parser.add_argument('-g', '--gpad_file', help="Filepath of GPAD source with annotations to model", required=True)
parser.add_argument('-s', '--specific_gene', help="If specified, will only translate model for annotations "
                                                  "to this specific gene")
parser.add_argument('-n', '--max_model_limit', help="Only translate specified number of models. Mainly for testing.")
parser.add_argument('-m', '--mod', help="MOD rules to follow for filtering and translating.")
parser.add_argument('-d', '--output_directory', help="Directory to output model ttl files to")


class FilterRule():
    # Default filter - Remove IEA, IBA, as well as IKRs originating from PAINT.
    def __init__(self, unwanted_evidence_codes=['IEA', 'IBA'],
                 unwanted_evi_code_ref_combos=[('IKR', 'PMID:21873635')],
                 required_attributes=[]):
        self.unwanted_evidence_codes = unwanted_evidence_codes
        self.unwanted_evi_code_ref_combos = unwanted_evi_code_ref_combos
        self.required_attributes = required_attributes
        self.unwanted_properties = []


class WBFilterRule(FilterRule):
    # So far same as default FilterRule.
    def __init__(self):
        FilterRule.__init__(self)


class MGIFilterRule(FilterRule):
    # Only provided_by=MGI lines are valid, also filtering out ISOs along with default FilterRule filters.
    def __init__(self):
        FilterRule.__init__(self, required_attributes=[{"provided_by": ["MGI"]}])
        self.unwanted_evidence_codes.append('ISO')
        self.unwanted_properties = ["noctua-model-id"]

mod_filter_map = {
    'WB': WBFilterRule,
    'MGI': MGIFilterRule
}

def extract_properties(assocs):
    new_assoc_list = []
    for a in assocs:
        cols = a["source_line"].rstrip().split("\t")
        if len(cols) >= 12:
            prop_col = cols[11]
            props = prop_col.split("|")
            props_dict = {}
            for p in props:
                k, v = p.split("=")
                if k in props_dict:
                    props_dict[k].append(v)
                else:
                    props_dict[k] = [v]
            a["annotation_properties"] = props_dict
        new_assoc_list.append(a)
    return new_assoc_list

if __name__ == "__main__":
    args = parser.parse_args()

    if args.mod in mod_filter_map:
        filter_rule = mod_filter_map[args.mod]()
    else:
        filter_rule = FilterRule()
    ecomap = EcoMap()
    ext_mapper = ExtensionsMapper()

    def translate_to_model(gene, assocs):
        model = AssocGoCamModel(gene, assocs)
        model.extensions_mapper = ext_mapper

        model.translate()

        out_filename = "{}.ttl".format(gene.replace(":", "_"))
        if args.output_directory:
            out_filename = path.join(args.output_directory, out_filename)
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
        if len(filter_rule.unwanted_properties) > 0 and "annotation_properties" in assoc:
            for up in filter_rule.unwanted_properties:
                if up in assoc["annotation_properties"]:
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
    assocs = extract_properties(assocs)
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

    # ont = OntologyFactory().create("go")

    if args.specific_gene:
        # python3 gen_models_by_gene.py --gpad_file wb.gpad --specific_gene WB:WBGene00003609
        for specific_gene in args.specific_gene.split(","):
            if specific_gene not in assocs_by_gene:
                print("ERROR: specific gene {} not found in filtered annotation list".format(specific_gene))
            else:
                model = translate_to_model(specific_gene, assocs_by_gene[specific_gene])
    else:
        count = 0
        for gene in assocs_by_gene:
            model = translate_to_model(gene, assocs_by_gene[gene])
            count += 1
            if args.max_model_limit and count == args.max_model_limit:
                break

