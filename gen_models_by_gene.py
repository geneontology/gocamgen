from gocamgen.gocamgen import AssocGoCamModel
from ontobio.io.gpadparser import GpadParser
from ontobio.ontol_factory import OntologyFactory
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('-g', '--gpad_file', help="Filepath of GPAD source with annotations to model", required=True)
parser.add_argument('-s', '--specific_gene', help="If specified, will only translate model for annotations "
                                                  "to this specific gene")
parser.add_argument('-n', '--max_model_limit', help="Only translate specified number of models. Mainly for testing.")

def translate_to_model(gene, assocs, ont):
    model = AssocGoCamModel(gene, assocs, ont)

    model.translate()

    out_filename = "{}.ttl".format(gene.replace(":", "_"))
    model.write(out_filename)

    print("Model for {} written to {}".format(gene, out_filename))
    return model

args = parser.parse_args()

gpad_parser = GpadParser()
assocs = gpad_parser.parse(args.gpad_file, skipheader=True)

assocs_by_gene = {}
for a in assocs:
    subject_id = a["subject"]["id"]
    if subject_id in assocs_by_gene:
        assocs_by_gene[subject_id].append(a)
    else:
        assocs_by_gene[subject_id] = [a]
print("{} distinct genes".format(len(assocs_by_gene)))

ont = OntologyFactory().create("go")

if args.specific_gene:
    # python3 gen_models_by_gene.py --gpad_file wb.gpad --specific_gene WB:WBGene00003609
    model = translate_to_model(args.specific_gene, assocs_by_gene[args.specific_gene], ont)
else:
    count = 0
    for gene in assocs_by_gene:
        model = translate_to_model(gene, assocs_by_gene[gene], ont)
        count += 1
        if args.max_model_limit and count == args.max_model_limit:
            break

