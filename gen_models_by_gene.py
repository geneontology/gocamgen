from gocamgen.gocamgen import AssocGoCamModel
from gpad_extensions_mapper import ExtensionsMapper
from filter_rule import AssocFilter, FilterRule, get_filter_rule
from ontobio.io.gpadparser import GpadParser
from ontobio.ontol_factory import OntologyFactory
# from ontobio.ecomap import EcoMap
import argparse
import logging
import requests
import gzip
import time
from os import path
# from abc import ABC, abstractmethod

logging.basicConfig()
logger = logging.getLogger(__name__)
logger.setLevel("DEBUG")

parser = argparse.ArgumentParser()
parser.add_argument('-g', '--gpad_file', help="Filepath of GPAD source with annotations to model", required=True)
parser.add_argument('-s', '--specific_gene', help="If specified, will only translate model for annotations "
                                                  "to this specific gene")
parser.add_argument('-n', '--max_model_limit', help="Only translate specified number of models. Mainly for testing.")
parser.add_argument('-m', '--mod', help="MOD rules to follow for filtering and translating.")
parser.add_argument('-d', '--output_directory', help="Directory to output model ttl files to")
parser.add_argument('-r', '--report', help="Generate report", action="store_const", const=True)


class GoCamBuilder:
    def __init__(self):
        self.ext_mapper = ExtensionsMapper()
        # self.ro_ontology = OntologyFactory().create("http://purl.obolibrary.org/obo/ro.owl")
        self.go_ontology = OntologyFactory().create("go")

    def translate_to_model(self, gene, assocs):
        model = AssocGoCamModel(gene, assocs)
        model.extensions_mapper = self.ext_mapper
        model.ontology = self.go_ontology
        model.translate()

        return model


class AssocExtractor:
    def __init__(self, gpad_file, filter_rule : FilterRule):
        gpad_parser = GpadParser()
        assocs = gpad_parser.parse(gpad_file, skipheader=True)
        self.assocs = extract_properties(assocs)
        self.assoc_filter = AssocFilter(filter_rule)

    def group_assocs(self):
        assocs_by_gene = {}
        for a in self.assocs:
            # validation function
            if not self.assoc_filter.validate_line(a):
                continue
            subject_id = a["subject"]["id"]
            if subject_id in assocs_by_gene:
                assocs_by_gene[subject_id].append(a)
            else:
                assocs_by_gene[subject_id] = [a]
        return assocs_by_gene


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


def unzip(filepath):
    input_file = gzip.GzipFile(filepath, "rb")
    s = input_file.read()
    input_file.close()

    target = path.splitext(filepath)[0]
    logger.info("Gunzipping file: {}".format(filepath))
    with open(target, "wb") as output:
        output.write(s)
    return target


def handle_gpad_file(args_gpad_file):
    if args_gpad_file.startswith("http://"):
        gpad_file_target = args_gpad_file.split("/")[-1]
        logger.info("Downloading GPAD from {} and saving to {}".format(args_gpad_file, gpad_file_target))
        response = requests.get(args_gpad_file)
        with open(gpad_file_target, "wb") as gft:
            gft.write(response.content)
        if gpad_file_target.endswith(".gz"):
            gpad_file = unzip(gpad_file_target)
        else:
            gpad_file = gpad_file_target
    else:
        gpad_file = args.gpad_file
    return gpad_file


def parse_header(gpad_file):
    header_data = {
        "date": ""
    }
    with open(gpad_file) as gf:
        for l in gf.readlines():
            if l.startswith("!"):
                date_key = "!date: "
                if l.startswith(date_key):
                    header_data["date"] = l.split(date_key)[1].split("$")[0].strip()
            else:
                break
    return header_data


if __name__ == "__main__":
    args = parser.parse_args()

    filter_rule = get_filter_rule(args.mod)

    gpad_file = handle_gpad_file(args.gpad_file)
    relevant_header_data = parse_header(gpad_file)
    gpad_file_metadata = {
        "source_path": args.gpad_file,
        # TODO: Figure out how to get real creation date from file
        "download_date": time.ctime(path.getmtime(gpad_file)),
        "header_date": relevant_header_data["date"]
    }

    extractor = AssocExtractor(gpad_file, filter_rule)
    assocs_by_gene = extractor.group_assocs()
    logger.debug("{} distinct genes".format(len(assocs_by_gene)))

    builder = GoCamBuilder()

    model_count = 0
    if args.specific_gene:
        for specific_gene in args.specific_gene.split(","):
            if specific_gene not in assocs_by_gene:
                logger.error("ERROR: specific gene {} not found in filtered annotation list".format(specific_gene))
            else:
                logger.debug("{} filtered annotations to translate for {}".format(len(assocs_by_gene[specific_gene]), specific_gene))
                model = builder.translate_to_model(specific_gene, assocs_by_gene[specific_gene])
                out_filename = "{}.ttl".format(specific_gene.replace(":", "_"))
                if args.output_directory:
                    out_filename = path.join(args.output_directory, out_filename)
                model.write(out_filename)
                logger.info("Model for {} written to {}".format(specific_gene, out_filename))
                model_count += 1
    else:
        for gene in assocs_by_gene:
            model = builder.translate_to_model(gene, assocs_by_gene[gene])
            out_filename = "{}.ttl".format(gene.replace(":", "_"))
            if args.output_directory:
                out_filename = path.join(args.output_directory, out_filename)
            model.write(out_filename)
            logger.info("Model for {} written to {}".format(gene, out_filename))
            model_count += 1
            if args.max_model_limit and model_count == args.max_model_limit:
                break

    if args.report:
        report_file_path = "{}.report".format(gpad_file)
        with open(report_file_path, "w+") as reportf:
            for k in gpad_file_metadata:
                reportf.write("{}: {}\n".format(k, gpad_file_metadata[k]))
            # TODO FilterRule().__str__() to display filters
            reportf.write("# of models generated: {}\n".format(model_count))
        logger.info("Report file generated at {}".format(report_file_path))
