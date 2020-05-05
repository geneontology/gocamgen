from gocamgen.gocamgen import AssocGoCamModel
from gocamgen.gpad_extensions_mapper import ExtensionsMapper
from gocamgen.filter_rule import AssocFilter, FilterRule, get_filter_rule
from gocamgen.collapsed_assoc import extract_properties
from gocamgen.errors import GocamgenException, GeneErrorSet
from gocamgen.utils import ShexException
from ontobio.io.gpadparser import GpadParser
from ontobio.ontol_factory import OntologyFactory
# from ontobio.ecomap import EcoMap
import argparse
import logging
import requests
from requests.exceptions import ConnectionError
import gzip
import time
from os import path
# from abc import ABC, abstractmethod
from rdflib.graph import ConjunctiveGraph
from rdflib.store import Store
from rdflib import plugin

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
parser.add_argument('-N', '--nquads', help="Filepath to write model file in N-Quads format")

# GoCamInputHandler


class GoCamBuilder:
    def __init__(self):
        self.ro_ontology = OntologyFactory().create("http://purl.obolibrary.org/obo/ro.owl")
        self.gorel_ontology = OntologyFactory().create("http://release.geneontology.org/2019-03-18/ontology/extensions/gorel.obo")
        # Can't get logical_definitions w/ ont.create("go"), need to load ontology via PURL
        self.go_ontology = OntologyFactory().create("http://purl.obolibrary.org/obo/go.owl")
        self.ext_mapper = ExtensionsMapper(go_ontology=self.go_ontology, ro_ontology=self.ro_ontology)
        self.store = plugin.get('IOMemory', Store)()

    def translate_to_model(self, gene, assocs):
        model = AssocGoCamModel(gene, assocs, store=self.store)
        model.extensions_mapper = self.ext_mapper
        model.ontology = self.go_ontology
        model.ro_ontology = self.ro_ontology
        model.gorel_ontology = self.gorel_ontology
        model.translate()

        return model


class AssocExtractor:
    def __init__(self, gpad_file, filter_rule : FilterRule):
        gpad_parser = GpadParser()
        assocs = gpad_parser.parse(gpad_file, skipheader=True)
        self.assocs = extract_properties_from_assocs(assocs)
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


def extract_properties_from_assocs(assocs):
    new_assoc_list = []
    for a in assocs:
        new_assoc_list.append(extract_properties(a))
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
    errors = GeneErrorSet()  # Errors by gene ID

    def make_model_and_write_out(gene, nquads=False):
        # All these shenanigans are to prevent mid-run crashes due to an external resource simply blipping
        # out for a second.
        retry_count = 0
        retry_limit = 5
        while True:
            try:
                start_time = time.time()
                model = builder.translate_to_model(gene, assocs_by_gene[gene])
                # add_to_conjunctive_graph(model, conjunctive_graph)
                if nquads:
                    logger.info(
                        "Model for {} added to graphstore in {} sec".format(gene, (time.time() - start_time)))
                else:
                    out_filename = "{}.ttl".format(gene.replace(":", "_"))
                    if args.output_directory:
                        out_filename = path.join(args.output_directory, out_filename)
                    model.write(out_filename)
                    logger.info("Model for {} written to {} in {} sec".format(gene, out_filename, (time.time() - start_time)))
            except GocamgenException as ex:
                errors.add_error(gene, ex)
            except (TimeoutError, ConnectionError) as ex:
                # This has been happening randomly and breaking full runs
                errors.add_error(gene, ex)
                if retry_count < retry_limit:
                    retry_count += 1
                    continue  # retry
                errors.add_error(gene, GocamgenException(f"Bailing on model for {gene} after {retry_count} retries"))
            break  # Done with this model. Move on to the next one.


    model_count = 0
    if args.specific_gene:
        for specific_gene in args.specific_gene.split(","):
            if specific_gene not in assocs_by_gene:
                logger.error("ERROR: specific gene {} not found in filtered annotation list".format(specific_gene))
            else:
                logger.debug("{} filtered annotations to translate for {}".format(len(assocs_by_gene[specific_gene]), specific_gene))
                make_model_and_write_out(specific_gene, args.nquads)
                model_count += 1
    else:
        for gene in assocs_by_gene:
            make_model_and_write_out(gene, args.nquads)
            model_count += 1
            if args.max_model_limit and model_count == int(args.max_model_limit):
                break

    if args.nquads:
        cg = ConjunctiveGraph(builder.store)
        cg.serialize(destination=args.nquads, format="nquads")
        logger.info(f"Full model graphstore written out in N-Quads format to {args.nquads}")

    if args.report:
        report_file_path = "{}.report".format(gpad_file)
        with open(report_file_path, "w+") as reportf:
            for k in gpad_file_metadata:
                reportf.write("{}: {}\n".format(k, gpad_file_metadata[k]))
            # TODO FilterRule().__str__() to display filters
            reportf.write("# of models generated: {}\n".format(model_count))
            for gene, errs in errors.errors.items():
                for ex in errs:
                    reportf.write(f"{type(ex).__name__} - {gene}: {ex}\n")
        logger.info("Report file generated at {}".format(report_file_path))
