from ontobio.rdfgen.assoc_rdfgen import CamRdfTransform, TurtleRdfWriter, genid, prefix_context
from ontobio.vocabulary.relations import OboRO, Evidence
from ontobio.vocabulary.upper import UpperLevel
from ontobio.util.go_utils import GoAspector
from ontobio.ontol_factory import OntologyFactory
from prefixcommons.curie_util import expand_uri, contract_uri
from rdflib.namespace import OWL, RDF
from rdflib import Literal
from rdflib.term import URIRef
from rdflib.namespace import Namespace
import rdflib
import networkx
import logging
import argparse
import datetime
import os.path as path
import logging
from triple_pattern_finder import TriplePattern, TriplePatternFinder, TriplePair, TriplePairCollection

# logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
logger.setLevel("INFO")

ro = OboRO()
evt = Evidence()
upt = UpperLevel()
LEGO = Namespace("http://geneontology.org/lego/")
LAYOUT = Namespace("http://geneontology.org/lego/hint/layout/")
PAV = Namespace('http://purl.org/pav/')
DC = Namespace("http://purl.org/dc/elements/1.1/")
RDFS = Namespace("http://www.w3.org/2000/01/rdf-schema#")

# Stealing a lot of code for this from ontobio.rdfgen:
# https://github.com/biolink/ontobio


def expand_uri_wrapper(id):
    uri = expand_uri(id, cmaps=[prefix_context])
    return uri

def contract_uri_wrapper(id):
    uri = contract_uri(id, cmaps=[prefix_context])
    return uri

HAS_SUPPORTING_REFERENCE = URIRef(expand_uri_wrapper("dc:source"))
ENABLED_BY = URIRef(expand_uri_wrapper(ro.enabled_by))
ENABLES = URIRef(expand_uri_wrapper(ro.enables))
INVOLVED_IN = URIRef(expand_uri_wrapper(ro.involved_in))
PART_OF = URIRef(expand_uri_wrapper(ro.part_of))
OCCURS_IN = URIRef(expand_uri_wrapper(ro.occurs_in))
COLOCALIZES_WITH = URIRef(expand_uri_wrapper(ro.colocalizes_with))
CONTRIBUTES_TO = URIRef(expand_uri_wrapper("RO:0002326"))
MOLECULAR_FUNCTION = URIRef(expand_uri_wrapper(upt.molecular_function))
REGULATES = URIRef(expand_uri_wrapper("RO:0002211"))

# RO_ONTOLOGY = OntologyFactory().create("http://purl.obolibrary.org/obo/ro.owl")  # Need propertyChainAxioms to parse (https://github.com/biolink/ontobio/issues/312)
INPUT_RELATIONS = {
    #TODO Create rule for deciding (MOD-specific?) whether to convert has_direct_input to has input
    # "has_direct_input": "RO:0002400",
    "has_direct_input": "RO:0002233",
    "has input": "RO:0002233"
}
ACTS_UPSTREAM_OF_RELATIONS = {
    "acts_upstream_of": "RO:0002263",
    "acts_upstream_of_or_within": "RO:0002264",
    "acts upstream of or within, positive effect": "RO:0004032",
    "acts upstream of or within, negative effect": "RO:0004033",
    "acts_upstream_of_positive_effect": "RO:0004034",
    "acts_upstream_of_negative_effect": "RO:0004035",
}
HAS_REGULATION_TARGET_RELATIONS = {
    # WB:WBGene00013591 involved_in GO:0042594
    "has_regulation_target": ""
}
# TODO: Grab these from RO property chain axioms
ENABLES_O_RELATION_LOOKUP = {
    # "acts upstream of": "causally upstream of",
    "RO:0002263": "RO:0002411",
    # "acts upstream of or within": "causally upstream of or within",
    "RO:0002264": "RO:0002418",
    # "acts upstream of or within, positive effect": "causally upstream of or within, positive effect",
    "RO:0004032": "RO:0004047",
    # "acts upstream of or within, negative effect": "causally upstream of or within, negative effect",
    "RO:0004033": "RO:0004046",
    "RO:0004034": "RO:0002304",
    "RO:0004035": "RO:0002305"
}


def has_regulation_target_bucket(ontology, term):
    ancestors = ontology.ancestors(term, reflexive=True)
    buckets = []
    if "GO:0065009" in ancestors:
        buckets.append("a")
    if "GO:0010468" in ancestors:
        buckets.append("b")
    if "GO:0002092" in ancestors or "GO:0042176" in ancestors:
        buckets.append("c")
    if "GO:0019538" in ancestors or "GO:0032880" in ancestors:
        buckets.append("d")
    return buckets


now = datetime.datetime.now()


class Annoton():
    def __init__(self, subject_id, assocs, connections=None):
        self.enabled_by = subject_id
        self.annotations = assocs
        self.connections = connections
        self.individuals = {}

class GoCamModel():
    relations_dict = {
        "has_direct_input": "RO:0002400",
        "has input": "RO:0002233",
        "has_regulation_target": "RO:0002211",  # regulates
        "regulates_activity_of": "RO:0002578",  # directly regulates
        "with_support_from": "RO:0002233",  # has input
        "directly_regulates": "RO:0002578",
        "directly_positively_regulates": "RO:0002629",
        "directly_negatively_regulates": "RO:0002630",
        "colocalizes_with": "RO:0002325",
        "contributes_to": "RO:0002326",
        "part_of": "BFO:0000050",
        "acts_upstream_of": "RO:0002263",
        "acts_upstream_of_negative_effect": "RO:0004035",
        "acts_upstream_of_or_within": "RO:0002264",
        "acts_upstream_of_positive_effect": "RO:0004034",
    }

    def __init__(self, modeltitle, connection_relations=None):
        cam_writer = CamTurtleRdfWriter(modeltitle)
        self.writer = AnnotonCamRdfTransform(cam_writer)
        self.modeltitle = modeltitle
        self.classes = []
        self.individuals = {}   # Maintain entity-to-IRI dictionary. Prevents dup individuals but we may want dups?
        self.graph = networkx.MultiDiGraph()  # networkx graph of individuals and relations? Could this replace self.individuals? Will this conflict with self.writer.writer.graph?
        # Each node:
        ## node_id
        ## class
        ## attributes
        # Each edge:
        ## source
        ## target
        ## relation
        ## other attributes?
        if connection_relations is None:
            self.connection_relations = GoCamModel.relations_dict
        else:
            self.connection_relations = connection_relations
        self.declare_properties()

    def write(self, filename):
        if path.splitext(filename)[1] != ".ttl":
            filename += ".ttl"
        with open(filename, 'wb') as f:
            self.writer.writer.serialize(destination=f)

    def declare_properties(self):
        # AnnotionProperty
        self.writer.emit_type(URIRef("http://geneontology.org/lego/evidence"), OWL.AnnotationProperty)
        self.writer.emit_type(URIRef("http://geneontology.org/lego/hint/layout/x"), OWL.AnnotationProperty)
        self.writer.emit_type(URIRef("http://geneontology.org/lego/hint/layout/y"), OWL.AnnotationProperty)
        self.writer.emit_type(URIRef("http://purl.org/pav/providedBy"), OWL.AnnotationProperty)
        self.writer.emit_type(URIRef("http://purl.org/dc/elements/1.1/contributor"), OWL.AnnotationProperty)
        self.writer.emit_type(URIRef("http://purl.org/dc/elements/1.1/date"), OWL.AnnotationProperty)
        self.writer.emit_type(URIRef("http://purl.org/dc/elements/1.1/source"), OWL.AnnotationProperty)

    def declare_class(self, class_id):
        if class_id not in self.classes:
            self.writer.emit_type(URIRef(expand_uri_wrapper(class_id)), OWL.Class)
            self.classes.append(class_id)

    def declare_individual(self, entity_id):
        entity = genid(base=self.writer.writer.base + '/')
        self.writer.emit_type(entity, self.writer.uri(entity_id))
        self.writer.emit_type(entity, OWL.NamedIndividual)
        self.individuals[entity_id] = entity
        self.graph.add_node(entity, **{"label": entity_id})
        return entity

    def add_axiom(self, statement, evidence=None):
        (source_id, property_id, target_id) = statement
        stmt_id = self.find_bnode(statement)
        if stmt_id is None:
            stmt_id = self.writer.blanknode()
            self.writer.emit_type(stmt_id, OWL.Axiom)
        self.writer.emit(stmt_id, OWL.annotatedSource, source_id)
        self.writer.emit(stmt_id, OWL.annotatedProperty, property_id)
        self.writer.emit(stmt_id, OWL.annotatedTarget, target_id)

        if evidence:
            self.add_evidence(stmt_id, evidence.evidence_code, evidence.references)

        return stmt_id

    def create_axiom(self, subject_id, relation_uri, object_id):
        subject_uri = subject_id if subject_id.__class__.__name__ == "URIRef" else self.declare_individual(subject_id)
        object_uri = object_id if object_id.__class__.__name__ == "URIRef" else self.declare_individual(object_id)
        axiom_id = self.add_axiom(self.writer.emit(subject_uri, relation_uri, object_uri))
        return axiom_id

    def find_or_create_axiom(self, subject_id : str, relation_uri : URIRef, object_id : str, annoton=None,
                             exact_length=False):
        # Maybe overkill but gonna try using find_pattern_recursive to find only one triple
        pattern = TriplePattern([(subject_id, relation_uri, object_id)])
        found_triples = TriplePatternFinder().find_pattern_recursive(self, pattern, exact_length=True)
        # found_triples = self.triples_by_ids(subject_id, relation_uri, object_id)
        if len(found_triples) > 0:
            # Gonna be a list of "triple-chains", itself a list of triples, and each triple is sort like a 3-index list.
            # So we just want the first triple from the first chain:
            found_triple = found_triples[0][0]
            subject_uri = found_triple[0]
            object_uri = found_triple[2]
            axiom_id = self.find_bnode(found_triple)
        else:
            # subject_uri = self.declare_individual(subject_id)
            subject_uri = subject_id if subject_id.__class__.__name__ == "URIRef" else self.declare_individual(subject_id)
            object_uri = object_id if object_id.__class__.__name__ == "URIRef" else self.declare_individual(object_id)
            # TODO Can emit() be changed to emit_axiom()?
            axiom_id = self.add_axiom(self.writer.emit(subject_uri, relation_uri, object_uri))
        if annoton and relation_uri == ENABLED_BY:
            annoton.individuals[subject_id] = subject_uri
            annoton.individuals[object_id] = object_uri
        return axiom_id

    def add_evidence(self, axiom, evidence_code, references, contributors=[], date="", comment=""):
        ev = GoCamEvidence(evidence_code, references, contributors=contributors, date=date, comment=comment)
        # Try finding existing evidence object containing same type and references
        # ev_id = self.writer.find_or_create_evidence_id(ev)
        ev_id = self.writer.create_evidence(ev)
        self.writer.emit(axiom, URIRef("http://geneontology.org/lego/evidence"), ev_id)
        ### Emit ev fields to axiom here TODO: Couple evidence and axiom emitting together
        self.writer.emit(axiom, DC.date, Literal(ev.date))
        self.writer.emit(axiom, RDFS.comment, Literal(ev.comment))
        # self.writer.emit(axiom, RDFS.comment, Literal(""))
        for c in ev.contributors:
            self.writer.emit(axiom, DC.contributor, Literal(c))


    def add_connection(self, gene_connection, source_annoton):
        # Switching from reusing existing activity node from annoton to creating new one for each connection - Maybe SPARQL first to check if annoton activity already used for connection?
        # Check annoton for existing activity.
        # if gene_connection.object_id in source_annoton.individuals:
        #     # If exists and activity has connection relation,
        #     # Look for two triples: (gene_connection.object_id, ENABLED_BY, source_annoton.enabled_by) and (gene_connection.object_id, connection_relations, anything)
        # Annot MF should be declared by now - don't declare object_id if object_id == annot MF?
        if gene_connection.gp_b not in self.individuals:
            return
        source_id = None
        uri_list = self.uri_list_for_individual(gene_connection.object_id)
        for u in uri_list:
            if gene_connection.relation in self.connection_relations:
                rel = URIRef(expand_uri_wrapper(self.connection_relations[gene_connection.relation]))
                # Annot MF should be declared by now - don't declare object_id if object_id == annot MF?
                try:
                    annot_mf = source_annoton.molecular_function["object"]["id"]
                except:
                    annot_mf = ""
                if self.writer.writer.graph.__contains__((u,rel,None)) and gene_connection.object_id != annot_mf:
                    source_id = self.declare_individual(gene_connection.object_id)
                    source_annoton.individuals[gene_connection.object_id] = source_id
                    break

        if source_id is None:
            try:
                source_id = source_annoton.individuals[gene_connection.object_id]
            except KeyError:
                source_id = self.declare_individual(gene_connection.object_id)
                source_annoton.individuals[gene_connection.object_id] = source_id
        # Add enabled by stmt for object_id - this is essentially adding another annoton connecting gene-to-extension/with-MF to the model
        self.writer.emit(source_id, ENABLED_BY, source_annoton.individuals[source_annoton.enabled_by])
        self.writer.emit_axiom(source_id, ENABLED_BY, source_annoton.individuals[source_annoton.enabled_by])
        property_id = URIRef(expand_uri_wrapper(self.connection_relations[gene_connection.relation]))
        target_id = self.individuals[gene_connection.gp_b]
        # Annotate source MF GO term NamedIndividual with relation code-target MF term URI
        self.writer.emit(source_id, property_id, target_id)
        # Add axiom (Source=MF term URI, Property=relation code, Target=MF term URI)
        self.writer.emit_axiom(source_id, property_id, target_id)

    def uri_list_for_individual(self, individual):
        uri_list = []
        graph = self.writer.writer.graph
        for t in graph.triples((None,None,self.writer.uri(individual))):
            uri_list.append(t[0])
        return uri_list

    def triples_by_ids(self, subject, relation_uri, object_id):
        graph = self.writer.writer.graph

        triples = []
        if subject.__class__.__name__ == "URIRef" or subject is None:
            subjects = [subject]
        else:
            subjects = self.uri_list_for_individual(subject)
        if object_id.__class__.__name__ == "URIRef" or object_id is None:
            objects = [object_id]
        else:
            objects = self.uri_list_for_individual(object_id)
        for object_uri in objects:
            for subject_uri in subjects:
                # if (subject_uri, relation_uri, object_uri) in graph:
                #     triples.append((subject_uri, relation_uri, object_uri))
                for t in graph.triples((subject_uri, relation_uri, object_uri)):
                    triples.append(t)
        return triples

    def individual_label_for_uri(self, uri):
        ind_list = []
        graph = self.writer.writer.graph
        for t in graph.triples((uri, RDF.type, None)):
            if t[2] != OWL.NamedIndividual: # We know OWL.NamedIndividual triple does't contain the label so don't return it
                ind_list.append(t[2])
        return ind_list

    def class_for_uri(self, uri):
        try:
            class_curie = contract_uri_wrapper(self.individual_label_for_uri(uri)[0])[0]
            return class_curie
        except:
            return None

    def axioms_for_source(self, source, property_uri=None):
        if property_uri is None:
            property_uri = OWL.annotatedSource
        axiom_list = []
        graph = self.writer.writer.graph
        for uri in self.uri_list_for_individual(source):
            for t in graph.triples((None, property_uri, uri)):
                axiom_list.append(t[0])
        return axiom_list

    def find_bnode(self, triple):
        (subject,predicate,object_id) = triple
        s_triples = self.writer.writer.graph.triples((None, OWL.annotatedSource, subject))
        s_bnodes = [s for s,p,o in s_triples]
        p_triples = self.writer.writer.graph.triples((None, OWL.annotatedProperty, predicate))
        p_bnodes = [s for s,p,o in p_triples]
        o_triples = self.writer.writer.graph.triples((None, OWL.annotatedTarget, object_id))
        o_bnodes = [s for s,p,o in o_triples]
        bnodes = set(s_bnodes) & set(p_bnodes) & set(o_bnodes)
        if len(bnodes) > 0:
            return list(bnodes)[0]

    def triples_involving_individual(self, ind_id, relation=None):
        # "involving" meaning individual (URI) is either subject or object
        graph = self.writer.writer.graph
        found_triples = list(graph.triples((ind_id, relation, None)))
        for t in graph.triples((None, relation, ind_id)):
            if t not in found_triples:
                found_triples.append(t)
        return found_triples

class AssocGoCamModel(GoCamModel):

    def __init__(self, modeltitle, assocs, connection_relations=None):
        GoCamModel.__init__(self, modeltitle, connection_relations)
        self.associations = assocs
        self.ontology = None
        # self.ro_ontology = ro_ontology
        self.extensions_mapper = None
        self.default_contributor = "http://orcid.org/0000-0002-6659-0416"

    def translate(self):
        # input_relations = {
        #     #TODO Create rule for deciding (MOD-specific?) whether to convert has_direct_input to has input
        #     # "has_direct_input": "RO:0002400",
        #     "has_direct_input": "RO:0002233",
        #     "has input": "RO:0002233"
        # }

        for a in self.associations:

            annoton = Annoton(a["subject"]["id"], [a])
            term = a["object"]["id"]

            # Paul's current rules are based on aspect, similar to rdfgen's current state, which may change
            # since relation is explicitly stated in GPAD
            # Standardize aspect using GPAD relations?

            # TODO stuff annot_date and contributors into annot data structure for reuse
            # Add evidence tied to axiom_ids
            annot_date = "{0:%Y-%m-%d}".format(datetime.datetime.strptime(a["date"], "%Y%m%d"))
            source_line = a["source_line"].rstrip().replace("\t", " ")
            # contributors = handle_annot_properties() # Need annot_properties to be parsed w/ GpadParser first
            contributors = []
            if "annotation_properties" in a and "contributor" in a["annotation_properties"]:
                contributors = a["annotation_properties"]["contributor"]
            if len(contributors) == 0:
                contributors = [self.default_contributor]

            annotation_extensions = get_annot_extensions(a)

            # Translate extension - maybe add function argument for custom translations?
            if len(annotation_extensions) == 0:
                self.translate_primary_annotation(a, annoton)
            else:
                aspect = self.extensions_mapper.go_aspector.go_aspect(term)

                for uo in annotation_extensions['union_of']:
                    int_bits = []
                    for rel in uo["intersection_of"]:
                        int_bits.append("{}({})".format(rel["property"], rel["filler"]))
                    ext_str = ",".join(int_bits)

                    anchor_uri = self.translate_primary_annotation(a, annoton)
                    intersection_extensions = self.extensions_mapper.dedupe_extensions(uo['intersection_of'])
                    is_cool = self.extensions_mapper.annot_following_rules(intersection_extensions, aspect)
                    if is_cool:
                        logger.debug("GOOD: {}".format(ext_str))
                        for rel in intersection_extensions:
                            ext_relation = rel["property"]
                            ext_target = rel["filler"]
                            if ext_relation in INPUT_RELATIONS:
                                logger.debug("Adding connection {} {} {}".format(annoton.enabled_by, ext_relation, ext_target))
                                target_gene_id = self.declare_individual(ext_target)
                                annoton.individuals[ext_target] = target_gene_id
                                axiom_id = self.find_or_create_axiom(anchor_uri, URIRef(expand_uri_wrapper(INPUT_RELATIONS[ext_relation])), target_gene_id)
                                self.add_evidence(axiom_id, a["evidence"]["type"],
                                                  a["evidence"]["has_supporting_reference"],
                                                  contributors=contributors,
                                                  date=annot_date,
                                                  comment=source_line)
                                # Nice-to-have functions:
                                # add_axiom(triple, evidence=[])
                                # add_evidence_to_axiom(axiom_id, evidence)
                            elif ext_relation in HAS_REGULATION_TARGET_RELATIONS:
                                buckets = has_regulation_target_bucket(self.ontology, term)
                                print("{} {}: {}".format(annoton.enabled_by, term, buckets))
                    else:
                        logger.debug("BAD: {}".format(ext_str))
        self.extensions_mapper.go_aspector.write_cache()

    def translate_primary_annotation(self, annotation, annoton):
        term = annotation["object"]["id"]

        anchor_uri = None
        axiom_ids = []
        for q in annotation["qualifiers"]:
            if q == "enables":
                axiom_id = self.find_or_create_axiom(term, ENABLED_BY, annoton.enabled_by, annoton=annoton)
                # Get enabled_by URI (owl:annotatedTarget) using axiom_id (a hack because I'm still using Annoton object with gene_connections)
                enabled_by_uri = list(self.writer.writer.graph.triples((axiom_id, OWL.annotatedTarget, None)))[0][2]
                anchor_uri = list(self.writer.writer.graph.triples((axiom_id, OWL.annotatedSource, None)))[0][2]
                annoton.individuals[annoton.enabled_by] = enabled_by_uri
                axiom_ids.append(axiom_id)
            elif q == "involved_in":
                # Try to find chain of two connected triples # TODO: Write function to find chain of any length
                query_pair = TriplePair((upt.molecular_function, ENABLED_BY, annoton.enabled_by), (upt.molecular_function, PART_OF, term), connecting_entity=upt.molecular_function)
                query_pair_collection = TriplePairCollection()
                query_pair_collection.chain_collection.append(query_pair)
                result_pair_collection = TriplePatternFinder().find_connected_pattern(self, query_pair_collection)
                # From result_pair_collection, we need axiom_ids of both triples (self.find_bnode()) and anchor_uri (pair.connecting_entity_uri())
                if result_pair_collection is not None:
                    # Only need one
                    triple_pair = result_pair_collection.chain_collection[0]
                    # we need axiom_ids of both triples
                    for tr in triple_pair.triples:
                        axiom_ids.append(self.find_bnode(tr))
                    anchor_uri = triple_pair.connecting_entity_uri(self)  # Root MF will be the connecting entity
                else:
                    mf_root_uri = self.declare_individual(upt.molecular_function)
                    gp_uri = self.declare_individual(annoton.enabled_by)
                    term_uri = self.declare_individual(term)
                    axiom_id = self.add_axiom(self.writer.emit(mf_root_uri, ENABLED_BY, gp_uri))
                    axiom_ids.append(axiom_id)
                    # Get enabled_by URI (owl:annotatedTarget) using axiom_id
                    enabled_by_uri = list(self.writer.writer.graph.triples((axiom_id, OWL.annotatedTarget, None)))[0][2]
                    annoton.individuals[annoton.enabled_by] = enabled_by_uri
                    anchor_uri = mf_root_uri
                    axiom_ids.append(self.add_axiom(self.writer.emit(mf_root_uri, PART_OF, term_uri)))
            elif q in ACTS_UPSTREAM_OF_RELATIONS:
                # Look for existing GP <- enabled_by [root MF] -> causally_upstream_of BP
                causally_relation = ENABLES_O_RELATION_LOOKUP[ACTS_UPSTREAM_OF_RELATIONS[q]]
                causally_relation_uri = URIRef(expand_uri_wrapper(causally_relation))
                query_pair = TriplePair((upt.molecular_function, ENABLED_BY, annoton.enabled_by),
                                        (upt.molecular_function, causally_relation_uri, term),
                                        connecting_entity=upt.molecular_function)
                query_pair_collection = TriplePairCollection()
                query_pair_collection.chain_collection.append(query_pair)
                result_pair_collection = TriplePatternFinder().find_connected_pattern(self, query_pair_collection)
                if result_pair_collection is not None:
                    triple_pair = result_pair_collection.chain_collection[0]
                    for tr in triple_pair.triples:
                        axiom_ids.append(self.find_bnode(tr))
                    anchor_uri = triple_pair.connecting_entity_uri(self)
                else:
                    # If no chain found create one.
                    mf_root_uri = self.declare_individual(upt.molecular_function)
                    gp_uri = self.declare_individual(annoton.enabled_by)
                    term_uri = self.declare_individual(term)
                    axiom_id = self.add_axiom(self.writer.emit(mf_root_uri, ENABLED_BY, gp_uri))
                    axiom_ids.append(axiom_id)
                    anchor_uri = mf_root_uri
                    axiom_ids.append(self.add_axiom(self.writer.emit(mf_root_uri, causally_relation_uri, term_uri)))
            elif q == "NOT":
                # Try it in UI and look at OWL
                do_stuff = 1
            else:
                relation_uri = URIRef(expand_uri_wrapper(self.relations_dict[q]))
                # TODO: should check that existing axiom/triple isn't connected to anything else; length matches exactly
                axiom_id = self.find_or_create_axiom(annoton.enabled_by, relation_uri, term)
                # Get enabled_by URI (owl:annotatedSource) using axiom_id
                enabled_by_uri = list(self.writer.writer.graph.triples((axiom_id, OWL.annotatedSource, None)))[0][2]
                annoton.individuals[annoton.enabled_by] = enabled_by_uri
                term_uri = list(self.writer.writer.graph.triples((axiom_id, OWL.annotatedTarget, None)))[0][2]
                # anchor_uri = enabled_by_uri
                anchor_uri = term_uri
                axiom_ids.append(axiom_id)

        # Add evidence tied to axiom_ids
        annot_date = "{0:%Y-%m-%d}".format(datetime.datetime.strptime(annotation["date"], "%Y%m%d"))
        source_line = annotation["source_line"].rstrip().replace("\t", " ")
        # contributors = handle_annot_properties() # Need annot_properties to be parsed w/ GpadParser first
        contributors = []
        if "annotation_properties" in annotation and "contributor" in annotation["annotation_properties"]:
            contributors = annotation["annotation_properties"]["contributor"]
        if len(contributors) == 0:
            # contributors = [""]
            contributors = [self.default_contributor]
        for a_id in axiom_ids:
            self.add_evidence(a_id, annotation["evidence"]["type"],
                              annotation["evidence"]["has_supporting_reference"],
                              contributors=contributors,
                              date=annot_date,
                              comment=source_line)

        return anchor_uri


def get_annot_extensions(annot):
    if "object_extensions" in annot:
        return annot["object_extensions"]
    elif "extensions" in annot["object"]:
        return annot["object"]["extensions"]
    return {}


class ReferencePreference():
    def __init__(self):
        # List order in python should be persistent
        self.order_of_prefix_preference = [
            "PMID",
            "GO_REF",
            "doi"
        ]

    def pick(self, references):
        for pfx in self.order_of_prefix_preference:
            for ref in references:
                if ref.startswith(pfx):
                    return ref


class GoCamEvidence():
    def __init__(self, code, references, contributors=[], date="", comment=""):
        self.evidence_code = code
        self.references = references
        self.date = date
        self.contributors = contributors
        self.comment = comment
        self.id = None


class CamTurtleRdfWriter(TurtleRdfWriter):
    def __init__(self, modeltitle):
        self.base = genid(base="http://model.geneontology.org")
        self.graph = rdflib.Graph(identifier=self.base)
        self.graph.bind("owl", OWL)
        self.graph.bind("obo", "http://purl.obolibrary.org/obo/")
        self.graph.bind("dc", DC)
        self.graph.bind("rdfs", RDFS)

        self.graph.add((self.base, RDF.type, OWL.Ontology))

        # Model attributes TODO: Should move outside init
        self.graph.add((self.base, URIRef("http://purl.org/pav/providedBy"), Literal("http://geneontology.org")))
        self.graph.add((self.base, DC.date, Literal(str(now.year) + "-" + str(now.month) + "-" + str(now.day))))
        self.graph.add((self.base, DC.title, Literal(modeltitle)))
        self.graph.add((self.base, DC.contributor, Literal("http://orcid.org/0000-0002-6659-0416"))) #TODO
        self.graph.add((self.base, URIRef("http://geneontology.org/lego/modelstate"), Literal("development")))
        self.graph.add((self.base, OWL.versionIRI, self.base))
        self.graph.add((self.base, OWL.imports, URIRef("http://purl.obolibrary.org/obo/go/extensions/go-lego.owl")))


class AnnotonCamRdfTransform(CamRdfTransform):
    def __init__(self, writer=None):
        CamRdfTransform.__init__(self, writer)
        self.annotons = []
        self.classes = []
        self.evidences = []
        self.ev_ids = []
        self.bp_id = None

    # TODO Remove "find" feature
    def find_or_create_evidence_id(self, evidence):
        for existing_evidence in self.evidences:
            if evidence.evidence_code == existing_evidence.evidence_code and set(evidence.references) == set(existing_evidence.references):
                if existing_evidence.id is None:
                    existing_evidence.id = genid(base=self.writer.base + '/')
                    self.ev_ids.append(existing_evidence.id)
                return existing_evidence.id
        return self.create_evidence(evidence)

    def create_evidence(self, evidence):
        # Use/figure out standard for creating URIs
        # Find minerva code to generate URI, add to Noctua doc
        ev_id = genid(base=self.writer.base + '/')
        evidence.id = ev_id
        # ev_cls = self.eco_class(self.uri(evidence.evidence_code))
        # ev_cls = self.eco_class(evidence.evidence_code) # This is already ECO:##### due to a GPAD being used
        ev_cls = self.uri(evidence.evidence_code)
        self.emit_type(ev_id, OWL.NamedIndividual)
        self.emit_type(ev_id, ev_cls)
        self.emit(ev_id, DC.date, Literal(evidence.date))
        for c in evidence.contributors:
            self.emit(ev_id, DC.contributor, Literal(c))
        ref_to_emit = ReferencePreference().pick(evidence.references)
        o = Literal(ref_to_emit)  # Needs to go into Noctua like 'PMID:####' rather than full URL
        self.emit(ev_id, HAS_SUPPORTING_REFERENCE, o)
        self.evidences.append(evidence)
        return evidence.id

    # Use only for OWLAxioms
    # There are two of these methods. AnnotonCamRdfTransform.find_bnode and GoCamModel.find_bnode. Which one is used?
    def find_bnode(self, triple):
        (subject,predicate,object_id) = triple
        s_triples = self.writer.graph.triples((None, OWL.annotatedSource, subject))
        s_bnodes = [s for s,p,o in s_triples]
        p_triples = self.writer.graph.triples((None, OWL.annotatedProperty, predicate))
        p_bnodes = [s for s,p,o in p_triples]
        o_triples = self.writer.graph.triples((None, OWL.annotatedTarget, object_id))
        o_bnodes = [s for s,p,o in o_triples]
        bnodes = set(s_bnodes) & set(p_bnodes) & set(o_bnodes)
        if len(bnodes) > 0:
            return list(bnodes)[0]

    def emit_axiom(self, source_id, property_id, target_id):
        stmt_id = self.blanknode()
        self.emit_type(stmt_id, OWL.Axiom)
        self.emit(stmt_id, OWL.annotatedSource, source_id)
        self.emit(stmt_id, OWL.annotatedProperty, property_id)
        self.emit(stmt_id, OWL.annotatedTarget, target_id)
        return stmt_id

    def find_annotons(self, enabled_by, annotons_list=None):
        found_annotons = []
        if annotons_list is not None:
            annotons = annotons_list
        else:
            annotons = self.annotons
        for annoton in annotons:
            if annoton.enabled_by == enabled_by:
                found_annotons.append(annoton)
        return found_annotons

    def add_individual(self, individual_id, annoton):
        obj_uri = self.uri(individual_id)
        if individual_id not in annoton.individuals:
            tgt_id = genid(base=self.writer.base + '/')
            annoton.individuals[individual_id] = tgt_id
            self.emit_type(tgt_id, obj_uri)
            self.emit_type(tgt_id, OWL.NamedIndividual)
        else:
            tgt_id = annoton.individuals[individual_id]
