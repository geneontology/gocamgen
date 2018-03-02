from ontobio.rdfgen.assoc_rdfgen import CamRdfTransform, TurtleRdfWriter, RdfTransform, genid
from ontobio.vocabulary.relations import OboRO, Evidence
from ontobio.vocabulary.upper import UpperLevel
from ontobio.ontol_factory import OntologyFactory
from ontobio.assoc_factory import AssociationSetFactory
from prefixcommons.curie_util import expand_uri
from rdflib.namespace import OWL, RDF
from rdflib import Literal
from rdflib.term import URIRef
from rdflib.namespace import Namespace
import rdflib
import logging
import argparse
import datetime
import os.path as path

# logging.basicConfig(level=logging.INFO)

ro = OboRO()
evt = Evidence()
upt = UpperLevel()
LEGO = Namespace("http://geneontology.org/lego/")
LAYOUT = Namespace("http://geneontology.org/lego/hint/layout/")
PAV = Namespace('http://purl.org/pav/')
DC = Namespace("http://purl.org/dc/elements/1.1/")

# Stealing a lot of code for this from ontobio.rdfgen:
# https://github.com/biolink/ontobio

HAS_SUPPORTING_REFERENCE = URIRef(expand_uri(evt.has_supporting_reference, cmaps=[evt._prefixmap]))
ENABLED_BY = URIRef(expand_uri(ro.enabled_by))
ENABLES = URIRef(expand_uri(ro.enables))
INVOLVED_IN = URIRef(expand_uri(ro.involved_in))
PART_OF = URIRef(expand_uri(ro.part_of))
OCCURS_IN = URIRef(expand_uri(ro.occurs_in))
COLOCALIZES_WITH = URIRef(expand_uri(ro.colocalizes_with))
MOLECULAR_FUNCTION = URIRef(expand_uri(upt.molecular_function))
REGULATES = URIRef("http://purl.obolibrary.org/obo/RO_0002211")

now = datetime.datetime.now()

class Annoton():
    def __init__(self, gene_info, subject_id):
        self.enabled_by = subject_id
        self.molecular_function = self.get_aspect_object(gene_info, "molecular_function")
        self.cellular_component = self.get_aspect_object(gene_info, "cellular_component")
        self.biological_process = self.get_aspect_object(gene_info, "bp")
        self.connections = gene_info["connections"]
        self.individuals = {}

    def get_aspect_object(self, gene_info, aspect):
        if aspect in gene_info:
            return gene_info[aspect]

class GoCamModel():
    relations_dict = {
        "has_direct_input" : "RO:0002400",
        "has input" : "RO:0002233",
        "has_regulation_target" : "RO:0002211", # regulates
        "regulates_activity_of" : "RO:0002578", # directly regulates
        "with_support_from" : "RO:0002233" # has input
    }

    def __init__(self, filepath, connection_relations=None):
        self.filepath = filepath
        self.modeltitle = path.basename(self.filepath)
        # if self.modeltitle.endswith(".ttl"):
        if path.splitext(self.filepath)[1] != ".ttl":
            self.filepath += ".ttl"
        else:
            self.modeltitle = self.modeltitle[:-4]
        cam_writer = CamTurtleRdfWriter(self.modeltitle)
        self.writer = AnnotonCamRdfTransform(cam_writer)
        self.classes = []
        self.individuals = {}   # Maintain entity-to-IRI dictionary. Prevents dup individuals but we may want dups?
        if connection_relations is None:
            self.connection_relations = GoCamModel.relations_dict
        else:
            self.connection_relations = connection_relations
        self.declare_properties()

    def declare_properties(self):
        # AnnotionProperty
        self.writer.emit_type(URIRef("http://geneontology.org/lego/evidence"), OWL.AnnotationProperty)
        self.writer.emit_type(URIRef("http://geneontology.org/lego/hint/layout/x"), OWL.AnnotationProperty)
        self.writer.emit_type(URIRef("http://geneontology.org/lego/hint/layout/y"), OWL.AnnotationProperty)
        self.writer.emit_type(URIRef("http://purl.org/pav/providedBy"), OWL.AnnotationProperty)

    def declare_class(self, class_id):
        self.writer.emit_type(URIRef("http://identifiers.org/" + class_id), OWL.Class)
        self.classes.append(class_id)

    def declare_individual(self, entity_id):
        entity = genid(base=self.writer.writer.base + '/')
        self.writer.emit_type(entity, self.writer.uri(entity_id))
        self.writer.emit_type(entity, OWL.NamedIndividual)
        self.individuals[entity_id] = entity
        return entity

    def add_axiom(self, statement):
        (source_id, property_id, target_id) = statement
        stmt_id = self.find_bnode(statement)
        if stmt_id is None:
            stmt_id = self.writer.blanknode()
            self.writer.emit_type(stmt_id, OWL.Axiom)
        self.writer.emit(stmt_id, OWL.annotatedSource, source_id)
        self.writer.emit(stmt_id, OWL.annotatedProperty, property_id)
        self.writer.emit(stmt_id, OWL.annotatedTarget, target_id)
        return stmt_id

    def add_evidence(self, axiom, evidence_code, references):
        ev = {'type' : evidence_code,
            'has_supporting_reference' : references}
        # Try finding existing evidence object containing same type and references
        ev_id = self.writer.find_or_create_evidence_id(ev)
        self.writer.emit(axiom, URIRef("http://geneontology.org/lego/evidence"), ev_id)

    def add_connection(self, gene_connection, source_annoton):
        # Switching from reusing existing activity node from annoton to creating new one for each connection - Maybe SPARQL first to check if annoton activity already used for connection?
        # Check annoton for existing activity.
        # if gene_connection.object_id in source_annoton.individuals:
        #     # If exists and activity has connection relation,
        #     # Look for two triples: (gene_connection.object_id, ENABLED_BY, source_annoton.enabled_by) and (gene_connection.object_id, connection_relations, anything)
        source_id = None
        uri_list = self.uri_list_for_individual(gene_connection.object_id)
        for u in uri_list:
            if gene_connection.relation in self.connection_relations:
                rel = URIRef(expand_uri(self.connection_relations[gene_connection.relation]))
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
        property_id = URIRef(expand_uri(self.connection_relations[gene_connection.relation]))
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

class CamTurtleRdfWriter(TurtleRdfWriter):
    def __init__(self, modeltitle):
        self.base = genid(base="http://model.geneontology.org")
        self.graph = rdflib.Graph(identifier=self.base)
        self.graph.bind("owl", OWL)
        self.graph.bind("obo", "http://purl.obolibrary.org/obo/")
        self.graph.bind("dc", DC)

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
    
    def find_or_create_evidence_id(self, evidence):
        for existing_evidence in self.evidences:
            if evidence["type"] == existing_evidence["type"] and set(evidence["has_supporting_reference"]) == set(existing_evidence["has_supporting_reference"]):
                if "id" not in existing_evidence:
                    existing_evidence["id"] = genid(base=self.writer.base + '/')
                    self.ev_ids.append(existing_evidence["id"])
                return existing_evidence["id"]
        ev_id = genid(base=self.writer.base + '/')
        evidence["id"] = ev_id
        ev_cls = self.eco_class(self.uri(evidence['type']))
        self.emit_type(ev_id, OWL.NamedIndividual)
        self.emit_type(ev_id, ev_cls)
        if 'with_support_from' in evidence:
            for w in evidence['with_support_from']:
                self.emit(ev_id, self.uri(evt.evidence_with_support_from), self.uri(w))
        for ref in evidence['has_supporting_reference']:
            o = Literal(ref) # Needs to go into Noctua like 'PMID:####' rather than full URL
            self.emit(ev_id, HAS_SUPPORTING_REFERENCE, o)
        if 'with_support_from' in evidence:
            for ref in evidence['with_support_from']:
                self.emit(ev_id, self.uri(evt.evidence_with_support_from), self.uri(ref))
        self.evidences.append(evidence)
        return evidence["id"]

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
