from ontobio.rdfgen.assoc_rdfgen import prefix_context
from prefixcommons.curie_util import expand_uri, contract_uri


def expand_uri_wrapper(id):
    uri = expand_uri(id, cmaps=[prefix_context])
    return uri


def contract_uri_wrapper(id):
    uri = contract_uri(id, cmaps=[prefix_context])
    return uri


def sort_terms_by_ontology_specificity(terms):
    # Used primarily for sorting occurs_in annotation extensions
    # What's first? EMAPA or UBERON? Shouldn't matter for extensions since assertion
    # should chain occurs_in to both EMAPA and UBERON (They would be split into separate assertions or thrown out)
    ONTOLOGY_ORDER = {'GO': 1, 'CL': 2, 'WBbt': 3, 'EMAPA': 4, 'UBERON': 5}  # From most specific to most general
    terms.sort(key=lambda t: ONTOLOGY_ORDER[t.split(":")[0]])

    return terms

