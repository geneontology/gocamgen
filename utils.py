from ontobio.rdfgen.assoc_rdfgen import prefix_context
from prefixcommons.curie_util import expand_uri, contract_uri

def expand_uri_wrapper(id):
    uri = expand_uri(id, cmaps=[prefix_context])
    return uri

def contract_uri_wrapper(id):
    uri = contract_uri(id, cmaps=[prefix_context])
    return uri