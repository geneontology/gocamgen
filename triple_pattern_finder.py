#from gocamgen.gocamgen import AssocGoCamModel

### Simple example: 'GP --enabled_by--> MF --part_of--> BP'
### First find all 'GP --enabled_by--> MF' triples
### Then find a 'MF --part_of--> BP' where URI of MF equals URI of MF of any of triples in first set.
### So subsequent links (triples) in the chain will be passed the whole set of candidate preceding chains
### As long as current triple matches next link in pattern and URI of subject matches URI of previous triple's object, keep going


class TriplePattern:
    def __init__(self, ordered_triples):
        self.ordered_triples = ordered_triples


class ConnectedTriplePattern(TriplePattern):
    def __init__(self, ordered_triples, connecting_entities):
        TriplePattern.__init__(self, ordered_triples)
        self.connecting_entities = connecting_entities


class TriplePair:
    def __init__(self, triple_a, triple_b, connecting_entity):
        self.triples = (triple_a, triple_b)
        # Check to make sure connecting_entity is in both triples
        self.connecting_entity = connecting_entity

    def is_connected_by_uri(self, model):
        #TODO Check that subjects and objects are URIs
        for uri in [self.triples[0][0], self.triples[0][2]]:
            class_curie = model.class_for_uri(uri)
            if class_curie == self.connecting_entity and uri in self.triples[1]:
                return True
        return False

    def connecting_entity_uri(self, model):
        # entity is either [0] or [2] of either self.triples[0] or self.triples[1]
        # Scan these, checking class_for_uri == connecting_entity
        for uri in [self.triples[0][0], self.triples[0][2], self.triples[1][0], self.triples[1][2]]:
            class_curie = model.class_for_uri(uri)
            if class_curie == self.connecting_entity:
                return uri
        return None

# AKA pattern? Holds all pairs and will be query input for TriplePatternFinder
class TriplePairCollection:
    def __init__(self):
        self.chain_collection = []


class TriplePatternFinder:

    def __init__(self):
        self.thing = 1

    # TODO: Add 'exact' arg for requiring that found pattern is of exact length as query (nothing continuing on either side)
    def find_pattern_recursive(self, model, pattern: TriplePattern, candidate_chains=[]):
        # break down pattern into component triples
        p = pattern.ordered_triples[0]
        found_triples = model.triples_by_ids(*p)
        # print(found_triples)

        if len(found_triples) > 0:
            if len(candidate_chains) == 0:
                candidate_chains = [[t] for t in found_triples]
            else:
                candidate_chains_local = []
                for chain in candidate_chains:
                    for triple in found_triples:
                        #TODO: Make more selective (e.g. subject or object must be in previous triple)
                        candidate_chains_local.append(chain + [triple])
                candidate_chains = candidate_chains_local

            if len(pattern.ordered_triples) > 1:
                rest_of_pattern = TriplePattern(pattern.ordered_triples[1:len(pattern.ordered_triples)])
                candidate_chains = self.find_pattern_recursive(model, rest_of_pattern, candidate_chains)
        else:
            candidate_chains = []

        return candidate_chains

    # def find_connected_pattern(self, model: AssocGoCamModel, pair_collection: TriplePairCollection, candidate_chains=[]):
    def find_connected_pattern(self, model, pair_collection: TriplePairCollection):
        # pair_collection is a TriplePairCollection w/ chain_collection of TriplePair[]. Each of these TriplePairs
        # consists of two triples, e.g. ("MF-term", relation_uri, "GP-class")
        # Returned output should be a TriplePairCollection w/ chain_collection of TriplePairs consisting of all-URI triples.
        connected_pair_collection = TriplePairCollection()  # will be triples of URIs
        for pair in pair_collection.chain_collection:
            pattern = TriplePattern([*pair.triples])
            found_pairs = self.find_pattern_recursive(model, pattern)  # only using this 2-deep
            connected_pair = None
            for fpair in found_pairs:
                # fpair is list of two triples
                #TODO put into TriplePairs
                f_triple_pair = TriplePair(*fpair, connecting_entity=pair.connecting_entity)
                if f_triple_pair.is_connected_by_uri(model):
                    connected_pair = f_triple_pair
                    break
            if connected_pair:
                connected_pair_collection.chain_collection.append(connected_pair)
            else:
                return None
        return connected_pair_collection
