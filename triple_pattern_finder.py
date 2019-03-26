from gocamgen.gocamgen import AssocGoCamModel

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

class TriplePatternFinder:

    def __init__(self):
        self.thing = 1

    # What should this return? Ordered list of triples filled with URIs?
    # def find_pattern(self, model, pattern : TriplePattern):
    #     # break down pattern into component triples
    #     connecting_uri = None
    #     candidate_chains = [] # [[chain1
    #     for p in pattern.ordered_triples:
    #         found_triples = model.triples_by_ids(*p)
    #         print(found_triples)
    #
    #         if len(found_triples) > 0:
    #             # if p is last in chain,
    #             for mf_triple in found_triples:
    #                 # If none of triple URIs in previous triple, throw out chain
    #                 found_mf_uri = mf_triple[0]
    #                 found_triples = self.triples_by_ids(found_mf_uri, PART_OF, term)
    #         #         if len(found_triples) > 0:
    #         #             # Found both triples. Add em and get out of here
    #         #             bp_triple = found_triples[0]
    #         #             anchor_uri = found_mf_uri
    #         #             make_new = False
    #         #             break
    #     # Can we support non-linear (branching) patterns?
    #     # Can we accomplish same thing by constructing sub-model from query pattern and checking against model?
    #     do_stuff = 1

    def find_pattern_recursive(self, model, pattern: TriplePattern, candidate_chains=[]):
        # break down pattern into component triples
        connecting_uri = None
        # for p in pattern.ordered_triples:
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

        return candidate_chains

            # candidate_chains_local = []
            # for triple in found_triples:
            #     for c in candidate_chains:
            #         if triple[0] in c[-1] or triple[1] in c[-1] or triple[2] in c[-1]:
            #             candidate_chains_local.append(c + triple)
            #         else:
            #             candidate_chains_local.append(c)
            #     # If none of triple URIs in previous triple, throw out chain
            #     found_mf_uri = mf_triple[0]
            #     found_triples = self.triples_by_ids(found_mf_uri, PART_OF, term)
