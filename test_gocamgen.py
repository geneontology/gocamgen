import gocamgen
from gocamgen.gocamgen import expand_uri_wrapper, contract_uri_wrapper
import unittest
import logging
from gen_models_by_gene import WBFilterRule, MGIFilterRule, AssocExtractor, GoCamBuilder
from triple_pattern_finder import TriplePattern, TriplePatternFinder, TriplePair, TriplePairCollection
from rdflib.term import URIRef

# logging.basicConfig(level=logging.DEBUG)
# logger = logging.getLogger("gen_models_by_gene")
# logger.setLevel(logging.DEBUG)


class TestReferencePreference(unittest.TestCase):

    def test_ref_picker(self):
        test_refs = [
            "GO_REF:0000483",
            "doi:485930",
            "WB_REF:WBPaper00003384",
            "PMID:9834189",
        ]
        ref_picker = gocamgen.gocamgen.ReferencePreference()

        result = ref_picker.pick(test_refs)
        self.assertEqual(result, "PMID:9834189")

class TestGoCamModel(unittest.TestCase):

    def test_triple_finder(self):
        gpad_file = "resources/test/wb_6498.gpad"
        test_gene = "WB:WBGene00006498"
        filter_rule = WBFilterRule()

        # gpad_file = "mgi_87859.gpad"
        # test_gene = "MGI:MGI:87859"
        # filter_rule = MGIFilterRule()

        extractor = AssocExtractor(gpad_file, filter_rule)
        assocs_by_gene = extractor.group_assocs()
        print("{} distinct genes".format(len(assocs_by_gene)))

        builder = GoCamBuilder()
        if test_gene not in assocs_by_gene:
            print("ERROR: specific gene {} not found in filtered annotation list".format(test_gene))
        else:
            model = builder.translate_to_model(test_gene, assocs_by_gene[test_gene])

            # Get model.writer.graph whatever and check for loose evidence (not attached to axioms)
            # Orphaned evidence - how to I find these in debugger? They're in rdflib writer somewhere

            # Can this work?
            #   (MF, ENABLED_BY, GP) & (Same MF, has input, GP)
            pattern_a = TriplePattern([("GO:0003674", gocamgen.gocamgen.ENABLED_BY, test_gene)])
            pattern_b = TriplePattern([("GO:0003674", URIRef(expand_uri_wrapper("BFO:0000050")), "GO:0019953")])
            whole_pattern = TriplePattern([("GO:0003674", gocamgen.gocamgen.ENABLED_BY, test_gene),
                                     ("GO:0003674", URIRef(expand_uri_wrapper("BFO:0000050")), "GO:0019953")])

            triple_finder = TriplePatternFinder()
            a_triples = triple_finder.find_pattern_recursive(model, pattern_a)
            print("A count: {}".format(len(a_triples)))
            b_triples = triple_finder.find_pattern_recursive(model, pattern_b)
            print("B count: {}".format(len(b_triples)))
            found_chains = triple_finder.find_pattern_recursive(model, whole_pattern)
            # print(found_chains)
            print("Chain count: {}".format(len(found_chains)))
            for fc in found_chains:
                print(contract_uri_wrapper(model.individual_label_for_uri(fc[1][2])[0])[0])

            triple_pair = TriplePair(pattern_a.ordered_triples[0], pattern_b.ordered_triples[0], connecting_entity="GO:0003674")
            tp_collection = TriplePairCollection()
            tp_collection.chain_collection.append(triple_pair)
            uri_tp_collection = triple_finder.find_connected_pattern(model, tp_collection)
            for pair in uri_tp_collection.chain_collection:
                # print(model.class_for_uri(pair.triples[0][0]))
                # print(model.class_for_uri(pair.triples[0][2]))
                # print(model.class_for_uri(pair.triples[1][0]))
                # print(model.class_for_uri(pair.triples[1][2]))
                print(pair.triples[0][0])
                print(pair.triples[0][2])
                print(pair.triples[1][0])
                print(pair.triples[1][2])

        self.assertEqual(1, 1)

    def test_evidence(self):
        # gpad_file = "wb_903.gpad"
        # test_gene = "WB:WBGene00000903"
        # filter_rule = WBFilterRule()

        self.assertEqual(1, 1)

if __name__ == '__main__':
    unittest.main()