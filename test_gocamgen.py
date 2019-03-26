import gocamgen
from gocamgen.gocamgen import expand_uri_wrapper
import unittest
import logging
from gen_models_by_gene import WBFilterRule, MGIFilterRule, AssocExtractor, GoCamBuilder
from triple_pattern_finder import TriplePattern, TriplePatternFinder
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
        # gpad_file = "wb_903.gpad"
        # test_gene = "WB:WBGene00000903"
        # filter_rule = WBFilterRule()


        gpad_file = "mgi_87859.gpad"
        test_gene = "MGI:MGI:87859"
        filter_rule = MGIFilterRule()

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
            # pattern = TriplePattern([("WB:WBGene00000903", gocamgen.gocamgen.PART_OF, "GO:0005615")])
            # pattern = TriplePattern([("GO:0004672", gocamgen.gocamgen.ENABLED_BY, "MGI:MGI:87859")])

            pattern = TriplePattern([("GO:0004672", gocamgen.gocamgen.ENABLED_BY, "MGI:MGI:87859"), ("GO:0004672", URIRef(expand_uri_wrapper("RO:0002233")), "MGI:MGI:88508")])
            # input_pattern = TriplePattern([("GO:0004672", URIRef(expand_uri_wrapper("RO:0002233")), "MGI:MGI:88508")])
            # Can this work?
            #   (MF, ENABLED_BY, GP) & (Same MF, has input, GP)

            triple_finder = TriplePatternFinder()
            # triple_finder.find_pattern(model, pattern)
            found_chains = triple_finder.find_pattern_recursive(model, pattern)
            print(found_chains)
            # triple_finder.find_pattern(model, input_pattern)

        self.assertEqual(1, 1)

    def test_evidence(self):
        # gpad_file = "wb_903.gpad"
        # test_gene = "WB:WBGene00000903"
        # filter_rule = WBFilterRule()

        self.assertEqual(1, 1)

if __name__ == '__main__':
    unittest.main()