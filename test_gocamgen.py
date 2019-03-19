import gocamgen
import unittest
import logging
from gen_models_by_gene import WBFilterRule, AssocExtractor, GoCamBuilder

# logging.basicConfig(level=logging.DEBUG)
# logger = logging.getLogger("gen_models_by_gene")
# logger.setLevel(logging.INFO)

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

    def test_evidence(self):
        gpad_file = "wb.gpad"
        # gpad_file = "wb.gpad.WBGene00000903"
        test_gene = "WB:WBGene00000903"

        filter_rule = WBFilterRule()
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

        self.assertEqual(1, 1)

if __name__ == '__main__':
    unittest.main()