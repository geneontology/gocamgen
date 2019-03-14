import gocamgen
import unittest

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
        # model = AssocGoCamModel(gene, assocs)
        # model.extensions_mapper = ext_mapper
        # model.translate()

        # Get model.writer.graph whatever and check for loose evidence (not attached to axioms)
        # Orphaned evidence - how to I find these in debugger? They're in rdflib writer somewhere

        self.assertEqual(1, 1)

if __name__ == '__main__':
    unittest.main()