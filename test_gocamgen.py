import gocamgen
from gocamgen.gocamgen import expand_uri_wrapper, contract_uri_wrapper, ACTS_UPSTREAM_OF_RELATIONS, ENABLES_O_RELATION_LOOKUP
import unittest
import logging
from filter_rule import WBFilterRule, MGIFilterRule
from gen_models_by_gene import AssocExtractor, GoCamBuilder
from triple_pattern_finder import TriplePattern, TriplePatternFinder, TriplePair, TriplePairCollection
from rdflib.term import URIRef
from rdflib_sparql_wrapper import RdflibSparqlWrapper

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
    BUILDER = GoCamBuilder()  # Takes a sec to init so only make once

    def setUp(self):
        pass

    def gen_model(self, gpad_file, test_gene, filter_rule):
        extractor = AssocExtractor(gpad_file, filter_rule)
        assocs_by_gene = extractor.group_assocs()

        if test_gene not in assocs_by_gene:
            # print("ERROR: specific gene {} not found in filtered annotation list".format(test_gene))
            return None
        else:
            model = TestGoCamModel.BUILDER.translate_to_model(test_gene, assocs_by_gene[test_gene])
            return model

    def test_triple_finder(self):
        test_gene = "WB:WBGene00006498"
        model = self.gen_model(gpad_file="resources/test/wb_6498.gpad", test_gene=test_gene,
                          filter_rule=WBFilterRule())
        if model:
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
            # print("A count: {}".format(len(a_triples)))
            b_triples = triple_finder.find_pattern_recursive(model, pattern_b)
            # print("B count: {}".format(len(b_triples)))
            found_chains = triple_finder.find_pattern_recursive(model, whole_pattern)
            # print(found_chains)
            # print("Chain count: {}".format(len(found_chains)))
            # for fc in found_chains:
            #     print(contract_uri_wrapper(model.individual_label_for_uri(fc[1][2])[0])[0])

            triple_pair = TriplePair(pattern_a.ordered_triples[0], pattern_b.ordered_triples[0],
                                     connecting_entity="GO:0003674")
            tp_collection = TriplePairCollection()
            tp_collection.chain_collection.append(triple_pair)
            uri_tp_collection = triple_finder.find_connected_pattern(model, tp_collection)

            self.assertGreaterEqual(len(uri_tp_collection.chain_collection), 1)
        else:
            self.fail("Couldn't generate model for WB:WBGene00006498")

    def test_evidence(self):
        # gpad_file = "wb_903.gpad"
        # test_gene = "WB:WBGene00000903"
        # filter_rule = WBFilterRule()

        self.assertEqual(1, 1)

    def test_has_input(self):
        # See https://github.com/geneontology/gocamgen/issues/39#issuecomment-479988904 for background
        model = self.gen_model(gpad_file="resources/test/wb.gpad.WBGene00003167", test_gene="WB:WBGene00003167",
                          filter_rule=WBFilterRule())
        if model:
            # Look for translation of 'GO:0000977 has_direct_input(WB:WBGene00036254)'
            found_triples = model.triples_by_ids("GO:0000977", URIRef(expand_uri_wrapper("RO:0002233")),
                                                 "WB:WBGene00036254")
            self.assertGreaterEqual(len(found_triples), 1, "No has_input extensions translated")
        else:
            self.fail("Couldn't generate model for WB:WBGene00003167")

    def test_extension_pipe_separation(self):
        # See https://github.com/geneontology/gocamgen/issues/40
        model = self.gen_model(gpad_file="resources/test/wb.gpad.WBGene00003167", test_gene="WB:WBGene00003167",
                          filter_rule=WBFilterRule())

        if model:
            # Look for count of 'WB:WBGene00003167 contributes_to GO:0000977'
            found_triples = model.triples_by_ids("WB:WBGene00003167", gocamgen.gocamgen.CONTRIBUTES_TO,
                                                 "GO:0000977")
            self.assertGreaterEqual(len(found_triples), 3,
                                    "Less than 3 annotations for WB:WBGene00003167 contributes_to GO:0000977")
        else:
            self.fail("Couldn't generate model for WB:WBGene00003167")

    def test_no_dup_individuals(self):
        # See https://github.com/geneontology/gocamgen/issues/40
        model = self.gen_model(gpad_file="resources/test/mgi.gpa.MGI_2159711", test_gene="MGI:MGI:2159711",
                          filter_rule=MGIFilterRule())

        if model:
            # Look for 'MGI:MGI:2159711 PART_OF GO:0044297'. Should only be one.
            found_triples = model.triples_by_ids("MGI:MGI:2159711", gocamgen.gocamgen.PART_OF,
                                                 "GO:0044297")
            self.assertEqual(len(found_triples), 1)
        else:
            self.fail("Couldn't generate model for MGI:MGI:2159711")

        # Case of multiple experimental annotations to same GP and term.
        model = self.gen_model(gpad_file="resources/test/mgi.gpa.MGI_98956", test_gene="MGI:MGI:98956",
                               filter_rule=MGIFilterRule())

        if model:
            # This example has like 10-20 exp assertions for MGI:MGI:98956 to GO:0060070
            sparql_wrapper = RdflibSparqlWrapper()
            gp = "MGI:MGI:98956"
            term = "GO:0060070"
            qres = sparql_wrapper.find_involved_in(model.graph, gp, term)
            self.assertEqual(len(qres), 1)
        else:
            self.fail("Couldn't generate model for MGI:MGI:98956")


    def test_has_regulation_target(self):
        # Examples:
        # F - MGI:MGI:107771 GO:0005096 'has_regulation_target(MGI:MGI:97846)|has_regulation_target(MGI:MGI:2180784)'
        # P - WB:WBGene00013591 GO:0042594 'causally_upstream_of(GO:0001934),has_regulation_target(WB:WBGene00008480)'
        # Which has_regulation_target bucket does this fall into? None so far (GO:0042594 is "response to starvation")
        # bucket = gocamgen.gocamgen.has_regulation_target_bucket(ont, "GO:0001934")

        # model = gen_model(gpad_file="resources/test/wb.gpad.WBGene00013591", test_gene="WB:WBGene00013591",
        #                   filter_rule=WBFilterRule())
        model = self.gen_model(gpad_file="resources/test/wb.gpad.WBGene00003167", test_gene="WB:WBGene00003167",
                          filter_rule=WBFilterRule())

        self.assertEqual(1, 1)

    def test_acts_upstream_of(self):
        # TODO: Test for MGI:MGI:1206591
        self.assertEqual(1, 1)

        # what about MGI:MGI:1914305?
        # resources/test/mgi.gpa.MGI_1914305
        model = self.gen_model(gpad_file="resources/test/mgi.gpa.MGI_1914305", test_gene="MGI:MGI:1914305",
                               filter_rule=MGIFilterRule())

        if model:
            sparql_wrapper = RdflibSparqlWrapper()
            gp = "MGI:MGI:1914305"
            term = "GO:0007416"
            causally_relation = ENABLES_O_RELATION_LOOKUP[ACTS_UPSTREAM_OF_RELATIONS["acts_upstream_of_or_within"]]
            qres = sparql_wrapper.find_acts_upstream_of_translated(model.graph, gp, causally_relation, term)
            self.assertEqual(len(qres), 1)
        else:
            self.fail("Couldn't generate model for MGI:MGI:1914305")

    def test_sparql(self):
        # Just gonna see what we can sparql out of this guy.
        model = self.gen_model(gpad_file="resources/test/wb.gpad.WBGene00003167", test_gene="WB:WBGene00003167",
                               filter_rule=WBFilterRule())

        if model:
            sparql_wrapper = RdflibSparqlWrapper()
            gp = "WB:WBGene00003167"
            term = "GO:0007638"
            qres = sparql_wrapper.find_involved_in(model.graph, gp, term)

            self.assertEqual(len(qres), 1)
        else:
            self.fail("Couldn't generate model for WB:WBGene00003167")


if __name__ == '__main__':
    unittest.main()