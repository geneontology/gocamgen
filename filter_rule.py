from ontobio.ecomap import EcoMap
from abc import ABC, abstractmethod

class FilterRule(ABC):
    # Default filter - Remove IEA, IBA, as well as IKRs originating from PAINT.
    def __init__(self, unwanted_evidence_codes=['IEA', 'IBA'],
                 unwanted_evi_code_ref_combos=[('IKR', 'PMID:21873635')],
                 required_attributes=None):
        self.unwanted_evidence_codes = unwanted_evidence_codes
        self.unwanted_evi_code_ref_combos = unwanted_evi_code_ref_combos
        if required_attributes is None:
            self.required_attributes = [{"provided_by": [self.mod_id()]}]
        else:
            self.required_attributes = required_attributes
        self.unwanted_properties = []

    @abstractmethod
    def mod_id(self):
        pass


class WBFilterRule(FilterRule):
    # So far same as default FilterRule.

    def mod_id(self):
        return "WB"


class MGIFilterRule(FilterRule):
    # Only provided_by=MGI lines are valid, also filtering out ISOs along with default FilterRule filters.
    def __init__(self):
        FilterRule.__init__(self)
        self.unwanted_evidence_codes.append('ISO')
        self.unwanted_properties = ["noctua-model-id"]

    def mod_id(self):
        return "MGI"

mod_filter_map = {
    'WB': WBFilterRule,
    'MGI': MGIFilterRule
}


class AssocFilter:
    def __init__(self, filter_rule : FilterRule):
        self.filter_rule = filter_rule
        self.ecomap = EcoMap()

    def validate_line(self, assoc):
        evi_code = self.ecomap.ecoclass_to_coderef(assoc["evidence"]["type"])[0]
        if evi_code in self.filter_rule.unwanted_evidence_codes:
            return False
        if len(self.filter_rule.unwanted_evi_code_ref_combos) > 0:
            references = assoc["evidence"]["has_supporting_reference"]
            for evi_ref_combo in self.filter_rule.unwanted_evi_code_ref_combos:
                if evi_ref_combo[0] == evi_code and evi_ref_combo[1] in references:
                    return False
        if len(self.filter_rule.unwanted_properties) > 0 and "annotation_properties" in assoc:
            for up in self.filter_rule.unwanted_properties:
                if up in assoc["annotation_properties"]:
                    return False
        if len(self.filter_rule.required_attributes) > 0:
            meets_requirement = False
            for attr in self.filter_rule.required_attributes:
                # a[attr] is dict
                for k in attr.keys():
                    if assoc[k] in attr[k]:
                        meets_requirement = True
            if not meets_requirement:
                return False
        return True