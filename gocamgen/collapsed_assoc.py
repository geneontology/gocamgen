
class CollapsedAssociationSet:
    def __init__(self, associations):
        self.associations = associations
        self.collapsed_associations = []
        self.assoc_dict = {}

    def collapse_annotations(self):
        # Here we shall decide the distinct assertion instances going into the model
        # This will reduce/eliminate need to SPARQL model graph
        # Group by:
        # 		1. ID
        # 		2. qualifiers (normalize order; any array gotta do this)
        # 		3. primary term
        # 		4. With/From (if primary term is GO:0005515 or descendant)
        # 		5. Extensions
        # 	Collapse multiple:
        # 		1. Reference
        # 		2. Evidence Code
        #       3. With/From (if primary term is not GO:0005515 or descendant)
        # 		4. Source line
        # 		5. Date
        # 		6. Assigned by
        # 		7. Properties
        for a in self.associations:
            # Header
            subj_id = a["subject"]["id"]
            qualifiers = a["qualifiers"]
            term = a["object"]["id"]
            with_from = a["evidence"]["with_support_from"]
            extensions = get_annot_extensions(a)
            ca = self.find_or_create_collapsed_association(subj_id, qualifiers, term, with_from, extensions)

            # Line
            source_line = a["source_line"]
            references = a["evidence"]["has_supporting_reference"]
            evidence_code = a["evidence"]["type"]
            date = a["date"]
            assigned_by = a["provided_by"]
            association_line = {
                'source_line': source_line,
                'evidence': {
                    'type': evidence_code,
                    'has_supporting_reference': sorted(references)
                },
                'date': date,
                'provided_by': assigned_by
            }
            if "annotation_properties" in a:
                association_line["annotation_properties"] = a["annotation_properties"]
            if term != "GO:0005515":
                association_line['evidence']['with_support_from'] = sorted(with_from)
            ca.lines.append(association_line)

    def find_or_create_collapsed_association(self, subj_id, qualifiers, term, with_from, extensions):
        query_header = {
            'subject': {
                'id': subj_id
            },
            'qualifiers': sorted(qualifiers),
            'object': {
                'id': term
            },
            'object_extensions': extensions
        }
        # TODO: Get "GO:0005515" in onto.ancestors(term, reflexive=True) hooked up
        if term == "GO:0005515":
            query_header['evidence'] = {'with_support_from': sorted(with_from)}
        for ca in self.collapsed_associations:
            if ca.header == query_header:
                return ca
        new_ca = CollapsedAssociation(query_header)
        self.collapsed_associations.append(new_ca)
        return new_ca

    def __iter__(self):
        return iter(self.associations)


class CollapsedAssociation:
    def __init__(self, header):
        self.header = header
        self.lines = []

    def subject_id(self):
        if "subject" in self.header and "id" in self.header["subject"]:
            return self.header["subject"]["id"]

    def object_id(self):
        if "object" in self.header and "id" in self.header["object"]:
            return self.header["object"]["id"]

    def __str__(self):
        return "{} - {}".format(self.subject_id(), self.object_id())


def get_annot_extensions(annot):
    if "object_extensions" in annot:
        return annot["object_extensions"]
    elif "extensions" in annot["object"]:
        return annot["object"]["extensions"]
    return {}


def extract_properties(annot):
    cols = annot["source_line"].rstrip().split("\t")
    if len(cols) >= 12:
        prop_col = cols[11]
        props = prop_col.split("|")
        props_dict = {}
        for p in props:
            k, v = p.split("=")
            if k in props_dict:
                props_dict[k].append(v)
            else:
                props_dict[k] = [v]
        annot["annotation_properties"] = props_dict
    return annot
