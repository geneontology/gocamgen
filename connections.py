

class GeneConnectionSet():
    def __init__(self):
        self.gene_connections = []

    def append(self, gene_connection):
        self.gene_connections.append(gene_connection)

    def contains(self, gene_connection):
        for gc in self.gene_connections:
            if gc.equals(gene_connection):
                return True
            if gene_connection.relation == "with_support_from" and (self.find(gene_connection.gp_a, gene_connection.gp_b) or self.find(gene_connection.gp_b, gene_connection.gp_a)):
                return True
        return False

    def merge(self, other_set):
        for connection in other_set.gene_connections:
            if connection.annotation["object"]["id"] == "GO:0005515": # Check if triple through extension already exists
                if self.find(connection.gp_a, connection.gp_b, "has_direct_input"):
                    continue
            if not self.contains(connection):
                self.append(connection)

    def find(self, gp_a, gp_b, relation=None):
        connections = []
        for connection in self.gene_connections:
            if connection.gp_a == gp_a and connection.gp_b == gp_b:
                if relation is None or connection.relation == relation:
                    connections.append(connection)
        return connections

class GeneConnection():
    def __init__(self, gene_a, gene_b, object_id, relation, annotation=None):
        self.gp_a = gene_a
        self.gp_b = gene_b
        self.object_id = object_id
        self.relation = relation
        self.annotation = annotation

    def equals(self, gene_connection):
        if self.gp_a == gene_connection.gp_a and self.gp_b == gene_connection.gp_b and self.relation == gene_connection.relation:
            return True
        else:
            return False

    def print_connection(self, a_set):
        # Need a better way to pass in GP info to get labels
        # a_set = afactory.create(ontology, subject_category='gene', object_category='function')  # Labels for debugging
        return self.gp_a + " (" + a_set.label(self.gp_a) + ") " + self.relation + " " + self.gp_b + " (" + a_set.label(self.gp_b) + ") through " + self.object_id + " (" + a_set.label(self.object_id) + ")"