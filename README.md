# gocamgen
Base repo for constructing GO-CAM model RDF

## Installation
```
pip install gocamgen
```

## Usage
```
from gocamgen.gocamgen import GoCamModel

model = GoCamModel("output_file.ttl")
model.declare_class("PomBase:SPBC12C2.02c")
uri_a = model.declare_individual("GO:0016757")
uri_b = model.declare_individual("PomBase:SPBC12C2.02c")
axiom = model.add_axiom(uri_a, URIRef(expand_uri("RO:0002333")), uri_b)
model.add_evidence(axiom, "EXP", "PMID:1234567")

with open(model.filepath, 'wb') as f:
    model.writer.writer.serialize(destination=f)
```
