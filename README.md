# gocamgen
Base repo for constructing GO-CAM model RDF

## Installation
```
pip install gocamgen
```

## Usage
```
from gocamgen.gocamgen import GoCamModel

model = GoCamModel("model title")
model.declare_class("PomBase:SPBC12C2.02c")
uri_a = model.declare_individual("GO:0016757")
uri_b = model.declare_individual("PomBase:SPBC12C2.02c")
axiom = model.add_axiom(uri_a, URIRef(expand_uri("RO:0002333")), uri_b)
model.add_evidence(axiom, "EXP", "PMID:1234567")

model.write("output_file.ttl")
```

## Quick generation of models from GPAD
Specify source GPAD file. All possible models will be generated and exported to `.ttl`.
```
python3 gen_models_by_gene.py --gpad_file wb.gpad
```
Additionally, a gene product identifier can be specified to only translate and export that GP's model.
```
python3 gen_models_by_gene.py --gpad_file wb.gpad --specific_gene WB:WBGene00004055
```
In general, annotation lines will be grouped by gene product identifier (col 2) with some lines filtered out due to various evidence code/reference rules.
