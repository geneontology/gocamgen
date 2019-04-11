# gocamgen
Base repo for constructing GO-CAM model RDF

## Installation
```
pip install gocamgen
```

## Tests
```
python3 test_gocamgen.py
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

## Generating annotation extensions usage spreadsheet
To be clear, this creates a spreadsheet aggregating all extensions "patterns" (e.g. `happens_during(GO:P)`, `positively_regulates(GO:P)`) that are decided to be invalid according to the rules currently encoded in this same script. These same rules are used by the model generation import.

Example to generate spreadsheet from one GPAD file:
```
python3 gpad_extensions_mapper.py --filename wb.gpad --out_file bad_extensions.tsv
```
Example to generate spreadsheet (results compiled into one .tsv) from directory of GPADs:
```
python3 gpad_extensions_mapper.py --dir gpad_files/ --out_file bad_extensions.tsv
```
Snippet from created spreadsheet:

| Aspect | Total count | Extension                                 | UniProt | WB | ParkinsonsUK-UCL |
|--------|-------------|-------------------------------------------|---------|----|------------------|
| F      | 2           | happens_during(GO:P),happens_during(WBls) | 2       | 0  | 0                |
| F      | 2           | happens_during(WBls)                      | 0       | 2  | 0                |
| F      | 4           | regulates_activity_of(geneID)             | 0       | 4  | 0                |
| F      | 1           | occurs_in(WBbt),part_of(GO:P)             | 1       | 0  | 0                |
| F      | 3           | occurs_in(WBbt)                           | 2       | 1  | 0                |
| F      | 1           | activated_by(geneID)                      | 0       | 1  | 0                |

#### Drilling down to offending GPAD lines "by extension pattern"
If you'd like to get the actual GPAD containing the invalid extension pattern, just copy-paste the pattern into the `--pattern` option and rerun on the same GPAD(s):
```
python3 gpad_extensions_mapper.py --filename wb.gpad --pattern "regulates_activity_of(geneID)"
```
By default this will output the lines to a `[pattern].gpad` file like `regulates_activity_of(geneID).gpad` but this can be overidden with `--pattern_outfile`. There's also an option for running multiple patterns at once by plugging the path to a newline-separated list of patterns into the `--pattern_sourcefile`, which then creates multiple GPAD files named after each pattern.
