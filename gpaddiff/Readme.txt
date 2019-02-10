GpadDiffer compares to gpad files, 'source' and 'target'.
It looks through each line in source and tries to find a match in target.
Run by
> -g1 source.gpad -g2 target.gpad

It will generate an output file called 'compare.txt' in the current directory.
The output file will contains the tab delimited gpad from source with a new column first column
The first column indicates whether the line in the gpad had a match in the target file.
Matches range from 'no match' to 'ont qualifier eco with ref match' depending on how many parameters matched
'ont' means the gene id and the ontology id matched.  etc.  If there is no ont match, no match will be reported.