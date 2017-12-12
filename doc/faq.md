# Frequently asked questions

## What if I have a lot of small scaffolds in my genome file?
As Proteny performs a pairwise analysis between each sequence, it does not handle such a computational explosion very well.
If you remove the smaller scaffolds from your assembly, You should also remove references to these scaffolds in your GFF file.

## What if I want to compare more than 2 genomes?
Proteny currently only supports two genomes at a time.
To compare more genomes, you could perform all pairwise comparissons to other genomes, but then you must integrate those results yourself somehow.
