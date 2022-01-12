This is the `Resources/metacyc` folder.

This is where you should extract all the metacyc database files.

You can request your academic license here:
https://biocyc.org/download.shtml

When you obtain your license, download the `Metacyc flat files`. The file is usually named `meta.tar.gz`.

Extract all the files in the folder `data` to `DRAX/Resources/metacyc/`.
Your need the following files:
- compounds.dat
- proteins.dat
- reactions.dat
- gene-links.dat
- genes.dat

**When you execute DRAX, all other files in this folder are removed!**

So basically, this folder will only contain the previous files plus this `README.md`.

**If you don't have these files, DRAX will not integrate Metacyc**

If everything is working correctly, a `metacyc.db` file should be created in the `Resources` folder


If using Metacyc date please cite:

Ron Caspi, Richard Billington, Ingrid M Keseler, Anamika Kothari, Markus Krummenacker,
Peter E Midford, Wai Kit Ong, Suzanne Paley, Pallavi Subhraveti, Peter D Karp
The MetaCyc database of metabolic pathways and enzymes-a 2019 update
Nucleic Acids Research, Volume 48, Issue D1, 08 January 2020, Pages D445–D453,
https://doi.org/10.1093/nar/gkz862
