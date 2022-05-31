# What is this workflow for?

This workflow is used to create an universal UniFuncNet input tsv file.

It collects all available KOs,enzyme ECs, Metacyc protein ids and Rhea reactions, and saves a UniFuncNet input tsv with
those IDs with the 'global' search mode.

You need the Rhea and metacyc SQLite dbs
Also, place the KEGG ECs json [ko01000.json](https://www.kegg.jp/brite/ko01000) and KEGG KOs
json [ko00001.json](https://www.genome.jp/kegg-bin/show_brite?ko00001.keg) in this folder.

## Executing workflow

- `o`: output folder
- `ec`: EC json path, downloaded from https://www.kegg.jp/brite/ko01000
- `ko`: KO json path, downloaded from https://www.genome.jp/kegg-bin/show_brite?ko00001.keg'
- `db`: database to generate input, `all` (default) or `ko`
