# DRAX


DRAX is a biological data collector, fetching information on gene-protein-reaction-compound relationship across multiple databases.
It was built to aid in the expansion of metabolic networks, automating the often monotonous, manual data collection.

## Installation 


1. `git clone git@github.com:PedroMTQ/DRAX.git`  
2. Go to the cloned DRAX folder and run `conda env create -f drax_env.yml`
3. Run `conda activate drax_env`
4. Download the browser driver for your current OS and current browser version (e.g., [geckodriver](https://github.com/mozilla/geckodriver/releases) for Mozilla Firefox)
5. Move downloaded browser driver to `DRAX/External_Tools/Browser_Drivers/`

#### Browser drivers

DRAX supports two browsers, Mozilla Firefox and Google Chrome, please download the driver for your current browser version and add it to `DRAX/Browser_Drivers`:
- Mozilla Firefox please go to https://github.com/mozilla/geckodriver/releases
- Google Chrome please go to https://chromedriver.chromium.org/downloads 


## Using DRAX


Running DRAX within a personal computer should be straightforward, however keep in mind that in order to scrape websites with Javascript, DRAX (more specifically the package `Selenium`) needs to use a browser. To install browsers in a server you may need admin access. 

You can run the code below to test the execution:

    python DRAX --example

A typical run would look like:

    python DRAX pr protein_search -i input.tsv 




Data collection is possible for multiple databases:
- KEGG
- HMDB
- Biocyc

To avoid overloading these database websites, a 10 seconds pause between requests was added.

DRAX accepts the following parameters:


    python DRAX search_direction search_type -i input_path -o output_folder -rm -db db1,db2

    Mandatory arguments:  search_direction:
                                accepted: 'global', 'rpg', 'pg', 'gpr', 'pr','gp','na'
                          search_type
                                accepted: 'gene_search', 'protein_search', 'reaction_search', 'compound_search'
                          --input_path / -i
    Optional arguments:  --output_folder / -o
                         --reaction_metabolites / -rm
                         --databases / -db

Where each parameter corresponds to the following:

- `search_direction`  - direction of the search, `global` for searching in both directions, `na` for searching only the provided IDs, without connecting to other instances, `rpg` for searching from reaction->protein->gene, and `pg`, `gpr`, `pr`,`gp`
- `search_type` - starting point of the search, if the user aims to provide IDs for proteins, then it would be `protein_search`, and the same for the other types of biological instances
- `input_path` - the input tsv file path. 
- `output_folder` - the output folder where the spreadsheets are stored
- `reaction_metabolites` - searches for data on each reaction's compound, which improves reaction matching across different databases
- `databases ` - databases that DRAX can search in, by default HMDB,Biocyc, and KEGG
- `politeness_timer` - time (seconds) between requests. Default is 10. Please be careful not to overleaf the corresponding databases, you might get blocked from doing future requests.


All search types except `compound_search` provide network-based scraping, we have decided to keep `compound_search` as an option since the user may also use it to search for more information on compounds (e.g., start with a synonym or KEGG ID and find more information on these)




Data is  retrieved according to the information provided, for example, if the user provides the KEGG gene ID hsa:150763, then, given that the gpr_convergence is used, DRAX would fetch information on this gene, the KEGG protein entries connected to this gene (i.e., 2.3.1.15) and by extent the reactions these protein(s) catalyze (i.e., R00851,R02617,R09380).




### Formatting input file

The input file should be a tab separated file that looks something like this:
 
  |      |      |
  | ---  | ---  |
  | enzyme_ec  | 2.7.1.1  |
  | kegg_ko    | K00844   |


The first column should correspond to an ID type, where following nomenclature is accepted:
- biocyc
- kegg
- hmdb
- enzyme_ec
- uniprot
- inchi_key
- chemspider
- synonym

The second column would then correspond to the respective ID that the user would like to gather information for.

### Supported input IDs

Several IDs are allowed per biological instance:

- Reaction:
    - Biocyc (e.g., "RXN66-521")
    - KEGG (e.g., "R02848")
    - HMDB (e.g., "14073")
    
- Protein:
    - enzyme EC number (e.g., "2.7.1.1")
    - KEGG (e.g., "2.7.1.1")
    - KEGG KO (e.g., "K00844")
    - Biocyc (e.g., "2.7.1.1")
    - Uniprot (e.g., "P10632")
    - HMDB (e.g., "HMDBP00609")
    

- Gene:
    - KEGG (e.g., "hsa:150763")
    - Uniprot (e.g., "P10632")
    - Biocyc (e.g., "HS08548")
    - HMDB (e.g., "HMDBP00087")

- Compound:
    - Biocyc (e.g., "WATER")
    - KEGG (e.g., "C00001")
    - HMDB (e.g., "HMDB02111")
    - InChI key (e.g., "XLYOFNOQVPJJNP-UHFFFAOYSA-N")
    - Chemspider (e.g., "937")
    - Synonym  (e.g., "water")

### Output

DRAX outputs 4 tsv files: `Compounds.tsv`,`Reactions.tsv`,`Proteins.tsv`,`Genes.tsv`
Each file contains multiple instances (e.g., compound) with a tab-separated list of identifiers or other relevant information.
Specifically, all instances contain an `internal_id` which can be used for graph-based approaches cross-linking (e.g., `manuscript.py`), and often a list of identifiers and synonyms.
In the case of reactions, proteins and genes, cross-linking is available in the form of `<instance>_connected`. For example, if the user searches for all reactions of a set of proteins, then the retrieved proteins would have a list of `reactions_connected:<reaction internal_id>` depicting which reactions this protein is connected to. The same would apply for other search directions or search starting points.
Reactions also contain the list of compounds involved in this reaction, e.g.: `reaction_compounds:<C1> + <C2> <=> <C3> + <C4> + <C5>` where `CX` corresponds to a compound's `internal_id`.

Using the example above as an example (with input the enzyme EC 2.7.8.26), the output for each instance would look somewhat like:
- `Proteins.tsv`
  - internal_id:**270**	cas:DA85_04120	enzyme_ec:2.7.8.26	kegg_ko:K02233	reactions_connected:**25550**	genes_connected:**9035**	synonyms:adenosylcobinamide-gdp ribazoletransferase	synonyms:cobalamin (5'-phosphate) synthase	synonyms:cobalamin (5′-phosphate) synthase	synonyms:cobalamin synthase	synonyms:cobalamin-5'-phosphate synthase	synonyms:cobalamin-5′-phosphate synthase	synonyms:cobs	synonyms:α-ribazole ribazoletransferase
- `Reactions.tsv`
  - internal_id:**25550**	biocyc:RXN-19297	pathways:4-methylphenyl adenosylcobamide biosynthesis from adenosylcobinamide-GDP	reaction_compounds:**10310** + 6731 => 21252 + 24415 + 8385	reaction_str:adenosylcobinamide-GDP + 4-methylphenyl ribotide phosphate => 4-methylphenyl-Coβ-adenosylcobamide 5'-phosphate + GMP + H+
- `Compounds.tsv`
  - internal_id:**10310**	bigg:agdpcbi	biocyc:ADENOSYLCOBINAMIDE-GDP	chebi:60487	hmdb:HMDB12185	kegg:C06510	pubchem_cid:135398566	seed:cpd03920	synonyms:adenosine-gdp-cobinamide	synonyms:adenosylcobinamide-gdp

As can be seen, the protein (i.e., `internal_id:270`) shown above is connected to the reaction `25550` which in turn is described as the following interaction between compounds: 10310 + 6731 => 21252 + 24415 + 8385. These compounds are then listed in the `Compounds.tsv` as shown above. For visualization purposess only a small transcript is shown above.


### On search directions

The search direction  `na` removes any type of extra search besides the initial IDs provided, meaning that if the user provides the KEGG gene ID hsa:150763, we would still search for it (and retrieve gene-related data) but we would not search for its respective proteins.

The search direction `global` is an extensive search, where we can search in both directions (i.e. g->p->r and r->p->g), meaning that if the user provides the KEGG gene ID hsa:150763 we would retrieve the respective proteins 2.3.1.15, but we would then also search for all the gene IDs for the protein 2.3.1.15 (e.g., hsa:150763, ptr:107971389, csab:103215676). The same would apply to the protein and reaction search.

Compound search is limited to retrieve only compound-related information. Should the user provide a compound name (e.g. "water") the compound search may also retrieve related compounds (since DRAX uses the each website's search bar to retrieve the most likely compound entry). However, if an ID is provided, DRAX will first search for the ID, if information is not found, then the synonyms are used as a search method.
This also applies to when the option `reaction_metabolites` is enabled and the reaction does not contain any compound ID, in that case the reaction string (e.g. "sn-Glycerol 3-phosphate + Acyl-CoA <=> 1-Acyl-sn-glycerol 3-phosphate + CoA") is parsed  and its compounds are searched using the method previously described.


# License and copyright

This project is available under the [MIT license](https://github.com/PedroMTQ/DRAX/LICENSE).

# References and acknowledgements

> Minoru Kanehisa, Susumu Goto, KEGG: Kyoto Encyclopedia of Genes and Genomes, Nucleic Acids Research, Volume 28, Issue 1, 1 January 2000, Pages 27–30, https://doi.org/10.1093/nar/28.1.27
> 
> Caspi R, Billington R, Keseler IM, Kothari A, Krummenacker M, Midford PE, Ong WK, Paley S, Subhraveti P, Karp PD. The MetaCyc database of metabolic pathways and enzymes - a 2019 update. Nucleic Acids Res. 2020 Jan 8;48(D1):D445-D453. doi: 10.1093/nar/gkz862. PMID: 31586394; PMCID: PMC6943030.
>
> Wishart DS, Feunang YD, Marcu A, Guo AC, Liang K, Vázquez-Fresno R, Sajed T, Johnson D, Li C, Karu N, Sayeeda Z, Lo E, Assempour N, Berjanskii M, Singhal S, Arndt D, Liang Y, Badran H, Grant J, Serra-Cayuela A, Liu Y, Mandal R, Neveu V, Pon A, Knox C, Wilson M, Manach C, Scalbert A. HMDB 4.0: the human metabolome database for 2018. Nucleic Acids Res. 2018 Jan 4;46(D1):D608-D617. doi: 10.1093/nar/gkx1089. PMID: 29140435; PMCID: PMC5753273.