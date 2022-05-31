# Running the UniFuncNet_Neo4j_Connector workflow

This workflow can be used to export UniFuncNet's data to neo4j.
This database can then be queried by the user using the browser, python, etc.

I've also included an option to extract SIF-formatted functional networks from this database. To do so, first export the
UniFuncNet to Neo4j. Then run the workflow again, but now use the `consensus_annotation.tsv` coming
from [Mantis](https://github.com/PedroMTQ/mantis).

Keep in mind each time you run the UniFuncNet->Neo4j part of this workflow, the database is reset!
Running this workflow to extract the SIF networks, does not reset the database.

## Required installation

This workflow requires the [Neo4j graph database](https://neo4j.com/product/neo4j-graph-database/), either
the [enterprise](https://neo4j.com/download-center/#enterprise)
or [community](https://neo4j.com/download-center/#community) edition.
The community edition is more limited, however, it will fit within the scope of this workflow, i.e., exporting
UniFuncNet's data to Neo4j.

Since only one database can be created per Neo4j community edition instance, this workflow always deletes previous
databases in the Neo4j database instance.

To run Neo4j in UniFuncNet do the following:

1. Go to the url provided above for enterprise or community editions
2. Download your specific OS version and follow the instructions
3. Go to your Neo4j folder and run  `./bin/neo4j console`. Now your Neo4j database should be running in the background.
4. Install the Neo4j driver in your environment with `conda install -c conda-forge neo4j-python-driver`
5. Go to http://localhost:7474/browser/ and set up your database (basically, set a username and password). The initial
   username and password are `neo4j`

## Executing workflow

6. run the workflow
   with `unifuncnet neo4j -i unifuncnet_input.tsv -o output_folder -db kegg,rhea,biocyc,hmdb -username <neo4j_username> --password <neo4j_password>`

- `input`: input folder if creating the Neo4j database or consensus_annotation.tsv if exporting network
- `output_folder`: output directory
- `database_name`: database name, `neo4j` by default
- `username`: database username, `neo4j` by default
- `password`: database password, `unifuncnet_neo4j` by default
- `port`: port used by No4j, 7687 by default