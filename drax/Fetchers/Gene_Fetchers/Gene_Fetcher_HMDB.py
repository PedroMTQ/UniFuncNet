
from drax.Fetchers.Gene_Fetchers.Gene_Fetcher import *


#hmdb doesnt separate genes and enzymes so we just use the same id for both. If we come form the enzyme fetcher we can just reuse the enzyme soup
#otherwise we fetch a new soup
class Gene_Fetcher_HMDB(Gene_Fetcher):
    def __init__(self,gene_id=None,extra_args={},memory_storage=None):
        #hmdb is a database for human data
        Gene_Fetcher.__init__( self, gene_id=gene_id,extra_args=extra_args,memory_storage=memory_storage)
        self.db='hmdb'
        self.set_convergence_args(extra_args)
        self.gene=self.get_gene_hmdb()
        self.add_gene()

    def set_convergence_args(self,extra_args):
        #in hmdb, genes and proteins are in the same page, so effectively we are just parsing the same page
        #from PG we would parse the webpage coming from the protein fethcher, in GP we just send the wepbage to the protein fetcher instead
        #args for getting gene
        self.convergence_args['protein_soup']=None
        #args for convergence
        if 'gene_soup' in extra_args: self.convergence_args['gene_soup']=extra_args['gene_soup']
        else: self.convergence_args['gene_soup']=None


    def get_gene_hmdb(self):

        gene_soup=self.convergence_args['gene_soup']
        if not gene_soup:
            url = f'http://www.hmdb.ca/proteins/{self.gene_id}'
            webpage = self.get_with_fetcher(url)
            if not webpage: return None
            gene_soup = BeautifulSoup(webpage, 'lxml')
            self.convergence_args['gene_soup']=gene_soup
        wanted_gene_details=['Gene Name',
                 'Chromosome Location',
                 'Locus',
                 'Gene Sequence',
                 'Protein Sequence',
                 'UniProtKB/Swiss-Prot ID',
                 'UniProtKB/Swiss-Prot Entry Name',
                 'GenBank Gene ID',
                 'GeneCard ID',
                 'GenAtlas ID',
                 'HGNC ID'
                 ]
        details={}
        for i in wanted_gene_details:
            wanted_detail= gene_soup.find('th', text=re.compile(i))
            if wanted_detail:
                wanted_detail=wanted_detail.findNext().text
                if 'Not Available' not in wanted_detail:
                    details[i]=wanted_detail.strip()
                else:
                    details[i]=None
            else:
                details[i]=None

        gene_instance = Gene({'synonyms':details['Gene Name'],
                              'ncbi_genbank_gene':details['GenBank Gene ID'],
                              'genecard':details['GeneCard ID'],
                              'genatlas':details['GenAtlas ID'],
                              'hgnc':details['HGNC ID'],
                              'uniprot':details['UniProtKB/Swiss-Prot ID'],
                              'uniprot_name':details['UniProtKB/Swiss-Prot Entry Name']
                              })
        return gene_instance



    def converge_gene_gpr(self):
        self.converge_gene_to_protein()

    def converge_gene_to_protein(self):
        #this function will merely send the gene soup to the protein fetcher
        print(f'Linking from gene {self.gene_id} in {self.db} to protein {self.gene_id}')

        protein_instance = self.find_protein(query_id=self.gene_id,
                                             extra_args={'protein_soup': self.convergence_args['gene_soup']},
                                             )
        if protein_instance:
            self.get_gene().set_detail('protein_instances',protein_instance,converged_in='hmdb')


if __name__ == '__main__':
    from drax.Biological_Components.Protein import Protein
    from drax.Fetchers.Reaction_Fetchers.Reaction_Fetcher import *
    gene=Gene_Fetcher_HMDB('HMDBP00087').get_gene()


