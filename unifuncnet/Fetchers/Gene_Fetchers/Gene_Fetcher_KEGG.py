
from unifuncnet.Fetchers.Gene_Fetchers.Gene_Fetcher import *

class Gene_Fetcher_KEGG(Gene_Fetcher):
    def __init__(self,gene_id,extra_args={},memory_storage=None):
        Gene_Fetcher.__init__( self, gene_id=gene_id,extra_args=extra_args,memory_storage=memory_storage)
        self.db='kegg'
        self.set_convergence_args(extra_args)
        self.gene=self.get_gene_kegg()
        self.add_gene()

    def set_convergence_args(self,extra_args):
        #args for convergence
        self.convergence_args['proteins_list']=[]

    def get_gene_kegg(self):
        #match wasnt found so we keep going
        url = f'https://www.kegg.jp/dbget-bin/www_bget?{self.gene_id}'
        webpage = self.get_with_fetcher(url)
        if not webpage: return None
        gene_soup = BeautifulSoup(webpage, 'lxml')
        gene_name = self.textbox_KEGG(gene_soup, 'Gene name')
        if gene_name:
            gene_name=gene_name.split(',')
            gene_name=[i.strip() for i in gene_name]
        ec_pattern = re.compile('\[EC:(\d+\.){2,3}((\d+)|-)(\s(\d+\.){2,3}((\d+)|-))*\]')
        ko_and_ec = self.textbox_KEGG(gene_soup, 'KO')
        if ko_and_ec:
            ecs= re.search(ec_pattern,ko_and_ec)
            if ecs:
                ecs=ecs.group()[4:-1].split(' ')
                ecs=[i for i in ecs if not i.endswith('-')]
                self.convergence_args['proteins_list'].extend(ecs)
            ko = ko_and_ec.split()[0].strip()
        else: ko=None
        other_dbs = self.textbox_KEGG(gene_soup, 'Other DBs', to_split=True)
        ncbi_genbank_gene,ncbi_genbank_protein,hgnc,uniprot_id,ensembl,vgnc=None,None,None,None,None,None
        if other_dbs:
            for i in other_dbs:
                if 'ncbi-geneid' in i.lower():
                    ncbi_genbank_gene = i.split(':')[-1].strip()
                    ncbi_genbank_gene= ncbi_genbank_gene.split(' ')
                if 'ncbi-proteinid' in i.lower():
                    ncbi_genbank_protein = i.split(':')[-1].strip()
                    ncbi_genbank_protein= ncbi_genbank_protein.split(' ')
                if 'hgnc' in i.lower():
                    hgnc= i.split(':')[-1].strip()
                    hgnc= hgnc.split(' ')
                if 'ensembl' in i.lower():
                    ensembl= i.split(':')[-1].strip()
                    ensembl= ensembl.split(' ')
                if 'vgnc' in i.lower():
                    vgnc= i.split(':')[-1].strip()
                    vgnc= vgnc.split(' ')
                if 'uniprot' in i.lower():
                    uniprot_id = i.split(':')[-1].strip()
                    uniprot_id= uniprot_id.split(' ')
        res = {'synonyms':gene_name,
                              'kegg_ko':ko,
                              'kegg':self.gene_id,
                              'hgnc':hgnc,
                              'ensembl':ensembl,
                              #'vgnc':vgnc,
                              'ncbi_genbank_gene':ncbi_genbank_gene,
                              'ncbi_genbank_protein':ncbi_genbank_protein,
                              'uniprot': uniprot_id}
        return Gene(res)

    def converge_gene_gpr(self):
        self.converge_gene_to_protein()

    def converge_gene_to_protein(self):
        if self.convergence_args['proteins_list']:
            #this function will merely send the gene soup to the protein fetcher
            #the same gene can code different proteins
            for enzyme_ec in self.convergence_args['proteins_list']:
                print(f'Linking from gene {self.gene_id} in {self.db} to protein {enzyme_ec}')

                protein_instance = self.find_protein(query_id=enzyme_ec,
                                                     extra_args={},
                                                     )
                if protein_instance:
                    self.get_gene().set_detail('protein_instances',protein_instance,converged_in=self.db)



if __name__ == '__main__':
    from unifuncnet.Biological_Components.Protein import Protein
    from unifuncnet.Fetchers.Reaction_Fetchers.Reaction_Fetcher import *
    gene=Gene_Fetcher_KEGG('hsa:150763')
    print(gene.get_gene())