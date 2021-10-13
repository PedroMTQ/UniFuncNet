

from source.Fetchers.Protein_Fetchers.Protein_Fetcher import *

class Protein_Fetcher_KEGG(Protein_Fetcher):
    def __init__(self,protein_id,extra_args={},memory_storage=None,init_Fetcher=True):
        Protein_Fetcher.__init__( self, protein_id=protein_id,memory_storage=memory_storage)
        self.db='kegg'
        self.set_convergence_args(extra_args)
        if init_Fetcher:
            self.protein=self.get_protein_KEGG()
            self.add_protein()

    def set_convergence_args(self, extra_args):
        # args for convergence
        self.convergence_args['reactions_list'] = []
        self.convergence_args['genes_list'] = []



    def get_genes(self, ec_soup):
        genes=self.textbox_KEGG(ec_soup, 'Genes',to_split=True)
        if not genes: return
        res=[]
        for g in genes:
            organism=g.split(':')[0].lower()
            gene_ids=g.split(':')[-1].split()
            for g_id in gene_ids:
                temp_g_id=g_id.split('(')[0]
                res.append(f'{organism}:{temp_g_id}')
        self.convergence_args['genes_list']=res

    def get_reactions(self, ec_soup):
        reactions=self.textbox_KEGG(ec_soup, 'Reaction(KEGG)',to_split=False)
        if not reactions: return
        self.convergence_args['reactions_list'] = re.findall('R\d+',reactions)

    def get_protein_KEGG(self):

        url = 'https://www.genome.jp/dbget-bin/www_bget?ec:' + self.protein_id
        webpage = self.get_with_fetcher(url)
        if not webpage: return None
        ec_soup = BeautifulSoup(webpage, 'lxml')
        names = self.textbox_KEGG(ec_soup, 'Name')
        if names:
            names=names.split(';')
            names=[i.strip('\n').lower() for i in names]
        orthology = self.textbox_KEGG(ec_soup, 'Orthology',to_split=True)
        if orthology: orthology=[i.split(' : ')[0] for i in orthology]
        external_dbs =[
                       'BRENDA, the Enzyme Database:',
                       'CAS:'
                       ]
        external_dbs_ids={}
        for db in external_dbs:
            ext_db_id = ec_soup.find('td', text=re.compile('.*'+db+'.*'))
            if ext_db_id:
                external_dbs_ids[db]=ext_db_id.findNextSibling().text.split()
            else:   external_dbs_ids[db]=None
        #PG
        self.get_genes(ec_soup)
        #PR
        self.get_reactions(ec_soup)
        res = {'synonyms': names,
             'kegg_ko':orthology,
             'kegg': self.protein_id,
             'cas':external_dbs_ids['CAS:'],
             'brenda':external_dbs_ids['BRENDA, the Enzyme Database:'],
             }
        return Protein(res)

    def converge_protein_global(self):
        self.converge_protein_to_gene()
        self.converge_protein_to_reaction()

    def converge_protein_gpr(self):
        self.converge_protein_to_reaction()

    def converge_protein_rpg(self):
        #PG part of the CRPG pipeline
        self.converge_protein_to_gene()


    def converge_protein_to_gene(self):
        for gene_id in self.convergence_args['genes_list']:
            #so we dont do unnecessary queries for genes from organisms we dont want
            if self.get_wanted_org_kegg_codes():
                org_code=gene_id.split(':')[0]
                if org_code in self.get_wanted_org_kegg_codes():
                    gene_instance= self.find_gene(query_id=gene_id)
                    if gene_instance:
                        self.get_protein().set_detail('gene_instances',gene_instance,converged_in=self.db)
            else:
                gene_instance= self.find_gene(query_id=gene_id)
                if gene_instance:
                    self.get_protein().set_detail('gene_instances',gene_instance,converged_in=self.db)


    def converge_protein_to_reaction(self):
        for reaction_id in self.convergence_args['reactions_list']:
            reaction_instance= self.find_reaction(query_id=reaction_id)

            if reaction_instance:
                self.get_protein().set_detail('reaction_instances',reaction_instance,converged_in=self.db)


if __name__ == '__main__':
    from source.Biological_Components.Protein import Protein
    from source.Fetchers.Reaction_Fetchers.Reaction_Fetcher import *
    import re
    #p0=Protein_Fetcher_KEGG('2.3.1.15')
    p0=Protein_Fetcher_KEGG('2.3.1.15')
    #p0.converge_protein_gpr()
    #p1=Protein_Fetcher_KEGG('1.14.14.73').get_protein()
