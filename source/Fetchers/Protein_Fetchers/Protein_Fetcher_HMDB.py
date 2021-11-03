
from source.Fetchers.Protein_Fetchers.Protein_Fetcher import *

class Protein_Fetcher_HMDB(Protein_Fetcher):
    def __init__(self,protein_id,extra_args={},memory_storage=None,init_Fetcher=True):
        Protein_Fetcher.__init__( self, protein_id=protein_id,memory_storage=memory_storage)
        self.db='hmdb'
        self.set_convergence_args(extra_args)
        if init_Fetcher:
            self.protein = self.get_protein_HMDB()
            self.add_protein()

    def set_convergence_args(self,extra_args):
        #args for getting protein
        if 'protein_soup' in extra_args: self.convergence_args['protein_soup']=extra_args['protein_soup']
        else: self.convergence_args['protein_soup']=[]
        #args for convergence
        #PR
        self.convergence_args['reactions_list']=[]

    def get_protein_HMDB(self):

        protein_soup=self.convergence_args['protein_soup']
        if not protein_soup:
            url = 'http://www.hmdb.ca/proteins/' + self.protein_id
            webpage = self.get_with_fetcher(url)
            if not webpage: return None
            protein_soup = BeautifulSoup(webpage, 'lxml')
            self.convergence_args['protein_soup']=protein_soup
        enz_syns = [protein_soup.find('th', text=re.compile('Name')).findNext().text]
        wanted_enz_details=[#'General Function',
                            # 'Specific Function',
                             'UniProtKB/Swiss-Prot ID',
                             'UniProtKB/Swiss-Prot Entry Name',
                             'PDB IDs',
                             'Pathways'
                 ]
        details={}
        for i in wanted_enz_details:
            wanted_detail= protein_soup.find('th', text=re.compile(i)).findNext().text
            if 'Not Available' not in wanted_detail:
                fixed_detail=wanted_detail.strip()
                fixed_detail=fixed_detail.replace('\n',',')
                if i!='Pathways':
                    fixed_detail=fixed_detail.replace('"','')
                    fixed_detail=fixed_detail.replace(']','')
                    fixed_detail=fixed_detail.replace('[','')
                fixed_detail=fixed_detail.split(',')
                fixed_detail=[j for j in fixed_detail if j]
                fixed_detail=[j.strip() for j in fixed_detail]
                details[i]=fixed_detail
            else:
                details[i]=None


        protein_instance = Protein(
                                {'synonyms':enz_syns,
                                   'hmdb':self.protein_id,
                                   'uniprot':details['UniProtKB/Swiss-Prot ID'],
                                   'uniprot_name':details['UniProtKB/Swiss-Prot Entry Name'],
                                   'pdb':details['PDB IDs'],
                                 })
        self.get_reactions_from_soup(protein_soup)
        return protein_instance

    def get_reactions_from_soup(self,soup):
        reactions_ids=soup.find_all(href=re.compile('/reactions/\d+'))
        self.convergence_args['reactions_list'] = [re.search('\d+', str(i)).group() for i in reactions_ids]

    def converge_protein_global(self):
        self.converge_protein_to_gene()
        self.converge_protein_to_reaction()

    def converge_protein_gpr(self):
        self.converge_protein_to_reaction()

    def converge_protein_rpg(self):
        #PG part of the CRPG pipeline
        self.converge_protein_to_gene()

    def converge_protein_to_gene(self):
        gene_instance= self.find_gene(query_id=self.protein_id,
                                     extra_args={'gene_soup':self.convergence_args['protein_soup']},
                                      )
        if gene_instance:
            self.get_protein().set_detail('gene_instances',gene_instance,converged_in='hmdb')


    def converge_protein_to_reaction(self):
        #here all we do is fetch each reaction in the reactions_list
        for reaction_id in self.convergence_args['reactions_list']:
            reaction_instance= self.find_reaction(query_id=reaction_id,
                                                  extra_args={},
                                                  )
            if reaction_instance:
                self.get_protein().set_detail('reaction_instances',reaction_instance,converged_in='hmdb')


if __name__ == '__main__':
    from source.Biological_Components.Protein import Protein
    from source.Fetchers.Reaction_Fetchers.Reaction_Fetcher import *
    import re
    protein=Protein_Fetcher_HMDB('HMDBP00609').get_protein()

