
from unifuncnet.Fetchers.Protein_Fetchers.Protein_Fetcher import *

class Protein_Fetcher_HMDB(Protein_Fetcher):
    def __init__(self,protein_id,extra_args={},memory_storage=None):
        Protein_Fetcher.__init__( self, protein_id=protein_id,memory_storage=memory_storage)
        self.db='hmdb'
        self.set_convergence_args(extra_args)
        self.protein = self.get_protein_HMDB()
        self.add_protein()

    def set_convergence_args(self,extra_args):
        #args for getting protein
        if 'protein_soup' in extra_args: self.convergence_args['protein_soup']=extra_args['protein_soup']
        else: self.convergence_args['protein_soup']=[]
        if 'cpd_id' in extra_args: self.convergence_args['cpd_id']=extra_args['cpd_id']
        else: self.convergence_args['cpd_id']=None
        #args for convergence
        #PR
        self.convergence_args['reactions_list']=[]

    def get_protein_HMDB(self):

        protein_soup=self.convergence_args['protein_soup']
        if not protein_soup:
            url = f'http://www.hmdb.ca/proteins/{self.protein_id}'
            webpage = self.get_with_fetcher(url)
            if not webpage: return None
            protein_soup = BeautifulSoup(webpage, 'lxml')
            self.convergence_args['protein_soup']=protein_soup
        enz_syns = protein_soup.find('th', text=re.compile('Name'))
        if enz_syns:
            enz_syns=[enz_syns.findNext().text]
        wanted_enz_details=[#'General Function',
                            # 'Specific Function',
                             'UniProtKB/Swiss-Prot ID',
                             'UniProtKB/Swiss-Prot Entry Name',
                             'PDB IDs',
                             'Pathways'
                 ]
        details={}
        for i in wanted_enz_details:
            try:
                wanted_detail= protein_soup.find('th', text=re.compile(i)).findNext().text
                if 'Not Available' in wanted_detail: wanted_detail=None
            except:
                wanted_detail=None
            if wanted_detail:
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

    def get_reactions_for_cpd(self,reactions_dict):
        res=set()
        for reaction_id in reactions_dict:
            reaction_str=reactions_dict[reaction_id].replace('+',' + ').replace('=',' = ')
            reaction_str = fix_html_sign(reaction_str)
            for syn in self.convergence_args['main_cpd'].get_detail('synonyms',all_possible=True):
                if syn in reaction_str:
                    res.add(reaction_id)
        return res

    def get_reactions_from_soup(self,soup):
        reactions_ids=soup.find_all(href=re.compile('/reactions/\d+'))
        self.convergence_args['reactions_list'] = [re.search('\d+', str(i)).group() for i in reactions_ids]


    def get_reactions_from_soup_legacy(self,soup):
        reactions_dict={}
        reactions_box= soup.find('th',text=re.compile('Reactions'))
        if reactions_box:
            reactions_box=reactions_box.findNext()
            for elements in reactions_box:
                elements = str(elements).strip('\n').split('\n')
                for element in elements:
                    if element:
                        if element.startswith('<td>') and element.endswith('</td>'):
                            reaction_str=element.replace('</td>','').replace('<td>','')
                        if element.startswith('<a class="btn '):
                            reaction_id=re.search('/reactions/\d+',element)
                            if reaction_id:
                                reaction_id=reaction_id.group().split('/')[-1]
                                reactions_dict[reaction_id]=reaction_str
        #when coming from compound->protein
        if self.convergence_args['main_cpd']:
            self.convergence_args['reactions_list'] = self.get_reactions_for_cpd(reactions_dict)
        else:
            self.convergence_args['reactions_list'] = list(reactions_dict.keys())

    def converge_protein_global(self):
        self.converge_protein_to_gene()
        self.converge_protein_to_reaction()

    def converge_protein_gpr(self):
        self.converge_protein_to_reaction()

    def converge_protein_rpg(self):
        #PG part of the CRPG pipeline
        self.converge_protein_to_gene()

    def converge_protein_to_gene(self):
        print(f'Linking from protein {self.protein_id} in {self.db} to gene {self.protein_id}')

        gene_instance= self.find_gene(query_id=self.protein_id,
                                     extra_args={'gene_soup':self.convergence_args['protein_soup']},
                                      )
        if gene_instance:
            self.get_protein().set_detail('gene_instances',gene_instance,converged_in='hmdb')


    def converge_protein_to_reaction(self):
        #here all we do is fetch each reaction in the reactions_list
        for reaction_id in self.convergence_args['reactions_list']:
            print(f'Linking from protein {self.protein_id} in {self.db} to reaction {reaction_id}')

            reaction_instance= self.find_reaction(query_id=reaction_id,
                                                  extra_args={'cpd_id':self.convergence_args['cpd_id']},
                                                  )
            if reaction_instance:
                self.get_protein().set_detail('reaction_instances',reaction_instance,converged_in='hmdb')


if __name__ == '__main__':
    from unifuncnet.Biological_Components.Protein import Protein
    from unifuncnet.Fetchers.Reaction_Fetchers.Reaction_Fetcher import *
    import re
    cpd=Compound({
        'metacyc':'SER',
        'kegg':'C00065',
        'hmdb':'HMDB0000187',
        'chemspider':'5736',
        'inchi_key':{'MTCFGRXMJLQNBG-REOHCLBHSA-N','MTCFGRXMJLQNBG-REOHCLBHSA-N'},
        'bigg':'33717',
        'synonyms':'l-serine',
        'chebi':'17115',
        'drugbank':'DB00133',
    })
    protein=Protein_Fetcher_HMDB('HMDBP00629',extra_args={'cpd_id':'HMDB0000187'})
