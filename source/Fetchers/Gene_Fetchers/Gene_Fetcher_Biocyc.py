
from source.Fetchers.Gene_Fetchers.Gene_Fetcher import *

class Gene_Fetcher_Biocyc(Gene_Fetcher):
    def __init__(self,gene_id=None,extra_args={},memory_storage=None,init_Fetcher=True):
        Gene_Fetcher.__init__( self, gene_id=gene_id,extra_args=extra_args,memory_storage=memory_storage)
        self.db='biocyc'
        self.set_convergence_args(extra_args)
        if init_Fetcher:
            self.gene=self.get_gene_biocyc()
            self.tidy_up_convergence_args()
            self.add_gene()

    #PG
    def get_gene_biocyc(self):
        gene_url='https://biocyc.org/gene?orgid=META&id='+self.gene_id+'#tab=showAll'
        webpage = self.get_with_fetcher(gene_url, selenium=True)
        if not webpage: return None
        soup = BeautifulSoup(webpage, 'lxml')
        headers = soup.find_all('font', class_='header')
        gene_info = {'biocyc':None}
        for h in headers:
            previous = h.parent
            if 'gene' in previous.text:
                gene_info['synonyms'] = h.text
            elif 'enzyme' in previous.text:
                enzyme_name = h.text
                self.convergence_args['protein_name'] = enzyme_name
        self.get_ids_soup(soup,gene_info)
        if gene_info['biocyc']:
            self.gene_id = gene_info['biocyc'][0]
        if not gene_info['biocyc']: gene_info['biocyc']=self.gene_id
        self.get_reactions_convergence(soup)
        return Gene(gene_info)

    def get_ids_soup(self,soup,gene_info):
        corresponding_dbs = {
            'Entrez': 'ncbi_genbank_protein',
            'Entrez-Nucleotide': 'ncbi_genbank_nucleotide',
            'Entrez-gene': 'ncbi_genbank_gene',
            'GeneCards': 'genecard',
            'UniProt': 'uniprot',
            'MetaCyc': 'biocyc',
            'CAZy In-Family': 'cazy',
            'InterPro In-Family': 'interpro',
            'Panther In-Family': 'panther',
            'Pfam In-Family': 'pfam',
        }
        links = soup.find('table', class_='linkBox')
        if links:
            links = links.find_all('tr')
            for l in links:
                search = l.find_all('td')
                if len(search) > 1:
                    db, db_id = search
                    if db.text in corresponding_dbs:
                        corresponding_db_info= db_id.text.split(',')
                        corresponding_db_info=[i.strip() for i in corresponding_db_info]
                        gene_info[corresponding_dbs[db.text]]=corresponding_db_info
        accession_ids = soup.find(text=re.compile('Accession IDs'))
        if accession_ids:
            db_ids = accession_ids.findNext('td')
            temp=[]
            for ele in db_ids:
                if ele.name =='br':
                    if len(temp)<2: temp=temp[0].split(' ')
                    if len(temp)<2: temp=temp[0].split(' ')
                    for j in temp:
                        if j in corresponding_dbs:
                            temp.remove(j)
                            gene_info[corresponding_dbs[j]]=temp
                    temp=[]
                else:
                    ele = strip_tags(str(ele)).replace('\xa0', ' ').replace('\n', '').replace('(', '').replace(')', '').strip()
                    temp.append(ele)
            if len(temp) < 2: temp = temp[0].split(' ')
            for j in temp:
                if j in corresponding_dbs:
                    temp.remove(j)
                    gene_info[j] = temp
        #adding family ids to placeholder protein
        if 'protein_gene_instance' in self.convergence_args:
            if 'pfam' in gene_info:
                self.convergence_args['protein_gene_instance'].set_detail('pfam',gene_info['pfam'])
                gene_info.pop('pfam')
            if 'uniprot' in gene_info:
                self.convergence_args['protein_gene_instance'].set_detail('uniprot',gene_info['uniprot'])
            if 'panther' in gene_info:
                self.convergence_args['protein_gene_instance'].set_detail('panther',gene_info['panther'])
                gene_info.pop('panther')
            if 'cazy' in gene_info:
                self.convergence_args['protein_gene_instance'].set_detail('cazy',gene_info['cazy'])
                gene_info.pop('cazy')
            if 'interpro' in gene_info:
                self.convergence_args['protein_gene_instance'].set_detail('interpro',gene_info['interpro'])
                gene_info.pop('interpro')


    def get_reactions_convergence(self,soup):
        reaction_boxes_soup=soup.find_all(text=re.compile('Enzymatic activity:.*'))
        if reaction_boxes_soup:
            c=0
            for box_soup in reaction_boxes_soup:
                box_soup=box_soup.parent.findNext('table',class_='infoBox2')
                protein_id= box_soup.find('a',class_='EC-NUMBER')
                reaction_soup = box_soup.find('a', class_='REACTION')
                if reaction_soup:
                    r_link = reaction_soup['href']
                    if protein_id:
                        protein_id=protein_id.text
                    else:
                        protein_id=find_ecs(str(box_soup))
                        if protein_id: protein_id=protein_id[0]
                        else: protein_id=c
                    if protein_id not in self.convergence_args['reactions']:    self.convergence_args['reactions'][protein_id]= []
                    reaction_id = r_link[re.search('object=', r_link).span()[1]:]
                    #removing scripts from the text
                    unwanted= reaction_soup.find_all('script')
                    for u in unwanted: u.extract()
                    reaction_str=reaction_soup.text
                    if not reaction_str:
                        reaction_str = strip_tags(reaction_soup['onmouseover'].split('<br>')[1])

                    #sometimes these boxes dont have the reaction str so we ignore them
                    if reaction_str:
                        reaction_str=fix_html_sign(reaction_str)
                        self.convergence_args['reactions'][protein_id].append([reaction_id, reaction_str])
                        c+=1


    def set_convergence_args(self,extra_args):
        self.convergence_args['reactions_complete_ec'] = []
        self.convergence_args['reactions_incomplete_ec'] = []
        self.convergence_args['reactions_not_ec'] = []
        self.convergence_args['reactions_list'] = []
        self.convergence_args['protein_name'] = None
        #args for convergence
        self.convergence_args['reactions']={}
        #args for getting gene and to give to protein

    def tidy_up_convergence_args(self):
        if self.convergence_args['reactions']:
            for p_id in self.convergence_args['reactions']:
                #gene to reaction
                if isinstance(p_id, int):
                    for r_id_str in self.convergence_args['reactions'][p_id]:
                        self.convergence_args['reactions_list'].append(r_id_str)
                else:
                    if is_ec(p_id,4):
                        self.convergence_args['reactions_complete_ec'].append([p_id,self.convergence_args['reactions'][p_id]])
                    elif is_ec(p_id,3):
                        pass
                    else:
                        self.convergence_args['reactions_not_ec'].append([p_id,self.convergence_args['reactions'][p_id]])




    def converge_gene_gpr(self):
        if self.get_gene():
            self.converge_gene_to_reaction_biocyc()
            self.converge_gene_to_protein_biocyc_complete_ec()
            #self.converge_gene_to_protein_biocyc_incomplete_ec_temp_Protein()
            self.converge_gene_to_protein_biocyc_not_ec()


    # GR
    def converge_gene_to_reaction_biocyc(self):
        #sometimes there wont be an enzyme EC in the reaction of the gene page.
        #In this case we just assume the protein has no EC and merely go by protein name
        if self.convergence_args['protein_name']:
            prot_info = {'synonyms': self.convergence_args['protein_name']}
            protein_instance = Protein(prot_info)
            match= self.get_protein_match(protein_instance)
            if match: protein_instance=match
        else: protein_instance=None
        for r_id_str in self.convergence_args['reactions_list']:
            reaction_id, reaction_str = r_id_str
            print(f'Linking from gene {self.gene_id} in {self.db} to reaction {reaction_id}')

            # convergence penalty will be lower because we want to connect the gene to the reaction as the protien is merely a placeholder
            reaction_instance = self.find_reaction(query_id=reaction_id,
                                                   extra_args={'reaction_str': reaction_str},
                                                   )
            if protein_instance:
                reaction_instance.set_detail('protein_instances',protein_instance,converged_in=self.db)
                protein_instance.set_detail('reaction_instances',reaction_instance,converged_in=self.db)
                self.get_gene().set_detail('protein_instances',protein_instance,converged_in=self.db)

            self.get_gene().set_detail('reaction_instances',reaction_instance,converged_in=self.db)

    # GP
    def converge_gene_to_protein_biocyc_complete_ec(self):
        for p in self.convergence_args['reactions_complete_ec']:

            p_id,reactions_list= p
            print(f'Linking from gene {self.gene_id} in {self.db} to protein {p_id}')

            protein_instance = self.find_protein(query_id=p_id,
                                                 extra_args={'reactions_list': reactions_list},
                                                 )
            if protein_instance:
                self.get_gene().set_detail('protein_instances',protein_instance,converged_in=self.db)

    # GP
    def converge_gene_to_protein_biocyc_incomplete_ec(self):
        #when we want full convergence, we go to the incomplete EC page and get the IDS for the specific enzyme instances,reactions and genes
        #this most likely wont happen in any pipeline, only if we want to populate the DB
        for p in self.convergence_args['reactions_incomplete_ec']:
            p_id,reactions_list= p
            print(f'Linking from gene {self.gene_id} in {self.db} to protein {p_id}')

            protein_instance = self.find_protein(query_id=p_id,
                                                 extra_args={'reactions_list': reactions_list},
                                                 )
            if protein_instance:
                self.get_gene().set_detail('protein_instances',protein_instance,converged_in=self.db)

    # GP
    def converge_gene_to_protein_biocyc_incomplete_ec_temp_Protein(self):
        # when we come from GPR pipeline we want a proper GPR network, so we connect the gene to an incomplete EC protein (which will only have basic info)
        # nonetheless this protein will then connected to a set of reaction instances, thus bridging the connection gene-general protein- reaction
        for p in self.convergence_args['reactions_incomplete_ec']:
            p_id,reactions_list= p
            prot_info = {'enzyme_ec':p_id}
            protein_instance = Protein(prot_info)
            match= self.get_protein_match(protein_instance)
            if match: protein_instance=match
            for r_id_str in reactions_list:
                reaction_id, reaction_str = r_id_str
                #convergence penalty will be lower because we want to connect the gene to the reaction as the protien is merely a placeholder
                reaction_instance = self.find_reaction(query_id=reaction_id,
                                                       extra_args={'reaction_str': reaction_str},
                                                       )
                if reaction_instance:
                    reaction_instance.set_detail('protein_instances',protein_instance,converged_in=self.db)
                    protein_instance.set_detail('reaction_instances',reaction_instance,converged_in=self.db)
                    self.get_gene().set_detail('reaction_instances',reaction_instance,converged_in=self.db)
            self.get_gene().set_detail('protein_instances',protein_instance,converged_in=self.db)


    #GP
    def converge_gene_to_protein_biocyc_not_ec(self):
        for p in self.convergence_args['reactions_not_ec']:
            p_id,reactions_list= p
            print(f'Linking from gene {self.gene_id} in {self.db} to protein {p_id}')

            protein_instance = self.find_protein(query_id=p_id,
                                                 extra_args={'reactions_list': reactions_list},
                                                 )
            if protein_instance:
                protein_instance.set_detail('gene_instances',self.get_gene(),converged_in=self.db)
                self.get_gene().set_detail('protein_instances',protein_instance,converged_in=self.db)





if __name__ == '__main__':
    from source.Biological_Components.Protein import Protein
    from source.Fetchers.Reaction_Fetchers.Reaction_Fetcher import *
    gene=Gene_Fetcher_Biocyc('HS08548')
