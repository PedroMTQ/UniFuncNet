
from source.Fetchers.Protein_Fetchers.Protein_Fetcher import *

class Protein_Fetcher_Biocyc(Protein_Fetcher):
    def __init__(self,protein_id=None,extra_args={},memory_storage=None,init_Fetcher=True):
        Protein_Fetcher.__init__( self, protein_id= protein_id, memory_storage=memory_storage)
        self.db='biocyc'
        self.set_convergence_args(extra_args)
        if init_Fetcher:
            self.protein=self.get_protein_biocyc()
            self.add_protein()

    def set_convergence_args(self,extra_args):
        #args for convergence
        if self.get_protein():
            if 'enzyme_ec' in self.convergence_args: self.protein.set_detail('enzyme_ec', self.convergence_args['enzyme_ec'])
        if 'reactions_list' not in self.convergence_args: self.convergence_args['reactions_list'] = []
        if 'reactions_list' in extra_args: self.convergence_args['reactions_list'].extend(extra_args['reactions_list'])
        #when we get an EC from a family rather than the specific enzyme (we cant just ignore these as they contain reacitons which may not have been attributed to a specific enzyme!)
        if 'enzyme_ecs' not in self.convergence_args: self.convergence_args['enzyme_ecs']=[]
        if 'protein_and_genes' not in self.convergence_args: self.convergence_args['protein_and_genes']=[]
        #args for getting protein
        if 'protein_names' in extra_args:  self.convergence_args['protein_names'] = extra_args['protein_names']
        else:                              self.convergence_args['protein_names'] = []
        if 'genes' in extra_args:          self.convergence_args['genes'] = extra_args['genes']
        else:                              self.convergence_args['genes'] = []


    def get_protein_biocyc(self):
        protein_names=self.convergence_args['protein_names']
        genes=self.convergence_args['genes']
        if genes:
            return self.get_protein_genes(genes,protein_names)
        else:
            if is_ec(self.protein_id):
                self.convergence_args['enzyme_ec']=self.protein_id
                url = 'https://biocyc.org/META/NEW-IMAGE?type=EC-NUMBER&object=EC-' + clean_ec(self.protein_id)
                return self.get_protein_biocyc_ec_number(url)
            else:
                url='https://biocyc.org/META/NEW-IMAGE?type=ENZYME&object='+self.protein_id
                return self.get_protein_biocyc_enzyme(url)

    def get_protein_genes(self,genes,protein_names):
        protein_info={'synonyms': protein_names,'biocyc':self.protein_id,}
        protein_gene_instance = Protein(protein_info)
        self.protein= protein_gene_instance
        self.add_protein()
        for gene in genes:
            fetched_gene = self.find_gene(query_id=gene,
                                          extra_args={'protein_gene_instance':protein_gene_instance},
                                          )
            if fetched_gene:
                protein_gene_instance.set_detail('gene_instances',fetched_gene,converged_in=self.db)
        # this will be a placeholder protein that establishes the connection protein-gene.
        # It can generate from reaction fetcher or ec_protein fetcher
        return protein_gene_instance

    def get_protein_biocyc_ec_number(self,enzyme_url):
        #this is more of a general page, here we can find several proteins too
        webpage = self.get_with_fetcher(enzyme_url)
        if not webpage: return None
        soup = BeautifulSoup(webpage, 'lxml')
        s1 = 'Systematic Name: \n'
        s2 = 'Synonyms: \n'
        s3 = 'Summary: \n'
        s4 = 'Unification Links: \n'
        enz_syns, extra_syns, enz_summary = [], [], None
        paragraphs = soup.find_all('p', class_='ecoparagraph')
        unification_links = {}
        for i in paragraphs:
            if s1 in i.text:
                enz_syns = [i.text[len(s1):].strip().lower()]
            elif s2 in i.text:
                extra_syns = i.text[len(s2):].strip().lower()
                if ', ' in extra_syns:
                    extra_syns = extra_syns.split(', ')
                else:
                    extra_syns = [extra_syns]
            elif s1 in i.text:
                enz_summary = i.text[len(s3):].strip().replace('\n', '')
            elif s4 in i.text:
                unification_links = [j.strip().split(':') for j in
                                     i.text[len(s4):].strip().replace('\n', '').split(',')]
                unification_links = {i[0].lower(): i[1] for i in unification_links}

        if extra_syns:
            enz_syns.extend(extra_syns)  # so I get all the possible names in one list
        for i_syn in range(len(enz_syns)):
            if ':' in enz_syns[i_syn]:
                enz_syns[i_syn] = enz_syns[i_syn].split(':')[1]
        enz_syns = [i.strip('\n') for i in enz_syns]
        #this is the main EC Protein, in truth it will probably be merged with others further along
        protein_info={'synonyms': enz_syns,
             'biocyc': self.protein_id,
             'enzyme_ec': self.protein_id,
             'brenda': unification_links['brenda'] if 'brenda' in unification_links else None,
             }
        self.get_other_info_ec_page(soup)
        return Protein(protein_info)

    def get_other_info_ec_page(self,soup):
        self.get_reactions_info(soup)
        self.get_transferred_ecs(soup)
        self.get_enzyme_ecs_info(soup)
        self.get_enzymes_and_genes_info(soup)

    def get_reactions_info(self,soup):
        reactions_soup = soup.find_all('a', class_='REACTION')
        if reactions_soup:
            for r in reactions_soup:
                r_link = r['href']
                reaction_id = r_link[re.search('object=', r_link).span()[1]:]
                reaction_str = r.text.replace('\n', '')
                reaction_str=fix_html_sign(reaction_str)
                self.convergence_args['reactions_list'].append([reaction_id, reaction_str])

    def get_transferred_ecs(self,soup):
        enzyme_ecs = soup.find(text=re.compile('has been transferred by the Enzyme Commission\. It is now listed as'))
        if not enzyme_ecs: return None
        print('TRANSFERRED ECS',self.protein_id)
        if enzyme_ecs: enzyme_ecs = enzyme_ecs.parent
        enzyme_ecs=enzyme_ecs.find_all('a')
        for i in enzyme_ecs:
            protein_id,protein_name=i.text.split('--')
            protein_name=protein_name.strip()
            protein_id=protein_id.strip()
            if protein_id!=self.protein_id:
                self.convergence_args['enzyme_ecs'].append({'protein_id': protein_id, 'protein_name': protein_name})

    def get_enzyme_ecs_info(self,soup):
        enzyme_ecs = soup.find(text=re.compile('.*Instances:.*'))
        if not enzyme_ecs: return None
        if enzyme_ecs: enzyme_ecs = enzyme_ecs.parent
        temp = None
        for i in enzyme_ecs:
            if i.name == 'br':
                if temp: self.convergence_args['enzyme_ecs'].append(temp)
                temp = {'protein_id': None, 'protein_name': None}
            else:
                if not isinstance(i, str):
                    if i['href']:
                        if 'EC-NUMBER' in i['class']:
                            if 'type=ENZYME&object=' in i['href']:
                                protein_id = re.search('\?type=ENZYME&object=.*', i['href']).group()[20:]
                            elif 'gene?orgid=META&id=' in i['href']:
                                protein_id = re.search('gene\?orgid=META&id=.*', i['href']).group()[19:]
                            else:
                                protein_id = re.search('\?type=EC-NUMBER&object=EC-.*', i['href']).group()[26:]
                            temp['protein_id'] = protein_id
                            protein_name = i.text.split('--')
                            if protein_name:
                                protein_name = protein_name[1].strip()
                                temp['protein_name'] = protein_name
            self.convergence_args['enzyme_ecs'].append(temp)


    def get_enzymes_and_genes_info(self,soup):
        enzymes_and_genes = soup.find(text=re.compile('.*Enzymes and Genes:.*'))
        if not enzymes_and_genes: return None
        if enzymes_and_genes: enzymes_and_genes = enzymes_and_genes.parent
        temp = None
        for i in enzymes_and_genes:
            if i.name == 'br':
                if temp: self.convergence_args['protein_and_genes'].append(temp)
                temp = {'protein_id': None,'protein_name':None, 'genes': []}
            else:
                if not isinstance(i, str):
                    if i['href']:
                        if 'ENZYME' in i['class']:
                            if 'type=ENZYME&object=' in i['href']:
                                protein_id = re.search('\?type=ENZYME&object=.*', i['href']).group().split('object=')[-1]
                            elif 'gene?orgid=' in i['href']:
                                protein_id = re.search('gene\?orgid=.*', i['href']).group().split('id=')[-1]
                            else:
                                protein_id = re.search('\?type=EC-NUMBER&object=EC-.*', i['href']).group().split('object=EC-')[-1]
                            temp['protein_id'] = protein_id
                            temp['protein_name'] = i.text
                        elif 'GENE' in i['class']:
                            if 'gene?orgid=' in i['href']:
                                gene_id = re.search('gene\?orgid=.*', i['href']).group().split('id=')[-1]
                                temp['genes'].append(gene_id)

        self.convergence_args['protein_and_genes'].append(temp)

    def get_protein_biocyc_enzyme(self,enzyme_url):
        webpage = self.get_with_fetcher(enzyme_url, selenium=True,original_response=True)
        if not webpage: return None
        page_source=webpage.page_source
        current_url=webpage.current_url
        if 'selenium.webdriver' in repr(webpage): webpage.quit()
        if current_url!= enzyme_url:
            #sometimes biocyc will redirect to gene page, we want to extract info from this page too
            #this will only happen when coming from the protein_searcher method
            return self.redirected_protein_page(page_source)
        else:
            return self.NOT_redirected_protein_page(page_source)

    def redirected_protein_page(self,webpage):
        soup = BeautifulSoup(webpage, 'lxml')
        headers = soup.find_all('font', class_='header')
        gene, enzyme_name = None, None
        protein_info = {'biocyc': self.protein_id}
        for h in headers:
            previous = h.parent
            if 'enzyme' in previous.text:
                enzyme_name = h.text
        protein_info['synonyms']=enzyme_name
        self.get_ids_soup(soup,protein_info)
        return Protein(protein_info)

    def get_ids_soup(self,soup,protein_info):
        corresponding_dbs = {
            'UniProt': 'uniprot',
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
                        protein_info[corresponding_dbs[db.text]]=corresponding_db_info



    def NOT_redirected_protein_page(self,webpage):
        soup = BeautifulSoup(webpage, 'lxml')
        paragraphs = soup.find_all('p', class_='ecoparagraph')
        if not paragraphs: return None
        enz_syns=soup.find('font',class_='header pageTitle')
        if enz_syns: enz_syns=enz_syns.text.split(':')[1]
        protein_info = {'database': 'biocyc',
             'synonyms': enz_syns if enz_syns else [],
             'biocyc': self.protein_id,
             }
        #this page also has reactions
        reaction_boxes_soup = soup.find_all(text=re.compile('Enzymatic reaction:.*'))
        if reaction_boxes_soup:
            for box_soup in reaction_boxes_soup:
                reaction_soup = box_soup.parent.findNext('a', class_='REACTION')
                r_link = reaction_soup['href']
                reaction_id = r_link[re.search('object=', r_link).span()[1]:]
                reaction_str = reaction_soup.text.replace('\n', '')
                reaction_str=fix_html_sign(reaction_str)
                self.convergence_args['reactions_list'].append([reaction_id,reaction_str])
        return Protein(protein_info)



    def converge_protein_global(self):
        self.converge_protein_to_gene()
        self.converge_protein_to_reaction_biocyc()

    def converge_protein_gpr(self):
        self.converge_protein_to_reaction_biocyc()

    def converge_protein_rpg(self):
        #PG part of the CRPG pipeline
        self.converge_protein_to_gene()

    #PR
    def converge_protein_to_reaction_biocyc(self):
        if self.convergence_args['reactions_list']:
            for reaction_id,reaction_str in self.convergence_args['reactions_list']:
                reaction_instance= self.find_reaction(reaction_id,
                                                      extra_args={'reaction_str':reaction_str},
                                                      )
                if reaction_instance:
                    self.get_protein().set_detail('reaction_instances',reaction_instance,converged_in=self.db)

    #PP (when we start from an enzyme family)
    def converge_protein_to_ec(self):
        if self.convergence_args['enzyme_ecs']:
            for prt in self.convergence_args['enzyme_ecs']:
                prt_id=prt['protein_id']
                protein_name=prt['protein_name']
                #this will be a protein with both protein and gene info
                #exceptionally here we recursively call upon the class type
                protein_instance  = Protein_Fetcher_Biocyc (protein_id=prt_id,
                                        extra_args={'protein_names':protein_name},
                                        memory_storage=self.memory_storage,
                                        ).get_protein()
                #biocyc has info on enzyme families,we only set two protein instances with the same id if they actually have info on all levels
                if is_ec(self.protein_id,4):
                    if protein_instance:
                        protein_instance.set_detail('enzyme_ec',self.protein_id)
                        self.get_protein().unite_instances(protein_instance)

    #PG
    def converge_protein_to_gene(self):
        if self.convergence_args['protein_and_genes']:
            #this will retrieve the genes and also create separate protein instances that will be merged with the current protein
            for prt_gene in self.convergence_args['protein_and_genes']:
                genes= prt_gene['genes']
                protein_id= prt_gene['protein_id']
                protein_name= prt_gene['protein_name']
                #this will be a protein with only genes info
                protein_instance  = Protein_Fetcher_Biocyc (protein_id,extra_args={'protein_names':protein_name,
                                                                        'genes':genes},
                                                                        memory_storage=self.memory_storage,
                                                                        ).get_protein()
                if is_ec(self.protein_id,4):
                    if protein_instance:
                        protein_instance.set_detail('enzyme_ec',self.protein_id)
                        self.get_protein().unite_instances(protein_instance)








if __name__ == '__main__':
    from source.Biological_Components.Protein import Protein
    from source.Fetchers.Reaction_Fetchers.Reaction_Fetcher import *
    protein=Protein_Fetcher_Biocyc('2.7.1.1')




