
from source.Fetchers.Reaction_Fetchers.Reaction_Fetcher import *

class Reaction_Fetcher_Biocyc(Reaction_Fetcher):
    def __init__(self,reaction_id,extra_args={},memory_storage=None,init_Fetcher=True):
        Reaction_Fetcher.__init__(self,reaction_id=reaction_id,memory_storage=memory_storage)
        self.db = 'biocyc'
        self.set_convergence_args(extra_args)
        if init_Fetcher:
            self.reaction= self.get_reaction_biocyc()
            self.add_reaction()

    def reaction_IDs_biocyc(self, reaction_id):
        url = f'https://biocyc.org/rxn-svg?r={reaction_id}&o=META&n=n&c=y&ai=0&h=y'
        res = None
        webpage = self.get_with_fetcher(url)
        if webpage:
            reaction_soup = BeautifulSoup(webpage, 'lxml')
            r_ids = reaction_soup.find_all('a', class_='cpd_a')
            res = []
            for r_id in r_ids:
                if 'xlink:href' in r_id.attrs:
                    search_pattern = re.compile('META&id=|&object=')
                    search = re.search(search_pattern, r_id['xlink:href'])
                    search_delimiter = search.span()[1]
                    res.append(r_id['xlink:href'][search_delimiter:])
        return res

    def clean_ec(self,original_ec):
        temp_ec = str(original_ec)
        if '-' in temp_ec[-1]:
            temp_ec = temp_ec[:-2]
        return temp_ec

    def get_rn_str(self,enzyme_ec):
        temp_ec=self.clean_ec(str(enzyme_ec))
        webpage = self.get_with_fetcher(f'https://biocyc.org/META/NEW-IMAGE?type=EC-NUMBER&object=EC-{temp_ec}', selenium=True)
        if not webpage: return None
        soup = BeautifulSoup(webpage, 'lxml')
        for extra_rn_id in [self.reaction_id]+self.convergence_args['extra_rn_ids']:
            print('trying',extra_rn_id)
            reaction_str = soup.find('a',{'href':re.compile('.*'+extra_rn_id)})
            if reaction_str:
                return fix_html_sign(reaction_str.text)

    def get_genes_protein_orgs(self):
        url=f'https://biocyc.org/META/reaction-genes?object={self.reaction_id}'
        webpage = self.get_with_fetcher(url)
        if webpage:
            get_r_str=False
            start_recording=False
            reaction_str=None
            to_parse= webpage.split('\n')
            protein_id=None
            for line in to_parse:
                if start_recording:
                    line=line.split('\t')
                    if len(line)>1:
                        gene_id,gene_accession,gene_name,enzymatic_activity,evidence,organism = line
                        self.convergence_args['genes'].append([protein_id,enzymatic_activity, gene_id])
                if not start_recording:
                    if get_r_str:
                        reaction_str= line
                        reaction_str = fix_html_sign(reaction_str)
                        get_r_str=False
                    if line.strip() == self.reaction_id: get_r_str=True
                    if 'EC-' in line:
                        protein_id=line.replace('EC-','').strip()
                    if 'Gene ID' in line:
                        start_recording=True
            if reaction_str:
                if not self.convergence_args['reaction_str']:
                   self.convergence_args['reaction_str'] = reaction_str

    def clean_reaction_id(self,xpath_target):
        current_reaction_id = xpath_target.replace('id="', '').replace('-toggle"', '')
        search=re.search('_(TAX|ORG)',current_reaction_id)
        if search:
            current_reaction_id=current_reaction_id[0:search.span()[0]]
            return current_reaction_id

    def get_protein_orgs(self,soup):
        xpaths=[]
        xpaths_ids=[]
        xpaths_ids_to_load=[]
        show_mores= soup.find_all(text=re.compile(' \[Show \d+ more.*\]'))
        for s in show_mores:
            xpath_target= re.search('id=\"'+self.reaction_id+'.*-toggle"',str(s.parent))
            if xpath_target:
                xpath_target=xpath_target.group()
                current_reaction_id=self.reaction_id
            else:
                xpath_target = re.search('id=\".*RXN.*-toggle"', str(s.parent))
                xpath_target=xpath_target.group()
                current_reaction_id=self.clean_reaction_id(xpath_target)
            xpath_id= re.search('id=\"'+regex_escape(current_reaction_id)+'.*-toggle',str(s.parent)).group()[4:-9]


            #somtimes we can find alternative ids for the same reaction
            if current_reaction_id not in self.convergence_args['extra_rn_ids']: self.convergence_args['extra_rn_ids'].append(current_reaction_id)
            xpaths_ids.append(xpath_id)
            xpath_target='//*[@'+xpath_target+']'
            xpaths.append(xpath_target)
        for x_id in xpaths_ids:
            ids_to_load = soup.find_all(id=re.compile('.*'+x_id+'.*'))
            for i in ids_to_load:
                toggled_id = re.search('.*-\d+-content',i['id'])
                if toggled_id:
                    found=toggled_id.group()
                    if '-1-content' not in found:
                        xpaths_ids_to_load.append(found)
        if xpaths:
            url = f'https://biocyc.org/META/NEW-IMAGE?type=REACTION&object={self.reaction_id}'
            webpage = self.get_with_fetcher(url,selenium=True,xpath=xpaths,ids_to_load=xpaths_ids_to_load)
            if not webpage: return None
            soup = BeautifulSoup(webpage, 'lxml')


        organisms = soup.find_all(class_='ORGANISM')
        if not organisms: return None
        organisms = [org.find_parent('tr') for org in organisms]
        for o in organisms:
            if o:
                temp_table = [[]]
                info=o.td.findNext('td')
                line = 0
                for i in info:
                    if '<br/>' in str(i):
                        line += 1
                        temp_table.append([])
                    else:
                        if '<script>' not in str(i):
                            temp_table[line].append(i)
                temp_table = [i for i in temp_table if i]
                for i in temp_table:
                    if not re.search('.*/gene\?orgid=.*',str(i)):
                        if 'href' in i[0]:
                            to_add = i[0]['href'].split('=')[-1]
                            if to_add not in self.convergence_args['proteins']:
                                self.convergence_args['proteins'].append(to_add)

    def get_reaction_biocyc(self):  # this basically gets all the details from biocyc, may be able to implement a more concise search, see later
        if not self.reaction_id: return None
        self.get_genes_protein_orgs()
        reaction_str=self.convergence_args['reaction_str']
        url = f'https://biocyc.org/META/NEW-IMAGE?type=REACTION&object={self.reaction_id}'
        enzyme_ec,Unification_links, Relationship_links, pathways = None, [], [], []
        webpage = self.get_with_fetcher(url)
        if not webpage: return None
        soup = BeautifulSoup(webpage, 'lxml')
        self.get_protein_orgs(soup)

        enzyme_ecs = soup.find('a', class_='EC-NUMBER')
        if enzyme_ecs:
            self.convergence_args['enzyme_ecs']=enzyme_ecs.text.split(',')
            if not reaction_str:
                for enzyme_ec in self.convergence_args['enzyme_ecs']:
                    if reaction_str: break
                    else: reaction_str=self.get_rn_str(enzyme_ec)
        #THIS MIGHT BE USEFUL FOR FASTER COMPOUND MATCHING
        rn_with_ids = self.reaction_IDs_biocyc(self.reaction_id)
        try:
            rn_with_ids, complete_l, len_sub = get_stoichiometry(reaction_str, rn_with_ids)
        except:
            try:
                rn_with_ids, complete_l, len_sub  = get_stoichiometry(reaction_str, reaction_str)
            except:
                return None
        rn_with_instances = self.reaction_met_instances(reaction_str, rn_with_ids, 'biocyc')
        paragraphs = soup.find_all('p',class_='ecoparagraph')
        for i in paragraphs:
            s0 = ' Unification Links: \n'
            s1 = ' Relationship Links: \n'
            if s0 in i.text:
                Unification_links = [j.strip() for j in i.text[len(s0):].strip().replace('\n', '').split(',')]
            elif s1 in i.text:
                Relationship_links = [j.strip() for j in i.text[len(s1):].strip().replace('\n', '').split(',')]
        pathway_search = soup.find_all('a', class_='PATHWAY')
        pathways = [i.text for i in pathway_search]
        kegg_id,rhea_id=None,None
        for uni in Unification_links:
            if 'kegg' in uni.lower():
                kegg_id=uni.split(':')[1]
            elif 'rhea' in uni.lower():
                rhea_id=uni.split(':')[1]
        reaction_dict = {
            'biocyc': self.reaction_id,
            'reaction_str': reaction_str,
            'kegg':kegg_id,
            'rhea':rhea_id,
            'pathways': pathways
        }
        if rn_with_instances:  reaction_dict['reaction_with_instances']= rn_with_instances
        else:
            if rn_with_ids: reaction_dict['rn_with_ids']= [reaction_str, rn_with_ids, 'biocyc']
        reaction_instance = Reaction(reaction_dict)
        for extra_id in self.convergence_args['extra_rn_ids']: reaction_instance.set_detail('biocyc',extra_id)
        return reaction_instance


    def set_convergence_args(self,extra_args):
        #args for getting reaction
        if 'reaction_str' in extra_args: self.convergence_args['reaction_str'] =extra_args['reaction_str']
        else:self.convergence_args['reaction_str']= None
        #args for convergence
        if 'enzyme_ecs' not in self.convergence_args: self.convergence_args['enzyme_ecs'] = []
        if 'proteins' not in self.convergence_args: self.convergence_args['proteins'] = []
        if 'genes' not in self.convergence_args: self.convergence_args['genes'] = []
        if 'extra_rn_ids' not in self.convergence_args: self.convergence_args['extra_rn_ids']=[]


    def converge_reaction_rpg(self):
        #RP and RG part of the CRPG pipeline
        self.converge_reaction_to_protein_ec_biocyc()
        self.converge_reaction_to_protein_biocyc()
        self.converge_reaction_to_genes_biocyc()

    def merge_enzyme_ec_genes(self):
        self.enzyme_ec_to_gene={}
        for enzyme_ec in self.convergence_args['enzyme_ecs']:
            for gene in self.convergence_args['genes']:
                protein_id,enzyme_names,gene_id= gene
                if enzyme_ec==protein_id:
                    if enzyme_ec not in self.enzyme_ec_to_gene: self.enzyme_ec_to_gene[enzyme_ec]=[]
                    self.enzyme_ec_to_gene[enzyme_ec]


    #RP with enzyme ec
    def converge_reaction_to_protein_ec_biocyc(self):
        if self.convergence_args['enzyme_ecs']:
            for enzyme_ec in self.convergence_args['enzyme_ecs']:
                print(f'Linking from reaction {self.reaction_id} in {self.db} to protein {enzyme_ec}')

                protein_instance = self.find_protein(query_id=enzyme_ec)
                if protein_instance:
                    self.get_reaction().set_detail('protein_instances',protein_instance,converged_in=self.db)

    #RP with normal metacyc enzyme id
    def converge_reaction_to_protein_biocyc(self):
        if self.convergence_args['proteins']:
            for prt_id in self.convergence_args['proteins']:
                print(f'Linking from reaction {self.reaction_id} in {self.db} to protein {prt_id}')

                #this will be a protein with both protein and gene info
                protein_instance  = self.find_protein (query_id=prt_id,
                                                       extra_args={},
                                                       )
                if protein_instance:
                    self.reaction.set_detail('protein_instances',protein_instance,converged_in=self.db)

    #RG
    def converge_reaction_to_genes_biocyc(self):
        if self.convergence_args['genes']:
            for gene in self.convergence_args['genes']:
                protein_id,enzyme_names,gene_id= gene
                if protein_id:
                    print(f'Linking from reaction {self.reaction_id} in {self.db} to protein {protein_id} (gene {gene_id})')
                else:
                    print(f'Linking from reaction {self.reaction_id} in {self.db} gene {gene_id}')


                #this will be a protein with only genes info
                #convergence penalty will be lower because we want to connect the gene to the reaction as the protien is merely a placeholder
                protein_instance  = self.find_protein(query_id=protein_id,extra_args={'protein_names':enzyme_names,'genes':[gene_id]})
                if protein_instance:
                    self.reaction.set_detail('protein_instances',protein_instance,converged_in=self.db)





if __name__ == '__main__':
    #rn_search=Reaction_Fetcher_Biocyc('RXN-2043')
    rn_search=Reaction_Fetcher_Biocyc('RXN-20993')
    rn_search.reaction.get_all_info()


