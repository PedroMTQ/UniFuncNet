
from drax.Fetchers.Reaction_Fetchers.Reaction_Fetcher import *

class Reaction_Fetcher_HMDB(Reaction_Fetcher):
    def __init__(self,reaction_id,extra_args={},memory_storage=None):
        Reaction_Fetcher.__init__(self,reaction_id=reaction_id,memory_storage=memory_storage)
        self.db= 'hmdb'
        self.set_convergence_args(extra_args)
        self.reaction=self.get_reaction_HMDB()
        self.add_reaction()

    def set_convergence_args(self,extra_args):
        #args for convergence
        if 'proteins_list' in extra_args: self.convergence_args['proteins_list']=extra_args['proteins_list']
        else: self.convergence_args['proteins_list']=[]
        if 'protein_id' in extra_args: self.convergence_args['proteins_list']=self.convergence_args['proteins_list'].append(extra_args['protein_id'])
        if 'protein_soup' in extra_args: self.convergence_args['protein_soup']=extra_args['protein_soup']
        else: self.convergence_args['protein_soup']=None
        #args for getting reaction
        if 'cpd_id' in extra_args: self.convergence_args['cpd_id']=extra_args['cpd_id']
        else: self.convergence_args['cpd_id']=None


    #######HMDB
    def synonyms_compound_HMDB(self, compound_ID, compound):
        url = f'http://www.hmdb.ca/metabolites/{compound_ID}'
        webpage = self.get_with_fetcher(url)
        if not webpage: return None
        soup = BeautifulSoup(webpage, 'lxml')
        try:
            common_name = soup.find('th', text='Common Name').findNext().text
            if compound.lower() == common_name.lower(): return True
        except:
            pass
        try:
            chemical_formula = soup.find('th', string='Chemical Formula')
            if chemical_formula:
                chemical_formula = chemical_formula.findNext()
                if chemical_formula.text: chemical_formula = chemical_formula.text
                if chemical_formula.lower() == compound.lower(): return True

        except:
            pass
        try:
            synonyms = soup.find('th', text='Synonyms').findNext()
            if synonyms.text != 'Not Available':
                headings = [th.get_text() for th in synonyms.find("tr").find_all("th")]
                datasets = []
                for row in synonyms.find_all("tr")[1:]:
                    dataset = zip(headings, (td.get_text() for td in row.find_all("td")))
                    datasets.append(dataset)
                for dataset in datasets:
                    for field in dataset:
                        if field[0] == 'Value':
                            synonym = field[1]
                            if synonym.lower() == compound.lower(): return True
        except:
            pass
        return False

    def list_dict_keys(self, l1, dictio):
        for i in range(len(l1)):
            if l1[i].lower() in dictio:
                l1[i] = dictio[l1[i].lower()]
            else:
                l1[i] = self.find_compound_ID_HMDB(l1[i])
        return l1

    def find_compound_ID_HMDB(self, compound):
        # this will only check the first page, which should already contain the best results
        url = f'http://www.hmdb.ca/unearth/q?utf8=%E2%9C%93&query={compound}&searcher=metabolites&button='
        webpage = self.get_with_fetcher(url)
        if not webpage: return None
        soup = BeautifulSoup(webpage, 'lxml')
        metabolites = soup.find_all('div', class_='result-link')
        for result in metabolites:
            metabolite_ID = result.find('a', class_='btn-card').text
            found = self.synonyms_compound_HMDB(metabolite_ID, compound)
            if found:  return metabolite_ID

    def reaction_IDs_HMDB(self, soup):
        search = soup.find('div', class_='reaction-panel')
        if not search: return None
        search=search.find('div',class_='panel-heading').strong
        sub, prod = search.text.replace('\n', '').replace('\t', '').split('=')
        sub = [i for i in sub.split('+')]
        prod = [i for i in prod.split('+')]
        temp = search.find_all('a')
        compound_dic = {}
        for compound in temp:
            name = compound.text.lower()
            reaction_id = compound['href'][re.search('metabolites/', compound['href']).span()[1]:]
            compound_dic[name] = reaction_id
        sub, prod = self.list_dict_keys(sub, compound_dic), self.list_dict_keys(prod, compound_dic)
        return sub + prod

    def get_details_reaction_compound_page(self, details_source):
        details_list = [[th.get_text(), ''] for th in details_drax.find_all("dt")]
        details_text = details_drax.find_all("dd")
        for row in range(len(details_list)):
            details_list[row][1] = details_text[row].text.strip()
        details = {i[0]: i[1] for i in details_list}
        return details

    def get_table_detail(self, hmdb_id, database, detail):  # gets a certain detail from HMDB
        # database= proteins for enzymes, metabolites for compounds, only works for these 2
        url = f'http://www.hmdb.ca/{database}/{hmdb_id}'
        try:
            webpage = self.get_with_fetcher(url)
            if not webpage: return None
            soup = BeautifulSoup(webpage, 'lxml')
            search = soup.find('th', string=re.compile(detail))
            return search.findNext().text
        except:
            return None


    def get_reaction_HMDB(self):

        url = f'http://www.hmdb.ca/reactions/{self.reaction_id}'
        webpage = self.get_with_fetcher(url)
        if not webpage: return None
        reaction_soup = BeautifulSoup(webpage, 'lxml')
        #if self.check_correct_reaction_HMDB(reaction_soup,self.convergence_args['cpd_id']):
        rn_with_ids = self.reaction_IDs_HMDB(reaction_soup)
        reaction_str = reaction_soup.find('div', class_='reaction-panel').strong.text
        reaction_str=reaction_str.replace('+',' + ').replace('=',' = ')
        reaction_str = fix_html_sign(reaction_str)
        try:
            rn_with_ids, complete_l, len_sub = get_stoichiometry(reaction_str, rn_with_ids)
        except:
            rn_with_ids, complete_l, len_sub = get_stoichiometry(reaction_str, reaction_str)
        rn_with_instances = self.reaction_met_instances(reaction_str, rn_with_ids, 'hmdb')
        rn_info={
            'hmdb': self.reaction_id,
            'reaction_str': reaction_str,
        }
        external_links= reaction_soup.find('h4',text='External Links')
        if external_links:
            external_links=external_links.findNext('ol').findChildren('li')
            for ext_link in external_links:
                if 'kegg' in ext_link.text.lower():
                    kegg_id= ext_link.text.split(':')[1].strip()
                    rn_info['kegg']=kegg_id

        if rn_with_instances: rn_info['reaction_with_instances'] = rn_with_instances
        else:                 rn_info['rn_with_ids'] = [reaction_str, rn_with_ids, 'hmdb']
        reaction_instance = Reaction(rn_info)
        self.get_proteins_from_soup(reaction_soup)
        return reaction_instance

    def get_proteins_from_soup(self,soup):
        protein_ids=soup.find_all(href=re.compile('/proteins/HMDBP\d+'))
        self.convergence_args['proteins_list'] = [re.search('HMDBP\d+', str(i)).group() for i in protein_ids]


    def extract_IDs_HMDB(self, soup):
        search = soup.find('div', class_='reaction-panel')
        res = []
        if search:
            search=search.strong
            temp = search.find_all('a')
            compound_dic = {}
            for compound in temp:
                name = compound.text
                compound_id = compound['href'][re.search('metabolites/', compound['href']).span()[1]:]
                compound_dic[name] = compound_id
                res.append(re.search('HMDB\d+', str(compound)).group())
        return res

    def check_correct_reaction_HMDB(self, reaction_soup,cpd_id):
        if not cpd_id: return True
        res = self.extract_IDs_HMDB(reaction_soup)
        for i in res:
            if i == cpd_id: return True
        return False


    def converge_reaction_rpg(self):
        #RP part of the CRPG pipeline
        if self.convergence_args['proteins_list']:
            self.converge_reaction_to_protein()

    def converge_reaction_to_protein(self):
        for protein_id in self.convergence_args['proteins_list']:
            print(f'Linking from reaction {self.reaction_id} in {self.db} to protein {protein_id}')

            protein_instance = self.find_protein(query_id=protein_id,
                                                 extra_args={'protein_soup': self.convergence_args['protein_soup']},
                                                 )
            if protein_instance:
                self.get_reaction().set_detail('protein_instances',protein_instance,converged_in='hmdb')


if __name__ == '__main__':
    r=Reaction_Fetcher_HMDB('1746')
    r.find_compound_ID_HMDB('A phenyl acetate')



