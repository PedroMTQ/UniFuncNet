import re

from source.Fetchers.Compound_Fetchers.Compound_Fetcher import *

class Compound_Fetcher_HMDB(Compound_Fetcher):
    def __init__(self,compound_id,memory_storage=None):
        Compound_Fetcher.__init__( self, compound_id,memory_storage)
        self.db='hmdb'
        self.set_convergence_args()

        self.compound=self.get_compound_hmdb()
        self.add_compound()

    def set_convergence_args(self):
        #args for convergence
        self.convergence_args['proteins'] = set()
        self.convergence_args['soup'] = set()

    def remove_wrong_protein_type(self, table_dic, desired_type):
        correct_type, incorrect_type=set(),set()
        for r in range(len(table_dic['Type'])):
            protein_id = table_dic['HMDBP ID'][r]
            if table_dic['Type'][r] == desired_type:
                correct_type.add(protein_id)
            else:
                incorrect_type.add(protein_id)

        return correct_type, incorrect_type

    def get_all_enzymes(self):
        url=f'https://hmdb.ca//metabolites/{self.compound_id}/metabolite_protein_links'
        webpage = self.get_with_fetcher(url)
        if not webpage: return None
        soup = BeautifulSoup(webpage, 'lxml')
        table = soup.find('table', class_='table table-condensed table-striped metabolite-protein-links proteins')
        if not table: return None
        table_head = table.thead.tr.find_all('th')
        table_dic = {}
        for i in table_head: table_dic[i.text.strip()] = []
        # to see how many pages
        number_pages = soup.find('li', class_='last next')
        if number_pages:
            number_pages = number_pages.a['href']
            number_pages = int(number_pages[re.search('page=', number_pages).span()[1]:])
        else:
            number_pages = 1
        i = 1
        while i <= number_pages:
            if i > 1:
                url = f'http://www.hmdb.ca/metabolites/{self.compound_id}/metabolite_protein_links?c=hmdb_id&d=up&page={i}'
                webpage = self.get_with_fetcher(url)
                soup = BeautifulSoup(webpage, 'lxml')
            table = soup.find('table', class_='table table-condensed table-striped metabolite-protein-links proteins')
            table_body = table.tbody.find_all('tr')
            for r in table_body:
                cols = r.find_all('td')
                for c in range(len(cols)):
                    table_dic[list(table_dic.keys())[c]].append(cols[c].text.strip())
            i += 1
        # not interested in transporters for example
        correct_type,incorrect_type = self.remove_wrong_protein_type(table_dic, 'Enzyme')
        return correct_type


    def get_enzymes(self,soup):
        redirect=soup.find(string='Show all enzymes and transporters')
        if redirect:
            return self.get_all_enzymes()
        enzymes = soup.find('h2', string='Enzymes')
        if enzymes: enzymes = enzymes.findNext()
        enzymes=re.findall('href="/proteins/HMDBP\d+',str(enzymes))
        res=set([i.replace('href="/proteins/','') for i in enzymes])
        return res






    def get_compound_hmdb(self):
        external_DBs = {'DrugBank ID': 'drugbank',
                        'Phenol Explorer Compound ID': 'phenol_explorer',
                        'FoodDB ID': 'food_db',
                        'KNApSAcK ID': 'knapsack',
                        'Chemspider ID': 'chemspider',
                        'KEGG Compound ID': 'kegg',
                        'BioCyc ID': 'metacyc',
                        'BiGG ID': 'bigg',
                        'METLIN ID': 'metlin',
                        'PubChem Compound': 'puchem_cid',
                        'PDB ID': 'pdb',
                        'ChEBI ID': 'chebi'}
        url = f'http://www.hmdb.ca/metabolites/{self.compound_id}'
        res = {}
        res['hmdb'] = self.compound_id
        webpage = self.get_with_fetcher(url)
        if not webpage: return Compound(res)
        soup = BeautifulSoup(webpage, 'lxml')
        location_search = soup.find('th', text='Biospecimen Locations')
        if location_search:
            location_search = location_search.findNext()
            location = location_search.text.split('\n')
            if 'Not Available' in location: location=None
            res['Location'] = location
        SMILES = soup.find('th', string='SMILES')
        if SMILES: SMILES = SMILES.findNext().text
        res['smiles'] = SMILES
        inchi_key = soup.find('th', string='InChI Key')
        if inchi_key: inchi_key = remove_inchi_key_equal(inchi_key.findNext().text)
        res['inchi_key'] = inchi_key
        inchi = soup.find('th', string='InChI Identifier')
        if inchi: inchi = remove_inchi_key_equal(inchi.findNext().text)
        res['inchi'] = inchi
        synonyms = []
        common_name = soup.find('th', string='Common Name')
        if common_name: common_name = common_name.findNext().text.lower()
        synonyms.insert(0, common_name)
        res['synonyms'] = synonyms
        chemical_formula = soup.find('th', string='Chemical Formula')
        if chemical_formula:
            chemical_formula = chemical_formula.findNext()
            if chemical_formula.text: res['chemical_formula'] = chemical_formula.text
        cas = soup.find('th', string='CAS Registry Number')
        if cas:
            cas=cas.findNext()
            if cas.text: res['cas'] = cas.text
        for i in external_DBs:
            external_link_search = soup.find('th', string=i)
            if external_link_search:
                external_link = external_link_search.findNext().text.strip()
                res[external_DBs[i]] = external_link
            res = self.remove_unwanted_info(res)
        res=self.remove_unwanted_info(res)
        if number_of_nones_dict(res)==len(res): return None
        compound_instance= Compound(res)
        self.convergence_args['soup']=soup
        return compound_instance

    def converge_compound_global(self):
        self.converge_compound_to_protein()





    #RP with enzyme ec
    def converge_compound_to_protein(self):
        self.convergence_args['proteins']=self.get_enzymes(self.convergence_args['soup'])
        if self.convergence_args['proteins']:
            for protein_id in self.convergence_args['proteins']:
                print(f'Linking from compound {self.compound_id} in {self.db} to protein {protein_id}')
                self.find_protein(query_id=protein_id,extra_args={'cpd_id':self.compound_id})
        self.convergence_args['soup']=None

if __name__ == '__main__':
    c=Compound_Fetcher_HMDB('HMDB0002111')
    c.get_compound().get_all_info()
    #c=Compound_Fetcher_HMDB('HMDB0005794')
    #print(c.get_compound())