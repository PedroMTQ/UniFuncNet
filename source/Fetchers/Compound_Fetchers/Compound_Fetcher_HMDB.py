
from source.Fetchers.Compound_Fetchers.Compound_Fetcher import *

class Compound_Fetcher_HMDB(Compound_Fetcher):
    def __init__(self,compound_id,memory_storage=None):
        Compound_Fetcher.__init__( self, compound_id,memory_storage)
        self.db='hmdb'
        self.compound=self.get_compound_hmdb()
        self.add_compound()

    def get_compound_hmdb(self):
        external_DBs = {'DrugBank ID': 'drugbank',
                        'Phenol Explorer Compound ID': 'phenol_explorer',
                        'FoodDB ID': 'food_db',
                        'KNApSAcK ID': 'knapsack',
                        'Chemspider ID': 'chemspider',
                        'KEGG Compound ID': 'kegg',
                        'BioCyc ID': 'biocyc',
                        'BiGG ID': 'bigg',
                        'METLIN ID': 'metlin',
                        'PubChem Compound': 'puchem_cid',
                        'PDB ID': 'pdb',
                        'ChEBI ID': 'chebi'}
        url = 'http://www.hmdb.ca/metabolites/' + self.compound_id
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
        SMILES = soup.find('th', sstirng='SMILES')
        if SMILES: SMILES = SMILES.findNext().text
        res['SMILES'] = SMILES
        inchi_key = soup.find('th', string='InChI Key')
        if inchi_key: inchi_key = remove_inchi_key_equal(inchi_key.findNext().text)
        res['inchi_key'] = inchi_key
        synonyms = []
        common_name = soup.find('th', string='Common Name')
        if common_name: common_name = common_name.findNext().text.lower()
        synonyms.insert(0, common_name)
        res['synonyms'] = synonyms
        chemical_formula = soup.find('th', string='Chemical Formula')
        if chemical_formula:
            chemical_formula = chemical_formula.findNext()
            if chemical_formula.text: res['Chemical_formula'] = chemical_formula.text
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
        return compound_instance

if __name__ == '__main__':
    c=Compound_Fetcher_HMDB('HMDB0000538')
    print(c.get_compound())