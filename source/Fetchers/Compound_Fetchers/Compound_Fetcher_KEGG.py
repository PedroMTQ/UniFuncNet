
from source.Fetchers.Compound_Fetchers.Compound_Fetcher import *

class Compound_Fetcher_KEGG(Compound_Fetcher):
    def __init__(self,compound_id,memory_storage=None):
        Compound_Fetcher.__init__( self, compound_id,memory_storage)
        self.db='kegg'
        self.compound=self.get_compound_KEGG()
        self.add_compound()

    #we could also do c->p or c->r but chose not to

    def get_compound_KEGG(self):
        wanted_DBs = {'CAS': 'cas',
                      'PubChem': 'pubchem_sid',
                      'KNApSAcK': 'knapsack',
                      'HMDB': 'hmdb',
                      'ChEBI': 'chebi',
                      'ChEMBL': 'chembl'
                      }
        res = {}
        if self.compound_id[0]=='D':
            url = 'https://www.genome.jp/dbget-bin/www_bget?dr:' + self.compound_id
        elif self.compound_id[0]=='G':
            url = 'https://www.genome.jp/dbget-bin/www_bget?gl:' + self.compound_id
        else:
            url = 'https://www.genome.jp/dbget-bin/www_bget?cpd:' + self.compound_id
        webpage = self.get_with_fetcher(url)
        if not webpage: return None
        soup = BeautifulSoup(webpage, 'lxml')
        if soup.find('pre', text='No such data was found.\n'): return None
        synonyms = self.textbox_KEGG(soup, 'Name')
        if synonyms:
            synonyms = synonyms.split(';')
            synonyms = [i.strip('\n').lower() for i in synonyms]
        else:
            synonyms = []
        formula = self.textbox_KEGG(soup, 'Formula')
        other_dbs = soup.find(string='Other DBs')
        if other_dbs:
            s = other_dbs.findNext().find_all('td', {'valign': "top"})
            for i in s:
                db = i.text.split(':')[0]
                db_id = i.findNext().findNext().text  # sometimes 2 ids appear, the first one is correct, the second is disambiguation
                if db in wanted_DBs.keys():
                    db_possible_ids = db_id.split()
                    if len(db_possible_ids) > 1:
                        res[wanted_DBs[db]] = {possible_id: 1 for possible_id in db_possible_ids}
                    else:
                        res[wanted_DBs[db]] = db_possible_ids[0]
        res['synonyms'] = synonyms
        res['kegg'] = self.compound_id
        res=self.remove_unwanted_info(res)
        if formula: res['Chemical_formula'] = formula
        if number_of_nones_dict(res)==len(res): return None
        compound_instance= Compound(res)
        return compound_instance

if __name__ == '__main__':
    search= Compound_Fetcher_KEGG('C00093')
