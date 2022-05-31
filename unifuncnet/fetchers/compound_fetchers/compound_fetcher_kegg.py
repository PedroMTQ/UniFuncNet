from unifuncnet.fetchers.compound_fetchers.compound_fetcher import *


class CompoundFetcherKegg(CompoundFetcher):
    def __init__(self, compound_id, memory_storage=None):
        CompoundFetcher.__init__(self, compound_id, memory_storage)
        self.db = 'kegg'
        self.set_convergence_args()
        self.compound = self.get_compound_kegg()
        self.add_compound()

    def set_convergence_args(self):
        # args for convergence
        self.convergence_args['reactions'] = set()
        self.convergence_args['soup'] = None

    def get_reactions(self, soup):
        res = self.textbox_KEGG(soup, 'Reaction')
        res = res.split()
        res = [i.strip('\n') for i in res]
        res = set(res)
        return res

    def get_compound_kegg(self):
        wanted_dbs = {'CAS': 'cas',
                      'PubChem': 'pubchem_sid',
                      'KNApSAcK': 'knapsack',
                      'HMDB': 'hmdb',
                      'ChEBI': 'chebi',
                      'ChEMBL': 'chembl'
                      }
        res = {}
        if self.compound_id[0] == 'D':
            url = f'https://www.genome.jp/dbget-bin/www_bget?dr:{self.compound_id}'
        elif self.compound_id[0] == 'G':
            url = f'https://www.genome.jp/dbget-bin/www_bget?gl:{self.compound_id}'
        else:
            url = f'https://www.genome.jp/dbget-bin/www_bget?cpd:{self.compound_id}'
        webpage = self.get_with_fetcher(url)
        if not webpage:
            return None
        soup = BeautifulSoup(webpage, 'lxml')
        if soup.find('pre', text='No such data was found.\n'):
            return None
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
                if db in wanted_dbs.keys():
                    db_possible_ids = db_id.split()
                    if len(db_possible_ids) > 1:
                        res[wanted_dbs[db]] = {possible_id: 1 for possible_id in db_possible_ids}
                    else:
                        res[wanted_dbs[db]] = db_possible_ids[0]
        res['synonyms'] = synonyms
        res['kegg'] = self.compound_id
        res = self.remove_unwanted_info(res)
        if formula: res['chemical_formula'] = formula
        if number_of_nones_dict(res) == len(res):
            return None
        compound_instance = Compound(res)
        self.convergence_args['soup'] = soup
        return compound_instance

    def converge_compound_global(self):
        self.converge_compound_to_reaction()

    # RP with enzyme ec
    def converge_compound_to_reaction(self):
        self.convergence_args['reactions'] = self.get_reactions(self.convergence_args['soup'])
        if self.convergence_args['reactions']:
            for reaction_id in self.convergence_args['reactions']:
                print(f'Linking from compound {self.compound_id} in {self.db} to reaction {reaction_id}')
                self.find_reaction(query_id=reaction_id)
        self.convergence_args['soup'] = None


if __name__ == '__main__':
    search = CompoundFetcherKegg('C00093')
