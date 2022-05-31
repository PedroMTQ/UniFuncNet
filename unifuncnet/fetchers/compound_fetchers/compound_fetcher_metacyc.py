from unifuncnet.fetchers.compound_fetchers.compound_fetcher import *


class CompoundFetcherMetacyc(CompoundFetcher):
    def __init__(self, compound_id, memory_storage=None):
        CompoundFetcher.__init__(self, compound_id=compound_id, memory_storage=memory_storage)
        self.db = 'metacyc'
        self.set_convergence_args()
        self.compound = self.get_compound_metacyc()
        self.add_compound()

    def set_convergence_args(self):
        # args for convergence
        self.convergence_args['reactions'] = set()

    def get_compound_metacyc(self):
        if not self.compound_id: return None

        fetched_info = self.fetch_metacyc_id_info(self.compound_id, 'compound')
        self.convergence_args['reactions'] = self.fetch_metacyc_rxn_from_cpd(self.compound_id)
        fetched_info['metacyc'] = self.compound_id
        return Compound(fetched_info)

    def converge_compound_global(self):
        self.converge_compound_to_reaction()

    # RP with enzyme ec
    def converge_compound_to_reaction(self):
        for reaction_id in self.convergence_args['reactions']:
            print(f'Linking from compound {self.compound_id} in {self.db} to reaction {reaction_id}')
            self.find_reaction(query_id=reaction_id)


if __name__ == '__main__':
    search = CompoundFetcherMetacyc('OXYGEN-MOLECULE')
    search.get_compound().get_all_info()
