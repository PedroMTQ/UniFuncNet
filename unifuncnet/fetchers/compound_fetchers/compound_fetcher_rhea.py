from unifuncnet.fetchers.compound_fetchers.compound_fetcher import *
from unifuncnet.utils.rhea_sqlite_connector import RheaSqliteConnector


class CompoundFetcherRhea(CompoundFetcher, RheaSqliteConnector):
    def __init__(self, compound_id, memory_storage=None):
        CompoundFetcher.__init__(self, compound_id=compound_id, memory_storage=memory_storage)
        self.db = 'rhea'
        self.set_convergence_args()
        self.compound = self.get_compound_rhea()
        self.add_compound()

    def set_convergence_args(self):
        # args for convergence
        self.convergence_args['reactions'] = set()

    def get_compound_rhea(self):
        compound_instance = Compound({'chebi': self.compound_id})
        self.convergence_args['reactions'] = self.fetch_reactions_rhea_from_chebi(self.compound_id)
        return compound_instance

    def converge_compound_global(self):
        self.converge_compound_to_reaction()

    def converge_compound_to_reaction(self):
        if self.convergence_args['reactions']:
            for reaction_id in self.convergence_args['reactions']:
                print(f'Linking from compound {self.compound_id} in {self.db} to reaction {reaction_id}')
                self.find_reaction(query_id=reaction_id)


if __name__ == '__main__':
    search = CompoundFetcherRhea('7580')
    search.compound.get_all_info()
