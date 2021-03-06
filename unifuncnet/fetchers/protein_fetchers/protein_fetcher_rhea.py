from unifuncnet.fetchers.protein_fetchers.protein_fetcher import *
from unifuncnet.utils.rhea_sqlite_connector import RheaSqliteConnector


class ProteinFetcherRhea(ProteinFetcher, RheaSqliteConnector):
    def __init__(self, protein_id, extra_args={}, memory_storage=None):
        ProteinFetcher.__init__(self, protein_id=protein_id, memory_storage=memory_storage)
        self.db = 'rhea'
        self.protein = self.get_protein_rhea()
        self.add_protein()

    def get_protein_rhea(self):
        if not self.protein_id: return None
        if is_ec(self.protein_id):
            id_type = 'enzyme_ec'
        else:
            id_type = 'uniprot'

        self.convergence_args['reactions_list'] = self.fetch_rhea_from_id(id_type, self.protein_id)

        res = {id_type: self.protein_id, 'rhea': self.protein_id,
               }
        return Protein(res)

    def converge_protein_global(self):
        self.converge_protein_to_reaction()

    def converge_protein_gpr(self):
        self.converge_protein_to_reaction()

    def converge_protein_rpg(self):
        pass

    def converge_protein_to_reaction(self):
        for reaction_id in self.convergence_args['reactions_list']:
            print(f'Linking from protein {self.protein_id} in {self.db} to reaction {reaction_id}')

            reaction_instance = self.find_reaction(query_id=reaction_id)

            if reaction_instance:
                self.get_protein().set_detail('reaction_instances', reaction_instance, converged_in=self.db)


if __name__ == '__main__':
    from unifuncnet.biological_components.Protein import Protein
    from unifuncnet.fetchers.Reaction_Fetchers.Reaction_Fetcher import *
    import re

    # p0=Protein_Fetcher_Rhea('2.3.1.15')
    p0 = ProteinFetcherRhea('2.3.1.15')
    print(p0.get_protein().get_details_list())
    # p0.converge_protein_gpr()
    # p1=Protein_Fetcher_KEGG('1.14.14.73').get_protein()
