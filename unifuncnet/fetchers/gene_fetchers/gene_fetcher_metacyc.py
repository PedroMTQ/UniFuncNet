from unifuncnet.fetchers.gene_fetchers.gene_fetcher import *


class GeneFetcherMetacyc(GeneFetcher):
    def __init__(self, gene_id=None, extra_args={}, memory_storage=None):
        GeneFetcher.__init__(self, gene_id=gene_id, extra_args=extra_args, memory_storage=memory_storage)
        self.db = 'metacyc'
        self.set_convergence_args()
        self.gene = self.get_gene_metacyc()
        self.add_gene()

    def get_gene_metacyc(self):
        if not self.gene_id: return None

        fetched_info = self.fetch_metacyc_id_info(self.gene_id, 'gene')
        if 'protein_ids' in fetched_info:
            self.convergence_args['proteins'] = fetched_info.pop('protein_ids')
        fetched_info['metacyc'] = self.gene_id
        return Gene(fetched_info)

    def set_convergence_args(self):
        self.convergence_args['proteins'] = set()

    def converge_gene_gpr(self):
        if self.get_gene():
            self.converge_gene_to_protein()

    # GP
    def converge_gene_to_protein(self):
        for protein_id in self.convergence_args['proteins']:
            print(f'Linking from gene {self.gene_id} in {self.db} to protein {protein_id}')
            protein_instance = self.find_protein(query_id=protein_id)
            if protein_instance:
                self.get_gene().set_detail('protein_instances', protein_instance, converged_in=self.db)


if __name__ == '__main__':
    from unifuncnet.biological_components.Protein import Protein
    from unifuncnet.fetchers.Reaction_Fetchers.Reaction_Fetcher import *

    gene = GeneFetcherMetacyc('G1G01-1865')
    gene.get_gene().get_all_info()
