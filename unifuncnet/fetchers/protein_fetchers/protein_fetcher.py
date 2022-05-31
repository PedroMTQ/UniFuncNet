from unifuncnet.fetchers.global_fetcher import *
from unifuncnet.searchers.protein_searcher import ProteinSearcher


class ProteinFetcher(GlobalFetcher):
    def __init__(self, protein_id, memory_storage=None):
        GlobalFetcher.__init__(self, memory_storage=memory_storage)
        self.protein_id = protein_id
        self.db = None
        self.protein = None
        self.convergence_args = {}
        if not self.memory_storage:
            self.memory_storage = ProteinSearcher()

    def add_protein(self):
        if self.get_protein():
            match = self.get_protein_match()
            if match:
                if match is not self.get_protein():
                    # since enzyme ecs are used as identifiers we need to match them
                    match.unite_instances(self.get_protein(), always_unite=True)
                    self.set_protein(match)
            else:
                self.memory_storage.add_protein(self.get_protein())

    def get_protein(self):
        return self.protein

    def set_protein(self, match_instance):
        self.protein = match_instance

    def get_protein_match(self):
        if self.get_protein():
            return self.memory_storage.get_biological_instance('proteins', self.get_protein())
        else:
            return self.memory_storage.get_biological_instance('proteins', self.protein_id, self.db)

    def find_gene(self, query_id=None, extra_args={}):
        memory_type = get_instance_type(self.memory_storage)
        if 'GeneSearcher' in memory_type:
            return self.memory_storage.run_searcher(bio_db=self.db, bio_query=query_id, extra_args=extra_args)
        else:
            return self.memory_storage.gene_searcher.run_searcher(bio_db=self.db, query_id=query_id,
                                                                  extra_args=extra_args)

    def find_protein(self, query_id=None, extra_args={}, convergence_search=False):
        memory_type = get_instance_type(self.memory_storage)
        if 'ProteinSearcher' in memory_type:
            return self.memory_storage.run_searcher(bio_db=self.db, bio_query=query_id, extra_args=extra_args,
                                                    convergence_search=convergence_search)
        else:
            return self.memory_storage.protein_searcher.run_searcher(bio_db=self.db, bio_query=query_id,
                                                                     extra_args=extra_args,
                                                                     convergence_search=convergence_search)

    def find_reaction(self, query_id=None, extra_args={}):
        memory_type = get_instance_type(self.memory_storage)
        if memory_type == 'ReactionSearcher':
            return self.memory_storage.run_searcher(bio_db=self.db, bio_query=query_id, extra_args=extra_args)
        else:
            return self.memory_storage.reaction_searcher.run_searcher(bio_db=self.db, bio_query=query_id,
                                                                      extra_args=extra_args)

    def get_wanted_org_kegg_codes(self):
        memory_type = get_instance_type(self.memory_storage)
        if memory_type == 'ProteinSearcher':
            return self.memory_storage.wanted_org_kegg_codes
        else:
            return self.memory_storage.protein_searcher.wanted_org_kegg_codes
