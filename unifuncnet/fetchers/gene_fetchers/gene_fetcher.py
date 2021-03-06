from unifuncnet.fetchers.global_fetcher import *
from unifuncnet.searchers.gene_searcher import GeneSearcher


class GeneFetcher(GlobalFetcher):
    def __init__(self, gene_id, extra_args={}, memory_storage=None):
        GlobalFetcher.__init__(self, memory_storage=memory_storage)
        self.gene_id = gene_id
        self.db = None
        self.gene = None
        self.convergence_args = {}
        # if no memory_storage is present from one of the pipelines or previous fetchers we assign it one and initialize the memory
        if not self.memory_storage:
            self.memory_storage = GeneSearcher()

    def add_gene(self):
        if self.get_gene():
            match = self.get_gene_match()
            if match:
                if match is not self.get_gene():
                    match.unite_instances(self.get_gene(), always_unite=True)
                    self.set_gene(match)
            else:
                self.memory_storage.add_gene(self.get_gene())

    def get_gene(self):
        return self.gene

    def set_gene(self, match_instance):
        self.gene = match_instance

    def get_gene_match(self):
        if self.get_gene():
            return self.memory_storage.get_biological_instance('genes', self.get_gene())
        else:
            return self.memory_storage.get_biological_instance('genes', self.gene_id, self.db)

    def find_protein(self, query_id=None, extra_args={}, convergence_search=False):
        memory_type = get_instance_type(self.memory_storage)
        if memory_type == 'ProteinSearcher':
            return self.memory_storage.run_searcher(bio_db=self.db, bio_query=query_id, extra_args=extra_args,
                                                    convergence_search=convergence_search)
        else:
            return self.memory_storage.protein_searcher.run_searcher(bio_db=self.db, bio_query=query_id,
                                                                     extra_args=extra_args,
                                                                     convergence_search=convergence_search)

    def find_reaction(self, query_id=None, extra_args={}):
        memory_type = get_instance_type(self.memory_storage)
        if memory_type == 'Reaction_Searcher':
            return self.memory_storage.run_searcher(bio_db=self.db, bio_query=query_id, extra_args=extra_args)
        else:
            return self.memory_storage.reaction_searcher.run_searcher(bio_db=self.db, bio_query=query_id,
                                                                      extra_args=extra_args)
