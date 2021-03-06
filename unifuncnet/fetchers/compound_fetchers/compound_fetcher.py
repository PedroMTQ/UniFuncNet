from unifuncnet.fetchers.global_fetcher import *
from types import GeneratorType as generator
from unifuncnet.utils.util import remove_inchi_key_equal
from unifuncnet.Searchers.compound_searcher import CompoundSearcher


class CompoundFetcher(GlobalFetcher):
    def __init__(self, compound_id, memory_storage=None):
        GlobalFetcher.__init__(self, memory_storage=memory_storage)
        self.compound_id = compound_id
        self.db = None
        self.compound = None
        self.convergence_args = {}
        # if no memory_storage is present from one of the pipelines or previous fetchers we assign it one and initialize the memory
        if not self.memory_storage:
            self.memory_storage = CompoundSearcher()

    def add_compound(self):
        if self.get_compound():
            match = self.get_compound_match()
            if match:
                if match is not self.get_compound():
                    match.unite_instances(self.get_compound(), always_unite=True)
                    self.set_compound(match)
            else:
                self.memory_storage.add_compound(self.get_compound())

    def get_compound(self):
        return self.compound

    def set_compound(self, match_instance):
        self.compound = match_instance

    def get_compound_match(self):
        if self.get_compound():
            return self.memory_storage.get_biological_instance('compounds', self.get_compound())
        else:
            return self.memory_storage.get_biological_instance('compounds', self.compound_id, self.db)

    def replace_reaction_met_instances(self, new_instance, old_instance):
        self.memory_storage.replace_reaction_met_instances(new_instance, old_instance)

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
        if memory_type == 'ReactionSearcher':
            return self.memory_storage.run_searcher(bio_db=self.db, bio_query=query_id, extra_args=extra_args)
        else:
            return self.memory_storage.reaction_searcher.run_searcher(bio_db=self.db, bio_query=query_id,
                                                                      extra_args=extra_args)

    def remove_unwanted_info(self, dict_to_change):
        res = {}
        for detail in dict_to_change:
            if dict_to_change[detail]:
                res[detail] = dict_to_change[detail]
        for detail in res:
            if isinstance(res[detail], str):
                res[detail] = {res[detail]}
        for detail in res:
            without_nones = {i for i in res[detail] if i}
            res[detail] = {i for i in without_nones if not re.search('not available', i.lower())}
        return res
