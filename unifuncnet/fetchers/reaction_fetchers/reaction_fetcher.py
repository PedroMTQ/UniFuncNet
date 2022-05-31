# UniFuncNet modules
from unifuncnet.fetchers.global_Fetcher import *
from unifuncnet.utils.util import get_stoichiometry, sub_prod_to_reaction
from unifuncnet.searchers.reaction_searcher import ReactionSearcher


class ReactionFetcher(GlobalFetcher):
    def __init__(self, reaction_id, memory_storage=None):
        GlobalFetcher.__init__(self, memory_storage=memory_storage)
        self.reaction_id = reaction_id
        self.db = None
        self.reaction = None
        self.convergence_args = {}
        if not self.memory_storage:
            self.memory_storage = ReactionSearcher()

    def get_reaction(self):
        return self.reaction

    def set_reaction(self, match_instance):
        self.reaction = match_instance

    def reaction_met_instances(self, rn, rn_with_ids, db):
        # this wont run with the fetchers. this was intentionally designed this way
        return self.memory_storage.reaction_met_instances(rn, rn_with_ids, db)

    def reaction_met_instances_simple(self, reaction_stoichiometry, db):
        # this wont run with the fetchers. this was intentionally designed this way
        return self.memory_storage.reaction_met_instances_simple(reaction_stoichiometry, db)

    def add_reaction(self):
        if self.get_reaction():
            match = self.get_reaction_match()
            if match:
                # we would only unite instances when we find a match accross databases, so in truth this rarely happens
                # when we converge though we need to add new info from the convergence args
                if match is not self.get_reaction():
                    match.unite_instances(self.get_reaction(), always_unite=True)
                    self.set_reaction(match)
            else:
                self.memory_storage.add_reaction(self.get_reaction())

    def find_protein(self, query_id=None, extra_args={}, convergence_search=False):
        memory_type = get_instance_type(self.memory_storage)
        if memory_type == 'ProteinSearcher':
            return self.memory_storage.run_searcher(bio_db=self.db, bio_query=query_id, extra_args=extra_args,
                                                    convergence_search=convergence_search)
        else:
            return self.memory_storage.protein_searcher.run_searcher(bio_db=self.db, bio_query=query_id,
                                                                     extra_args=extra_args,
                                                                     convergence_search=convergence_search)

    def find_gene(self, query_id=None, extra_args={}):
        memory_type = get_instance_type(self.memory_storage)
        if memory_type == 'GeneSearcher':
            return self.memory_storage.run_searcher(bio_db=self.db, bio_query=query_id, extra_args=extra_args)
        else:
            return self.memory_storage.gene_searcher.run_searcher(bio_db=self.db, bio_query=query_id,
                                                                  extra_args=extra_args)

    def get_reaction_match(self):
        if self.get_reaction():
            return self.memory_storage.get_biological_instance('reactions', self.get_reaction())
        else:
            return self.memory_storage.get_biological_instance('reactions', self.reaction_id, self.db)


if __name__ == '__main__':
    import __main__ as main

    print(main.__file__)
