#DRAX modules
from source.Fetchers.Fetchers_Utils.Global_Fetcher import *
from source.Utils.util import get_stoichiometry,sub_prod_to_reaction


class Reaction_Fetcher(Global_Fetcher):
    def __init__(self,reaction_id,memory_storage=None):
        Global_Fetcher.__init__(self)
        self.reaction_id=reaction_id
        self.memory_storage=memory_storage
        self.db= None
        self.reaction=None
        self.convergence_args={}
        if not self.memory_storage:
            from source.Pipelines.Searchers.Reaction_Searcher import Reaction_Searcher
            self.memory_storage = Reaction_Searcher()

    def get_reaction(self):
        return self.reaction

    def set_reaction(self,match_instance):
        self.reaction=match_instance

    def reaction_met_instances(self, rn, rn_with_ids, db):
        #this wont run with the fetchers. this was intentionally designed this way
        return self.memory_storage.reaction_met_instances(rn,rn_with_ids,db)

    def add_reaction(self):
        if self.get_reaction():
            match = self.get_reaction_match()
            if match:
                #we would only unite instances when we find a match accross databases, so in truth this rarely happens
                #when we converge though we need to add new info from the convergence args
                if match is not self.get_reaction():
                    match.unite_instances(self.get_reaction(),always_unite=True)
                    self.set_reaction(match)
            else:
                self.memory_storage.add_reaction(self.get_reaction())





    def find_protein(self,query_id=None,extra_args={}):
        memory_type=get_instance_type(self.memory_storage)
        if memory_type=='Protein_Searcher':
            return self.memory_storage.find_protein(db=self.db,query_id=query_id,extra_args=extra_args)
        else:
            return self.memory_storage.protein_searcher.find_protein(db=self.db,query_id=query_id,extra_args=extra_args)

    def find_gene(self,query_id=None,extra_args={}):
        memory_type=get_instance_type(self.memory_storage)
        if memory_type=='Gene_Searcher':
            return self.memory_storage.find_gene(db=self.db,query_id=query_id,extra_args=extra_args)
        else:
            return self.memory_storage.gene_searcher.find_gene(db=self.db,query_id=query_id,extra_args=extra_args)


    def get_reaction_match(self):
        if self.get_reaction():
            return self.memory_storage.get_biological_instance('reactions',self.get_reaction())
        else:
            return self.memory_storage.get_biological_instance('reactions',self.reaction_id,self.db)

if __name__ == '__main__':
    import __main__ as main
    print(main.__file__)