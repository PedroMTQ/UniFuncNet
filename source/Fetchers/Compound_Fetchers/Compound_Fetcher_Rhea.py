
from source.Fetchers.Compound_Fetchers.Compound_Fetcher import *
from source.Utils.Rhea_SQLITE_Connector import Rhea_SQLITE_Connector


class Compound_Fetcher_Rhea(Compound_Fetcher,Rhea_SQLITE_Connector):
    def __init__(self,compound_id,memory_storage=None):
        Compound_Fetcher.__init__( self, compound_id,memory_storage)
        Rhea_SQLITE_Connector.__init__(self)
        self.db='rhea'
        self.set_convergence_args()
        self.compound=self.get_compound_rhea()
        self.add_compound()
        self.close_sql_connection()


    def set_convergence_args(self):
        #args for convergence
        self.convergence_args['reactions'] = set()


    def get_compound_rhea(self):
        compound_instance= Compound({'chebi':self.compound_id})
        self.convergence_args['reactions']=self.find_reactions_chebi(self.compound_id)
        return compound_instance


    def converge_compound_global(self):
        self.converge_compound_to_reaction()

    #RP with enzyme ec
    def converge_compound_to_reaction(self):
        if self.convergence_args['reactions']:
            for reaction_id in self.convergence_args['reactions']:
                print(f'Linking from compound {self.compound_id} in {self.db} to reaction {reaction_id}')
                self.find_reaction(query_id=reaction_id)



if __name__ == '__main__':
    search= Compound_Fetcher_Rhea('15377')
