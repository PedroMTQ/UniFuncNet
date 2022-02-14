
from drax.Fetchers.Reaction_Fetchers.Reaction_Fetcher import *

class Reaction_Fetcher_Metacyc(Reaction_Fetcher):
    def __init__(self,reaction_id,extra_args={},memory_storage=None):
        Reaction_Fetcher.__init__(self,reaction_id=reaction_id,memory_storage=memory_storage)
        self.db = 'metacyc'
        self.set_convergence_args()
        self.reaction= self.get_reaction_metacyc()
        self.add_reaction()


    def generate_reaction_str(self,reaction_stoichiometry):
        reaction_str=[]
        for i in range(len(reaction_stoichiometry)):
            if isinstance(reaction_stoichiometry[i],list):
                stoi,cpd_id=reaction_stoichiometry[i]
                reaction_str.append(f'{stoi} {cpd_id} +')
            else:
                reaction_str[i-1]=reaction_str[i-1].strip('+')
                reaction_str.append(reaction_stoichiometry[i])
        reaction_str[-1] = reaction_str[-1].strip('+')
        reaction_str=' '.join(reaction_str)
        reaction_str=reaction_str.strip()
        return reaction_str

    def get_reaction_metacyc(self):
        if not self.reaction_id: return None
        reaction_info=self.fetch_metacyc_id_info(self.reaction_id,'reaction')
        if 'enzyme_ec' in reaction_info:
            self.convergence_args['enzyme_ecs'].update(reaction_info.pop('enzyme_ec'))
        intermediate_reaction_ids=set()
        if 'intermediate_reaction_ids' in reaction_info:
            intermediate_reaction_ids.update(reaction_info.pop('intermediate_reaction_ids'))
        for int_reaction_id in intermediate_reaction_ids:
            fetched_info = self.fetch_metacyc_intermediate_rxn_ids(int_reaction_id)
            if 'protein_ids' in fetched_info:
                self.convergence_args['proteins'].update(fetched_info['protein_ids'])
        reaction_info['metacyc']=self.reaction_id
        if 'reaction_stoichiometry' in reaction_info:
            reaction_stoichiometry=reaction_info.pop('reaction_stoichiometry')
            reaction_info['reaction_with_instances'] = self.reaction_met_instances_simple(reaction_stoichiometry, 'metacyc')
            reaction_info['reaction_str'] = self.generate_reaction_str(reaction_stoichiometry)
        return Reaction(reaction_info)


    def set_convergence_args(self):
        self.convergence_args['enzyme_ecs']=set()
        self.convergence_args['proteins']=set()

    def converge_reaction_rpg(self):
        #RP and RG part of the CRPG pipeline
        self.converge_reaction_to_protein()

    def converge_reaction_to_protein(self):
        for protein_id in self.convergence_args['proteins']:
            print(f'Linking from gene {self.reaction_id} in {self.db} to protein {protein_id}')
            protein_instance = self.find_protein(query_id=protein_id)
            if protein_instance:
                self.get_reaction().set_detail('protein_instances',protein_instance,converged_in=self.db)





if __name__ == '__main__':
    #rn_search=Reaction_Fetcher_Metacyc('RXN-2043')
    rn_search=Reaction_Fetcher_Metacyc('RXN-14064')
    rn_search.reaction.get_all_info()


