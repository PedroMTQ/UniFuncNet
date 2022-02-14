
from drax.Fetchers.Reaction_Fetchers.Reaction_Fetcher import *
from drax.Utils.Rhea_SQLITE_Connector import Rhea_SQLITE_Connector



class Reaction_Fetcher_Rhea(Reaction_Fetcher,Rhea_SQLITE_Connector):
    def __init__(self,reaction_id,extra_args={},memory_storage=None):
        Reaction_Fetcher.__init__(self,reaction_id=reaction_id,memory_storage=memory_storage)
        self.db = 'rhea'
        self.set_convergence_args(extra_args)
        self.reaction= self.get_reaction_rhea()
        self.add_reaction()


    def get_reaction_rhea(self):  # this basically gets all the details from rhea, may be able to implement a more concise search, see later
        if not self.reaction_id: return None
        reaction_info=self.fetch_rhea_id_info(self.reaction_id)
        if reaction_info:
            rhea_ids=reaction_info['alt_ids']
            rhea_ids.append(self.reaction_id)
            kegg_ids=reaction_info['kegg']
            metacyc_ids=reaction_info['metacyc']
            uniprot_ids=reaction_info['uniprot']
            enzyme_ec_ids=reaction_info['enzyme_ec']
            reaction_str=reaction_info['reaction_str']
            rn_with_ids=reaction_info['chebi_equation']
            self.convergence_args['uniprot_ids']=uniprot_ids
            self.convergence_args['enzyme_ec_ids']=enzyme_ec_ids


            try:
                rn_with_ids, complete_l, len_sub = get_stoichiometry(reaction_str, rn_with_ids)
            except:
                try:
                    rn_with_ids, complete_l, len_sub  = get_stoichiometry(reaction_str, reaction_str)
                except:
                    return None
            rn_with_instances = self.reaction_met_instances(reaction_str, rn_with_ids, 'chebi')
            reaction_dict = {
                'rhea': rhea_ids,
                'reaction_str': reaction_str,
                'kegg':kegg_ids,
                'metacyc':metacyc_ids,
            }
            if rn_with_instances:  reaction_dict['reaction_with_instances']= rn_with_instances
            else:
                if rn_with_ids: reaction_dict['rn_with_ids']= [reaction_str, rn_with_ids, 'rhea']
            reaction_instance = Reaction(reaction_dict)

            return reaction_instance


    def set_convergence_args(self,extra_args):
        #args for convergence
        if 'enzyme_ecs' not in self.convergence_args: self.convergence_args['enzyme_ec_ids'] = []
        if 'uniprot_ids' not in self.convergence_args: self.convergence_args['uniprot_ids'] = []



    def converge_reaction_rpg(self):
        #RP and RG part of the CRPG pipeline
        self.converge_reaction_to_protein()

    #RP with enzyme ec
    def converge_reaction_to_protein(self):
        for enzyme_ec in self.convergence_args['enzyme_ec_ids']:
            print(f'Linking from reaction {self.reaction_id} in {self.db} to protein {enzyme_ec}')
            protein_instance = self.find_protein(query_id=enzyme_ec)
            if protein_instance:
                self.get_reaction().set_detail('protein_instances',protein_instance,converged_in=self.db)

        for uniprot_id in self.convergence_args['uniprot_ids']:
            print(f'Linking from reaction {self.reaction_id} in {self.db} to protein {uniprot_id}')
            protein_instance = self.find_protein(query_id=uniprot_id)
            if protein_instance:
                self.get_reaction().set_detail('protein_instances',protein_instance,converged_in=self.db)







if __name__ == '__main__':
    rn_search=Reaction_Fetcher_Rhea('10000')


