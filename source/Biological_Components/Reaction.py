
#DRAX modules
from source.Utils.util import       l_rn_ids_to_str,\
                                    unite_possible_ids,\
                                    list_has_common_items,\
                                    score_match_possible_ids
from source.Biological_Components.Biological_Components_Utils.Base_Component import *






class Reaction(Base_Component):
    def __init__(self,init_dictionary):
        Base_Component.__init__(self,init_dictionary)

    def __str__(self):
        if self.get_reaction():
            res='Reaction\n'
            return res+self.get_reaction()
        return ''

    def get_reaction_internal_id(self):
        if not self.reaction_with_instances: return None
        reaction_dict=self.get_detail('reaction_with_instances')
        res=[]
        for side in reaction_dict:
            for stoi,compound in reaction_dict[side]:
                if compound:
                    compound_id = str(compound.internal_id)
                else: compound_id=str(compound)
                res.append([stoi,compound_id])
        return l_rn_ids_to_str(self.get_reaction(),res,without_stoichiometry=True)

    ###MATCHING AND UNITING###




    def is_match_reaction_instances(self, reaction_instance_2):
        rn_inst_1=[]
        rn_inst_1_dict=self.get_detail('reaction_with_instances')
        for rn_key in rn_inst_1_dict:   rn_inst_1.extend(rn_inst_1_dict[rn_key])

        rn_inst_2=[]
        rn_inst_2_dict=reaction_instance_2.get_detail('reaction_with_instances')
        for rn_key in rn_inst_2_dict:   rn_inst_2.extend(rn_inst_2_dict[rn_key])

        if not rn_inst_1 or not rn_inst_2: return False
        if len(rn_inst_1) != len(rn_inst_2): return False
        rn_inst_1_list = []
        rn_inst_2_list = []
        for i in range(len(rn_inst_1)):
            rn_inst_1_list.append(rn_inst_1[i][1])
            rn_inst_2_list.append(rn_inst_2[i][1])
        test=0
        for r1 in rn_inst_1_list:
            #reaction may have Nones, we dont want these to be taken into account
            if not r1: return False
            for r2 in rn_inst_2_list:
                if not r2: return False
                if r1.is_match_instances(r2):
                    test+=1
        if test==len(rn_inst_1): return True
        return False



    def score_match_instances(self,instance_2):
        c = 0
        #reaction with instances and reaction str hold more weight than IDs
        if self.is_match_reaction_instances(instance_2):
            c+=2.25
        if list_has_common_items(self.get_strs_reaction(), instance_2.get_strs_reaction()):
            c+=2
        for unique_detail in self.get_unique_details(remove_from_list=['reaction_with_instances','reaction_str']):
            c+= score_match_possible_ids(self.get_detail(unique_detail), instance_2.get_detail(unique_detail))
        return c




    def unite_instances_bio_specific(self, instance_2):
        for detail_type in self.get_details_list(extra_instances=instance_2):
            if detail_type == 'reaction_str':                  self.set_detail(detail_type,instance_2.get_detail(detail_type,all_possible=True))
            elif detail_type == 'pathways':                     self.set_detail(detail_type,instance_2.get_detail(detail_type,all_possible=True))
            elif detail_type == 'protein_instances':            self.replace_instances(detail_type,instance_2)
            elif detail_type == 'gene_instances':               self.replace_instances(detail_type,instance_2)
            elif detail_type == 'rn_with_ids':
                for i2_rn_with_ids in instance_2.get_rn_with_ids():   self.set_detail(detail_type,i2_rn_with_ids)
            elif detail_type == 'reaction_with_instances':      self.set_detail(detail_type,instance_2.get_detail(detail_type,all_possible=True))
            else:                                               unite_possible_ids(self, instance_2, detail_type)






def test_instance_creator(test_string='test'):
    details=['reaction_str','kegg','rhea']
    d={}
    p=Reaction(d)

    for d in details:
        if d=='reaction_str':
            p.set_detail(d,test_string+' = '+test_string)
        else:
            p.set_detail(d,test_string)
    return p

if __name__ == '__main__':
    test_instance_creator()
