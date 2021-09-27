# DRAX modules
from source.Utils.util import  unite_possible_ids,\
                            test_match_possible_ids,\
                            score_match_possible_ids
from source.Biological_Components.Biological_Components_Utils.Base_Component import *

# External modules
import re


class Compound(Base_Component):
    def __init__(self, init_dictionary={}):
        Base_Component.__init__(self,init_dictionary)

    def __str__(self):
        res='Compound\n'
        if self.get_detail('synonyms'):
            res+='Name: '+self.get_most_common_synonym()
        return res

    def get_most_common_synonym(self):
        try:    return next(self.get_detail('synonyms'))
        except: return ''

    ###MATCHING AND UNITING###


#Unites two compounds. The main one will be self, the secondary will be provided as an argument.
#Main compound will inherint missing information from the argument compound
    def unite_instances_bio_specific(self, instance_2):

        for detail_type in self.get_details_list(extra_instances=instance_2):
            unite_possible_ids(self, instance_2, detail_type)



    def find_match_synonyms(self, synList):
        if not self.get_detail('synonyms') or not synList:            return False
        if isinstance(synList, str): synList = [synList]
        for syn1 in self.get_detail('synonyms'):
            for syn2 in synList:
                # add plurals also, which can happen with more general compounds
                if syn1.lower() == syn2.lower() or \
                        syn1.lower() == syn2.lower() + 's' or \
                        syn1.lower() + 's' == syn2.lower():
                    return True
        return False




    def score_match_instances(self,instance_2):
        c=0
        for detail_type in  self.get_unique_details(append_to_list=['synonyms']):
            c+=score_match_possible_ids(self.get_detail(detail_type,all_possible=True),instance_2.get_detail(detail_type,all_possible=True))
        return c





    def is_match(self,detail_type,detail_id):
        if not detail_id or not self.get_detail(detail_type): return False
        if self.get_detail(detail_type):
            if test_match_possible_ids(self.get_detail(detail_type),detail_id):
                return True
        else: return False


    def is_empty_metabolite(self):
        c=0
        for d in self.get_details_list():
            if self.get_detail(d): c+=1
        if c==0:    return True
        else:       return False


def test_instance_creator(test_string='test'):
    d= {i:test_string for i in ['kegg','enzyme_ec']}
    p=Compound(d)
    p.get_details_list()
    return p


if __name__ == '__main__':
    cpd=test_instance_creator()
