# DRAX modules
from source.Utils.util import  unite_possible_ids,\
                            test_match_possible_ids,\
                            score_match_possible_ids

from source.Biological_Components.Biological_Components_Utils.Base_Component import *
from source.Utils.CHEBI_SQLITE_Connector import CHEBI_SQLITE_Connector

# External modules
import re


class Compound(Base_Component,CHEBI_SQLITE_Connector):
    def __init__(self, init_dictionary={}):
        Base_Component.__init__(self,init_dictionary)
        CHEBI_SQLITE_Connector.__init__(self)
        self.add_chebi_ids()
        self.close_sql_connection()


    def __str__(self):
        res='Compound\n'
        if self.get_detail('synonyms'):
            res+='Name: '+self.get_most_common_synonym()
        return res

    def get_most_common_synonym(self):
        try:    return next(self.get_detail('synonyms'))
        except: return ''

    def add_chebi_ids(self):
        for chebi_id in self.get_detail('chebi',all_possible=True):
            chebi_mapping=self.fetch_chebi_id_info(chebi_id)
            for db in chebi_mapping:
                self.set_detail(db,chebi_mapping[db],converged_in='chebi')


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

    def add_artificial_hmdb_ids(self, dict_hmdb):
        # hmdb sometimes has ids with 0s in between, sometimes it doesnt. this can lead to discrepancies
        # in total there should be 4 letter- HMDB followed
        for current_id in dict(dict_hmdb):
            fixed_id=None
            if current_id.startswith('HMDB'):
                zero_search = re.search('HMDB0+', current_id)
                if zero_search:
                    fixed_id = f'HMDB{current_id[zero_search.span()[1]:]}'
            if fixed_id:
                dict_hmdb[fixed_id]=dict_hmdb[current_id]



    def score_match_instances(self,instance_2):
        c=0
        for detail_type in  self.get_unique_details(append_to_list=['synonyms']):
            inst1_ids=self.get_detail(detail_type,all_possible=True)
            inst2_ids=instance_2.get_detail(detail_type,all_possible=True)
            if detail_type=='hmdb':
                self.add_artificial_hmdb_ids(inst1_ids)
                self.add_artificial_hmdb_ids(inst2_ids)
            score=score_match_possible_ids(inst1_ids,inst2_ids)
            c+=score
        return c





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
    #cpd=test_instance_creator()
    cpd1=Compound({'chebi':'17968'})
    cpd1.get_all_info()