# UniFuncNet modules
from unifuncnet.Utils.util import (unite_possible_ids,
                                   test_match_possible_ids,
                                   score_match_possible_ids)

from unifuncnet.biological_components.utils.base_component import *
from unifuncnet.utils.chebi_sqlite_Connector import ChebiSqliteConnector

# External modules
import re


class Compound(BaseComponent, ChebiSqliteConnector):
    def __init__(self, init_dictionary={}):
        BaseComponent.__init__(self, init_dictionary)
        ChebiSqliteConnector.__init__(self)
        self.add_chebi_ids()
        self.chebi_close_sql_connection()

    def __str__(self):
        res = 'Compound\n'
        if self.get_detail('synonyms'):
            res += 'Name: ' + self.get_most_common_synonym()
        return res

    def get_most_common_synonym(self):
        try:
            return next(self.get_detail('synonyms'))
        except:
            return ''

    def get_top_n_synonyms(self, top_n=5):
        res = []
        try:
            synonyms = self.get_detail('synonyms', all_possible=True)
            for i in range(top_n):
                try:
                    res.append(synonyms.pop(0))
                except:
                    break
            print(res)
            return res
        except:
            return res

    def add_chebi_ids(self):
        for chebi_id in self.get_detail('chebi', all_possible=True):
            main_chebi_id, chebi_mapping = self.fetch_chebi_id_info(chebi_id)
            for db in chebi_mapping:
                self.set_detail(db, chebi_mapping[db], converged_in='chebi')
            self.set_detail('chebi', main_chebi_id, converged_in='chebi')

    ###MATCHING AND UNITING###

    # Unites two compounds. The main one will be self, the secondary will be provided as an argument.
    # Main compound will inherint missing information from the argument compound
    def unite_instances_bio_specific(self, instance_2):
        for detail_type in self.get_details_list(extra_instances=instance_2):
            unite_possible_ids(self, instance_2, detail_type)

    def find_match_synonyms(self, synList):
        if not self.get_detail('synonyms') or not synList:
            return False
        if isinstance(synList, str):
            synList = [synList]
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
            fixed_id = None
            if current_id.startswith('HMDB'):
                zero_search = re.search('HMDB0+', current_id)
                if zero_search:
                    fixed_id = f'HMDB{current_id[zero_search.span()[1]:]}'
            if fixed_id:
                dict_hmdb[fixed_id] = dict_hmdb[current_id]

    def score_match_instances(self, instance_2):
        c = 0
        for detail_type in self.get_unique_details(append_to_list=['synonyms']):
            inst1_ids = self.get_detail(detail_type, all_possible=True)
            inst2_ids = instance_2.get_detail(detail_type, all_possible=True)
            if detail_type == 'hmdb':
                self.add_artificial_hmdb_ids(inst1_ids)
                self.add_artificial_hmdb_ids(inst2_ids)
            score = score_match_possible_ids(inst1_ids, inst2_ids)
            c += score
        return c


def test_instance_creator(test_string='test'):
    d = {i: test_string for i in ['kegg', 'enzyme_ec']}
    p = Compound(d)
    p.get_details_list()
    return p


if __name__ == '__main__':
    # cpd=test_instance_creator()
    cpd1 = Compound({'chebi': {'141550', 'a', 'b'}, 'synonyms': {'a', 'b'}})
    print(cpd1.get_detail('synonyms', all_possible=True))
