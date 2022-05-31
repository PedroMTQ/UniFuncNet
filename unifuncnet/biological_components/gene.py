# UniFuncNet modules
from unifuncnet.Utils.util import (list_has_common_items,
                                   unite_possible_ids,
                                   test_match_possible_ids,
                                   score_match_possible_ids,
                                   get_instance_type)
from unifuncnet.biological_components.utils.base_component import *


class Gene(BaseComponent):
    def __init__(self, init_dictionary={}):
        BaseComponent.__init__(self, init_dictionary)

    def __str__(self):
        res = ''
        if self.get_detail('synonyms'):
            res += 'Gene name:\t' + self.get_most_common_synonym()

        if res:
            return f'Gene\n{res}'
        return res

    ###MATCHING AND UNITING###

    def score_match_instances(self, instance_2):
        c = 0
        for detail_type in self.get_unique_details():
            if detail_type == 'uniprot':
                match_c = 2
            else:
                match_c = 0.5
            c += score_match_possible_ids(self.get_detail(detail_type,
                                                          all_possible=True),
                                          instance_2.get_detail(detail_type,
                                                                all_possible=True),
                                          match_c=match_c)
        return c

    def unite_instances_bio_specific(self, instance_2):
        for detail_type in self.get_details_list(extra_instances=instance_2):
            if detail_type.endswith('_instances'):
                self.replace_instances(detail_type, instance_2)
            else:
                unite_possible_ids(self, instance_2, detail_type)


def test_instance_creator(test_string='test'):
    d = {}
    p = Gene(d)
    for d in ['kegg', 'uniprot']:
        p.set_detail(d, test_string)
    return p


if __name__ == '__main__':
    # test_instance_creator()
    g = Gene()
    g.set_detail('nt_seq', 'aaa')
    g.set_detail('nt_seq', 'bbb')
    a = test_match_possible_ids(g.get_detail('nt_seq', all_possible=True), g.get_detail('nt_seq', all_possible=True))
    test = test_instance_creator()
