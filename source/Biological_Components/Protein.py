

# DRAX modules
from source.Utils.util import is_ec,xstr,list_has_common_items,test_match_possible_ids, unite_possible_ids,score_match_possible_ids
from source.Biological_Components.Biological_Components_Utils.Base_Component import *





class Protein(Base_Component):
    def __init__(self, init_dictionary):
        Base_Component.__init__(self,init_dictionary)
        self.get_ec_from_all_ids()



    def __str__(self):
        try:
            return    '\nProtein ec: '          + xstr(self.get_detail('enzyme_ec')) \
                + '\nProtein Metacyc:'         + xstr(self.get_detail('metacyc')) \
                + '\nProtein Uniprot ID:'         + xstr(self.get_detail('uniprot'))
        except:
            return ''




    def get_ec_from_all_ids(self):
        for detail in self.get_unique_details():
            enz_id= self.get_detail(detail,all_possible=True,ignore_detail='enzyme_ec')
            if enz_id:
                for i in enz_id:
                    if is_ec(i):
                        self.set_detail('enzyme_ec',i)
                        return


    ###MATCHING AND UNITING###



    def score_match_instances(self,instance_2):
        c = 0
        for detail_type in self.get_unique_details(append_to_list=['synonyms']):
            if detail_type == 'synonyms':
                if self.get_detail(detail_type,all_possible=True):
                    if list_has_common_items(self.get_detail(detail_type,all_possible=True), instance_2.get_detail(detail_type,all_possible=True)):
                        c += 0.5
            else:
                c+= score_match_possible_ids(self.get_detail(detail_type,all_possible=True),instance_2.get_detail(detail_type,all_possible=True))
        if test_match_possible_ids(self.get_detail('enzyme_ec',all_possible=True),instance_2.get_detail('enzyme_ec',all_possible=True)):
            if len(self.get_detail('enzyme_ec',all_possible=True))==1 and len(instance_2.get_detail('enzyme_ec',all_possible=True))==1:
                c +=1
        return c



    def unite_instances_bio_specific(self, instance_2):
        for detail_type in self.get_details_list(extra_instances=instance_2):
            if detail_type == 'pathways':               self.set_detail(detail_type,instance_2.get_detail(detail_type,all_possible=True))
            elif detail_type.endswith('_instances'):       self.replace_instances(detail_type,instance_2)
            else:                                       unite_possible_ids(self, instance_2, detail_type)




def test_instance_creator(test_string='test'):
    details=['enzyme_ec','kegg']
    d={}
    p=Protein(d)
    for d in details:
        p.set_detail(d,test_string)
    return p

if __name__ == '__main__':
    test_instance_creator()