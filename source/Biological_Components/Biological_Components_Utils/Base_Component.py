from types import GeneratorType as generator
from source.Biological_Components.Info_Keepers.Counter import Counter
from source.Biological_Components.Info_Keepers.Item_Set import Item_Set
from source.Biological_Components.Info_Keepers.Synonyms import Synonyms
from source.Biological_Components.Biological_Components_Utils.Unique_details import get_details_list,generate_details_dict,get_unique_details
from source.Utils.util import       find_sign, \
                                    uniform_sign, \
                                    get_instance_type,\
                                    fix_html_sign,SCRAPPABLE_DBS,is_ec
from math import ceil


class Base_Component():
    def __init__(self,init_dictionary):
        self.identifiers={}
        self.instances={}
        for d_key in init_dictionary.keys():
            self.set_detail(d_key,init_dictionary[d_key])


    def __hash__(self): return hash(id(self))
    def __eq__(self, x): return x is self
    def __ne__(self, x): return x is not self


    def get_unique_details(self,remove_from_list=[],append_to_list=[]):
        return get_unique_details(bio_instance=self,remove_from_list=remove_from_list,append_to_list=append_to_list)

    def get_details_list(self,extra_instances=[],remove_from_list=[],append_to_list=[]):
        if not isinstance(extra_instances, list):
            extra_instances=[extra_instances]
        extra_instances.append(self)
        return get_details_list(bio_instances=extra_instances,remove_from_list=remove_from_list,append_to_list=append_to_list)

    def generate_details_dict(self):
        return generate_details_dict(bio_instance=self)

    def get_most_common_synonym(self):
        return self.synonyms.get_most_common_synonym()

    def saved(self):
        self.need_to_save = False


    def needs_saving(self):
        if self.need_to_save: return True
        return False



    def get_all_info(self):
        for d_key in self.get_details_list():
            if self.get_detail(d_key, all_possible=True):
                print(d_key, self.get_detail(d_key, all_possible=True))

    def export_all_info(self):
        line=[]
        for d_key in self.get_details_list(remove_from_list=[ 'internal_id',
                                                             'compound_instances',
                                                             ]):
            if d_key=='reaction_with_instances':
                reaction_compounds=self.get_reaction_internal_id()
                line.append(f'reaction_compounds:{reaction_compounds}')
            elif d_key in ['protein_instances',
                           'reaction_instances',
                           'gene_instances',]:
                details=self.get_detail(d_key, all_possible=True)
                if details:
                    key_to_write=d_key.replace('_instances','s')+'_connected'
                    for d in details:
                        line.append(f'{key_to_write}:{d.internal_id}')

            else:
                details=self.get_detail(d_key, all_possible=True)
                if details:
                    for d in details:
                        if is_ec(d,4) and d_key!='enzyme_ec':pass
                        else:
                            line.append(f'{d_key}:{d}')
        return '\t'.join(sorted(line))

    def set_detail(self,detail_type,detail,converged_in=None):
        if detail is not None:
            if detail_type == 'reaction_with_instances':       self.set_reaction_with_instances(detail)
            elif detail_type == 'rn_with_ids':                 self.set_rn_with_ids(detail)
            else:
                res = {}

                if not isinstance(detail, dict):
                    if isinstance(detail, list) or isinstance(detail,generator) or isinstance(detail,set):
                        for i in detail: res[i] = 1
                    else:
                        res[detail] = 1
                else:
                    res = dict(detail)
                for to_add in res:
                    count = res[to_add]
                    if detail_type=='synonyms':                       self.set_synonyms(detail)
                    elif detail_type.endswith('_instances'):            self.set_instances(instances_type=detail_type,instances_to_add=to_add,converged_in=converged_in)
                    elif detail_type == 'reaction_str':                 self.set_reaction(to_add)
                    elif detail_type == 'pathways':                     self.set_pathways(to_add)
                    elif detail_type=='enzyme_ec':
                        self.set_id(db=detail_type,id_to_add=to_add,count= count)
                        if get_instance_type(self)=='Protein':
                            # these DBs also use enzyme ec as ids
                            self.set_id(db='biocyc',id_to_add=to_add,count= 1)
                            self.set_id(db='kegg',id_to_add=to_add,count= 1)
                    else:                                               self.set_id(db=detail_type, id_to_add=to_add, count=count)


    def get_detail(self, detail_type,all_possible=False,ignore_detail=None,return_instance=False,return_convergence=False):
        if detail_type == ignore_detail:                return None
        if detail_type == 'reaction_with_instances':  return self.get_reaction_with_instances()
        elif detail_type == 'rn_with_ids':              return self.get_rn_with_ids()
        elif detail_type == 'synonyms':
            if not all_possible:                        return self.get_synonyms()
            else:                                       return self.get_possible_synonyms()
        elif detail_type == 'reaction_str':
            if all_possible:                            return self.get_possible_str_reaction()
            else:                                       return self.get_reaction()
        elif detail_type == 'pathways':                 return self.get_pathways()
        elif detail_type.endswith('_instances'):        return self.get_instances(instances_type=detail_type, return_instance=return_instance,return_convergence=return_convergence)
        else:
            if all_possible:                            return self.get_possible_ids(detail_type)
            else:                                       return self.get_id(detail_type)

    # Deleting information
    def remove_detail(self,detail_type,detail_to_remove):
        if detail_type == 'rn_with_ids':                   self.remove_rn_with_ids(detail_to_remove)
        elif detail_type == 'synonyms':                    self.remove_synonym(detail_to_remove)
        elif detail_type == 'reaction_str':                self.remove_reaction(detail_to_remove)
        elif detail_type == 'pathways':                    self.remove_pathways(detail_to_remove)
        elif detail_type.endswith('instances'):            return self.remove_detail_type_instances(self,detail_type,detail_to_remove)
        else:                                              self.remove_id(db=detail_type, id_to_remove=detail_to_remove)

    def get_id(self, db):

        if not db.endswith('_id'): db_str=db+'_id'
        else: db_str=db
        if db_str in self.identifiers:
            return self.identifiers[db_str].get_most_common_string()
        else:   return []
        #kegg= this id correponds to the bio componenets id in kegg
        #kegg_<others> = these others are non  bio componenents ids in kegg

    def get_IDs(self, db):
        if not db.endswith('_id'):          db_str = db + '_id'
        else:                               db_str = db
        if db_str in self.identifiers:
            return self.identifiers[db_str].get_strings()
        else:
            return []

    # For a dictionary {ID:counter} of all possible IDs
    def get_possible_ids(self, db):
        if not db.endswith('_id'):          db_str = db + '_id'
        else:                               db_str = db
        if db_str in self.identifiers:
            return self.identifiers[db_str].get_possible_strings()
        else:
            return []

    def remove_id(self, db, id_to_remove):
        if not db.endswith('_id'):          db_str = db + '_id'
        else:                               db_str = db
        if db_str in self.identifiers:
            self.identifiers[db_str].remove_string(id_to_remove)


# Setting information
# Databases IDs are actually instances of class Counter, which will define a compound's ID as the most frequent one from  a list of redundant IDs
    def set_id(self, db, id_to_add, count=1):
        if not db.endswith('_id'):          db_str = db + '_id'
        else:                               db_str = db
        if db_str not in self.identifiers:
            self.identifiers[db_str]=Counter(id_to_add)
        else:
            self.identifiers[db_str].append(id_to_add,count)



#### GET MAIN

    def get_reaction(self):
        if hasattr(self,'reaction_str'):
            return self.reaction_str.get_most_common_string()
        else:
            return []



    def get_rn_with_ids(self):
        if hasattr(self,'rn_with_ids'):
            return self.rn_with_ids
        else:
            return []

    def get_reaction_with_instances(self):
        if hasattr(self,'reaction_with_instances'):
            if self.reaction_with_instances:
                return self.reaction_with_instances
            else: return []
        else:
            return []

    #### GET LIST

    def get_instances(self,instances_type,return_instance=False,return_convergence=False):
        if instances_type in self.instances:
            if return_convergence:
                return self.instances[instances_type].convergence
            elif return_instance:
                return self.instances[instances_type]
            else:
                return self.instances[instances_type].get_item_set()
        else:
            return set()

    def get_pathways(self,return_instance=False):
        if hasattr(self,'pathways'):
            if return_instance:
                return self.pathways
            else:
                return self.pathways.get_item_set()
        else:   return set()

    def get_strs_reaction(self,return_instance=False):
        if hasattr(self,'reaction_str'):
            if return_instance:
                return self.reaction_str
            else:
                return self.reaction_str.get_strings()
        else:   return set()

    def get_synonyms(self,return_instance=False):
        if hasattr(self,'synonyms'):
            if return_instance:
                return self.synonyms
            else:
                return self.synonyms.get_synonyms()
        else:   return set()

#### GET POSSIBLE ENTRIES

    def get_possible_synonyms(self):
        if hasattr(self,'synonyms'):
            return self.synonyms.get_possible_synonyms()
        else:   return set()



    def get_possible_str_reaction(self):
        if hasattr(self,'reaction_str'):
            return self.reaction_str.get_possible_strings()
        else:   return set()

#### DELETING

    def remove_synonym(self, syn_to_remove):
        self.synonyms.remove_synonym(syn_to_remove)

    def remove_reaction(self, detail_to_remove):
        self.reaction_str.remove_string(detail_to_remove)

    def remove_references(self, detail_to_remove):
        self.references.remove_item(detail_to_remove)

    def remove_pathways(self, detail_to_remove):
        self.pathways.remove_item(detail_to_remove)

    def remove_rn_with_ids(self, detail_to_remove):
        self.rn_with_ids.remove(detail_to_remove)

    def clear_rn_with_ids(self):
        self.rn_with_ids = []

##### SETTING

    def set_synonyms(self, syns,count=1):
        if hasattr(self, 'synonyms'):
            self.synonyms.append(syns,count)
        else:
            self.synonyms = Synonyms(syns,bio_type=self.__class__.__name__)


#Databases IDs are actually instances of class Counter, which will define a compound's ID as the most frequent one from  a list of redundant IDs

    def set_reaction(self,reaction):
        #sometimes reactions have extra blank spaces
        if reaction:
            reaction=reaction.replace('  ',' ').strip()
            sign = find_sign(reaction)
            if sign:
                sign=sign.group()
                reaction = reaction.replace(sign,uniform_sign(sign))
            reaction = fix_html_sign(reaction)

        if hasattr(self, 'reaction_str'):
            self.reaction_str.append(reaction)
        else:
            self.reaction_str = Counter(reaction)

    def set_reaction_with_instances(self,reaction_with_instances):
        if hasattr(self,'reaction_with_instances'):
            if self.reaction_with_instances:
                if self.reaction_with_instances!=reaction_with_instances:

                    #print('reaction with instances already set!',self.get_detail('rhea'))
                    #print('reaction with instances already set!',self.get_reaction(),self.reaction_with_instances)
                    #print('trying to change to :',reaction_with_instances)
                    self.reaction_with_instances= reaction_with_instances
        else:
            self.reaction_with_instances= reaction_with_instances

    def set_rn_with_ids(self,rn_with_ids):
        if hasattr(self,'rn_with_ids'):
            self.rn_with_ids.append(rn_with_ids)
        else:
            self.rn_with_ids = []



    def set_instances(self,instances_type,instances_to_add,converged_in=None):
        if instances_type in self.instances:
            self.instances[instances_type].set_item_set(instances_to_add,converged_in)
        else:
            self.instances[instances_type] = Item_Set(instances_to_add,converged_in)

    def set_pathways(self, pathways):
        if hasattr(self,'pathways'):
            self.pathways.set_item_set(pathways)
        else:
            self.pathways = Item_Set(pathways)

    def get_convergence_db(self,detail,current_detail):
        return self.get_detail(detail,return_instance=True).get_convergence_db(current_detail)

    def check_if_same_instances(self,instance_2):
        if self is instance_2:
            return True
        return False

    #when the instances have too little info, they wont match, despite literally being the same
    #it may also be the case that one instance has 1 or 2 details and another has around 5
    def is_match_absolute(self, instance_2):
        if self.check_if_same_instances(instance_2): return True
        counter_max=0
        counter_match=0
        for detail_type in  self.get_unique_details():
            if detail_type=='reaction_with_instances':
                self_l = []
                rn_inst_1_dict = self.get_detail(detail_type)
                for rn_key in rn_inst_1_dict:   self_l.extend(rn_inst_1_dict[rn_key])

                instance_2_l = []
                rn_inst_2_dict = instance_2.get_detail(detail_type)
                for rn_key in rn_inst_2_dict:   instance_2_l.extend(rn_inst_2_dict[rn_key])

                if self_l!=instance_2_l: return False


            else:
                self_l=set(self.get_detail(detail_type,all_possible=True))
                instance_2_l=set(instance_2.get_detail(detail_type,all_possible=True))
                if self_l and instance_2_l:
                    counter_max+=1
                    if self_l.intersection(instance_2_l):
                        counter_match+=1
        if counter_max==counter_match and counter_match and counter_max:
            return True
        return False

    def is_match_instances(self, instance_2,threshold_for_match=None,no_penalties=False):
        if self.is_match_absolute(instance_2):
            return True
        c=self.score_match_instances(instance_2)
        if no_penalties:         c=ceil(c)
        node_type=get_instance_type(self)
        if threshold_for_match: pass
        elif node_type=='Compound': threshold_for_match=3
        elif node_type=='Gene': threshold_for_match=2
        elif node_type=='Protein': threshold_for_match=2
        elif node_type=='Reaction': threshold_for_match=2
        if c >= threshold_for_match:        return True
        else:                               return False

    def unite_instances(self, instance_2, always_unite=False):
        if not self.check_if_same_instances(instance_2):
            if always_unite:
                test_passed = True
            else:
                test_passed = self.is_match_instances(instance_2)
            if test_passed:
                self.unite_instances_bio_specific(instance_2)


    #when uniting, instance 1 will be the receiver, the instance 2 the giver
    def replace_instances(self,detail_type,giver_instance):
        #receiver_instance is self
        receiver_detail_type = get_instance_type(self).lower() + '_instances'
        #instances that we want received to get from giver_instance
        donated_instances=giver_instance.get_detail(detail_type)
        for donated in donated_instances:
            before_donated=str(donated.get_detail(receiver_detail_type,return_convergence=True))
            before_receiving=str(self.get_detail(detail_type,return_convergence=True))

            donated_converged_in = donated.remove_detail(receiver_detail_type, giver_instance)
            giver_instance_converged_in = giver_instance.remove_detail(detail_type, donated)
            donated.set_detail(receiver_detail_type, self, converged_in=donated_converged_in)
            self.set_detail(detail_type, donated, converged_in=giver_instance_converged_in)

            current_donated = str(donated.get_detail(receiver_detail_type, return_convergence=True))
            current_receiving = str(self.get_detail(detail_type, return_convergence=True))



    def remove_detail_type_instances(self,self_instance,detail_type, to_remove):
        res=[]
        #we remove the to_remove from the item list
        if hasattr(self_instance,'instances'):
            if detail_type in self_instance.instances:
                #then we check each db(key) that is in the convergence dictionary {db:converged_instance}
                for db in self_instance.get_detail(detail_type,return_convergence=True):
                    #and then we also remove the to_remove from the convergence dictionary
                    if to_remove in self_instance.get_detail(detail_type,return_convergence=True)[db]:
                        self_instance.get_detail(detail_type,return_convergence=True)[db].pop(to_remove)
                        #then we attach the dbs that should now set as convergence in the other instance
                        res.append(db)
        return res


if __name__ == '__main__':
    bc=Base_Component({'kegg_ko':'KOOOO1'})
    bc.set_detail('kegg_ko','1')
    print(bc.get_detail('kegg_ko',all_possible=True))
