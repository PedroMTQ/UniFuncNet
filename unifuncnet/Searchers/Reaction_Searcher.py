
from unifuncnet.Searchers.Global_Searcher import *
from unifuncnet.Fetchers.Reaction_Fetchers.Reaction_Fetcher_Metacyc import Reaction_Fetcher_Metacyc
from unifuncnet.Fetchers.Reaction_Fetchers.Reaction_Fetcher_HMDB import Reaction_Fetcher_HMDB
from unifuncnet.Fetchers.Reaction_Fetchers.Reaction_Fetcher_KEGG import Reaction_Fetcher_KEGG
from unifuncnet.Fetchers.Reaction_Fetchers.Reaction_Fetcher_Rhea import Reaction_Fetcher_Rhea


class Reaction_Searcher(Global_Searcher):
    def __init__(self,  memory_storage=None,search_mode='',db_name=None,wanted_org_kegg_codes=[],output_folder=None,politeness_timer=10):
        Global_Searcher.__init__(self,memory_storage,search_mode,db_name=db_name,wanted_org_kegg_codes=wanted_org_kegg_codes,output_folder=output_folder,politeness_timer=politeness_timer)



    def find_reaction(self,db,query_id=None,extra_args={}):
        if db in SCRAPPABLE_DBS:
            fetcher_reaction,fetcher= self.find_info(db, query_id, extra_args)
            return fetcher_reaction


    def select_fetcher(self,db,query_id,extra_args):
        if db in SCRAPPABLE_DBS and not self.check_already_searched_memory(db,query_id):
            if db == 'metacyc':        return  Reaction_Fetcher_Metacyc(query_id, extra_args=extra_args, memory_storage=self.memory_storage)
            elif 'kegg' in db:        return  Reaction_Fetcher_KEGG(query_id, extra_args=extra_args, memory_storage=self.memory_storage)
            elif db == 'hmdb':        return  Reaction_Fetcher_HMDB(query_id, extra_args=extra_args, memory_storage=self.memory_storage)
            elif db == 'rhea':        return  Reaction_Fetcher_Rhea(query_id, extra_args=extra_args, memory_storage=self.memory_storage)
            else:                     return  Global_Fetcher()

    def find_info(self, db, query_id, extra_args={}):
        if not query_id: return None,None
        fetcher=self.select_fetcher(db=db,query_id=query_id,extra_args=extra_args)
        self.add_to_already_tried_to_search(db, query_id)
        if fetcher:
            fetcher_reaction=fetcher.get_reaction()
            if fetcher_reaction:
                #converge only occurs in the searchers- these are the global classes
                if self.is_valid_search_mode({'global','rpg','rp','crp','crpg'}):
                    fetcher.converge_reaction_rpg()
                return fetcher_reaction,fetcher
            else:
                return None, None
        else:
            if self.check_already_searched_memory(db, query_id):
                return self.get_reaction_match(bio_query=query_id, bio_db=db), None
            else:
                return None, None

    def run_searcher(self,bio_query,bio_db,extra_args={}):
        args_to_search=[[bio_db, bio_query]]
        temp_inst=None
        while args_to_search:
            current_arg = args_to_search.pop()
            db_arg= current_arg[0]
            id_arg= current_arg[1]
            if isinstance(id_arg,str) or isinstance(id_arg,int): id_arg={id_arg}
            for s_id_arg in id_arg:
                if s_id_arg and not self.check_already_searched_memory(dict_input_key=db_arg,dict_input_value=s_id_arg):
                    temp_inst = self.find_reaction(db_arg,s_id_arg,extra_args=extra_args)
                self.add_to_already_tried_to_search(db_arg, s_id_arg)
                if temp_inst:   self.add_to_args_to_search(temp_inst, args_to_search)
        return self.get_reaction_match(bio_query,bio_db)


    def add_to_args_to_search(self, reaction_instance_or_list, args_to_search):
        if not isinstance(reaction_instance_or_list, list): reaction_instance_or_list = [reaction_instance_or_list]
        for reaction_instance in reaction_instance_or_list:
            if reaction_instance:
                id_to_add = reaction_instance.get_detail('kegg')
                if id_to_add and not self.check_already_searched_memory(dict_input_key='kegg',dict_input_value=id_to_add):
                    args_to_search.append(['kegg', id_to_add])

                id_to_add = reaction_instance.get_detail('metacyc')
                if id_to_add and not self.check_already_searched_memory(dict_input_key='metacyc',dict_input_value=id_to_add):
                    args_to_search.append(['metacyc', id_to_add])


                id_to_add = reaction_instance.get_detail('hmdb')
                if id_to_add and not self.check_already_searched_memory(dict_input_key='hmdb',dict_input_value=id_to_add):
                    args_to_search.append(['hmdb', id_to_add])

                rhea_ids = reaction_instance.get_detail('rhea',all_possible=True)
                for id_to_add in rhea_ids:
                    if id_to_add and not self.check_already_searched_memory(dict_input_key='rhea',dict_input_value=id_to_add):
                        args_to_search.append(['rhea', id_to_add])




if __name__ == '__main__':
    searcher=Reaction_Searcher(search_mode={'rc'})



    #r1=searcher.run_searcher('RXN-14064','metacyc')
    r1=searcher.run_searcher('R05135','kegg')
    #r1=searcher.run_searcher('27860','rhea')
    r1.get_all_info()
    #for cpd in r1.get_detail('reaction_with_instances'):
    #    cpd[1].get_all_info()

    #r1=searcher.run_searcher('RXN-20993','metacyc')
    #r1.get_all_info()
    #for cpd in r1.get_detail('reaction_with_instances'):
    #    cpd[1].get_all_info()
    #r1=searcher.find_reaction('metacyc','3.4.21.92-RXN')
    #r1=searcher.find_reaction('kegg','R01665')

