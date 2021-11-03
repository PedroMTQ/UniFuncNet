
from source.Pipelines.Pipelines_Utils.Global_Searcher import *
from source.Fetchers.Reaction_Fetchers.Reaction_Fetcher_Biocyc import Reaction_Fetcher_Biocyc
from source.Fetchers.Reaction_Fetchers.Reaction_Fetcher_HMDB import Reaction_Fetcher_HMDB
from source.Fetchers.Reaction_Fetchers.Reaction_Fetcher_KEGG import Reaction_Fetcher_KEGG

class Reaction_Searcher(Global_Searcher):
    def __init__(self,  memory_storage=None,search_direction='',do_reaction_met_instances=False,db_name=None,wanted_org_kegg_codes=[],output_folder=None,politeness_timer=10):
        Global_Searcher.__init__(self,memory_storage,search_direction,do_reaction_met_instances,db_name=db_name,wanted_org_kegg_codes=wanted_org_kegg_codes,output_folder=output_folder,politeness_timer=politeness_timer)



    def find_reaction(self,db,query_id=None,extra_args={}):
        if db in SCRAPPABLE_DBS:
            fetcher_reaction,fetcher= self.find_info(db, query_id, extra_args)
            return fetcher_reaction


    def select_fetcher(self,db,query_id,extra_args,init_Fetcher=True):
        if db in SCRAPPABLE_DBS and not self.check_already_searched_memory(db,query_id):
            if db == 'biocyc':        return  Reaction_Fetcher_Biocyc(query_id, extra_args=extra_args, memory_storage=self.memory_storage,init_Fetcher=init_Fetcher)
            elif 'kegg' in db:        return  Reaction_Fetcher_KEGG(query_id, extra_args=extra_args, memory_storage=self.memory_storage,init_Fetcher=init_Fetcher)
            elif db == 'hmdb':        return  Reaction_Fetcher_HMDB(query_id, extra_args=extra_args, memory_storage=self.memory_storage,init_Fetcher=init_Fetcher)
            else:                     return  Global_Fetcher()

    def find_info(self, db, query_id, extra_args={}):
        if not query_id: return None,None
        fetcher=self.select_fetcher(db=db,query_id=query_id,extra_args=extra_args)
        if fetcher:
            self.add_to_already_tried_to_search(db, query_id)
            fetcher_reaction=fetcher.get_reaction()
            if fetcher_reaction:
                #converge only occurs in the searchers- these are the global classes
                if  {'global'}.intersection(self.search_direction):     fetcher.converge_reaction_global()
                elif {'rpg','rp'}.intersection(self.search_direction):  fetcher.converge_reaction_rpg()
                return fetcher_reaction,fetcher
            else:
                return None, None
        else:
            return None,None

    def run_searcher(self,bio_query,bio_db):
        print(f'STARTING REACTION SEARCHER {bio_query} in {bio_db}')
        args_to_search=[[bio_db, bio_query]]
        temp_inst=None
        while args_to_search:
            current_arg = args_to_search.pop()
            db_arg= current_arg[0]
            id_arg= current_arg[1]
            if isinstance(id_arg,str) or isinstance(id_arg,int): id_arg={id_arg}
            for s_id_arg in id_arg:
                if s_id_arg and not self.check_already_searched_memory(dict_input_key=db_arg,dict_input_value=s_id_arg):
                    temp_inst = self.find_reaction(db_arg,s_id_arg)
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

                id_to_add = reaction_instance.get_detail('biocyc')
                if id_to_add and not self.check_already_searched_memory(dict_input_key='biocyc',dict_input_value=id_to_add):
                    args_to_search.append(['biocyc', id_to_add])


                id_to_add = reaction_instance.get_detail('hmdb')
                if id_to_add and not self.check_already_searched_memory(dict_input_key='hmdb',dict_input_value=id_to_add):
                    args_to_search.append(['hmdb', id_to_add])




if __name__ == '__main__':
    searcher=Reaction_Searcher(search_direction='rp',do_reaction_met_instances=False)


    #r1=searcher.run_searcher('R05188','kegg')
    #r1.get_all_info()
    #for cpd in r1.get_detail('reaction_with_instances'):
    #    cpd[1].get_all_info()

    r1=searcher.run_searcher('R00091','kegg')
    #r1.get_all_info()
    #for cpd in r1.get_detail('reaction_with_instances'):
    #    cpd[1].get_all_info()
    searcher.output_results()
    #r1=searcher.find_reaction('biocyc','3.4.21.92-RXN')
    #r1=searcher.find_reaction('kegg','R01665')

