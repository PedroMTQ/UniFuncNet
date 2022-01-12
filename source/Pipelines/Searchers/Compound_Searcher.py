
# DRAX modules
from source.Utils.util import strip_tags
from source.Pipelines.Pipelines_Utils.Global_Searcher import *
from source.Fetchers.Compound_Fetchers.Compound_Fetcher_Metacyc import Compound_Fetcher_Metacyc
from source.Fetchers.Compound_Fetchers.Compound_Fetcher_KEGG import Compound_Fetcher_KEGG
from source.Fetchers.Compound_Fetchers.Compound_Fetcher_PubChem import Compound_Fetcher_PubChem
from source.Fetchers.Compound_Fetchers.Compound_Fetcher_HMDB import Compound_Fetcher_HMDB
from source.Fetchers.Compound_Fetchers.Compound_Fetcher_Rhea import Compound_Fetcher_Rhea
from source.Utils.util import get_stoichiometry

# External modules
from urllib.parse import quote_plus

class Compound_Searcher(Global_Searcher):
    def __init__(self,memory_storage=None,search_mode=None,db_name=None,wanted_org_kegg_codes=[],output_folder=None,politeness_timer=10):
        Global_Searcher.__init__(self,memory_storage,search_mode,db_name=db_name,wanted_org_kegg_codes=wanted_org_kegg_codes,output_folder=output_folder,politeness_timer=politeness_timer)



    def find_compound(self,db,query_id,already_found=set(),convergence_search=False):
        if self.check_already_searched_memory(db,query_id):
            return self.get_compound_match(query_id,db),None
        if db=='synonyms':
            print('Finding compound', db, query_id)
            return self.find_compound_string(query_id,already_found,convergence_search=convergence_search)
        else:
            return self.find_info(db,query_id,convergence_search=convergence_search)


    def select_fetcher(self,db,query_id):
        #pubchem is a special case
        if db in ['pubchem_cid', 'pubchem_sid', 'inchi', 'inchi_key'] and 'pubchem' in SCRAPPABLE_DBS:
            return Compound_Fetcher_PubChem(query_id, db=db, memory_storage=self.memory_storage)
        elif db in SCRAPPABLE_DBS and not self.check_already_searched_memory(db,query_id):
            print('Finding compound', db, query_id)
            if db == 'metacyc':         return  Compound_Fetcher_Metacyc(query_id, memory_storage=self.memory_storage)
            elif db == 'kegg':          return  Compound_Fetcher_KEGG(query_id, memory_storage=self.memory_storage)
            elif db == 'hmdb':          return  Compound_Fetcher_HMDB(query_id, memory_storage=self.memory_storage)
            elif db == 'rhea':          return  Compound_Fetcher_Rhea(query_id, memory_storage=self.memory_storage)
            else:                       return  Global_Fetcher()

    #############Searching the Databases#############
    def find_info(self,db,query_id,convergence_search=False):
        if not query_id: return None,None
        fetcher=self.select_fetcher(db=db,query_id=query_id)
        if fetcher:
            self.add_to_already_tried_to_search(db, query_id)
            fetcher_cpd=fetcher.get_compound()
            if fetcher_cpd:
                if convergence_search:
                    # converge only occurs in the searchers- these are the global classes
                    if self.is_valid_search_mode({'global','crpg','crp','cr'}):
                        fetcher.converge_compound_global()
                return fetcher_cpd,fetcher
        else:
            if self.check_already_searched_memory(db, query_id):
                print('Already searched', db, query_id)
                return self.get_compound_match(bio_query=query_id, bio_db=db), None
        return None, None


######Retrieval of compounds for reactions#######

    def get_compound_reaction_met_instances(self,db,met_id, met_name):
        #####First memory search#####
        match = self.get_compound_match(met_name, 'synonyms')
        if match: return match

        #####Second memory search#####
        match = self.get_compound_match(met_id, db)
        if match: return match
        #####Web query search##### --> if not found in memory
        if met_id:
            print(f'Will now search for {met_id} in {db}')
            self.run_searcher(met_id, db)
            #####Third memory search#####
            match = self.get_compound_match(met_id, db)
        else:
            print(f'Will now search for {met_name}')
            self.run_searcher(met_name, 'synonyms')
            #####Third memory search#####
            match = self.get_compound_match(met_name, 'synonyms')

        if match: return match
        match = Compound({db: met_id, 'synonyms': [met_name]})
        print('POSSIBLE BROKEN LINK HERE', db, met_id, met_name)
        self.add_compound(match)
        return match

    def reaction_met_instances_simple(self, reaction_stoichiometry,  db):
        print('Reaction from', db, ': ', reaction_stoichiometry)
        rn_with_instances = {'left':[],'right':[]}
        side='left'
        for i in reaction_stoichiometry:
            if isinstance(i,list):
                stoi,met_id=i
                match = self.get_compound_match(met_id, db)
                if not match:
                    print(f'Will now search for {met_id} in {db}')
                    self.run_searcher(met_id, db)
                    match = self.get_compound_match(met_id, db)
                rn_with_instances[side].append([stoi, match])
            else:
                side='right'

        return rn_with_instances

    def reaction_met_instances(self, rn, rn_with_ids, db):
        """
        :param rn: reaction string as seen on the website
        :param rn_with_ids: reaction as a list of [n molecules][compound]
        :param db: database
        :return: returns a list of [n molecules][compound instance]
        """
        print('Reaction from', db, ': ', rn)
        rn,complete_l,len_sub = get_stoichiometry(rn, rn)
        # when the reaction is incomplete (either no products or no substrates
        if not rn: return None
        rn_with_instances = {'left':[],'right':[]}
        for i in range(len(rn)):
            if i<len_sub: side='left'
            else: side='right'
            met_id, met_name= None, None
            whole_met_name = rn[i]
            met_name = whole_met_name[1].lower()

            if rn != rn_with_ids and rn_with_ids:
                whole_met_id = rn_with_ids[i]
                met_id = whole_met_id[1]

            match=self.get_compound_reaction_met_instances(db,met_id, met_name)


            if not match.get_detail(db):
                found_met, _ = self.find_info(db, met_id)
                if found_met:
                    match.unite_instances(found_met)

                match.set_detail(db, met_id)
            if db != 'synonyms': match.set_detail('synonyms', met_name)
            rn_with_instances[side].append([whole_met_name[0], match])
        return rn_with_instances

    #############Getting all derivatives#############

    def derivatives_kegg(self, compound):
        search = self.get_with_fetcher(compound,api_kegg=True,database='cpd')
        res = []
        for i in search:
            i = i.split('\t')
            kegg_id = i[0].replace('cpd:','').strip()
            synonyms= i[1].split(';')
            synonyms=[i.strip() for i in synonyms]
            for syn in synonyms:
                if syn.lower()==compound:
                    fetcher_cpd,fetcher_inst=self.find_info('kegg',kegg_id)
                    if fetcher_inst:
                        if fetcher_cpd:
                            res.append(fetcher_inst)
        search = self.get_with_fetcher(compound,api_kegg=True,database='compound',kegg_option='formula')
        for i in search:
            i = i.split('\t')
            kegg_id = i[0].replace('cpd:','').strip()
            synonyms= i[1].split(';')
            synonyms=[i.strip() for i in synonyms]
            for syn in synonyms:
                if syn.lower()==compound:
                    fetcher_cpd,fetcher_inst=self.find_info('kegg',kegg_id)
                    if fetcher_inst:
                        if fetcher_cpd:
                            res.append(fetcher_inst)
        return res



    def derivatives_HMDB_ids(self, compound,limit_mets=20,number_search_not_exact=200):
        url = f'http://www.hmdb.ca/unearth/q?query={compound}&searcher=metabolites&button='
        webpage = self.get_with_fetcher(url)
        if not webpage: return None
        soup = BeautifulSoup(webpage, 'lxml')
        if soup.find(class_='no-results'): return None
        pages = soup.find_all(class_='page_info')
        if not pages: return None
        pages=pages[0]
        try:
            search_number=[int(strip_tags(i)) for i in re.findall('\d+\s',strip_tags(pages.text))]
        except: return None
        #when we have several pages of results
        if len(search_number)==3:
            several_pages=True
            page_met_limit=search_number[1]
            total_mets=search_number[2]
        #when we only have one page of results
        else:
            several_pages=False
            page_met_limit=search_number[0]
            total_mets=search_number[0]
        # when there's no exact match and many matches are returned (non relevant matches)
        if soup.find_all('div', class_='search-info',id='didyoumean') or total_mets>number_search_not_exact: return None
        #limit the number of search results. If limit_mets is 0, there's no limitaton
        if limit_mets:
            if total_mets>limit_mets: total_mets=limit_mets
        page_number=1
        current_met=0
        res = []
        while current_met<total_mets:
            derivatives = soup.find_all(class_='btn-card', href=re.compile('/metabolites/HMDB\d+'))
            for der in derivatives:
                if der.text not in res: res.append(der.text)
                current_met+=1
            if current_met==page_met_limit and page_met_limit != total_mets and several_pages:
                page_number+=1
                url = f'http://www.hmdb.ca/unearth/q?&page={page_number}&query={compound}&searcher=metabolites&button='
                webpage = self.get_with_fetcher(url)
                soup = BeautifulSoup(webpage, 'lxml')
                pages = soup.find(class_='page_info')
                search_number = [int(strip_tags(i)) for i in re.findall('\d+\s', pages.text)]
                page_met_limit = search_number[1]
        return res

    def derivatives_HMDB(self, compound):
        der_ids = self.derivatives_HMDB_ids(compound)
        if not der_ids: return None
        res = []
        for i in range(len(der_ids)):
            fetcher_cpd,fetcher_inst = self.find_info('hmdb',der_ids[i])
            # hmdb query doesn't work very well since there's some misinformation in synonyms, but this can fixed by making sure they match
            if fetcher_inst:
                if fetcher_cpd:
                    #first we try to find exact matches, if we find one we just return immediately, otherwise we keep adding to res to gather more info
                    if self.find_match_synonyms([compound], fetcher_cpd.get_detail('synonyms'), equal=True):
                        res=[fetcher_inst]
                        return res
                    if self.find_match_synonyms([compound], fetcher_cpd.get_detail('synonyms'), equal=False):
                        res.append(fetcher_inst)
        if not res: return None
        return res


    def derivatives_metacyc(self, compound):
        derivatives= self.fetch_metacyc_derivatives(compound)
        if not derivatives: return None
        res = []
        for cpd_id in derivatives:
            fetcher_cpd,fetcher_inst = self.find_info('metacyc',cpd_id)
            if fetcher_inst:
                if fetcher_cpd:
                    if self.find_match_synonyms([compound], fetcher_cpd.get_detail('synonyms'), equal=True):
                        res=[fetcher_inst]
                        return res
                    if self.find_match_synonyms([compound], fetcher_cpd.get_detail('synonyms'), equal=False):
                        res.append(fetcher_inst)
        if not res: return None
        return res



#############Joining queried info together#############

    def find_match_synonyms(self, cpd_syns, query_syns, equal=True):
        if not cpd_syns or not query_syns: return False
        for cpd_syn in cpd_syns:
            for query_syn in query_syns:
                if equal:
                    if cpd_syn.lower() == query_syn.lower():
                        return True
                if not equal:
                    if cpd_syn.lower() in query_syn.lower():
                        return True
        return False



#############Managing functions#############

    def search_derivatives(self,compound,db):
        if not self.check_already_searched_memory(db,compound):
            res=None
            if db == 'hmdb':      res= self.derivatives_HMDB(compound)
            if db == 'metacyc':    res= self.derivatives_metacyc(compound)
            if db == 'kegg':      res= self.derivatives_kegg(compound)
            self.add_to_already_tried_to_search(db,compound)
            return res



#############Types of search#############

    def all_derivatives(self, search_cpd,dbs_to_use,convergence_search=False):
        fetcher_insts=[]
        if dbs_to_use:
            for db in dbs_to_use:
                derivatives=self.search_derivatives(search_cpd,db)
                if derivatives:
                    fetcher_insts.extend(derivatives)
            compound_match=self.get_compound_match(search_cpd,'synonyms')
            if compound_match:
                passed_check=len(dbs_to_use)
                for db in dbs_to_use:
                    if compound_match.get_detail(db): passed_check-=1
                for fetcher in fetcher_insts:
                    fetcher_compound=fetcher.get_compound()
                    if fetcher_compound is compound_match and convergence_search:
                        if self.is_valid_search_mode({'global','crpg','crp','cr'}):
                            fetcher.converge_compound_global()
                if not passed_check:
                    return compound_match

    #this will find a compound and retrieve information, based on the best match found accross several databases
    #works when user inputs a compound name
    def find_compound_string(self,search_cpd,already_found,convergence_search=False):
        """
        :param de_novo_search: False if we want to retrieve data from neo4j (if available), else True
        :param search_cpd: string word of a compound to find
        :return:    returns a single compound. This is the best match from all the search results across the several databases
        """
        dbs_to_exclude=set()
        #first we just get the best match from the db
        match=self.get_compound_match(search_cpd,'synonyms')
        wanted_details = {i:None for i in SCRAPPABLE_DBS}
        #if its found we then check if this instance has all the ids from the scrappable Dbs
        if match:
            for i in wanted_details:
                if match.get_detail(i):
                    dbs_to_exclude.add(i)
                    wanted_details[i] = match.get_detail(i)
            #if all scrappable dbs have their own ID we stop here
            if number_of_nones_dict(wanted_details)==0:
                return match,None
        dbs_to_exclude.update(already_found)
        dbs_to_use = SCRAPPABLE_DBS
        dbs_to_use=set([i for i in dbs_to_use if i not in dbs_to_exclude])
        if dbs_to_use:
            str_dbs_to_use=', '.join(dbs_to_use)
            print(f'Searching for derivates of {search_cpd} in {str_dbs_to_use}')
            found_cpd=self.all_derivatives(search_cpd,dbs_to_use,convergence_search=convergence_search)
            if found_cpd: return found_cpd,None
            else:
                return None,None
        return None,None

    def run_searcher(self,bio_query,bio_db,convergence_search=False):
        print(f'STARTING COMPOUND SEARCHER {bio_query} in {bio_db}')
        temp_args_to_search=[]
        args_to_search=[]
        temp_inst=None
        #we roll out the search by filling out our args_to_search varible with all the possible data from what is provided
        if bio_db=='synonyms':
            temp_args_to_search.append(['synonyms',bio_query])
        if bio_db in ['pubchem_cid','pubchem_sid','inchi','inchi_key'] and 'pubchem' in SCRAPPABLE_DBS:
            temp_args_to_search.append([bio_db,bio_query])
        if bio_db=='chebi':
            temp_args_to_search.append(['rhea', bio_query])
            chebi_to_others=self.fetch_chebi_id_info(bio_query)
            for chebi_db in chebi_to_others:
                if chebi_db not in ['smiles','chemical_formula']: #we also have this info in the sql database, wont really change much, just avoid function calls since there's no fetcher for these
                    for chebi_db_id in chebi_to_others[chebi_db]:
                        temp_args_to_search.append([chebi_db,chebi_db_id])
        if bio_db in SCRAPPABLE_DBS:
            temp_args_to_search.append([bio_db,bio_query])
        for combo in temp_args_to_search:
            if combo not in args_to_search: args_to_search.append(combo)
        #now we just keep searching for all the possible compounds
        self.run_searcher_ids(args_to_search,convergence_search=convergence_search)
        #current best match
        res= self.get_compound_match(bio_query,bio_db)
        already_found=self.add_to_args_to_search_synonyms(res,args_to_search)
        self.run_searcher_synonyms(args_to_search,already_found,convergence_search=convergence_search)
        res = self.get_compound_match(bio_query,bio_db)
        return res




    def run_searcher_ids(self,args_to_search,convergence_search=False):
        temp_inst = None
        while args_to_search:
            current_arg = args_to_search.pop()
            db_arg= current_arg[0]
            id_arg= current_arg[1]
            if isinstance(id_arg,str) or isinstance(id_arg,int): id_arg={id_arg}
            for s_id_arg in id_arg:
                if s_id_arg and not self.check_already_searched_memory(dict_input_key=db_arg,dict_input_value=s_id_arg):
                    temp_inst,_ = self.find_compound(db_arg,s_id_arg,convergence_search=convergence_search)
                self.add_to_already_tried_to_search(db_arg, s_id_arg)
                if not temp_inst:
                    temp_inst = self.get_compound_match(id_arg, db_arg)
                if temp_inst:   self.add_to_args_to_search_ids(temp_inst, args_to_search)

    def add_to_args_to_search_ids(self,compound_instance,args_to_search):
        temp_args_to_search=[]
        if compound_instance:
            id_to_add = compound_instance.get_detail('pubchem_cid')
            if id_to_add and not self.check_already_searched_memory(dict_input_key='pubchem_cid',dict_input_value=id_to_add):
                temp_args_to_search.append(['pubchem_cid', id_to_add])

            id_to_add = compound_instance.get_detail('pubchem_sid')
            if id_to_add and not self.check_already_searched_memory(dict_input_key='pubchem_sid',dict_input_value=id_to_add):
                temp_args_to_search.append(['pubchem_sid', id_to_add])

            id_to_add = compound_instance.get_detail('inchi')
            if id_to_add and not self.check_already_searched_memory(dict_input_key='inchi',dict_input_value=id_to_add):
                temp_args_to_search.append(['inchi', id_to_add])

            id_to_add = compound_instance.get_detail('inchi_key')
            if id_to_add and not self.check_already_searched_memory(dict_input_key='inchi_key',dict_input_value=id_to_add):
                temp_args_to_search.append(['inchi_key', id_to_add])

            for db in SCRAPPABLE_DBS:
                id_to_add = compound_instance.get_detail(db)
                if id_to_add and not self.check_already_searched_memory(dict_input_key=db,dict_input_value=id_to_add):
                    temp_args_to_search.append([db, id_to_add])
        for combo in temp_args_to_search:
            if combo not in args_to_search:
                args_to_search.append(combo)

    #we do this last to complement with dbs we haven't found the compound, we only do this using the initial compound as seed
    def add_to_args_to_search_synonyms(self,compound_instance,args_to_search):
        already_found=set()
        temp_args_to_search=[]
        if compound_instance:
            for db in SCRAPPABLE_DBS:
                current_id=compound_instance.get_detail(db)
                if current_id and self.check_already_searched_memory(dict_input_key=db, dict_input_value=current_id):
                    already_found.add(db)
            all_syns = compound_instance.get_detail('synonyms',all_possible=True)
            #longer names tend  to have fewer results
            all_syns=sorted(all_syns, key=len,reverse=True)
            for syn in all_syns:
                if syn and not self.check_already_searched_memory(dict_input_key='synonyms', dict_input_value=syn):
                    temp_args_to_search.append(['synonyms',syn])
        for combo in temp_args_to_search:
            if combo not in args_to_search:
                args_to_search.append(combo)
        return already_found

    def run_searcher_synonyms(self,args_to_search,already_found,convergence_search=False):
        while args_to_search:
            current_arg = args_to_search.pop(0)
            db_arg= current_arg[0]
            id_arg= current_arg[1]
            if isinstance(id_arg,str) or isinstance(id_arg,int): id_arg={id_arg}
            for s_id_arg in id_arg:
                if s_id_arg and not self.check_already_searched_memory(dict_input_key=db_arg,dict_input_value=s_id_arg):
                    temp_inst,_ = self.find_compound(db_arg,s_id_arg,already_found=already_found,convergence_search=convergence_search)
                    if not temp_inst:
                        temp_inst=self.get_compound_match(id_arg,db_arg)
                    if temp_inst:
                        syns=temp_inst.get_detail('synonyms',all_possible=True)
                        #to avoid very long searches
                        if id_arg.intersection(syns): return
                self.add_to_already_tried_to_search(db_arg, s_id_arg)

if __name__ == '__main__':
    searcher = Compound_Searcher(search_mode={''})
    #searcher.derivatives_kegg('pi')
    #searcher.search_mode='global'
    #searcher.reset_db(delete_all=True,force_reset=True)
    #res=searcher.find_compound_string('quercetin')
    #print(res)
    r1=searcher.run_searcher('1S/C5H7NO2/c6-5-3(7)1-2-4(5)8/h5H,1-2,6H2','inchi')
