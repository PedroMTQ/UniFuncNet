
# DRAX modules
from source.Utils.util import strip_tags
from source.Pipelines.Pipelines_Utils.Global_Searcher import *
from source.Fetchers.Compound_Fetchers.Compound_Fetcher_Biocyc import Compound_Fetcher_Biocyc
from source.Fetchers.Compound_Fetchers.Compound_Fetcher_KEGG import Compound_Fetcher_KEGG
from source.Fetchers.Compound_Fetchers.Compound_Fetcher_Chemspider import Compound_Fetcher_Chemspider
from source.Fetchers.Compound_Fetchers.Compound_Fetcher_HMDB import Compound_Fetcher_HMDB
from source.Utils.util import get_stoichiometry
from source.Utils.CHEBI_SQLITE_Connector import CHEBI_SQLITE_Connector

# External modules
from urllib.parse import quote_plus

'''
This class has different methods for searching through databases.
1-Method all_derivatives:
    Search compound name and retrieve all the derivatives
This step starts by searching a database and extracting all search results it finds.
It then cross-matches all these compounds into a final list of compounds.
 
2- Method find_compound:
    Run all_derivavatives method and find best match between the found compounds and the initial query word.

3-Method collect_all_info:
    From a single or several IDs, other external IDs are collect through several databases. (step 1)
This method starts by collecting all the information from the provided IDs.
Then it can go through 2 methods: (step 2)
        The first applies a search by synonyms by using the method find_compound.
        The second collects an ID that was collected during step 1 by using the method find_info
It does so until all available DBs have been visited.

'''

########INFO GETTER
class Compound_Searcher(Global_Searcher,CHEBI_SQLITE_Connector):
    def __init__(self,memory_storage=None,search_direction=None,db_name=None,wanted_org_kegg_codes=[],output_folder=None,politeness_timer=10):
        Global_Searcher.__init__(self,memory_storage,search_direction,db_name=db_name,wanted_org_kegg_codes=wanted_org_kegg_codes,output_folder=output_folder,politeness_timer=politeness_timer)
        CHEBI_SQLITE_Connector.__init__(self)
        self.close_sql_connection()


    def find_compound(self,db,query_id,already_found=set(),convergence_search=False):
        if self.check_already_searched_memory(db,query_id):
            return self.get_compound_match(query_id,db)
        print('Finding compound',db,query_id)
        if db=='synonyms':
            return self.find_compound_string(query_id,already_found,convergence_search=convergence_search)
        else:
            return self.find_info(db,query_id,convergence_search=convergence_search)


    def select_fetcher(self,db,query_id,init_Fetcher=True):
        if db in SCRAPPABLE_DBS and not self.check_already_searched_memory(db,query_id):
            if db == 'biocyc':          return  Compound_Fetcher_Biocyc(query_id, memory_storage=self.memory_storage)
            elif db == 'kegg':          return  Compound_Fetcher_KEGG(query_id, memory_storage=self.memory_storage)
            elif db == 'chemspider':    return  Compound_Fetcher_Chemspider(query_id, memory_storage=self.memory_storage)
            elif db == 'hmdb':          return  Compound_Fetcher_HMDB(query_id, memory_storage=self.memory_storage)
            elif db == 'inchi_key' and 'chemspider' in SCRAPPABLE_DBS:     return  Compound_Fetcher_Chemspider(query_id, memory_storage=self.memory_storage,search_by_inchi_key=True)
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
                    if self.is_valid_search_direction({'global','crpg','crp','cr'}):
                        fetcher.converge_compound_global()
                return fetcher_cpd,fetcher
            else:
                return None, None
        else:
            return None,None


######Retrieval of compounds for reactions#######


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
            met_id, met_name, match, temp_met = None, None, None, None
            whole_met_name = rn[i]
            met_name = whole_met_name[1].lower()
            #####First match search#####
            search_by_name = self.get_compound_match(met_name, 'synonyms')

            if not search_by_name and rn != rn_with_ids:
                if rn_with_ids:
                    whole_met_id = rn_with_ids[i]
                    met_id = whole_met_id[1]
                if met_id:
                    search_by_id = self.get_compound_match(met_id, db)
                else:
                    search_by_id = None
                if not search_by_id:
                    if met_id:
                        # This is the last resource for finding matches, it will get all the information and try to match with that
                        print('Will now search for ' + db + ' ' + met_id)
                        temp_met = self.run_searcher(met_id,db)
                    else:
                        # sometimes reactions dont have an id, when this happens we search by name (should be rare as searching by name is less precise)
                        temp_met = self.run_searcher(met_name,'synonyms')
                    if temp_met: print('Had to search for compound ' + met_name + ' but found it!')
                    #####Third match search######
                    match = self.get_compound_match(temp_met)
                    if match and match is not temp_met:
                        match.unite_instances(temp_met)
                else:
                    match = search_by_id
                    if not match.get_detail(db):
                        found_met,_ = self.find_info(db, met_id)
                        if found_met:
                            match.unite_instances(found_met)
                    match.set_detail(db, met_id)
                    if db != 'synonyms': match.set_detail('synonyms', met_name)
            else:
                # some reactions wont have ids, when that happens
                if not search_by_name:
                    pass
                match = search_by_name
                if met_id:
                    if not match.get_detail(db):
                        found_met,_ = self.find_info(db, met_id)
                        if found_met:
                            match.unite_instances(found_met)
                    match.set_detail(db, met_id)
                    if db != 'synonyms': match.set_detail('synonyms', met_name)

            # if the metabolite has not been seen before
            if match:
                rn_with_instances[side].append([whole_met_name[0], match])
            elif not match and temp_met:
                rn_with_instances[side].append([whole_met_name[0], temp_met])
                if temp_met:
                    if not temp_met.is_empty_metabolite(): self.add_compound(temp_met)
            else:
                #last case scenario whent he db has broken links e.g. kegg G13127
                temp_met = Compound({db: met_id, 'synonyms': [met_name]})
                rn_with_instances[side].append([whole_met_name[0], temp_met])
                self.add_compound(temp_met)
                print('POSSIBLE BROKEN LINK HERE')

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
        url = 'http://www.hmdb.ca/unearth/q?query=' + compound + '&searcher=metabolites&button='
        webpage = self.get_with_fetcher(url)
        if not webpage: return None
        soup = BeautifulSoup(webpage, 'lxml')
        if soup.find(class_='no-results'): return None
        pages = soup.find_all(class_='page_info')
        if not pages: return None
        pages=pages[0]
        search_number=[int(strip_tags(i)) for i in re.findall('\d+\s',strip_tags(pages.text))]
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
                url = 'http://www.hmdb.ca/unearth/q?&page=' + str(page_number) + '&query=' + compound + '&searcher=metabolites&button='
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

    # this searches for a certain compound and all of its derivatives, returning a list of biocyc ids
    def derivatives_biocyc_ids(self,compound):
        cpd_search = quote_plus(compound)
        url = 'https://biocyc.org/META/substring-search?type=NIL&object=' + cpd_search + '&quickSearch=Quick+Search'
        webpage = self.get_with_fetcher(url,selenium=True,original_response=True)
        if not webpage: return None
        current_url=webpage.current_url
        page_source = webpage.page_source
        if 'selenium.webdriver' in repr(webpage): webpage.quit()
        if current_url!= url:
            #sometimes biocyc will redirect to gene page, we want to extract info from this page too
            #this will only happen when coming from the protein_searcher method
            try:
                print('REDIRECTED','current',webpage.current_url,'previous',url)
            except: pass
            return [self.redirected_compound_page(current_url)]
        else:
            return self.NOT_redirected_compound_page(page_source,compound)

    def redirected_compound_page(self,current_url):
        temp_url=current_url.replace('&redirect=T','')
        #sometimes it redirects to a protein... which shouldnt happen
        if 'EC-' in temp_url: return None
        has_id=temp_url.split('&id=')
        if len(has_id)>1: return has_id[1]
        has_object= temp_url.split('&object=')
        if len(has_object)>1: return has_object[1]
        print('COULDNT REDIRECT', temp_url)
        return None



    def NOT_redirected_compound_page(self,webpage,compound_str):
        res = []
        soup = BeautifulSoup(webpage, 'lxml')
        if soup.find_all(text='No exact matches were found.  Showing matches for '): return None
        search = soup.find_all('a', {'href': re.compile('/compound')})
        for i in search:
            cpd_syn=i.text
            cpd_id=i['href'][re.search('META&id=', i['href']).span()[1]:]
            if cpd_syn.lower()==compound_str.lower():
                return [cpd_id]
            if cpd_id not in res: res.append(cpd_id)
        return res

    def derivatives_biocyc(self, compound,limit_mets=20):
        res = []
        der_ids = self.derivatives_biocyc_ids(compound)
        if not der_ids: return None
        #some queries
        for compound_id in der_ids[0:limit_mets]:
            fetcher_cpd,fetcher_inst= self.find_info('biocyc',compound_id)
            if fetcher_inst:
                if fetcher_cpd:
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
            if db == 'biocyc':    res= self.derivatives_biocyc(compound)
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
                        if self.is_valid_search_direction({'global','crpg','crp','cr'}):
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
                #print('Compound string ',search_cpd,' wasn\'t found anywhere so we just created a placeholder of this compound!')
                #return Compound({'synonyms':search_cpd}),None
                return None,None

    def run_searcher(self,bio_query,bio_db,convergence_search=False):
        print(f'STARTING COMPOUND SEARCHER {bio_query} in {bio_db}')
        args_to_search=[]
        temp_inst=None
        #we roll out the search by filling out our args_to_search varible with all the possible data from what is provided
        if bio_db=='synonyms':
            args_to_search.append(['synonyms',bio_query])
        if bio_db=='inchi_key':
            args_to_search.append(['inchi_key',bio_query])
        if bio_db=='chebi':
            self.start_sqlite_cursor()
            chebi_to_others=self.fetch_chebi_id_info(bio_query)
            self.close_sql_connection()
            for chebi_db in chebi_to_others:
                for chebi_db_id in chebi_to_others[chebi_db]:
                    args_to_search.append([chebi_db,chebi_db_id])
        if bio_db in SCRAPPABLE_DBS:
            args_to_search.append([bio_db,bio_query])
        #now we just keep searching for all the possible compounds
        self.run_searcher_ids(args_to_search,convergence_search=convergence_search)
        #current best match
        res= self.get_compound_match(bio_query,bio_db)
        already_found=self.add_to_args_to_search_synonyms(res,args_to_search)
        self.run_searcher_synonyms(args_to_search,already_found,convergence_search=convergence_search)
        res = self.get_compound_match(bio_query,bio_db)




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
                if temp_inst:   self.add_to_args_to_search_ids(temp_inst, args_to_search)

    def add_to_args_to_search_ids(self,compound_instance,args_to_search):
        if compound_instance:
            id_to_add = compound_instance.get_detail('chemspider')
            if id_to_add and not self.check_already_searched_memory(dict_input_key='chemspider',dict_input_value=id_to_add):
                args_to_search.append(['chemspider', id_to_add])
            id_to_add = compound_instance.get_detail('inchi_key')
            if id_to_add and not self.check_already_searched_memory(dict_input_key='inchi_key',dict_input_value=id_to_add):
                args_to_search.append(['inchi_key', id_to_add])
            for db in SCRAPPABLE_DBS:
                id_to_add = compound_instance.get_detail(db)
                if id_to_add and not self.check_already_searched_memory(dict_input_key=db,dict_input_value=id_to_add):
                    args_to_search.append([db, id_to_add])

    #we do this last to complement with dbs we haven't found the compound, we only do this using the initial compound as seed
    def add_to_args_to_search_synonyms(self,compound_instance,args_to_search):
        already_found=set()
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
                    args_to_search.append(['synonyms',syn])
        return already_found

    def run_searcher_synonyms(self,args_to_search,already_found,convergence_search=False):
        print(f'Searching by synonyms')
        while args_to_search:
            current_arg = args_to_search.pop(0)
            db_arg= current_arg[0]
            id_arg= current_arg[1]
            if isinstance(id_arg,str) or isinstance(id_arg,int): id_arg={id_arg}
            for s_id_arg in id_arg:
                if s_id_arg and not self.check_already_searched_memory(dict_input_key=db_arg,dict_input_value=s_id_arg):
                    temp_inst,_ = self.find_compound(db_arg,s_id_arg,already_found=already_found,convergence_search=convergence_search)
                    if temp_inst:
                        syns=temp_inst.get_detail('synonyms',all_possible=True)
                        #to avoid very long searches
                        if id_arg.intersection(syns): return
                self.add_to_already_tried_to_search(db_arg, s_id_arg)
                #if temp_inst:   self.add_to_args_to_search_ids(temp_inst, args_to_search)

if __name__ == '__main__':
    searcher = Compound_Searcher()
    #searcher.derivatives_kegg('pi')
    #searcher.search_direction='global'
    #searcher.reset_db(delete_all=True,force_reset=True)
    #res=searcher.find_compound_string('quercetin')
    #print(res)
    searcher.run_searcher('C00080','kegg')
