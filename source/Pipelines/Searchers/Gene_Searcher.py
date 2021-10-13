
from source.Pipelines.Pipelines_Utils.Global_Searcher import *
from source.Fetchers.Gene_Fetchers.Gene_Fetcher_Biocyc import Gene_Fetcher_Biocyc
from source.Fetchers.Gene_Fetchers.Gene_Fetcher_HMDB import Gene_Fetcher_HMDB
from source.Fetchers.Gene_Fetchers.Gene_Fetcher_KEGG import Gene_Fetcher_KEGG



class Gene_Searcher(Global_Searcher):
    def __init__(self,memory_storage=None,search_direction='',do_reaction_met_instances=False,db_name=None,wanted_org_kegg_codes=[],output_folder=None,politeness_timer=10):
        Global_Searcher.__init__(self,memory_storage,search_direction,do_reaction_met_instances=do_reaction_met_instances,
                                 db_name=db_name,wanted_org_kegg_codes=wanted_org_kegg_codes,output_folder=output_folder,politeness_timer=politeness_timer)

    #this is db specific
    def find_gene(self,db,query_id=None,extra_args={}):
        if db in SCRAPPABLE_DBS:
            return self.find_info(db, query_id, extra_args)

    def select_fetcher(self,db,query_id,extra_args,init_Fetcher=True):
        if db in SCRAPPABLE_DBS and not self.check_already_searched_memory(db,query_id):
            if db == 'biocyc':        return    Gene_Fetcher_Biocyc(query_id,extra_args=extra_args, memory_storage=self.memory_storage,init_Fetcher=init_Fetcher)
            elif 'kegg' in db:        return    Gene_Fetcher_KEGG(query_id,extra_args=extra_args,memory_storage=self.memory_storage,init_Fetcher=init_Fetcher)
            elif db == 'hmdb':        return    Gene_Fetcher_HMDB(query_id,extra_args=extra_args, memory_storage=self.memory_storage,init_Fetcher=init_Fetcher)
            else:                     return    Global_Fetcher()

    def find_info(self,db,query_id,extra_args={}):
        if not query_id: return None

        fetcher=self.select_fetcher(db=db,query_id=query_id,extra_args=extra_args)
        if fetcher:
            self.add_to_already_tried_to_search(db, query_id)
            if fetcher.get_gene():
                #converge only occurs in the searchers- these are the global classes
                if self.search_direction=='global':                 fetcher.converge_gene_global()
                elif self.search_direction=='gpr' or\
                     self.search_direction=='gp':       fetcher.converge_gene_gpr()
                if fetcher.get_gene():
                    return fetcher.get_gene()

    def add_to_args_to_search(self,gene_instance_or_list,args_to_search):
        if not isinstance(gene_instance_or_list, list): gene_instance_or_list = [gene_instance_or_list]
        for gene_instance in gene_instance_or_list:
            if gene_instance:
                for db in SCRAPPABLE_DBS:
                    id_to_add = gene_instance.get_detail(db)
                    if db=='uniprot':
                        id_to_add=self.get_db_id_from_uniprot(id_to_add,db)
                    if not self.check_already_searched_memory(dict_input_key=db, dict_input_value=id_to_add):
                        args_to_search.append([db, id_to_add])




    #global search
    def run_searcher(self,bio_query,bio_db):
        """
        The only common ID in all DBs is Uniprot ID so this is the central point
        If a uniprot Id is provided we start with that, otherwise we go to the other dbs and extract the uniprot Id and go from there
        """
        print(f'STARTING GENE SEARCHER {bio_query} in {bio_db}')
        args_to_search=[]
        temp_inst=None

        if bio_db in SCRAPPABLE_DBS:
            if bio_db=='uniprot':
                for db in ['biocyc','hmdb','kegg']:
                    if db in SCRAPPABLE_DBS:
                        db_id = self.get_db_id_from_uniprot(bio_query,db)
                        if db_id:args_to_search.append([db,db_id])
            else:
                args_to_search.append([bio_db,bio_query])


        while args_to_search:
            current_arg = args_to_search.pop()
            db_arg= current_arg[0]
            id_arg= current_arg[1]
            if isinstance(id_arg,str) or isinstance(id_arg,int): id_arg={id_arg}
            for s_id_arg in id_arg:
                if s_id_arg and not self.check_already_searched_memory(dict_input_key=db_arg,dict_input_value=s_id_arg):
                    temp_inst = self.find_gene(db_arg,s_id_arg)
                self.add_to_already_tried_to_search(db_arg, s_id_arg)
                if temp_inst:   self.add_to_args_to_search(temp_inst, args_to_search)
        return self.get_gene_match(bio_query,bio_db)



    def get_db_id_from_uniprot(self,uniprot_id,db):
        if db=='biocyc': return self.get_db_id_from_uniprot_biocyc(uniprot_id)
        if db=='hmdb': return self.get_db_id_from_uniprot_hmdb(uniprot_id)
        if db=='kegg': return self.get_db_id_from_uniprot_kegg(uniprot_id)


    def get_db_id_from_uniprot_biocyc(self,uniprot_id):
        url='https://biocyc.org/META/search-query?type=GENE&pname='+uniprot_id
        webpage=self.get_with_fetcher(url)
        if not webpage: return None
        soup = BeautifulSoup(webpage, 'lxml')
        result_table= soup.find(class_='sortableSAQPoutputTable')
        if result_table:
            res=result_table.find_all('tr')
            res= res[1:]
            if len(res)>1:
                print(uniprot_id,url)
                print('TOO MANY RESULTS BIOCYC')
                return
            gene_id,protein_id,organism= res[0].find_all('a')
            gene_id = re.search('&id=.*',gene_id['href']).group()[4:]
            return gene_id

    def get_db_id_from_uniprot_hmdb(self,uniprot_id):
        url='http://www.hmdb.ca/unearth/q?utf8=%E2%9C%93&query='+uniprot_id+'&searcher=proteins&button='
        webpage=self.get_with_fetcher(url)
        if not webpage: return None
        soup = BeautifulSoup(webpage, 'lxml')
        res= soup.find('a',{'href':re.compile('.*uniprot/'+uniprot_id)})
        if res:
            hmdb_id=res.parent.parent.a.text
            return hmdb_id


    def get_db_id_from_uniprot_kegg(self,uniprot_id):
        if uniprot_id:
            #using kegg api
            webpage = self.get_with_fetcher(api_kegg=True, url='uniprot:' + uniprot_id, type_search='conv',
                                            database='genes')
            if webpage:
                kegg_id = webpage[0].split('\t')[1]
                return kegg_id
            else:
                #using uniprot api
                params = {'from': 'ACC', 'to': 'KEGG_ID', 'format': 'tab', 'query': uniprot_id}
                url = 'https://www.uniprot.org/uploadlists/'
                webpage = self.get_with_fetcher(url,data=params)
                split_lines = webpage.split('\n')
                if len(split_lines)>1:
                    kegg_id=split_lines[1]
                    if kegg_id:
                        kegg_id=kegg_id.split('\t')[1]
                        return kegg_id



if __name__ == '__main__':
    searcher=Gene_Searcher(search_direction='',do_reaction_met_instances=True,wanted_org_kegg_codes=[])
    gene=searcher.run_searcher('hsa:3098','kegg')
    gene.get_all_info()
    print('#####')
    gene=searcher.run_searcher('P19367','uniprot')
    gene.get_all_info()