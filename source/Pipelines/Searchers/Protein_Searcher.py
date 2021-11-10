
from source.Pipelines.Searchers.Gene_Searcher import Gene_Searcher
from source.Pipelines.Pipelines_Utils.Global_Searcher import *
from source.Biological_Components.Protein import Protein
from source.Fetchers.Protein_Fetchers.Protein_Fetcher_Biocyc import Protein_Fetcher_Biocyc
from source.Fetchers.Protein_Fetchers.Protein_Fetcher_KEGG import Protein_Fetcher_KEGG
from source.Fetchers.Protein_Fetchers.Protein_Fetcher_HMDB import Protein_Fetcher_HMDB
from source.Utils.util import is_ec


from types import GeneratorType as generator

class Protein_Searcher(Global_Searcher):
    def __init__(self,memory_storage=None,search_direction='',db_name=None,wanted_org_kegg_codes=[],output_folder=None,politeness_timer=10):
        Global_Searcher.__init__(self,memory_storage,search_direction,
                                 db_name=db_name,wanted_org_kegg_codes=wanted_org_kegg_codes,output_folder=output_folder,politeness_timer=politeness_timer)
        self.original_searcher=get_instance_type(self)

    def find_protein(self,db,query_id=None,extra_args={}):
        print(1,db,query_id,extra_args)
        if db in SCRAPPABLE_DBS:
            fetcher_protein,fetcher= self.find_info(db, query_id, extra_args)
            return fetcher_protein

    def select_fetcher(self,db,query_id,extra_args,init_Fetcher=True):
        print(2,db,query_id,extra_args)

        if db in SCRAPPABLE_DBS and not self.check_already_searched_memory(db,query_id):
            print(3,db, query_id, extra_args)

            if db == 'biocyc':        return    Protein_Fetcher_Biocyc(query_id,extra_args=extra_args,memory_storage=self.memory_storage,init_Fetcher=init_Fetcher)
            elif 'kegg' in db:        return    Protein_Fetcher_KEGG(query_id, extra_args=extra_args,memory_storage=self.memory_storage,init_Fetcher=init_Fetcher)
            elif db == 'hmdb':        return    Protein_Fetcher_HMDB(query_id, extra_args=extra_args,memory_storage= self.memory_storage,init_Fetcher=init_Fetcher)
            else:                     return    Global_Fetcher()


    def find_info(self, db,query_id,extra_args={}):
        if not query_id: return None,None
        fetcher=self.select_fetcher(db=db,query_id=query_id,extra_args=extra_args)
        if fetcher:
            self.add_to_already_tried_to_search(db,query_id)
            fetcher_protein=fetcher.get_protein()
            if fetcher_protein:
                #converge only occurs in the searchers- these are the global classes
                if self.is_valid_search_direction({'global'}):
                    fetcher.converge_protein_global()
                if self.is_valid_search_direction({'rpg','pg','crpg'}):
                    fetcher.converge_protein_rpg()
                if self.is_valid_search_direction({'gpr','pr','gprc','prc'}):
                    fetcher.converge_protein_gpr()

                return fetcher_protein,fetcher
            else:
                return None, None
        else:
            return None,None


    def run_searcher(self,bio_query,bio_db):
        """
        The only common ID in all DBs is Uniprot ID so this is the central point
        If a uniprot Id is provided we start with that, otherwise we go to the other dbs and extract the uniprot Id and go from there
        accepted dbs:
        enzyme_ec
        ko
        uniprot
        biocyc - uniprot id or enzyme ec or biocyc id
        kegg  - enzyme ec
        hmdb - uniprot id or hmdb id

        Ko and uniprot has 1:n but others 1:1 so, we assume that the 1:1 dbs input will be correct and will unite instances accordingly

        if only uniprot and/or ko we return several instances
        otherwise we find best match and return one or (if no unique match, several instances too)

        """
        print(f'STARTING PROTEIN SEARCHER {bio_query} in {bio_db}')
        args_to_search=[]
        if bio_db=='enzyme_ec' or is_ec(bio_query):
            if is_ec(bio_query,4):
                args_to_search.append(['enzyme_ec', bio_query])
        else:
            args_to_search.append([bio_db,bio_query])
        temp_inst=None
        #uniprot is only used in the begginning, since it has 1:n connections
        while args_to_search:
            #if args_to_search:
            current_arg = args_to_search.pop()
            db_arg= current_arg[0]
            id_arg= current_arg[1]
            if isinstance(id_arg,str) or isinstance(id_arg,int): id_arg={id_arg}
            for s_id_arg in id_arg:
                if db_arg == 'enzyme_ec':   temp_inst= self.get_protein_from_ec(s_id_arg)
                elif db_arg == 'uniprot':   temp_inst = self.get_proteins_from_uniprot(s_id_arg)
                elif db_arg=='kegg_ko':
                    ko_ecs = self.get_ec_from_ko(s_id_arg)
                    for ec in ko_ecs:
                        temp_inst = self.get_protein_from_ec(ec)
                        if temp_inst:
                            if not isinstance(temp_inst,list): temp_inst=[temp_inst]
                            for p in temp_inst:
                                if p:  p.set_detail('kegg_ko',s_id_arg)
                else:
                    if s_id_arg and not self.check_already_searched_memory(dict_input_key=db_arg, dict_input_value=s_id_arg):
                        temp_inst = self.find_protein(db=db_arg,query_id=s_id_arg)
                self.add_to_already_tried_to_search(db_arg, s_id_arg)
                if temp_inst:  self.add_to_args_to_search(temp_inst,args_to_search)
        return self.get_protein_match(bio_query,bio_db)



    def add_to_args_to_search(self,protein_instance_or_list,args_to_search):
        if not isinstance(protein_instance_or_list,list): protein_instance_or_list=[protein_instance_or_list]
        for protein_instance in protein_instance_or_list:
            if protein_instance:
                id_to_add = protein_instance.get_detail('uniprot')
                if id_to_add and not self.check_already_searched_memory(dict_input_key='uniprot', dict_input_value=id_to_add):
                    args_to_search.append(['uniprot', id_to_add])

                id_to_add = protein_instance.get_detail('enzyme_ec')
                if id_to_add and not self.check_already_searched_memory(dict_input_key='enzyme_ec', dict_input_value=id_to_add):
                    args_to_search.append(['enzyme_ec', id_to_add])

                id_to_add = protein_instance.get_detail('enzyme_ec')
                if id_to_add and not self.check_already_searched_memory(dict_input_key='kegg', dict_input_value=id_to_add):
                    args_to_search.append(['kegg', id_to_add])

                id_to_add = protein_instance.get_detail('biocyc')
                if id_to_add and not self.check_already_searched_memory(dict_input_key='biocyc', dict_input_value=id_to_add):
                    args_to_search.append(['biocyc', id_to_add])

                id_to_add = protein_instance.get_detail('hmdb')
                if id_to_add and not self.check_already_searched_memory(dict_input_key='hmdb', dict_input_value=id_to_add):
                    args_to_search.append(['hmdb', id_to_add])





    def get_proteins_from_uniprot(self,uniprot_id):
        """
        hmdb has 1:1 relation with uniprot (since its only for humans)
        biocyc and kegg have a 1:n relationship (n for families)
        brenda has no connection
        """
        res=[]
        for db in SCRAPPABLE_DBS:
            if uniprot_id and not self.check_already_searched_memory(dict_input_key=db+'_'+'uniprot',dict_input_value=uniprot_id):
                db_id = self.get_db_id_from_uniprot(db, uniprot_id=uniprot_id)
                if not (isinstance(db_id, list) and not isinstance(db_id,generator) and not isinstance(db_id,dict)): db_id = [db_id]
                for p_id in db_id:
                    temp_protein = self.find_protein(db=db,query_id= p_id)
                    #for kegg
                    if temp_protein:
                        if not isinstance(temp_protein,list): temp_protein=[temp_protein]
                        for p in temp_protein:
                            if p:
                                if not p.get_detail('uniprot_id'): p.set_detail('uniprot',uniprot_id)
                                res.append(p)
                    self.add_to_already_tried_to_search(db+'_'+'uniprot', uniprot_id)
        return res


    def get_protein_from_ec(self,enzyme_ec):
        """
        hmdb doesnt have EC info
        brenda and kegg have a 1:1 relationship
        biocyc has 1:n
        """
        res=[]
        protein_instance,kegg_prot,brenda_prot,biocyc_prot=None,None,None,None
        if is_ec(enzyme_ec,4):
            if enzyme_ec and not self.check_already_searched_memory(dict_input_key='kegg',dict_input_value=enzyme_ec):
                kegg_prot= self.find_protein(db='kegg',query_id=enzyme_ec)
                self.add_to_already_tried_to_search('kegg', enzyme_ec)
            if enzyme_ec and not self.check_already_searched_memory(dict_input_key='biocyc',dict_input_value=enzyme_ec):
                biocyc_prot=self.find_protein(db='biocyc',query_id=enzyme_ec)
                self.add_to_already_tried_to_search('biocyc', enzyme_ec)
        #depending which threadpool we are using (if we are using it), the pool will can be in different variables
        if isinstance(kegg_prot,list): res.extend(kegg_prot)
        else:
            if kegg_prot: res.append(kegg_prot)
        if isinstance(biocyc_prot,list): res.extend(biocyc_prot)
        else:
            if biocyc_prot: res.append(biocyc_prot)
        if isinstance(biocyc_prot,list): res.extend(biocyc_prot)
        else:
            if biocyc_prot: res.append(biocyc_prot)
        for p in res:
            if not protein_instance: protein_instance=p
            if p:
                protein_instance.unite_instances(p,always_unite=True)

        #if protein_instance:print('getting protein ec',enzyme_ec,protein_instance.get_possible_ids('enzyme_ec'))
        return protein_instance


    def get_db_id_from_uniprot(self,db,uniprot_id=None):
        if not uniprot_id: return None
        if db=='biocyc':
            db_id = self.get_db_id_from_uniprot_api_biocyc(uniprot_id)
            if not db_id:
                #goes directly to biocyc website
                db_id=self.get_db_id_from_uniprot_biocyc(uniprot_id)
            return db_id
        elif db=='kegg':
            return self.get_ecs_from_uniprot_kegg(uniprot_id)
        elif db=='hmdb':
            return self.get_db_id_from_uniprot_api_hmdb(uniprot_id)


    def get_db_id_from_uniprot_api_biocyc(self,uniprot_id):
        if uniprot_id:
            # using uniprot api
            params = {'from': 'ACC', 'to': 'BIOCYC_ID', 'format': 'tab', 'query': uniprot_id}
            url = 'https://www.uniprot.org/uploadlists/'
            webpage = self.get_with_fetcher(url, data=params)
            split_lines = webpage.split('\n')
            if len(split_lines)>1:
                res=[]
                split_lines=[i.split('\t') for i in split_lines[1:]]
                for line in split_lines:
                    if len(line)>1:
                        uniprot_id,biocyc_id=line[0],line[1].split('MetaCyc:')[1]
                        res.append(biocyc_id)
                return res

    def get_db_id_from_uniprot_api_hmdb(self,uniprot_id):
        url='http://www.hmdb.ca/unearth/q?utf8=%E2%9C%93&query='+uniprot_id+'&searcher=proteins&button='
        webpage=self.get_with_fetcher(url)
        if not webpage: return None
        soup = BeautifulSoup(webpage, 'lxml')
        res= soup.find('a',{'href':re.compile('.*uniprot/'+uniprot_id)})
        if res:
            hmdb_id=res.parent.parent.a.text
            return hmdb_id




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
                print('TOO MANY RESULTS BIOCYC')
                return
            gene_id,protein_id,organism= res[0].find_all('a')
            gene_id = re.search('&id=.*',gene_id['href']).group()[4:]
            return gene_id

    def get_proteins_from_ec_biocyc(self,enzyme_ec):
        fetcher=Protein_Fetcher_Biocyc(enzyme_ec)
        protein_and_genes= fetcher.convergence_args['protein_and_genes']
        res=[]
        for p in protein_and_genes:
            res.append(p['protein_id'])
        return res


    def get_ko_from_gene_kegg(self,gene_id):
        webpage = self.get_with_fetcher(api_kegg=True, url=gene_id, database='ko',type_search='link')
        if webpage:
            res = []
            for line in webpage:
                kegg_id, ko = line.split('\t')
                if kegg_id == gene_id:
                    ko = ko.split(':')[1]
                    res.append(ko)
            return res

    def get_ec_from_gene_kegg(self,gene_id):
        webpage = self.get_with_fetcher(api_kegg=True, url=gene_id, database='ec',type_search='link')
        if webpage:
            res=[]
            for line in webpage:
                kegg_id,ec= line.split('\t')
                if kegg_id==gene_id:
                    ec=ec.split(':')[1]
                    if is_ec(ec,4): res.append(ec)
            return res

    def get_ecs_from_uniprot_kegg(self,uniprot_id):
        if uniprot_id:
            #using uniprot api
            params = {'from': 'ACC', 'to': 'KEGG_ID', 'format': 'tab', 'query': uniprot_id}
            url = 'https://www.uniprot.org/uploadlists/'
            webpage = self.get_with_fetcher(url,data=params)
            split_lines = webpage.split('\n')
            if len(split_lines)>1:
                kegg_id=split_lines[1]
                if kegg_id:
                    kegg_id=kegg_id.split('\t')[1]
                    ecs=self.get_ec_from_gene_kegg(kegg_id)
                    if not ecs:
                        ko= self.get_ko_from_gene_kegg(kegg_id)
                        ecs= self.get_ec_from_ko(ko)
                    return ecs


    def get_ec_from_ko(self,ko):
        webpage = self.get_with_fetcher(api_kegg=True, url='ko:' + ko, database='ko')
        if webpage:
            kegg_id = webpage[0].split('\t')[1]
            ec_pattern = re.compile('\[EC:(\d+\.){2,3}((\d+)|-)(\s(\d+\.){2,3}((\d+)|-))*\]')
            enzyme_ec = re.search(ec_pattern, kegg_id)
            if enzyme_ec: return enzyme_ec.group()[4:-1].split()
            else: return []
        return []


if __name__ == '__main__':

    searcher=Protein_Searcher(search_direction='pr',output_folder='/home/pedroq/PycharmProjects/DRAX/test/test2/testout')
    p1=searcher.run_searcher('1.1.1.178','biocyc')
    #p1.get_all_info()
    searcher.output_results()