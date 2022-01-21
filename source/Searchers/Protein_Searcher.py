
from source.Searchers.Gene_Searcher import Gene_Searcher
from source.Searchers.Global_Searcher import *
from source.Biological_Components.Protein import Protein
from source.Fetchers.Protein_Fetchers.Protein_Fetcher_Metacyc import Protein_Fetcher_Metacyc
from source.Fetchers.Protein_Fetchers.Protein_Fetcher_KEGG import Protein_Fetcher_KEGG
from source.Fetchers.Protein_Fetchers.Protein_Fetcher_HMDB import Protein_Fetcher_HMDB
from source.Fetchers.Protein_Fetchers.Protein_Fetcher_Rhea import Protein_Fetcher_Rhea
from source.Utils.util import is_ec


from types import GeneratorType as generator

class Protein_Searcher(Global_Searcher):
    def __init__(self,memory_storage=None,search_mode='',db_name=None,wanted_org_kegg_codes=[],output_folder=None,politeness_timer=10):
        Global_Searcher.__init__(self,memory_storage,search_mode,
                                 db_name=db_name,wanted_org_kegg_codes=wanted_org_kegg_codes,output_folder=output_folder,politeness_timer=politeness_timer)
        self.original_searcher=get_instance_type(self)

    def find_protein(self,db,query_id=None,extra_args={},convergence_search=False):
        if db in SCRAPPABLE_DBS:
            fetcher_protein,fetcher= self.find_info(db, query_id, extra_args,convergence_search=convergence_search)
            return fetcher_protein

    def select_fetcher(self,db,query_id,extra_args):
        if db in SCRAPPABLE_DBS and not self.check_already_searched_memory(db,query_id):
            if db == 'metacyc':       return   Protein_Fetcher_Metacyc(query_id,extra_args=extra_args,memory_storage=self.memory_storage)
            elif 'kegg' in db:        return    Protein_Fetcher_KEGG(query_id, extra_args=extra_args,memory_storage=self.memory_storage)
            elif db == 'hmdb':        return    Protein_Fetcher_HMDB(query_id, extra_args=extra_args,memory_storage= self.memory_storage)
            elif db == 'rhea':        return    Protein_Fetcher_Rhea(query_id, extra_args=extra_args,memory_storage= self.memory_storage)
            else:                     return    Global_Fetcher()



    def find_info(self, db,query_id,extra_args={},convergence_search=False):
        if not query_id: return None,None
        fetcher=self.select_fetcher(db=db,query_id=query_id,extra_args=extra_args)
        self.add_to_already_tried_to_search(db, query_id)
        if fetcher:
            fetcher_protein=fetcher.get_protein()
            if fetcher_protein:
                if convergence_search:
                    if db=='metacyc': fetcher.converge_protein_to_protein()
                if self.is_valid_search_mode({'global'}):
                    fetcher.converge_protein_global()
                if self.is_valid_search_mode({'rpg','pg','crpg'}):
                    fetcher.converge_protein_rpg()
                if self.is_valid_search_mode({'gpr','pr','gprc','prc'}):
                    fetcher.converge_protein_gpr()
                return fetcher_protein,fetcher
            else:
                return None, None
        else:
            if self.check_already_searched_memory(db, query_id):
                print('Already searched', db, query_id)
                return self.get_protein_match(bio_query=query_id, bio_db=db),None
            else:
                return None,None


    def run_searcher(self,bio_query,bio_db,convergence_search=False,extra_args={}):
        """
        The only common ID in all DBs is Uniprot ID so this is the central point
        If a uniprot Id is provided we start with that, otherwise we go to the other dbs and extract the uniprot Id and go from there
        accepted dbs:
        enzyme_ec
        ko
        uniprot
        metacyc - uniprot id or enzyme ec or metacyc id
        kegg  - enzyme ec
        hmdb - uniprot id or hmdb id

        Ko and uniprot has 1:n but others 1:1 so, we assume that the 1:1 dbs input will be correct and will unite instances accordingly

        """
        print(f'Starting protein search {bio_query} in {bio_db}')
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
                if db_arg == 'enzyme_ec':   temp_inst= self.get_protein_from_ec(s_id_arg,convergence_search=convergence_search)
                elif db_arg == 'uniprot':   temp_inst = self.get_proteins_from_uniprot(s_id_arg,convergence_search=convergence_search)
                elif db_arg=='kegg_ko':
                    ko_ecs = self.get_ec_from_ko(s_id_arg)
                    for ec in ko_ecs:
                        temp_inst = self.get_protein_from_ec(ec,convergence_search=convergence_search)
                        if temp_inst:
                            if not isinstance(temp_inst,list): temp_inst=[temp_inst]
                            for p in temp_inst:
                                if p:  p.set_detail('kegg_ko',s_id_arg)
                else:
                    if s_id_arg and not self.check_already_searched_memory(dict_input_key=db_arg, dict_input_value=s_id_arg):
                        temp_inst = self.find_protein(db=db_arg,query_id=s_id_arg,convergence_search=convergence_search,extra_args=extra_args)
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

                id_to_add = protein_instance.get_detail('metacyc')
                if id_to_add and not self.check_already_searched_memory(dict_input_key='metacyc', dict_input_value=id_to_add):
                    args_to_search.append(['metacyc', id_to_add])

                id_to_add = protein_instance.get_detail('hmdb')
                if id_to_add and not self.check_already_searched_memory(dict_input_key='hmdb', dict_input_value=id_to_add):
                    args_to_search.append(['hmdb', id_to_add])





    def get_proteins_from_uniprot(self,uniprot_id,convergence_search=False,extra_args={}):
        """
        hmdb has 1:1 relation with uniprot (since its only for humans)
        metacyc and kegg have a 1:n relationship (n for families)
        brenda has no connection
        """
        res=[]
        for db in SCRAPPABLE_DBS:
            if uniprot_id and not self.check_already_searched_memory(dict_input_key=f'{db}_uniprot',dict_input_value=uniprot_id):
                db_id = self.get_db_id_from_uniprot(db, uniprot_id=uniprot_id)
                if not db_id: db_id=[]
                if isinstance(db_id,str) or isinstance(db_id,int) or isinstance(db_id,float): db_id = [db_id]
                for p_id in db_id:
                    temp_protein = self.find_protein(db=db,query_id= p_id,convergence_search=convergence_search,extra_args=extra_args)
                    #for kegg
                    if temp_protein:
                        if not isinstance(temp_protein,list): temp_protein=[temp_protein]
                        for p in temp_protein:
                            if p:
                                if not p.get_detail('uniprot_id'): p.set_detail('uniprot',uniprot_id)
                                res.append(p)
                    self.add_to_already_tried_to_search(f'{db}_uniprot', uniprot_id)
        return res


    def get_protein_from_ec(self,enzyme_ec,convergence_search=False,extra_args={}):
        """
        hmdb doesnt have EC info
        brenda and kegg have a 1:1 relationship
        metacyc has 1:n
        """
        res=[]
        protein_instance,kegg_prot,brenda_prot,metacyc_prot=None,None,None,None
        if is_ec(enzyme_ec,4):
            if enzyme_ec and not self.check_already_searched_memory(dict_input_key='kegg',dict_input_value=enzyme_ec):
                kegg_prot= self.find_protein(db='kegg',query_id=enzyme_ec,convergence_search=convergence_search,extra_args=extra_args)
                self.add_to_already_tried_to_search('kegg', enzyme_ec)
            if enzyme_ec and not self.check_already_searched_memory(dict_input_key='metacyc',dict_input_value=enzyme_ec):
                metacyc_prot=self.find_protein(db='metacyc',query_id=enzyme_ec,convergence_search=convergence_search,extra_args=extra_args)
                self.add_to_already_tried_to_search('metacyc', enzyme_ec)

        if isinstance(kegg_prot,list): res.extend(kegg_prot)
        else:
            if kegg_prot: res.append(kegg_prot)
        if isinstance(metacyc_prot,list): res.extend(metacyc_prot)
        else:
            if metacyc_prot: res.append(metacyc_prot)

        for p in res:
            if not protein_instance: protein_instance=p
            if p:
                protein_instance.unite_instances(p,always_unite=True)

        #if protein_instance:print('getting protein ec',enzyme_ec,protein_instance.get_possible_ids('enzyme_ec'))
        return protein_instance


    def get_db_id_from_uniprot(self,db,uniprot_id=None):
        res=[]
        if 'uniprot' not in SCRAPPABLE_DBS: return res
        if not uniprot_id: return res
        if db=='metacyc':
            return self.fetch_metacyc_from_uniprot(uniprot_id)
        elif db=='kegg':
            return self.get_ecs_from_uniprot_kegg(uniprot_id)
        elif db=='hmdb':
            return self.get_db_id_from_uniprot_api_hmdb(uniprot_id)


    def get_db_id_from_uniprot_api_hmdb(self,uniprot_id):
        url=f'http://www.hmdb.ca/unearth/q?utf8=%E2%9C%93&query={uniprot_id}&searcher=proteins&button='
        webpage=self.get_with_fetcher(url)
        if not webpage: return None
        soup = BeautifulSoup(webpage, 'lxml')
        res= soup.find('a',{'href':re.compile('.*uniprot/'+uniprot_id)})
        if res:
            hmdb_id=res.parent.parent.a.text
            return hmdb_id

    def get_ko_from_gene_kegg(self,gene_id):
        webpage = self.get_with_fetcher(api_kegg=True, url=gene_id, database='ko',type_search='link')
        res = []
        if webpage:
            for line in webpage:
                kegg_id, ko = line.split('\t')
                if kegg_id == gene_id:
                    ko = ko.split(':')[1]
                    res.append(ko)
        return res


    def get_ec_from_gene_kegg(self,gene_id):
        webpage = self.get_with_fetcher(api_kegg=True, url=gene_id, database='ec',type_search='link')
        res = []
        if webpage:
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
                        ecs=set()
                        kos= self.get_ko_from_gene_kegg(kegg_id)
                        for ko in kos:
                            ecs.update(self.get_ec_from_ko(ko))
                    return ecs


    def get_ec_from_ko(self,ko):
        if 'kegg' in SCRAPPABLE_DBS:
            webpage = self.get_with_fetcher(api_kegg=True, url='ko:' + ko, database='ko')
            if webpage:
                kegg_id = webpage[0].split('\t')[1]
                ec_pattern = re.compile('\[EC:(\d+\.){2,3}((\d+)|-)(\s(\d+\.){2,3}((\d+)|-))*\]')
                enzyme_ec = re.search(ec_pattern, kegg_id)
                if enzyme_ec: return enzyme_ec.group()[4:-1].split()
                else: return []
        return []


if __name__ == '__main__':

    searcher=Protein_Searcher(search_mode={''},politeness_timer=2)
    print(searcher.get_ko_from_gene_kegg('hsa:54165'))