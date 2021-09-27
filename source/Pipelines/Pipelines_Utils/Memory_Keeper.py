


from source.Utils.Web_Connector import Web_Connector
from source.Biological_Components.Info_Keepers.Item_Set import Item_Set
from source.Utils.util import number_of_nones_dict,get_instance_type
from source.Biological_Components.Compound import Compound
from source.Biological_Components.Gene import Gene
from source.Biological_Components.Protein import Protein
from source.Biological_Components.Reaction import Reaction

class Memory_Keeper():
    def __init__(self,db_name=None,politeness_timer=10):
        self.fetcher_test=False
        self.fetcher_omit_error_messages=False
        self.searched={}
        self.politeness_timer=politeness_timer


    def check_already_searched_memory(self,dict_input_key,dict_input_value):
        #so that we dont add empties
        if not dict_input_value: return True
        for searched_dict in self.already_tried_to_search:
            key_searched_dict=list(searched_dict.keys())[0]
            value_searched_dict=searched_dict[key_searched_dict]
            if key_searched_dict==dict_input_key:
                if value_searched_dict==dict_input_value:
                    return True
        return False

    def add_to_already_tried_to_search(self,dict_input_key,dict_input_value):
        if dict_input_value:
            if not self.check_already_searched_memory(dict_input_key,dict_input_value):
                self.already_tried_to_search.append({dict_input_key:dict_input_value})

    def setup_memory_lists(self):
        self.reactions=Item_Set()
        self.proteins=Item_Set()
        self.genes=Item_Set()
        self.compounds=Item_Set()
        self.already_tried_to_search=[]

    def setup_fetchers(self):
        #so that startup is quicker we dont init all the fetchers, we just do launch them the first time we need them
        self.fetcher_biocyc = None
        self.fetcher_kegg = None
        self.fetcher_ncbi = None
        self.fetcher_uniprot = None
        # since non-specific fetcher accesses several databases, it's not as important to rotate proxies that often.
        self.fetcher_others = None




    def initialize_fetcher(self,url_or_db):
        if not url_or_db:
            self.fetcher_others =Web_Connector(test=self.fetcher_test, omit_error_messages=self.fetcher_omit_error_message,politeness_timer=self.politeness_timer)
            return self.fetcher_others

        elif 'kegg' in url_or_db.lower() or 'genome.jp' in url_or_db.lower():
            self.fetcher_kegg =Web_Connector(test=self.fetcher_test, omit_error_messages=self.fetcher_omit_error_messages,politeness_timer=self.politeness_timer)
            return self.fetcher_kegg

        elif 'biocyc' in url_or_db.lower():
            self.fetcher_biocyc =Web_Connector(test=self.fetcher_test, omit_error_messages=self.fetcher_omit_error_messages,politeness_timer=self.politeness_timer)
            return self.fetcher_biocyc

        elif 'ncbi' in url_or_db.lower():
            self.fetcher_ncbi =Web_Connector(test=self.fetcher_test, omit_error_messages=self.fetcher_omit_error_messages,politeness_timer=self.politeness_timer)
            return self.fetcher_ncbi
        elif 'uniprot' in url_or_db.lower():
            self.fetcher_uniprot =Web_Connector(test=self.fetcher_test, omit_error_messages=self.fetcher_omit_error_messages,politeness_timer=self.politeness_timer)
            return self.fetcher_uniprot

        else:
            self.fetcher_others=Web_Connector(test=self.fetcher_test, omit_error_messages=self.fetcher_omit_error_messages,politeness_timer=self.politeness_timer)
            return self.fetcher_others

    def borrow_fetchers(self,memory_storage):
        #print('Borrowing Fetchers!')
        self.fetcher_biocyc = memory_storage.fetcher_biocyc
        self.fetcher_kegg = memory_storage.fetcher_kegg
        self.fetcher_ncbi = memory_storage.fetcher_ncbi
        self.fetcher_uniprot = memory_storage.fetcher_uniprot
        self.fetcher_others = memory_storage.fetcher_others

    def borrow_memory_lists(self,memory_storage):
        self.compounds= memory_storage.compounds
        self.genes = memory_storage.genes
        self.proteins = memory_storage.proteins
        self.reactions =memory_storage.reactions
        self.already_tried_to_search =memory_storage.already_tried_to_search
        self.wanted_org_kegg_codes=memory_storage.wanted_org_kegg_codes

    def get_db_fetcher(self,url_or_db=None):
        if not url_or_db:                                                       fetcher= self.fetcher_others
        elif 'kegg' in url_or_db.lower() or 'genome.jp' in url_or_db.lower():   fetcher= self.fetcher_kegg
        elif 'biocyc' in url_or_db.lower():                                     fetcher= self.fetcher_biocyc
        elif 'ncbi' in url_or_db.lower():                                       fetcher= self.fetcher_ncbi
        elif 'uniprot' in url_or_db.lower():                                    fetcher= self.fetcher_uniprot
        else:
            fetcher= self.fetcher_others
        if not fetcher: fetcher= self.initialize_fetcher(url_or_db)
        return fetcher


    def get_with_fetcher(self,url,selenium=False,api_kegg=False,scripts=[],data=None,original_response=False,xpath=[],timer=1,ids_to_load=[],database=None, type_search='find'):
        print(f'Connecting to {url}')
        if selenium:
            res= self.get_db_fetcher(url).get_driver_selenium(url, script=scripts,xpath=xpath,timer=timer,original_response=original_response,ids_to_load=ids_to_load)
        elif api_kegg:
            res= self.get_db_fetcher('kegg').api_KEGG(url, database, type_search)
        else:
            res= self.get_db_fetcher(url).try_until_catch(url,data=data,original_response=original_response)
            #sometimes the page wont open even though its working so we try to go there with selenium
            if not res and not data:
                res = self.get_db_fetcher(url).get_driver_selenium(url,original_response=original_response,exceptional_try_limit=10)
        return res

    #From fetchers
    #Adders

    def add_instance(self,instance_to_add):
        node_type = get_instance_type(instance_to_add)
        if node_type == 'Compound':
            self.add_compound(instance_to_add)
        elif node_type == 'Gene':
            self.add_gene(instance_to_add)
        elif node_type == 'Protein':
            self.add_protein(instance_to_add)
        elif node_type == 'Reaction':
            self.add_reaction(instance_to_add)


    def add_protein(self,Protein_instance):
        self.proteins.set_item_set(Protein_instance)

    def add_gene(self, Gene_instance):
        self.genes.set_item_set(Gene_instance)

    def add_reaction(self,Reaction_instance):
        self.reactions.set_item_set(Reaction_instance)

    def add_compound(self,Compound_instance):
        self.compounds.set_item_set(Compound_instance)




    def get_all_lists(self):
        yield   self.get_compounds_all()
        yield   self.get_genes_all()
        yield   self.get_proteins_all()
        yield   self.get_reactions_all()


    #Instance Getter
    def get_biological_instance(self,bio_type,bio_query,bio_db=None):
        if not bio_query: return None
        if      bio_type=='proteins' :
            bio_list=self.get_proteins_all()
        elif    bio_type=='genes' :
            bio_list=self.get_genes_all()
        elif    bio_type=='reactions' :
            bio_list=self.get_reactions_all()
        elif    bio_type=='compounds' :
            bio_list=self.get_compounds_all()
        else:
            bio_list=[]
        for bio_entity in bio_list:
            if bio_db:
                if bio_type=='proteins':        temp_bio_entity=Protein({bio_db:bio_query})
                elif bio_type=='genes':         temp_bio_entity=Gene({bio_db:bio_query})
                elif bio_type=='reactions':     temp_bio_entity=Reaction({bio_db:bio_query})
                elif bio_type=='compounds':     temp_bio_entity=Compound({bio_db:bio_query})
                if bio_entity.is_match_instances(temp_bio_entity,threshold_for_match=1):   return bio_entity
            else:
                if bio_entity.is_match_instances(bio_query):   return bio_entity
        return None

    def merge_all_instances(self):
        for yielder in [
            self.get_compounds_all(),
            self.get_genes_all(),
            self.get_proteins_all(),
            self.get_reactions_all(),
        ]:
            current_set=set()
            for current_inst in yielder:
                current_set.add(current_inst)
            current_length=len(current_set)
            print(current_set)



    #Getters all
    def get_proteins_all(self):
        return self.proteins.get_item_set()

    def get_genes_all(self):
        return self.genes.get_item_set()

    def get_reactions_all(self):
        return self.reactions.get_item_set()

    def get_compounds_all(self):
        return self.compounds.get_item_set()

    #Setters
    def set_proteins(self, list_proteins):
        self.proteins.set_item_set(list_proteins)

    def set_genes(self, list_genes):
        self.genes.set_item_set(list_genes)

    def set_reactions(self, list_reactions):
        self.reactions.set_item_set(list_reactions)

    def set_compounds(self, list_compounds):
        self.compounds.set_item_set(list_compounds)


if __name__ == '__main__':
    m=Memory_Keeper()
    m.setup_fetchers()
