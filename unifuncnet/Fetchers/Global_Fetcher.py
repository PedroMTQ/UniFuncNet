
from unifuncnet.Biological_Components.Compound import Compound
from unifuncnet.Biological_Components.Gene import Gene
from unifuncnet.Biological_Components.Protein import Protein
from unifuncnet.Biological_Components.Reaction import Reaction
from unifuncnet.Biological_Components.utils.Unique_details import generate_details_dict,get_details_list,get_unique_details
from unifuncnet.Utils.util import strip_tags,number_of_nones_dict,get_instance_type,find_ecs,is_ec,find_path,regex_escape,fix_html_sign


#External modules

import re
from unicodedata import normalize
from bs4 import BeautifulSoup

class Global_Fetcher():
    def __init__(self,memory_storage=None):
        self.memory_storage=memory_storage

    def get_with_fetcher(self,url,api_kegg=False,data=None,original_response=False,database=None, type_search='find',kegg_option=None):
        return self.memory_storage.get_with_fetcher(url=url,api_kegg=api_kegg,data=data,
                                                    original_response=original_response,database=database,
                                                    type_search=type_search,kegg_option=kegg_option)


    def fetch_metacyc_intermediate_rxn_ids(self,wanted_id):
        return self.memory_storage.fetch_metacyc_intermediate_rxn_ids(wanted_id)

    def fetch_metacyc_rxn_from_cpd(self,wanted_id):
        return self.memory_storage.fetch_metacyc_rxn_from_cpd(wanted_id)

    def fetch_metacyc_rxn_from_ec(self,wanted_id):
        return self.memory_storage.fetch_metacyc_rxn_from_ec(wanted_id)

    def fetch_metacyc_from_uniprot(self,wanted_id):
        return self.memory_storage.fetch_metacyc_from_uniprot(wanted_id)

    def fetch_metacyc_id_info(self,wanted_id,table_name):
        return self.memory_storage.fetch_metacyc_id_info(wanted_id,table_name)

    def fetch_metacyc_derivatives(self,wanted_id,table_name):
        return self.memory_storage.fetch_metacyc_derivatives(wanted_id,table_name)

    def fetch_rhea_id_info(self,wanted_id):
        return self.memory_storage.fetch_rhea_id_info(wanted_id)

    def fetch_rhea_from_id(self,id_type,input_id):
        return self.memory_storage.fetch_rhea_from_id(id_type,input_id)

    def fetch_reactions_rhea_from_chebi(self,wanted_id):
        return self.memory_storage.fetch_reactions_rhea_from_chebi(wanted_id)

    def fetch_chebi_id_info(self,wanted_id):
        return self.memory_storage.fetch_chebi_id_info(wanted_id)

    def get_reaction(self):
        return

    def get_compound(self):
        return

    def get_gene(self):
        return

    def get_protein(self):
        return


    def get_db_fetcher(self,url_or_db=None):
        return self.memory_storage.get_db_fetcher(url_or_db)

    def get_gene_match(self,bio_query,bio_db=None):
        return self.memory_storage.get_biological_instance('genes',bio_query,bio_db)

    def get_protein_match(self,bio_query,bio_db=None):
        return self.memory_storage.get_biological_instance('proteins',bio_query,bio_db)

    def get_reaction_match(self,bio_query,bio_db=None):
        return self.memory_storage.get_biological_instance('reactions',bio_query,bio_db)

    def get_compound_match(self,bio_query,bio_db=None):
        return self.memory_storage.get_biological_instance('compounds',bio_query,bio_db)



    def textbox_KEGG(self,ec_soup, box,to_split=False):
        to_search = ec_soup.find('nobr', text=re.compile(regex_escape(box)))
        if to_split: res = []
        else: res=''
        if not to_search: return res
        to_search = to_search.findNext()
        if to_split:
            c=0
            td_items=to_search.find_all('td')
            for i in range(len(td_items)):
                t=normalize('NFKD',td_items[i].text).strip()
                if  i%2 != 0:
                    res[c]+=f' : {t}'
                    c += 1
                else:
                    res.append([])
                    res[c]=t
        else:
            res=normalize('NFKD', to_search.text).strip('\n')
        return res



if __name__ == '__main__':
    pass