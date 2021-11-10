
from source.Biological_Components.Compound import Compound
from source.Biological_Components.Gene import Gene
from source.Biological_Components.Protein import Protein
from source.Biological_Components.Reaction import Reaction
from source.Biological_Components.Biological_Components_Utils.Unique_details import generate_details_dict,get_details_list,get_unique_details
from source.Utils.util import strip_tags,number_of_nones_dict,get_instance_type,find_ecs,is_ec,find_path,regex_escape,fix_html_sign



#External modules

import re
from unicodedata import normalize
from bs4 import BeautifulSoup

class Global_Fetcher():
    def __init__(self):
        pass

    def get_with_fetcher(self,url,selenium=False,api_kegg=False,scripts=[],data=None,original_response=False,xpath=[],timer=1,ids_to_load=[],database=None, type_search='find'):
        return self.memory_storage.get_with_fetcher(url=url,selenium=selenium,api_kegg=api_kegg,scripts=scripts,
                                                    data=data,original_response=original_response,xpath=xpath,
                                                    timer=timer,ids_to_load=ids_to_load,database=database,
                                                    type_search=type_search)

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
                    res[c]+=' : '+t
                    c += 1
                else:
                    res.append([])
                    res[c]=t
        else:
            res=normalize('NFKD', to_search.text).strip('\n')
        return res



if __name__ == '__main__':
    pass