from source.Fetchers.Fetchers_Utils.Global_Fetcher import *
from types import GeneratorType as generator
from source.Utils.util import remove_inchi_key_equal



class Compound_Fetcher(Global_Fetcher):
    def __init__(self,compound_id,memory_storage=None):
        Global_Fetcher.__init__(self)
        self.compound_id=compound_id
        self.memory_storage=memory_storage
        self.db= None
        self.compound=None
        #if no memory_storage is present from one of the pipelines or previous fetchers we assign it one and initialize the memory
        if not self.memory_storage:
            from source.Pipelines.Searchers.Compound_Searcher import Compound_Searcher
            self.memory_storage = Compound_Searcher()

    def add_compound(self):
        if self.get_compound():
            match= self.get_compound_match()
            if match:
                if match is not self.get_compound():
                    match.unite_instances(self.get_compound(),always_unite=True)
                    self.set_compound(match)
            else:
                self.memory_storage.add_compound(self.get_compound())

    def get_compound(self):
        return self.compound

    def set_compound(self,match_instance):
        self.compound=match_instance

    def get_compound_match(self):
        if self.get_compound():
            return self.memory_storage.get_biological_instance('compounds',self.get_compound())
        else:
            return self.memory_storage.get_biological_instance('compounds',self.compound_id,self.db)

    def get_wanted_info(self):
        return ['hmdb',
                'kegg',
                'chemspider',
                'chebi',
                'pubchem_cid',
                'bigg',
                'biocyc',
                'seed',
                'synonyms',
                'brenda',
                'inchi_key',
                'drugbank',
                'pdb',
                'pubchem_sid',
                ]

    def remove_unwanted_info(self,dict_to_change):
        wanted_info=self.get_wanted_info()
        res={i.lower():dict_to_change[i] for i in dict_to_change if dict_to_change[i] and i.lower() in wanted_info}
        for detail in res:
            if isinstance(res[detail],str):
                res[detail]={res[detail]}
        for detail in res:
            res[detail]={i for i in res[detail] if not re.search('not available',i.lower())}
        return res
