import re
import sys
import os
import pickle
import json

from source.Utils.Rhea_SQLITE_Connector import Rhea_SQLITE_Connector
from source.Utils.Metacyc_SQLITE_Connector import Metacyc_SQLITE_Connector


if sys.platform.startswith('win'):
    SPLITTER = '\\'
else:
    SPLITTER = '/'

DRAX_FOLDER = os.path.abspath(os.path.dirname(__file__)).split(SPLITTER)[0:-2]
DRAX_FOLDER = SPLITTER.join(DRAX_FOLDER) + SPLITTER
RESOURCES_FOLDER=f'{DRAX_FOLDER}Resources{SPLITTER}'


class Universal_input(Rhea_SQLITE_Connector,Metacyc_SQLITE_Connector):
    def __init__(self,output_path,pickle_path):
        Rhea_SQLITE_Connector.__init__(self)
        Metacyc_SQLITE_Connector.__init__(self)
        self.output_path=output_path
        self.ecs_kos_pickle=pickle_path
        if not os.path.exists(self.ecs_kos_pickle):
            self.generate_pickle_ecs_kos(pickle_path, ec_json, ko_json)

    def save_metrics(self,pickle_path, to_pickle):
        with open(pickle_path, 'wb') as handle:
            pickle.dump(to_pickle, handle)

    def load_metrics(self,pickle_path):
        if os.path.exists(pickle_path):
            with open(pickle_path, 'rb') as handle:
                pickled_results = pickle.load(handle)
                return pickled_results

    def get_ecs(self,file_path):
        res = set()
        with open(file_path) as file:
            json_modules = json.load(file)['children'][0]['children']
            for main_path in json_modules:
                sub_pathways = main_path['children']
                for sub_path in sub_pathways:
                    modules = sub_path['children']
                    for module in modules:
                        sub_sub_description = module['name']
                        sub_sub_ec = sub_sub_description.split()[0]
                        if '-' not in sub_sub_ec:
                            res.add(sub_sub_ec)
        return res

    def get_kos(self,file_path):
        res = set()
        with open(file_path) as file:
            json_modules = json.load(file)['children'][0]['children']
            for main_path in json_modules:
                sub_pathways = main_path['children']
                for sub_path in sub_pathways:
                    modules = sub_path['children']
                    for module in modules:
                        sub_sub_description = module['name']
                        sub_sub_ec = sub_sub_description.split()[0]
                        res.add(sub_sub_ec)
        return res

    def generate_pickle_ecs_kos(self,pickle_path, ec_json, ko_json):
        all_ecs = self.get_ecs(ec_json)
        all_kos = self.get_kos(ko_json)
        res = {'enzyme_ec': all_ecs, 'kegg_ko': all_kos}
        self.save_metrics(pickle_path, res)

    def yield_all_lines(self):
        for db_id in self.rhea_fetch_all_reactions():
            yield f'{db_id}\trhea\treaction\tglobal\n'
        for db_id in self.metacyc_fetch_all_proteins():
            yield f'{db_id}\tmetacyc\tprotein\tglobal\n'
        ecs_kos=self.load_metrics(self.ecs_kos_pickle)
        for db in ecs_kos:
            for db_id in ecs_kos[db]:
                yield f'{db_id}\t{db}\tprotein\tglobal\n'

    def generate_drax_input(self):
        with open(self.output_path,'w+') as file:
            for drax_line in self.yield_all_lines():
                file.write(drax_line)


if __name__ == '__main__':
    pickle_path='ecs_kos.pickle'
    # from https://www.kegg.jp/brite/ko01000
    ec_json = 'ko01000.json'
    # from https://www.genome.jp/kegg-bin/show_brite?ko00001.keg
    ko_json = 'ko00001.json'


    output_path=f'/home/pedroq/Desktop/test_drax/universal_input.tsv'
    s=Universal_input(output_path,pickle_path)
    s.generate_drax_input()


