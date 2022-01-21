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


class Input_Generator(Rhea_SQLITE_Connector,Metacyc_SQLITE_Connector):
    def __init__(self,output_path):
        Rhea_SQLITE_Connector.__init__(self)
        Metacyc_SQLITE_Connector.__init__(self)
        self.output_path=output_path

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

    def yield_all_lines(self,pickle_path):
        for db_id in self.rhea_fetch_all_reactions():
            yield f'{db_id}\trhea\treaction\tglobal\n'
        for db_id in self.metacyc_fetch_all_proteins():
            yield f'{db_id}\tmetacyc\tprotein\tglobal\n'
        ecs_kos=self.load_metrics(pickle_path)
        for db in ecs_kos:
            for db_id in ecs_kos[db]:
                yield f'{db_id}\t{db}\tprotein\tglobal\n'

    def generate_universal_input(self,pickle_path):
        if not os.path.exists(pickle_path):
            self.generate_pickle_ecs_kos(pickle_path, ec_json, ko_json)
        with open(self.output_path,'w+') as file:
            for drax_line in self.yield_all_lines(pickle_path):
                file.write(drax_line)

    def parse_tsv(self,input_file):
        res = {}
        with open(input_file) as file:
            file.readline()
            for line in file:
                line = line.strip('\n')
                line = line.split('\t')
                for line_tab in line:
                    if ':' in line_tab:
                        db = line_tab.split(':')[0]
                        if db not in ['description']:
                            if db not in res: res[db]=set()
                            # to avoid bad splitting when dealing with descriptions
                            annot = line_tab[len(db) + 1:]
                            res[db].add(annot)
        return res

    #this is not general enough
    def generate_input_from_annotations(self,input_annotations):
        parsed_annotations=self.parse_tsv(input_annotations)
        with open(self.output_path, 'w+') as file:
            for id_type in parsed_annotations:
                for annot in parsed_annotations[id_type]:
                    if id_type in ['kegg_ko','metacyc','cog','kegg_module']:
                        file.write(f'{annot}\t{id_type}\tprotein\tprc\n')
                    elif id_type in ['kegg_reaction']:
                        file.write(f'{annot}\t{id_type}\treaction\tcrp\n')
                    else:
                        pass



if __name__ == '__main__':
    pickle_path='ecs_kos.pickle'
    # from https://www.kegg.jp/brite/ko01000
    ec_json = 'ko01000.json'
    # from https://www.genome.jp/kegg-bin/show_brite?ko00001.keg
    ko_json = 'ko00001.json'


    output_path=f'universal_input.tsv'
    mantis_tsv=f'/home/pedroq/Desktop/test_mantis/test2/consensus_annotation.tsv'
    s=Input_Generator(output_path)
    s.generate_universal_input(pickle_path)


