import re
import sys
import os
import pickle
import json

from unifuncnet.utils.rhea_sqlite_connector import RheaSqliteConnector
from unifuncnet.utils.metacyc_sqlite_connector import MetacycSqliteConnector

if sys.platform.startswith('win'):
    SPLITTER = '\\'
else:
    SPLITTER = '/'


class InputGenerator(RheaSqliteConnector, MetacycSqliteConnector):
    def __init__(self, output_path):
        RheaSqliteConnector.__init__(self)
        MetacycSqliteConnector.__init__(self)
        self.output_path = output_path

    def save_metrics(self, pickle_path, to_pickle):
        with open(pickle_path, 'wb') as handle:
            pickle.dump(to_pickle, handle)

    def load_metrics(self, pickle_path):
        if os.path.exists(pickle_path):
            with open(pickle_path, 'rb') as handle:
                pickled_results = pickle.load(handle)
                return pickled_results

    def get_ecs(self, file_path):
        res = set()
        with open(file_path) as file:
            json_modules = json.load(file)['children']
            for level_1 in json_modules:
                level_1_name = level_1['name'].split()[0]
                level_1_description = level_1['name'].replace(level_1_name, '').strip()

                level_1_children = level_1['children']
                for level_2 in level_1_children:
                    level_2_name = level_2['name'].split()[0]
                    level_2_description = level_2['name'].replace(level_2_name, '').strip()
                    level_2_children = level_2['children']

                    for level_3 in level_2_children:
                        level_3_name = level_3['name'].split()[0]
                        level_3_description = level_3['name'].replace(level_3_name, '').strip()
                        level_3_children = level_3['children']

                        for level_4 in level_3_children:
                            level_4_name = level_4['name'].split()[0]
                            level_4_description = level_4['name'].replace(level_4_name, '').strip()
                            if not level_4_name.endswith('-'):
                                res.add(level_4_name)

        return res

    def get_kos(self, file_path):
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

    def generate_pickle_ecs_kos(self, pickle_path, ec_json, ko_json):
        all_ecs = self.get_ecs(ec_json)
        all_kos = self.get_kos(ko_json)
        res = {'enzyme_ec': all_ecs, 'kegg_ko': all_kos}
        self.save_metrics(pickle_path, res)

    def yield_all_lines(self, pickle_path):
        for db_id in self.rhea_fetch_all_reactions():
            yield f'{db_id}\trhea\treaction\tcrp\n'
        for db_id in self.metacyc_fetch_all_proteins():
            yield f'{db_id}\tmetacyc\tprotein\tprc\n'
        ecs_kos = self.load_metrics(pickle_path)
        for db in ecs_kos:
            for db_id in ecs_kos[db]:
                yield f'{db_id}\t{db}\tprotein\tprc\n'

    def yield_all_kos(self, pickle_path):
        ecs_kos = self.load_metrics(pickle_path)
        for db_id in ecs_kos['kegg_ko']:
            yield f'{db_id}\tkegg_ko\tprotein\tprc\n'

    def generate_universal_input(self, pickle_path, ec_json, ko_json):
        if not os.path.exists(pickle_path):
            self.generate_pickle_ecs_kos(pickle_path, ec_json, ko_json)
        with open(self.output_path, 'w+') as file:
            for unifuncnet_line in self.yield_all_lines(pickle_path):
                file.write(unifuncnet_line)

    def generate_ko_input(self, pickle_path, ec_json, ko_json):
        if not os.path.exists(pickle_path):
            self.generate_pickle_ecs_kos(pickle_path, ec_json, ko_json)
        with open(self.output_path, 'w+') as file:
            for unifuncnet_line in self.yield_all_kos(pickle_path):
                file.write(unifuncnet_line)

    def parse_tsv(self, input_file):
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
                            if db not in res: res[db] = set()
                            # to avoid bad splitting when dealing with descriptions
                            annot = line_tab[len(db) + 1:]
                            res[db].add(annot)
        return res


if __name__ == '__main__':
    pickle_path = 'ecs_kos.pickle'
    # from https://www.kegg.jp/brite/ko01000
    ec_json = 'ko01000.json'
    # from https://www.genome.jp/kegg-bin/show_brite?ko00001.keg
    ko_json = 'ko00001.json'
    # output_path=f'universal_input.tsv'
    output_path = f'kegg_ko_input.tsv'

    s = InputGenerator(output_path)
    s.generate_ko_input(pickle_path)
