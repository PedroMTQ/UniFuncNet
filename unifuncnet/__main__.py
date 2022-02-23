import argparse
import sys
import os
from datetime import datetime
from sys import platform
import uuid

from unifuncnet.Utils.util import SCRAPPABLE_DBS,set_scrappable_dbs,VALID_DIRECTIONS,print_version,UniFuncNet_FOLDER,RESOURCES_FOLDER,check_all_resources,check_all_resources
from unifuncnet.Searchers.Gene_Searcher import  Gene_Searcher
from unifuncnet.Searchers.Protein_Searcher import  Protein_Searcher
from unifuncnet.Searchers.Compound_Searcher import  Compound_Searcher
from unifuncnet.Searchers.Reaction_Searcher import  Reaction_Searcher
from time import time

if platform.startswith('win'):
    SPLITTER = '\\'
else:
    SPLITTER = '/'



def run_test():
    datetime_str = str(datetime.now().strftime("%Y-%m-%dT%H%M%S"))
    output_folder = os.getcwd() + SPLITTER + 'UniFuncNet_OUTPUT_' + datetime_str + SPLITTER
    os.mkdir(output_folder)
    print(f'No output folder provided! Saving data to: {output_folder}')
    searcher = Protein_Searcher(search_mode={'pr'}, output_folder=output_folder,politeness_timer=2)
    searcher.run_searcher('1.1.1.178', 'enzyme_ec')
    searcher.output_results()

def argv_gsmm_expansion_function():
    from Workflows.GSMM_Expansion.GSMM_Expansion import GSMM_Expansion
    print('Executing command:\n', ' '.join(sys.argv))
    parser = argparse.ArgumentParser(description='This workflow suggests new connections for Carveme metabolic models\n',formatter_class=argparse.RawTextHelpFormatter)
    #just a placeholder
    parser.add_argument('workflow')

    parser.add_argument('-i', '--input_folder', help='[required]\tInput folder with protein sequences fastas')
    parser.add_argument('-o', '--output_folder', help='[required]\tOutput directory')
    parser.add_argument('-db','--database', help='[optional]\tDatabases to be used in UniFuncNet')
    parser.add_argument('-mr', '--metacyc_ref', help='[optional]\tMetacyc reference folder.')
    parser.add_argument('-pt', '--politeness_timer', help='[optional]\tTime (seconds) between requests. Default is 10. Please be careful not to overleaf the corresponding databases, you might get blocked from doing future requests.')
    parser.add_argument('-oc','--only_connected', action='store_true', help='[optional]\tExpand network with only nodes that are connected to the original network (this is off by default)')
    args = parser.parse_args()
    input_folder = args.input_folder
    output_folder = args.output_folder
    database = args.database
    metacyc_ref = args.metacyc_ref
    politeness_timer = args.politeness_timer
    only_connected = args.only_connected
    if politeness_timer: politeness_timer=int(politeness_timer)
    else: politeness_timer=10
    if input_folder and output_folder:
        GSMM_Expansion(input_folder=input_folder,output_folder=output_folder,metacyc_ref=metacyc_ref,database=database,only_connected=only_connected,politeness_timer=politeness_timer)
    else:
        print('Missing input or output folders')

    pass

def argv_compounds_to_organism_function():
    from Workflows.Compounds_to_Organisms_Mapping.Compounds_to_Organisms_Mapping import Compounds_to_Organisms_Mapping

    print('Executing command:\n', ' '.join(sys.argv))
    parser = argparse.ArgumentParser(
        description='This workflow suggests new connections for Carveme metabolic models\n',
        formatter_class=argparse.RawTextHelpFormatter)
    #just a placeholder
    parser.add_argument('workflow')
    parser.add_argument('-i', '--input_samples', help='[required]\tInput folder with protein sequences fastas')
    parser.add_argument('-m', '--metabolites',
                        help='[required]\tMetabolites list (synonyms or ChEBI IDs), with each metabolite in a separate line')
    parser.add_argument('-mr', '--metacyc_ref', help='[required]\tMetacyc reference folder.')
    parser.add_argument('-o', '--output_folder', help='[required]\tOutput directory')
    parser.add_argument('-pt', '--politeness_timer',
                        help='[optional]\tTime (seconds) between requests. Default is 10. Please be careful not to overleaf the corresponding databases, you might get blocked from doing future requests.')
    parser.add_argument('-db', '--database', help='[optional]\tDatabases to be used in UniFuncNet')
    args = parser.parse_args()
    input_samples = args.input_samples
    metabolites = args.metabolites
    metacyc_ref = args.metacyc_ref
    output_folder = args.output_folder
    politeness_timer = args.politeness_timer
    database = args.database
    if politeness_timer:
        politeness_timer = int(politeness_timer)
    else:
        politeness_timer = 10
    if os.path.exists(input_samples) and os.path.exists(metabolites) and output_folder and metacyc_ref:
        Compounds_to_Organisms_Mapping(input_samples=input_samples, metabolites=metabolites,
                                       output_folder=output_folder, metacyc_ref=metacyc_ref, database=database,
                                       politeness_timer=politeness_timer)
    else:
        print(f'Missing a parameter:\n'
              f'input_samples {input_samples}'
              f'metabolites {metabolites}'
              f'output_folder {output_folder}'
              f'metacyc_ref {metacyc_ref}'
              )

def argv_neo4j_function():
    from Workflows.UniFuncNet_Neo4j_Connector.UniFuncNet_Neo4j_Connector import UniFuncNet_Neo4j_Connector
    print('Executing command:\n', ' '.join(sys.argv))
    parser = argparse.ArgumentParser(
        description='This workflow creates a Neo4j database and can export networks from mantis annotations\n',
        formatter_class=argparse.RawTextHelpFormatter)
    #just a placeholder
    parser.add_argument('workflow')
    parser.add_argument('-i', '--input', help='[required]\tInput folder if creating the Neo4j database or consensus_annotation.tsv if exporting network')
    parser.add_argument('-o', '--output_folder', help='[optional]\tOutput directory')
    parser.add_argument('-db', '--database_name', help='[optional]\tDatabase name, <neo4j> by default')
    parser.add_argument('-u', '--username', help='[optional]\tUsername, <neo4j> by default')
    parser.add_argument('-pw', '--password', help='[optional]\tPassword, <unifuncnet_neo4j> by default')
    parser.add_argument('-p', '--port', help='[optional]\tPort, <7687> by default')
    args = parser.parse_args()

    input_path = args.input
    output_folder = args.output_folder
    database_name = args.database_name
    username = args.username
    password = args.password
    port = args.port

    if not database_name:
        database_name='neo4j'
    if not username:
        username='neo4j'
    if not password:
        password='unifuncnet_neo4j'
    if not port:
        port='7687'



    uri = f'bolt://localhost:{port}'
    browser_link=f'http://localhost:7474/browser'
    print(f'Connecting to neo4j with:\n'
          f'Database name:\t{database_name}\n'
          f'Port:\t{port}\n'
          f'Username:\t{username}\n'
          f'Password:\t{password}\n'
          f'URI name:\t{uri}\n'
          f'Browser link:\t{browser_link}\n')
    neo4j_driver = UniFuncNet_Neo4j_Connector(username, password, db_name=database_name, uri=uri)
    print(f'Connected to neo4j.')


    if input_path.endswith('.tsv'):
        if output_folder:
            if not output_folder.endswith(SPLITTER): output_folder+=SPLITTER
        else:
            datetime_str = str(datetime.now().strftime("%Y-%m-%dT%H%M%S"))
            output_folder = os.getcwd()+SPLITTER+'UniFuncNet_OUTPUT_' + datetime_str+SPLITTER
            print(f'No output folder provided! Saving data to: {output_folder}')
        if os.path.exists(output_folder):
            if os.listdir(output_folder):
                datetime_str = str(datetime.now().strftime("%Y-%m-%dT%H%M%S"))
                hex_random = '_hex_' + uuid.uuid4().hex[:10]
                if not output_folder.endswith(SPLITTER): output_folder+=SPLITTER
                output_folder = f'{output_folder}{datetime_str}{hex_random}{SPLITTER}'
                print(f'The output folder already contains something! New output folder will be: {output_folder}')
                os.mkdir(output_folder)
        else:
            os.mkdir(output_folder)

        neo4j_driver.get_network_from_annotations(input_path,output_folder)

    else:
        neo4j_driver.export_unifuncnet_to_neo4j(input_path)

def argv_input_generator_function():
    from Workflows.Input_Generator.Input_Generator import Input_Generator
    print('Executing command:\n', ' '.join(sys.argv))
    parser = argparse.ArgumentParser(
        description='This workflow creates a UniFuncNet input tsv with Rhea reaction IDs, Metacyc protein IDs, KOs, and ECs\n',
        formatter_class=argparse.RawTextHelpFormatter)
    #just a placeholder
    parser.add_argument('workflow')
    parser.add_argument('-o', '--output_folder', help='[required]\tOutput folder')
    parser.add_argument('-ec', '--ec_json', help='[required]\tEC json from https://www.kegg.jp/brite/ko01000')
    parser.add_argument('-ko', '--ko_json', help='[required]\tKO json from https://www.genome.jp/kegg-bin/show_brite?ko00001.keg')
    parser.add_argument('-db', '--database', help='[optional]\tDatabase to generate input. <all> by default',choices = ['all', 'ko'])

    args = parser.parse_args()
    output_folder = args.output_folder
    ec_json = args.ec_json
    ko_json = args.ko_json
    if not ec_json or not ko_json:
        print('Missing json files')
        return
    database = args.database
    if not database: database='all'

    if output_folder:
        if not output_folder.endswith(SPLITTER): output_folder += SPLITTER
    else:
        datetime_str = str(datetime.now().strftime("%Y-%m-%dT%H%M%S"))
        output_folder = os.getcwd() + SPLITTER + 'UniFuncNet_OUTPUT_' + datetime_str + SPLITTER
        print(f'No output folder provided! Saving data to: {output_folder}')
    if os.path.exists(output_folder):
        if os.listdir(output_folder):
            datetime_str = str(datetime.now().strftime("%Y-%m-%dT%H%M%S"))
            hex_random = '_hex_' + uuid.uuid4().hex[:10]
            if not output_folder.endswith(SPLITTER): output_folder += SPLITTER
            output_folder = f'{output_folder}{datetime_str}{hex_random}{SPLITTER}'
            print(f'The output folder already contains something! New output folder will be: {output_folder}')
            os.mkdir(output_folder)
    else:
        os.mkdir(output_folder)



    pickle_path=f'{output_folder}ecs_kos.pickle'
    output_path=f'{output_folder}universal_input.tsv'
    s=Input_Generator(output_path=output_path)
    if database=='ko':
        s.generate_ko_input(output_path,pickle_path, ec_json, ko_json)
    else:
        s.generate_universal_input(pickle_path,ec_json, ko_json)

def set_search_mode(searchers_list,search_mode):
    for searcher in searchers_list:
        searcher.set_search_mode(search_mode)

def set_kegg_org_codes(searchers_list,kegg_org_codes):
    for searcher in searchers_list:
        searcher.set_kegg_org_codes(kegg_org_codes)

def check_validity_input(target_path):
    with open(target_path) as infile:
        line=infile.readline()
        while line:
            line=line.strip('\n')
            if line:
                line_split=line.split('\t')
                search_mode=''
                if len(line_split) == 5:
                    db_id, db_type, instance_type, search_mode, kegg_org_codes = line_split
                elif len(line_split) == 4:
                    db_id, db_type, instance_type, search_mode = line_split
                else:
                    db_id, db_type, instance_type = line_split
                search_mode=[i.lower() for i in search_mode.split(',')]
                search_mode=set(search_mode)
                if search_mode.difference(VALID_DIRECTIONS):
                    print(f'Invalid search direction, please fix this line:\n{line}')
                    raise Exception
                instance_type=instance_type.lower()
                if instance_type not in ['gene','protein','compound','reaction']:
                    print(f'Invalid instance type, please fix this line:\n{line}')
                    raise Exception
            line=infile.readline()

def run_searcher(target_path,output_folder,politeness_timer):
    start=time()
    print(f'Available databases:\n{",".join(SCRAPPABLE_DBS)}')
    current_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(f"UniFuncNet started running at  {current_time}")

    gene_searcher= Gene_Searcher(output_folder=output_folder,politeness_timer=politeness_timer)
    protein_searcher= gene_searcher.protein_searcher
    compound_searcher= gene_searcher.compound_searcher
    reaction_searcher= gene_searcher.reaction_searcher
    check_list=open(f'{output_folder}finished.tsv','w+')
    with open(target_path) as infile:
        line=infile.readline()
        while line:
            line=line.strip('\n')
            if line:
                line_split=line.split('\t')
                kegg_org_codes=[]
                search_mode=''
                if len(line_split)<3:
                    print('Invalid input line')

                else:
                    if len(line_split)==5:
                        db_id,db_type,instance_type,search_mode,kegg_org_codes = line_split
                    elif len(line_split)==4:
                        db_id,db_type,instance_type,search_mode = line_split
                    else:
                        db_id,db_type,instance_type= line_split

                    instance_type=instance_type.lower()

                    if db_type.startswith('syn'): db_type='synonyms'
                    search_mode=[i.lower() for i in search_mode.split(',')]
                    search_mode=set(search_mode)
                    instance_type=instance_type.lower()
                    if instance_type=='gene':       searcher=gene_searcher
                    elif instance_type=='protein':  searcher=protein_searcher
                    elif instance_type=='compound': searcher=compound_searcher
                    elif instance_type=='reaction': searcher=reaction_searcher

                    set_search_mode([gene_searcher,protein_searcher,compound_searcher,reaction_searcher],search_mode)
                    set_kegg_org_codes([gene_searcher,protein_searcher,compound_searcher,reaction_searcher],kegg_org_codes)
                    if instance_type == 'compound' or instance_type == 'protein':
                        searcher.run_searcher(db_id,db_type,convergence_search=True)
                    else:
                        searcher.run_searcher(db_id,db_type)
                print(f'{line}',file=check_list,flush=True)
            line=infile.readline()
        searcher.output_results()
    print(f'UniFuncNet took {time()-start} seconds to run')

def main():
    if '--example' in sys.argv:
        run_test()
    elif 'compounds_to_organism' in sys.argv:
        argv_compounds_to_organism_function()
    elif 'gsmm_expansion' in sys.argv:
        argv_gsmm_expansion_function()
    elif 'neo4j' in sys.argv:
        argv_neo4j_function()
    elif 'input_generator' in sys.argv:
        argv_input_generator_function()
    elif '--version' in sys.argv:
        print_version('pedromtq','unifuncnet')
    else:
        search_modes_str='; '.join(VALID_DIRECTIONS)
        print('Executing command:\n', ' '.join(sys.argv))
        #https://patorjk.com/software/taag/#p=display&f=Doom&t=UniFuncNet
        drax_styled=f' _   _       _______                _   _      _\n' \
                    '| | | |     (_)  ___|              | \ | |    | |\n' \
                    '| | | |_ __  _| |_ _   _ _ __   ___|  \| | ___| |_\n' \
                    '| | | | \'_ \| |  _| | | | \'_ \ / __| . ` |/ _ \ __|\n' \
                    '| |_| | | | | | | | |_| | | | | (__| |\  |  __/ |_\n' \
                    ' \___/|_| |_|_\_|  \__,_|_| |_|\___\_| \_/\___|\__|, a biological network annotator.'
        parser = argparse.ArgumentParser(description=drax_styled,
                                         formatter_class=argparse.RawTextHelpFormatter)
        parser.add_argument('-i', '--input_path', help='[required]\tTSV file path with a list of identifiers.')
        parser.add_argument('-o', '--output_folder', help='[required]\tOutput folder path for the information collected from the different sources')
        parser.add_argument('-db', '--databases', help='[optional]\tComma separated list of databases to search for information in.')
        parser.add_argument('-pt', '--politeness_timer', help='[optional]\tTime (seconds) between requests. Default is 10. Please be careful not to overleaf the corresponding databases, you might get blocked from doing future requests.')

        print(f'Resources folder is here:\n{RESOURCES_FOLDER}')

        args = parser.parse_args()
        target_path = args.input_path
        output_folder = args.output_folder
        politeness_timer = args.politeness_timer
        if not target_path:
            print('Missing input path, quitting!')
            return
        if politeness_timer:
            politeness_timer = int(politeness_timer)
        else: politeness_timer=10
        databases = args.databases
        if output_folder:
            if not output_folder.endswith(SPLITTER): output_folder+=SPLITTER
        else:
            datetime_str = str(datetime.now().strftime("%Y-%m-%dT%H%M%S"))
            output_folder = os.getcwd()+SPLITTER+'UniFuncNet_OUTPUT_' + datetime_str+SPLITTER
            print(f'No output folder provided! Saving data to: {output_folder}')
        if os.path.exists(output_folder):
            if os.listdir(output_folder):
                datetime_str = str(datetime.now().strftime("%Y-%m-%dT%H%M%S"))
                hex_random = '_hex_' + uuid.uuid4().hex[:10]
                if not output_folder.endswith(SPLITTER): output_folder+=SPLITTER
                output_folder = f'{output_folder}{datetime_str}{hex_random}{SPLITTER}'
                print(f'The output folder already contains something! New output folder will be: {output_folder}')
                os.mkdir(output_folder)
        else:
            os.mkdir(output_folder)

        if databases:
            databases=databases.split(',')
            set_scrappable_dbs(databases)
        print(f'Time between requests: {politeness_timer} seconds')
        if os.path.exists(target_path):
            passed_check=False
            check_validity_input(target_path)

            try:
                check_validity_input(target_path)
                passed_check=True
            except:
                print('Invalid input!')
            if passed_check:
                run_searcher(target_path=target_path,
                             output_folder=output_folder,
                             politeness_timer=politeness_timer)

        else:
            print('Input file not found!')

if __name__ == '__main__':
    main()
