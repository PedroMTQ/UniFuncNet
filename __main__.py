import argparse
import sys
import os
from datetime import datetime
from sys import platform
import uuid

from source.Utils.util import SCRAPPABLE_DBS,set_scrappable_dbs,VALID_DIRECTIONS,print_version
from source.Pipelines.Searchers.Gene_Searcher import  Gene_Searcher
from source.Pipelines.Searchers.Protein_Searcher import  Protein_Searcher
from source.Pipelines.Searchers.Compound_Searcher import  Compound_Searcher
from source.Pipelines.Searchers.Reaction_Searcher import  Reaction_Searcher
#https://patorjk.com/software/taag/#p=display&f=Alpha&t=DRAX
from time import time

if platform.startswith('win'):
    SPLITTER = '\\'
else:
    SPLITTER = '/'

def run_test():
    run_test_browser_drivers()
    set_scrappable_dbs(['kegg','biocyc'])
    datetime_str = str(datetime.now().strftime("%Y-%m-%dT%H%M%S"))
    output_folder = os.getcwd() + SPLITTER + 'DRAX_OUTPUT_' + datetime_str + SPLITTER
    os.mkdir(output_folder)
    print(f'No output folder provided! Saving data to: {output_folder}')
    searcher = Protein_Searcher(search_mode={'gpr'}, output_folder=output_folder,do_reaction_met_instances=True)
    searcher.run_searcher('1.1.1.178', 'enzyme_ec')
    searcher.output_results()

def run_test_web():
    from source.Utils.Web_Connector import Web_Connector
    print('Launching connector')
    f=Web_Connector()
    url='https://biocyc.org/META/substring-search?type=NIL&object=%28s%29-2%2C5-diaminopentanoate&quickSearch=Quick+Search'
    print('Trying to catch url')
    a=f.try_until_catch(url)
    print(a)
    print('Getting with selenium')
    a=f.get_driver_selenium(url)
    print(a)

def run_test_browser_drivers():
    from source.Utils.Web_Connector import Web_Connector
    f=Web_Connector()
    url='https://www.kegg.jp/'
    try:
        f.get_driver_selenium(url)
        print('Passed selenium check')
    except:
        print('Did not manage to connect using selenium! Did you download the drivers?')
        raise Exception

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
    print('change')
    #run_test_browser_drivers()
    print(f'Available databases:\n{",".join(SCRAPPABLE_DBS)}')
    current_time = datetime.now()
    print(f"DRAX started running at  {current_time}")

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
    print(f'DRAX took {time()-start} seconds to run')

def main():
    if '--example' in sys.argv:
        run_test()
    elif '--test_web' in sys.argv:
        run_test_web()
    elif '--version' in sys.argv:
        print_version('pedromtq','drax')
    else:
        search_modes_str='; '.join(VALID_DIRECTIONS)
        print('Executing command:\n', ' '.join(sys.argv))
        parser = argparse.ArgumentParser(description='____________  ___  __   __\n'+
                                                        '|  _  \ ___ \/ _ \ \ \ / /\n'+
                                                        '| | | | |_/ / /_\ \ \ V / \n'+
                                                        '| | | |    /|  _  | /   \ \n'+
                                                        '| |/ /| |\ \| | | |/ /^\ \\\n'+
                                                        f'|___/ \_| \_\_| |_/\/   \/, a biological database scraper.\nThese are the valid search directions:{search_modes_str}',
                                         formatter_class=argparse.RawTextHelpFormatter)
        parser.add_argument('-i', '--input_path', help='[required]\tTSV file path with a list of identifiers.')
        parser.add_argument('-o', '--output_folder', help='[required]\tOutput folder path for the information collected from the different sources')
        parser.add_argument('-db', '--databases', help='[optional]\tComma separated list of databases to search for information in.')
        parser.add_argument('-pt', '--politeness_timer', help='[optional]\tTime (seconds) between requests. Default is 10. Please be careful not to overleaf the corresponding databases, you might get blocked from doing future requests.')


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
            output_folder = os.getcwd()+SPLITTER+'DRAX_OUTPUT_' + datetime_str+SPLITTER
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
