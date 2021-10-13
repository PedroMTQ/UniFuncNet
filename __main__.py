
import argparse
import sys
import os
from datetime import datetime
from sys import platform
import uuid

from source.Utils.util import SCRAPPABLE_DBS,set_scrappable_dbs
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
    searcher = Protein_Searcher(search_direction='gpr', output_folder=output_folder,do_reaction_met_instances=True)
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

def run_searcher(target_path,output_folder,search_type,search_direction,reaction_metabolites,politeness_timer=10):
    start=time()
    run_test_browser_drivers()
    if search_direction=='na': search_direction=''
    if search_type=='gene_search': searcher= Gene_Searcher(search_direction=search_direction,output_folder=output_folder,do_reaction_met_instances=reaction_metabolites,politeness_timer=politeness_timer)
    elif search_type=='protein_search': searcher= Protein_Searcher(search_direction=search_direction,output_folder=output_folder,do_reaction_met_instances=reaction_metabolites,politeness_timer=politeness_timer)
    elif search_type=='compound_search': searcher= Compound_Searcher(search_direction=search_direction,output_folder=output_folder,do_reaction_met_instances=reaction_metabolites,politeness_timer=politeness_timer)
    elif search_type=='reaction_search': searcher= Reaction_Searcher(search_direction=search_direction,output_folder=output_folder,do_reaction_met_instances=reaction_metabolites,politeness_timer=politeness_timer)
    check_list=open(f'{output_folder}finished.txt','w+')
    with open(target_path) as infile:
        line=infile.readline()
        while line:
            line=line.strip('\n')
            if line:
                db_type,db_id=line.split('\t')
                searcher.run_searcher(db_id,db_type)
                print(f'{line}',file=check_list,flush=True)
            line=infile.readline()
        searcher.output_results()
    print(f'DRAX took {time()-start} seconds to run')


if __name__ == '__main__':
    if '--example' in sys.argv:
        run_test()
    elif '--test_web' in sys.argv:
        run_test_web()
    else:
        print('Executing command:\n', ' '.join(sys.argv))
        parser = argparse.ArgumentParser(description='____________  ___  __   __\n'+
                                                        '|  _  \ ___ \/ _ \ \ \ / /\n'+
                                                        '| | | | |_/ / /_\ \ \ V / \n'+
                                                        '| | | |    /|  _  | /   \ \n'+
                                                        '| |/ /| |\ \| | | |/ /^\ \\\n'+
                                                        '|___/ \_| \_\_| |_/\/   \/, a biological database scraper',
                                         formatter_class=argparse.RawTextHelpFormatter)
        #run mantis
        parser.add_argument('search_direction',
                            help='Please choose a search direction from :'
                                 'For example, if given a gene ID with gpr, the search would start in gene,it would then go to protein level and then to reaction. '
                                 'If starting from protein level, the gpr would only do protein to reaction.',
                            choices=['global', 'rpg', 'pg', 'gpr', 'pr','gp','na'])
        parser.add_argument('search_type',
                            help='Please choose from :\n\tgene_search\n\tprotein_search\n\tcompound_search\n\treaction_search\n\n'
                                 'For example, if the user has a list of KEGG gene IDs, the searcher type would be gene_searcher, '
                                 'if the gpr would also be applied, the search would start from these gene IDs, go to protein level and then to reaction.',
                            choices=['gene_search', 'protein_search', 'compound_search', 'reaction_search'])
        parser.add_argument('-i', '--input_path', help='[required]\tTSV file path with a list of identifiers.')
        parser.add_argument('-o', '--output_folder', help='[required]\tOutput folder path for the information collected from the different sources')
        parser.add_argument('-db', '--databases', help='[optional]\tComma separated list of databases to search for information in.')
        parser.add_argument('-pt', '--politeness_timer', help='[optional]\tTime (seconds) between requests. Default is 10. Please be careful not to overleaf the corresponding databases, you might get blocked from doing future requests.')
        parser.add_argument('-rm', '--reaction_metabolites', action='store_true',
                            help='[optional]\tskip memory management. No HMMER memory management (less stable), but may allow for runs to finish in low memory environments')

        args = parser.parse_args()
        target_path = args.input_path
        output_folder = args.output_folder
        search_type = args.search_type
        search_direction = args.search_direction
        politeness_timer = int(args.politeness_timer)
        reaction_metabolites = args.reaction_metabolites
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
            run_searcher(target_path=target_path,output_folder=output_folder,search_type=search_type,search_direction=search_direction,reaction_metabolites=reaction_metabolites,politeness_timer=politeness_timer)
        else:
            print('Input file not found!')
