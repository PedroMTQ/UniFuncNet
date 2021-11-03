
'''
entry points:
compound:
synonyms
brenda
biocyc
kegg
chemspider
inchi_key
drugbank
hmdb

reaction:
brenda
biocyc
kegg
drugbank
hmdb


protein:
enzyme_ec
kegg
biocyc
brenda
uniprot
kegg_ko

gene:
kegg
uniprot
biocyc


python DRAX/ gpr protein_search -i DRAX/test/test.tsv -o drax_test -db kegg -rm

scp -P 8022 -r /home/pedroq/PycharmProjects/DRAX/source pedro.queiros@10.240.6.223:/home/pedro.queiros/python_projects/DRAX/
'''


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
    searcher = Protein_Searcher(search_direction={'gpr'}, output_folder=output_folder,do_reaction_met_instances=True)
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

def set_search_direction(searchers_list,search_direction):
    for searcher in searchers_list:
        searcher.set_search_direction(search_direction)

def set_do_reaction_met_instances(searchers_list,do_reaction_met_instances):
    for searcher in searchers_list:
        searcher.set_do_reaction_met_instances(do_reaction_met_instances)


def run_searcher(target_path,output_folder,reaction_metabolites,politeness_timer):
    start=time()
    #run_test_browser_drivers()
    gene_searcher= Gene_Searcher(output_folder=output_folder,do_reaction_met_instances=reaction_metabolites,politeness_timer=politeness_timer)
    protein_searcher= gene_searcher.protein_searcher
    compound_searcher= gene_searcher.compound_searcher
    reaction_searcher= gene_searcher.reaction_searcher
    check_list=open(f'{output_folder}finished.txt','w+')
    with open(target_path) as infile:
        line=infile.readline()
        while line:
            line=line.strip('\n')
            if line:
                instance_type,search_direction,db_type,db_id=line.split('\t')

                if db_type.startswith('syn'): db_type='synonyms'
                search_direction=search_direction.split(',')
                search_direction=set(search_direction)
                instance_type=instance_type.lower()
                if instance_type=='gene':       searcher=gene_searcher
                elif instance_type=='protein':  searcher=protein_searcher
                elif instance_type=='compound': searcher=compound_searcher
                elif instance_type=='reaction': searcher=reaction_searcher

                set_search_direction([gene_searcher,protein_searcher,compound_searcher,reaction_searcher],search_direction)
                if instance_type == 'compound':
                    set_do_reaction_met_instances([gene_searcher,protein_searcher,compound_searcher,reaction_searcher],True)
                    searcher.run_searcher(db_id,db_type,convergence_search=True)

                else:
                    set_do_reaction_met_instances([gene_searcher,protein_searcher,compound_searcher,reaction_searcher],reaction_metabolites)
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
        #parser.add_argument('search_direction',
        #                    help='Please choose a search direction from :'
        #                         'For example, if given a gene ID with gpr, the search would start in gene,it would then go to protein level and then to reaction. '
        #                         'If starting from protein level, the gpr would only do protein to reaction.',
         #                   choices=['global', 'rpg', 'pg', 'gpr', 'pr','gp','na'])
        #parser.add_argument('search_type',
        #                    help='Please choose from :\n\tgene_search\n\tprotein_search\n\tcompound_search\n\treaction_search\n\n'
        #                         'For example, if the user has a list of KEGG gene IDs, the searcher type would be gene_searcher, '
        #                         'if the gpr would also be applied, the search would start from these gene IDs, go to protein level and then to reaction.',
        #                    choices=['gene_search', 'protein_search', 'compound_search', 'reaction_search'])
        parser.add_argument('-i', '--input_path', help='[required]\tTSV file path with a list of identifiers.')
        parser.add_argument('-o', '--output_folder', help='[required]\tOutput folder path for the information collected from the different sources')
        parser.add_argument('-db', '--databases', help='[optional]\tComma separated list of databases to search for information in.')
        parser.add_argument('-pt', '--politeness_timer', help='[optional]\tTime (seconds) between requests. Default is 10. Please be careful not to overleaf the corresponding databases, you might get blocked from doing future requests.')
        parser.add_argument('-rm', '--reaction_metabolites', action='store_true',
                            help='[optional]\tskip memory management. No HMMER memory management (less stable), but may allow for runs to finish in low memory environments')

        args = parser.parse_args()
        target_path = args.input_path
        output_folder = args.output_folder
        politeness_timer = args.politeness_timer
        if politeness_timer:
            politeness_timer = int(politeness_timer)
        else: politeness_timer=10
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
            run_searcher(target_path=target_path,output_folder=output_folder,reaction_metabolites=reaction_metabolites,politeness_timer=politeness_timer)
        else:
            print('Input file not found!')
