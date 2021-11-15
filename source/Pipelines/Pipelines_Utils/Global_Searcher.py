import os.path

from source.Utils.util import get_instance_type, check_if_any_none, number_of_nones_dict, unite_instance_list,number_of_nones_dict,find_path,SCRAPPABLE_DBS,download_file_ftp,DRAX_FOLDER,RESOURCES_FOLDER,gunzip
from source.Pipelines.Pipelines_Utils.Memory_Keeper import *
from source.Biological_Components.Compound import Compound
from source.Biological_Components.Gene import Gene
from source.Biological_Components.Protein import Protein
from source.Biological_Components.Reaction import Reaction
from source.Fetchers.Fetchers_Utils.Global_Fetcher import Global_Fetcher

import datetime
import re
import ast
from bs4 import BeautifulSoup



#This is just the general class for the Searchers, instead of replicating code in the several searchers

class Global_Searcher(Memory_Keeper):
    def __init__(self,memory_storage=None,search_direction='global',db_name=None,wanted_org_kegg_codes=[],output_folder=None,politeness_timer=10):
        """
        :param memory_storage: if one is given, memory location will be shared
        :param search_direction:
        types of convergence:
            'gpr'
            'rpg'
            'global'
            None
        """

        self.output_folder=output_folder
        Memory_Keeper.__init__(self,db_name=db_name,politeness_timer=politeness_timer)
        self.original_searcher=get_instance_type(self)
        self.search_direction=search_direction
        self.wanted_org_kegg_codes=wanted_org_kegg_codes
        self.db_name=None
        #if we receive a memory storage, self should become that memory, otherwise we proceed as a simple Memory_Keeper instance
        if memory_storage:
            self.borrow_memory_lists(memory_storage)
            self.borrow_fetchers(memory_storage)
            self.memory_storage = memory_storage
            self.search_direction=memory_storage.search_direction
            self.borrow_searchers()
        else:
            self.setup_memory_lists()
            self.setup_fetchers()
            self.memory_storage=self
            #starting the fetcher for the corresponding main searcher instance
            self.setup_new_searchers()
            self.search_direction=search_direction
        if not os.path.exists(RESOURCES_FOLDER):
            os.mkdir(RESOURCES_FOLDER)




    def set_search_direction(self,search_direction):
        self_searcher_type= get_instance_type(self)
        if self_searcher_type!='Gene_Searcher':
            self.gene_searcher.search_direction=search_direction

        if self_searcher_type != 'Protein_Searcher':
            self.protein_searcher.search_direction=search_direction

        if self_searcher_type != 'Reaction_Searcher':
            self.reaction_searcher.search_direction=search_direction

        if self_searcher_type != 'Compound_Searcher':
            self.compound_searcher.search_direction=search_direction

    def is_valid_search_direction(self,allowed_search_directions):
        if allowed_search_directions.intersection(self.search_direction):
            return True
        else:
            return False

    def set_kegg_org_codes(self,kegg_org_codes):
        self_searcher_type= get_instance_type(self)
        if self_searcher_type!='Gene_Searcher':
            self.gene_searcher.wanted_org_kegg_codes=kegg_org_codes

        if self_searcher_type != 'Protein_Searcher':
            self.protein_searcher.wanted_org_kegg_codes=kegg_org_codes

        if self_searcher_type != 'Reaction_Searcher':
            self.reaction_searcher.wanted_org_kegg_codes=kegg_org_codes

        if self_searcher_type != 'Compound_Searcher':
            self.compound_searcher.wanted_org_kegg_codes=kegg_org_codes



    def flush_memory(self):
        print('Flushing memory!')
        self.setup_memory_lists()
        self.flush_memory_searchers()

    def flush_memory_searchers(self):
        self_searcher_type= get_instance_type(self)
        if self_searcher_type!='Gene_Searcher':
            self.gene_searcher.reactions=self.reactions
            self.gene_searcher.proteins=self.proteins
            self.gene_searcher.genes=self.genes
            self.gene_searcher.compounds=self.compounds

        if self_searcher_type != 'Protein_Searcher':
            self.protein_searcher.reactions=self.reactions
            self.protein_searcher.proteins=self.proteins
            self.protein_searcher.genes=self.genes
            self.protein_searcher.compounds=self.compounds

        if self_searcher_type != 'Reaction_Searcher':
            self.reaction_searcher.reactions=self.reactions
            self.reaction_searcher.proteins=self.proteins
            self.reaction_searcher.genes=self.genes
            self.reaction_searcher.compounds=self.compounds

        if self_searcher_type != 'Compound_Searcher':
            self.compound_searcher.reactions=self.reactions
            self.compound_searcher.proteins=self.proteins
            self.compound_searcher.genes=self.genes
            self.compound_searcher.compounds=self.compounds


    #when we want to run the cpd separetely we run this afterwards
    def run_reaction_met_instances(self,match):
        if self.is_valid_search_direction({'gprc','prc','rc','cr','crp','crpg','global'}):
            rn_with_ids = match.get_detail('rn_with_ids')
            for rn_with_ids_str, rn_with_ids_ids, rn_with_ids_db in rn_with_ids:
                rn_with_ids_ids=ast.literal_eval(rn_with_ids_ids)
                reaction_with_instances= self.reaction_met_instances(rn_with_ids_str, rn_with_ids_ids, rn_with_ids_db)
                if check_if_any_none(reaction_with_instances,pos=1):
                    print('Could not perform reaction_met_instances for:\n',match)
                    return
                match.set_detail('reaction_with_instances',reaction_with_instances)
            match.clear_rn_with_ids()



    def setup_new_searchers(self):
        self.gene_searcher=None
        self.protein_searcher=None
        self.reaction_searcher=None
        self.compound_searcher=None

        self_searcher_type= get_instance_type(self)

        if self_searcher_type!='Gene_Searcher':
            from source.Pipelines.Searchers.Gene_Searcher import Gene_Searcher
            self.gene_searcher = Gene_Searcher(memory_storage=self,
                                               search_direction=self.search_direction,
                                               db_name=self.db_name,
                                               output_folder=self.output_folder,
                                               politeness_timer=self.politeness_timer)
        #if it is a gene searcher, find_gene should point to itself so that self.gene_searcher.find_gene == self.find_Gene. the same applies for the other searchers
        else: self.gene_searcher=self


        if self_searcher_type != 'Protein_Searcher':
            from source.Pipelines.Searchers.Protein_Searcher import Protein_Searcher
            self.protein_searcher = Protein_Searcher(memory_storage=self,
                                                     search_direction=self.search_direction,
                                                     db_name=self.db_name,
                                                     output_folder=self.output_folder,
                                                     politeness_timer=self.politeness_timer)
        else: self.protein_searcher=self

        if self_searcher_type != 'Reaction_Searcher':
            from source.Pipelines.Searchers.Reaction_Searcher import Reaction_Searcher
            self.reaction_searcher = Reaction_Searcher(memory_storage=self,
                                                       search_direction=self.search_direction,
                                                       db_name=self.db_name,
                                                       output_folder=self.output_folder,
                                                       politeness_timer=self.politeness_timer)
        else: self.reaction_searcher=self

        if self_searcher_type != 'Compound_Searcher':
            from source.Pipelines.Searchers.Compound_Searcher import Compound_Searcher
            self.compound_searcher = Compound_Searcher(memory_storage=self,
                                                       search_direction=self.search_direction,
                                                       db_name=self.db_name,
                                                       output_folder=self.output_folder,
                                                       politeness_timer=self.politeness_timer)
        else: self.compound_searcher=self

        #since before creating the searchers we borrow during initizaliation of sub-searchers, we need to launch borrow again to make sure everything is pointing at the same direction
        for inst in [self.gene_searcher,self.protein_searcher,self.reaction_searcher,self.compound_searcher]:
            inst.borrow_searchers()


    def borrow_searchers(self):
        if self is not self.memory_storage:
            self_searcher_type= get_instance_type(self)
            memory_searcher_type= get_instance_type(self.memory_storage)

            self.gene_searcher = self.memory_storage.gene_searcher
            self.protein_searcher = self.memory_storage.protein_searcher
            self.reaction_searcher = self.memory_storage.reaction_searcher
            self.compound_searcher = self.memory_storage.compound_searcher



            if memory_searcher_type=='Gene_Searcher':
                self.gene_searcher = self.memory_storage

            if memory_searcher_type == 'Protein_Searcher':
                self.protein_searcher = self.memory_storage

            if memory_searcher_type == 'Reaction_Searcher':
                self.reaction_searcher = self.memory_storage

            if memory_searcher_type == 'Compound_Searcher':
                self.compound_searcher = self.memory_storage


            if self_searcher_type=='Gene_Searcher':
                self.gene_searcher = self

            if self_searcher_type == 'Protein_Searcher':
                self.protein_searcher = self

            if self_searcher_type == 'Reaction_Searcher':
                self.reaction_searcher = self

            if self_searcher_type == 'Compound_Searcher':
                self.compound_searcher = self





    ###Function inheritance###

    def find_protein(self, db, query_id, extra_args={}):
        return self.protein_searcher.find_protein(db=db, query_id=query_id, extra_args=extra_args)

    def find_gene(self, db, query_id, extra_args={}):
        return self.gene_searcher.find_gene(db=db, query_id=query_id, extra_args=extra_args)

    def find_reaction(self, db, query_id, extra_args={}):
        return self.reaction_searcher.find_reaction(db=db, query_id=query_id, extra_args=extra_args)

    def find_compound(self,db,query_id):
        return self.compound_searcher.find_compound(db=db,query_id=query_id)

    def reaction_met_instances(self, rn, rn_with_ids, db):
        if self.is_valid_search_direction({'gprc','prc','rc','cr','crp','crpg','global'}):
            return self.compound_searcher.reaction_met_instances(rn,rn_with_ids,db)
        else: return None

    ##########################




    def get_protein_match(self,bio_query,bio_db=None):
        return self.get_biological_instance('proteins',bio_query,bio_db)

    def get_reaction_match(self,bio_query,bio_db=None):
        return self.get_biological_instance('reactions',bio_query,bio_db)

    def get_gene_match(self,bio_query,bio_db=None):
        return self.get_biological_instance('genes',bio_query,bio_db)

    def get_compound_match(self,bio_query,bio_db=None):
        return self.get_biological_instance('compounds',bio_query,bio_db)


    def add_counters(self):
        self.instance_counter=0
        for inst_type in [self.get_genes_all(),self.get_proteins_all(),self.get_compounds_all(),self.get_reactions_all()]:
            for inst in inst_type:
                inst.internal_id = self.instance_counter
                self.instance_counter += 1

    def write_to_output_file(self,inst,outfile):
        instance_info = inst.export_all_info()
        outfile.write(f'internal_id:{inst.internal_id}\t{instance_info}\n')

    def output_results(self):
        self.add_counters()
        print('Exporting to spreadsheets')
        with open(self.output_folder+'Genes.tsv','w+') as outfile:
            for inst in self.get_genes_all():
                self.write_to_output_file(inst,outfile)
        with open(self.output_folder+'Proteins.tsv','w+') as outfile:
            for inst in self.get_proteins_all():
                self.write_to_output_file(inst,outfile)

        with open(self.output_folder+'Compounds.tsv','w+') as outfile:
            for inst in self.get_compounds_all():
                print('output',id(inst))
                self.write_to_output_file(inst,outfile)

        with open(self.output_folder+'Reactions.tsv','w+') as outfile:
            for inst in self.get_reactions_all():
                self.write_to_output_file(inst,outfile)




if __name__ == '__main__':
    gs=Global_Searcher()
    print(SCRAPPABLE_DBS)
    #p=Protein({'kegg':'1'})
    #r=Reaction({'kegg':'2'})
    #p.set_detail('reaction_instances',r,converged_in='biocyc')
    #r.set_detail('protein_instances',p,converged_in='kegg')
    #print(r.get_detail('protein_instances',return_convergence=True))
    #print(p.get_detail('reaction_instances',return_convergence=True))
    #b.save_instance_neo4j(r)
    #b.save_instance_neo4j(p)


    #b.save_instance_neo4j(r)
    #p1=b.get_best_match_db({'enzyme_ec':'1.1.1.346'},node_type='Protein')
    #print(p1)