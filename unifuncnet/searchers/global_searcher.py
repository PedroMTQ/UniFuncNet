from unifuncnet.Utils.util import (get_instance_type,
                                   check_if_any_none,
                                   unite_instance_list,
                                   number_of_nones_dict,
                                   find_path,
                                   SCRAPPABLE_DBS,
                                   download_file_ftp,
                                   UniFuncNet_FOLDER,
                                   RESOURCES_FOLDER,
                                   gunzip)
from unifuncnet.searchers.memory_keeper import *
from unifuncnet.biological_components.compound import Compound
from unifuncnet.biological_components.gene import Gene
from unifuncnet.biological_components.protein import Protein
from unifuncnet.biological_components.reaction import Reaction
from unifuncnet.fetchers.global_fetcher import Global_Fetcher
from unifuncnet.utils.metacyc_sqlite_connector import MetacycSqliteConnector
from unifuncnet.utils.rhea_sqlite_connector import RheaSqliteConnector
from unifuncnet.utils.chebi_sqlite_connector import ChebiSqliteConnector
from unifuncnet.searchers.gene_searcher import GeneSearcher
from unifuncnet.searchers.protein_searcher import ProteinSearcher
from unifuncnet.searchers.reaction_searcher import ReactionSearcher
from unifuncnet.searchers.compound_searcher import CompoundSearcher

import os
import datetime
import re
import ast
from bs4 import BeautifulSoup


# This is just the general class for the Searchers, instead of replicating code in the several searchers

class Global_Searcher(MemoryKeeper,
                      MetacycSqliteConnector,
                      RheaSqliteConnector,
                      ChebiSqliteConnector):
    def __init__(self,
                 memory_storage=None,
                 search_mode='global',
                 db_name=None,
                 wanted_org_kegg_codes=[],
                 output_folder=None,
                 politeness_timer=10):
        self.output_folder = output_folder
        MemoryKeeper.__init__(self, db_name=db_name, politeness_timer=politeness_timer)
        MetacycSqliteConnector.__init__(self)
        RheaSqliteConnector.__init__(self)
        ChebiSqliteConnector.__init__(self)
        self.start_cursors()
        self.original_searcher = get_instance_type(self)
        self.search_mode = search_mode
        self.wanted_org_kegg_codes = wanted_org_kegg_codes
        self.db_name = None
        # if we receive a memory storage, self should become that memory, otherwise we proceed as a simple Memory_Keeper instance
        if memory_storage:
            self.borrow_memory_lists(memory_storage)
            self.borrow_fetchers(memory_storage)
            self.memory_storage = memory_storage
            self.search_mode = memory_storage.search_mode
            self.borrow_searchers()
        else:
            self.setup_memory_lists()
            self.setup_fetchers()
            self.memory_storage = self
            # starting the fetcher for the corresponding main searcher instance
            self.setup_new_searchers()
            self.search_mode = search_mode
        if not os.path.exists(RESOURCES_FOLDER):
            os.mkdir(RESOURCES_FOLDER)

    def start_cursors(self):
        self.metacyc_start_sqlite_cursor()
        self.rhea_start_sqlite_cursor()
        self.chebi_start_sqlite_cursor()

    def set_search_mode(self, search_mode):
        self_searcher_type = get_instance_type(self)
        if self_searcher_type != 'GeneSearcher':
            self.gene_searcher.search_mode = search_mode

        if self_searcher_type != 'ProteinSearcher':
            self.protein_searcher.search_mode = search_mode

        if self_searcher_type != 'ReactionSearcher':
            self.reaction_searcher.search_mode = search_mode

        if self_searcher_type != 'CompoundSearcher':
            self.compound_searcher.search_mode = search_mode

    def is_valid_search_mode(self, allowed_search_modes):
        if allowed_search_modes.intersection(self.search_mode):
            return True
        else:
            return False

    def set_kegg_org_codes(self, kegg_org_codes):
        self_searcher_type = get_instance_type(self)
        if self_searcher_type != 'GeneSearcher':
            self.gene_searcher.wanted_org_kegg_codes = kegg_org_codes

        if self_searcher_type != 'ProteinSearcher':
            self.protein_searcher.wanted_org_kegg_codes = kegg_org_codes

        if self_searcher_type != 'ReactionSearcher':
            self.reaction_searcher.wanted_org_kegg_codes = kegg_org_codes

        if self_searcher_type != 'CompoundSearcher':
            self.compound_searcher.wanted_org_kegg_codes = kegg_org_codes

    def flush_memory(self):
        print('Flushing memory!')
        self.setup_memory_lists()
        self.flush_memory_searchers()

    def flush_memory_searchers(self):
        self_searcher_type = get_instance_type(self)
        if self_searcher_type != 'GeneSearcher':
            self.gene_searcher.reactions = self.reactions
            self.gene_searcher.proteins = self.proteins
            self.gene_searcher.genes = self.genes
            self.gene_searcher.compounds = self.compounds

        if self_searcher_type != 'ProteinSearcher':
            self.protein_searcher.reactions = self.reactions
            self.protein_searcher.proteins = self.proteins
            self.protein_searcher.genes = self.genes
            self.protein_searcher.compounds = self.compounds

        if self_searcher_type != 'ReactionSearcher':
            self.reaction_searcher.reactions = self.reactions
            self.reaction_searcher.proteins = self.proteins
            self.reaction_searcher.genes = self.genes
            self.reaction_searcher.compounds = self.compounds

        if self_searcher_type != 'CompoundSearcher':
            self.compound_searcher.reactions = self.reactions
            self.compound_searcher.proteins = self.proteins
            self.compound_searcher.genes = self.genes
            self.compound_searcher.compounds = self.compounds

    def setup_new_searchers(self):
        self.gene_searcher = None
        self.protein_searcher = None
        self.reaction_searcher = None
        self.compound_searcher = None

        self_searcher_type = get_instance_type(self)

        if self_searcher_type != 'GeneSearcher':
            self.gene_searcher = GeneSearcher(memory_storage=self,
                                               search_mode=self.search_mode,
                                               db_name=self.db_name,
                                               output_folder=self.output_folder,
                                               politeness_timer=self.politeness_timer)
        # if it is a gene searcher, find_gene should point to itself so that self.gene_searcher.find_gene == self.find_Gene. the same applies for the other searchers
        else:
            self.gene_searcher = self

        if self_searcher_type != 'ProteinSearcher':
            self.protein_searcher = ProteinSearcher(memory_storage=self,
                                                     search_mode=self.search_mode,
                                                     db_name=self.db_name,
                                                     output_folder=self.output_folder,
                                                     politeness_timer=self.politeness_timer)
        else:
            self.protein_searcher = self

        if self_searcher_type != 'ReactionSearcher':
            self.reaction_searcher = ReactionSearcher(memory_storage=self,
                                                       search_mode=self.search_mode,
                                                       db_name=self.db_name,
                                                       output_folder=self.output_folder,
                                                       politeness_timer=self.politeness_timer)
        else:
            self.reaction_searcher = self

        if self_searcher_type != 'CompoundSearcher':
            self.compound_searcher = CompoundSearcher(memory_storage=self,
                                                       search_mode=self.search_mode,
                                                       db_name=self.db_name,
                                                       output_folder=self.output_folder,
                                                       politeness_timer=self.politeness_timer)
        else:
            self.compound_searcher = self

        # since before creating the searchers we borrow during initizaliation of sub-searchers, we need to launch borrow again to make sure everything is pointing at the same direction
        for inst in [self.gene_searcher, self.protein_searcher, self.reaction_searcher, self.compound_searcher]:
            inst.borrow_searchers()

    def borrow_searchers(self):
        if self is not self.memory_storage:
            self_searcher_type = get_instance_type(self)
            memory_searcher_type = get_instance_type(self.memory_storage)
            self.gene_searcher = self.memory_storage.gene_searcher
            self.protein_searcher = self.memory_storage.protein_searcher
            self.reaction_searcher = self.memory_storage.reaction_searcher
            self.compound_searcher = self.memory_storage.compound_searcher

            if memory_searcher_type == 'GeneSearcher':
                self.gene_searcher = self.memory_storage

            if memory_searcher_type == 'ProteinSearcher':
                self.protein_searcher = self.memory_storage

            if memory_searcher_type == 'ReactionSearcher':
                self.reaction_searcher = self.memory_storage

            if memory_searcher_type == 'CompoundSearcher':
                self.compound_searcher = self.memory_storage

            if self_searcher_type == 'GeneSearcher':
                self.gene_searcher = self

            if self_searcher_type == 'ProteinSearcher':
                self.protein_searcher = self

            if self_searcher_type == 'ReactionSearcher':
                self.reaction_searcher = self

            if self_searcher_type == 'CompoundSearcher':
                self.compound_searcher = self

    ###Function inheritance###

    def find_protein(self, db, query_id, extra_args={}, convergence_search=False):
        return self.protein_searcher.run_searcher(db=db, query_id=query_id, extra_args=extra_args,
                                                  convergence_search=convergence_search)

    def find_gene(self, db, query_id, extra_args={}):
        return self.gene_searcher.run_searcher(db=db, query_id=query_id, extra_args=extra_args)

    def find_reaction(self, db, query_id, extra_args={}):
        return self.reaction_searcher.run_searcher(db=db, query_id=query_id, extra_args=extra_args)

    def find_compound(self, db, query_id):
        return self.compound_searcher.run_searcher(db=db, query_id=query_id)

    def reaction_met_instances(self, rn, rn_with_ids, db):
        if self.is_valid_search_mode({'gprc', 'prc', 'rc', 'cr', 'crp', 'crpg', 'global'}):
            return self.compound_searcher.reaction_met_instances(rn, rn_with_ids, db)
        else:
            return None

    def reaction_met_instances_simple(self, reaction_stoichiometry, db):
        if self.is_valid_search_mode({'gprc', 'prc', 'rc', 'cr', 'crp', 'crpg', 'global'}):
            return self.compound_searcher.reaction_met_instances_simple(reaction_stoichiometry, db)
        else:
            return None

    ##########################

    def get_protein_match(self, bio_query, bio_db=None):
        return self.get_biological_instance('proteins', bio_query, bio_db)

    def get_reaction_match(self, bio_query, bio_db=None):
        return self.get_biological_instance('reactions', bio_query, bio_db)

    def get_gene_match(self, bio_query, bio_db=None):
        return self.get_biological_instance('genes', bio_query, bio_db)

    def get_compound_match(self, bio_query, bio_db=None):
        return self.get_biological_instance('compounds', bio_query, bio_db)

    def add_counters(self):
        self.instance_counter = 0
        for inst_type in [self.get_genes_all(), self.get_proteins_all(), self.get_compounds_all(),
                          self.get_reactions_all()]:
            for inst in inst_type:
                inst.internal_id = self.instance_counter
                self.instance_counter += 1

    def write_to_output_file(self, inst, outfile):
        instance_info = inst.export_all_info()
        instance_info = instance_info.replace('\n', '')
        if instance_info:
            outfile.write(f'internal_id:{inst.internal_id}\t{instance_info}\n')

    def output_results(self):
        self.add_counters()
        print('Exporting to spreadsheets')
        with open(self.output_folder + 'Genes.tsv', 'w+') as outfile:
            for inst in self.get_genes_all():
                self.write_to_output_file(inst, outfile)
        with open(self.output_folder + 'Proteins.tsv', 'w+') as outfile:
            for inst in self.get_proteins_all():
                self.write_to_output_file(inst, outfile)
        with open(self.output_folder + 'Compounds.tsv', 'w+') as outfile:
            for inst in self.get_compounds_all():
                self.write_to_output_file(inst, outfile)
        with open(self.output_folder + 'Reactions.tsv', 'w+') as outfile:
            for inst in self.get_reactions_all():
                self.write_to_output_file(inst, outfile)
        self.export_graph()

    def export_graph(self):
        print('Exporting graph.tsv')
        with open(self.output_folder + 'Graph.sif', 'w+') as outfile:
            outfile.write('SOURCE\tINTERACTION\tTARGET\n')
            for gene in self.get_genes_all():
                edges = gene.export_graph_edges()
                outfile.write(edges)

            for protein in self.get_proteins_all():
                edges = protein.export_graph_edges()
                outfile.write(edges)

            for compound in self.get_compounds_all():
                edges = compound.export_graph_edges()
                outfile.write(edges)
            for reaction in self.get_reactions_all():
                edges = reaction.export_graph_edges()
                outfile.write(edges)


if __name__ == '__main__':
    gs = Global_Searcher()
    # p=Protein({'kegg':'1'})
    # r=Reaction({'kegg':'2'})
    # p.set_detail('reaction_instances',r,converged_in='metacyc')
    # r.set_detail('protein_instances',p,converged_in='kegg')
    # print(r.get_detail('protein_instances',return_convergence=True))
    # print(p.get_detail('reaction_instances',return_convergence=True))
    # b.save_instance_neo4j(r)
    # b.save_instance_neo4j(p)

    # b.save_instance_neo4j(r)
    # p1=b.get_best_match_db({'enzyme_ec':'1.1.1.346'},node_type='Protein')
    # print(p1)
