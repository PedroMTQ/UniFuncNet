import networkx as nx
import sys
from pathlib import Path
import os
import re
import subprocess
from sys import platform

if platform.startswith('win'):
    SPLITTER = '\\'
else:
    SPLITTER = '/'


class GsmmExpansion:
    def __init__(self,
                 input_folder,
                 output_folder,
                 metacyc_ref=None,
                 database=None,
                 only_connected=False,
                 politeness_timer=10,
                 cpds_to_ignore={}):
        if not input_folder.endswith('/'):
            input_folder += '/'
        self.input_folder = input_folder
        if not output_folder.endswith('/'):
            output_folder += '/'
        self.output_folder = output_folder
        self.metacyc_ref = metacyc_ref
        self.database = database
        self.cpds_to_ignore = cpds_to_ignore
        self.mantis_cfg = f'{output_folder}mantis.cfg'
        if database:
            if 'metacyc' not in database:
                self.metacyc_ref = None
        self.carveme_models = f'{self.output_folder}models{SPLITTER}'
        self.unifuncnet_input = f'{self.output_folder}unifuncnet_input.tsv'
        self.unifuncnet_output = f'{self.output_folder}unifuncnet_output{SPLITTER}'
        self.mantis_input = f'{self.output_folder}mantis_input.tsv'
        self.mantis_output = f'{self.output_folder}mantis_output{SPLITTER}'
        self.workflow_output = f'{self.output_folder}workflow_output{SPLITTER}'
        self.workflow_console_out = f'{self.output_folder}console.out'
        self.output_report = f'{self.workflow_output}Report.txt'
        if os.path.exists(self.output_report): os.remove(self.output_report)

        self.unwanted_mantis_dbs = ['nog', 'ncbi', 'tcdb']
        self.carveme_env = 'carveme_env'
        self.conda_prefix = self.get_conda_prefix()
        self.compounds_match = {'searched': set(), 'matched': {}}
        # only expand network if some of the nodes are connected to the initial gsmm
        self.only_connected = only_connected
        self.politeness_timer = politeness_timer
        for p in [self.output_folder, self.carveme_models, self.unifuncnet_output, self.mantis_output,
                  self.workflow_output]:
            Path(p).mkdir(parents=True, exist_ok=True)
        self.workflow()

    ###### parsing model

    def extract_model_ids_proteins_and_reactions(self, model_file, verbose=True):
        '''
        listOfSpecies = compounds
        listOfReactions = reactions
        listOfGeneProducts = gene matches
        '''
        res = {}
        record = False
        counter = {}
        self.db_type_model_pr = set()

        with open(model_file) as file:
            line = file.readline()
            while line:
                line = line.strip('\n')
                if '<reaction metaid' in line:
                    reaction_id = line.split('metaid="')[1].split('"')[0]
                    res[reaction_id] = {}
                    record = True
                elif 'rdf:resource' in line and record:

                    db_type = line.split('/')[-3]
                    identifier = line.split('/')[-2].strip('"')

                    if db_type == 'ec-code':
                        db_type = 'enzyme_ec'
                    elif db_type == 'kegg.reaction':
                        db_type = 'kegg_reaction'
                    elif db_type == 'reactome':
                        db_type = 'reactome_reaction'
                    elif db_type == 'rhea':
                        db_type = 'rhea_reaction'
                        identifier = identifier.split('#')[0]
                    elif db_type == 'biocyc':
                        db_type = 'metacyc_reaction'
                    elif db_type == 'bigg.reaction' or db_type == 'seed.reaction' or db_type == 'metanetx.reaction' or db_type == 'SBO':
                        identifier = None
                    else:
                        identifier = None
                    self.db_type_model_pr.add(db_type)
                    if identifier:
                        if db_type not in res[reaction_id]:
                            res[reaction_id][db_type] = set()
                        res[reaction_id][db_type].add(identifier)
                elif '</reaction>' in line:
                    record = False

                line = file.readline()
        c = 0
        for reaction_id in res:
            if 'reactome_reaction' in res[reaction_id] and 'enzyme_ec' not in res[
                reaction_id] and 'kegg_reaction' not in res[reaction_id] and 'rhea_reaction' not in res[
                reaction_id] and 'metacyc_reaction' not in res[reaction_id]:
                c += 1

        all_ids = {}
        for reaction_id in res:
            for id_type in res[reaction_id]:
                if id_type not in all_ids:
                    all_ids[id_type] = set()
                all_ids[id_type].update(res[reaction_id][id_type])
        if verbose:
            with open(self.output_report, 'a+') as file:
                outline = f'####################################################################################################\nPre-processing {model_file}\n####################################################################################################\n'
                file.write(outline)
                outline = f'Excluding {c} reactions due to it only having reactome ids\n'
                file.write(outline)
                outline = f'This model has a total of {len(res)} reactions, ' \
                          f'{len([i for i in res if "kegg_reaction" in res[i]])} of these with a KEGG ID, ' \
                          f'{len([i for i in res if "enzyme_ec" in res[i]])} of these with an enzyme EC ID, ' \
                          f'{len([i for i in res if "rhea_reaction" in res[i]])} of these with a Rhea ID, ' \
                          f'{len([i for i in res if "metacyc_reaction" in res[i]])} of these with a Metacyc ID, and ' \
                          f'{len([i for i in res if "reactome_reaction" in res[i]])} of these with a Reactome ID.\n'
                file.write(outline)

        return all_ids

    def extract_model_ids_compounds(self, model_file):
        res = {}
        self.db_type_model_cpd = set()
        record = False
        with open(model_file) as file:
            line = file.readline()
            while line:
                line = line.strip('\n')
                if '<species metaid' in line:
                    compound_id = line.split('metaid="')[1].split('"')[0]
                    res[compound_id] = {}
                    record = True
                elif 'rdf:resource' in line and record:
                    db_type = line.split('/')[-3]
                    identifier = line.split('/')[-2].strip('"')
                    db_type = db_type.lower()
                    if 'kegg' in db_type:
                        db_type = 'kegg'
                    elif db_type == 'inchikey':
                        db_type = 'inchi_key'
                    elif 'seed' in db_type:
                        db_type = 'seed'
                    elif 'metanetx' in db_type:
                        db_type = 'metanetx'
                    if db_type not in res[compound_id]: res[compound_id][db_type] = set()
                    res[compound_id][db_type].add(identifier)
                    self.db_type_model_cpd.add(db_type)
                elif '</species>' in line:
                    record = False
                line = file.readline()
        compound_dict = {}
        for i in res:
            if res[i]:
                compound_dict[i] = res[i]
        with open(self.output_report, 'a+') as file:
            outline = f'This model has a total of {len(res)} compounds, with {len(res) - len(compound_dict)} not having any metadata to search with UniFuncNet\n'
            file.write(outline)
        return compound_dict

    ###### generating unifuncnet input

    def extract_mantis_ids(self, mantis_annotations):
        res = {}
        with open(mantis_annotations) as file:
            line = file.readline()
            while line:
                line = line.strip('\n')
                line = line.split('\t')
                protein_annotations = line[6:]
                for annot in protein_annotations:
                    if annot.startswith('enzyme_ec') or annot.startswith('metacyc'):
                        id_type = annot.split(':')[0]
                        annotation = annot[len(id_type) + 1:]
                        if id_type not in res: res[id_type] = set()
                        res[id_type].add(annotation)
                line = file.readline()
        return res

    def ids_to_run_unifuncnet(self, mantis_annotations, model_file):
        mantis_ids = self.extract_mantis_ids(mantis_annotations)
        model_ids = self.extract_model_ids_proteins_and_reactions(model_file)
        res = {}
        for id_type in mantis_ids:
            if id_type not in res: res[id_type] = set()
            if id_type in model_ids:
                for annot in mantis_ids[id_type]:
                    if annot not in model_ids[id_type]:
                        res[id_type].add(annot)
            else:
                res[id_type] = mantis_ids[id_type]
        with open(self.output_report, 'a+') as file:
            outline = 'Number of IDs unique to Mantis:\n'
            for id_type in res:
                outline += f'{id_type}:{len(res[id_type])}\n'
            file.write(outline)
        return res

    def merge_ids_to_run(self, ids_to_run_list):
        res = {}
        for ids_to_run in ids_to_run_list:
            for id_type in ids_to_run:
                if id_type not in res: res[id_type] = set()
                res[id_type].update(ids_to_run[id_type])
        return res

    def compile_input_unifuncnet(self):
        print('Compiling UniFuncNet input tsv')
        unifuncnet_ids = []
        for model in os.listdir(self.carveme_models):
            if model.endswith('.xml'):
                model_name = model.split('.')[0]
                model_path = f'{self.carveme_models}{model}'
                mantis_annotations = f'{self.mantis_output}{model_name}{SPLITTER}consensus_annotation.tsv'
                model_ids = self.ids_to_run_unifuncnet(mantis_annotations, model_path)
                unifuncnet_ids.append(model_ids)
        unifuncnet_ids = self.merge_ids_to_run(unifuncnet_ids)
        with open(self.unifuncnet_input, 'w+') as file:
            for id_type in unifuncnet_ids:
                if id_type in ['enzyme_ec', 'metacyc']:
                    for annot in unifuncnet_ids[id_type]:
                        file.write(f'{annot}\t{id_type}\tprotein\tprc\n')

    ###### network analysis

    def create_network_model(self, dict_reactions):
        G = nx.DiGraph()
        for reaction_id in dict_reactions:
            for r in dict_reactions[reaction_id]['reactants']:
                G.add_edge(r, reaction_id)
            for p in dict_reactions[reaction_id]['products']:
                G.add_edge(reaction_id, p)

        return G

    def create_network_unifuncnet(self, reactions_dict):
        compounds_match = self.compounds_match['matched']
        G = nx.DiGraph()
        for reaction_id in reactions_dict:
            unifuncnet_reaction_id = f'RD_{reaction_id}'
            if 'reaction_compounds' in reactions_dict[reaction_id]:
                reaction_cpds = reactions_dict[reaction_id]['reaction_compounds']
                if ' => ' in reaction_cpds: reaction_cpds = reaction_cpds.replace(' => ', ' <=> ')
                if ' <= ' in reaction_cpds: reaction_cpds = reaction_cpds.replace(' <= ', ' <=> ')
                unifuncnet_reactants, unifuncnet_products = reaction_cpds.split('<=>')
                unifuncnet_reactants, unifuncnet_products = [j.strip() for j in unifuncnet_reactants.split('+')], [
                    j.strip() for j in unifuncnet_products.split('+')]
                for reactant in unifuncnet_reactants:
                    if reactant in compounds_match:
                        reactant = compounds_match[reactant]
                    else:
                        reactant = {reactant}
                    for product in unifuncnet_products:
                        if product in compounds_match:
                            product = compounds_match[product]
                        else:
                            product = {product}

                        for r in reactant:
                            G.add_edge(r, unifuncnet_reaction_id)
                        for p in product:
                            G.add_edge(unifuncnet_reaction_id, p)

        return G

    def create_network_expanded(self, reactions_model, reactions_unifuncnet, cpd_ids_to_ignore):
        compounds_match = self.compounds_match['matched']
        rejected_reactions = 0
        G = nx.DiGraph()
        for reaction_id in reactions_model:
            for r in reactions_model[reaction_id]['reactants']:
                G.add_edge(r, reaction_id)
            for p in reactions_model[reaction_id]['products']:
                G.add_edge(reaction_id, p)

        for reaction_id in reactions_unifuncnet:
            unifuncnet_reaction_id = f'RD_{reaction_id}'
            if 'reaction_compounds' in reactions_unifuncnet[reaction_id]:

                reaction_cpds = reactions_unifuncnet[reaction_id]['reaction_compounds']
                if ' => ' in reaction_cpds: reaction_cpds = reaction_cpds.replace(' => ', ' <=> ')
                if ' <= ' in reaction_cpds: reaction_cpds = reaction_cpds.replace(' <= ', ' <=> ')
                unifuncnet_reactants, unifuncnet_products = reaction_cpds.split('<=>')
                unifuncnet_reactants, unifuncnet_products = [j.strip() for j in unifuncnet_reactants.split('+')], [
                    j.strip() for j in unifuncnet_products.split('+')]
                if not self.only_connected:
                    connected = True
                else:
                    connected = False
                    for unifuncnet_cpd_id in unifuncnet_reactants + unifuncnet_products:
                        if unifuncnet_cpd_id in compounds_match:
                            connected = True
                if connected:
                    for reactant in unifuncnet_reactants:
                        if reactant in compounds_match:
                            reactant = compounds_match[reactant]
                        else:
                            reactant = {reactant}

                    for product in unifuncnet_products:
                        if product in compounds_match:
                            product = compounds_match[product]
                        else:
                            product = {product}

                    for r in reactant:
                        if r not in cpd_ids_to_ignore:
                            G.add_edge(r, unifuncnet_reaction_id)
                    for p in product:
                        if p not in cpd_ids_to_ignore:
                            G.add_edge(unifuncnet_reaction_id, p)
                else:
                    rejected_reactions += 1
            else:
                rejected_reactions += 1

        return G, rejected_reactions

    def check_dead_end_metabolites(self, graph):
        products = set()
        reactants = set()
        for n1, n2 in graph.edges():
            # reaction -> product
            if (n1.startswith('R_') or n1.startswith('RD_')) and not n2.endswith('_e'):
                products.add(n2)
        for n1, n2 in graph.edges():
            # reactant -> reaction
            if (n2.startswith('R_') or n2.startswith('RD_')) and not n1.endswith('_e'):
                reactants.add(n1)
        transported = set()
        for n1, n2 in graph.edges():
            if n1.endswith('_e'):
                transported.add(n1)
            if n2.endswith('_e'):
                transported.add(n2)

        res = set()

        for product in products:
            if not n1.endswith('_e'):
                if product not in reactants:
                    res.add(product)
        for reactant in reactants:
            if not n1.endswith('_e'):
                if reactant not in products:
                    res.add(reactant)
        all_metabolites = set(
            [n for n in graph.nodes() if not n.startswith('R_') and not n.startswith('RD_') and not n.endswith('_e')])
        # print('Dead end metabolites', len(res))
        # print('Transported metabolites', len(transported))
        return res, len(all_metabolites)

    def compare_dead_end_metabolites(self, graph1, graph2):
        dead_end1, n_nodes1 = self.check_dead_end_metabolites(graph1)
        dead_end2, n_nodes2 = self.check_dead_end_metabolites(graph2)
        connected = set()
        # original
        for dead in dead_end1:
            # expanded
            if dead not in dead_end2:
                # connected
                connected.add(dead)

        dead_end1_str = ','.join(sorted(dead_end1))
        dead_end2_str = ','.join(sorted(dead_end2))
        with open(self.output_report, 'a+') as file:
            outline = f'Original network dead end metabolites: {len(dead_end1)} ({round(100 * len(dead_end1) / len(graph1.nodes), 3)}%)\n{dead_end1_str}\n' \
                      f'Expanded network dead end metabolites: {len(dead_end2)} ({round(100 * len(dead_end2) / len(graph2.nodes), 3)}%)\n{dead_end2_str}\n' \
                      f'Newly connected dead end metabolites: {len(connected)} ({round(100 * len(connected) / len(graph1.nodes), 3)}%)\n'
            file.write(outline)

    def evaluate_network(self, graph, evaluation_function, network_name):
        dead_end, n_nodes = self.check_dead_end_metabolites(graph)

        subnetworks = sorted(evaluation_function(graph), key=len, reverse=True)
        periplasmic = 0
        extracellular = 0
        cytosol = 0
        size_subnetworks_nodes = {}
        size_subnetworks_reactions = {}
        for node in graph.nodes():
            if node.endswith('_p'):
                periplasmic += 1
            elif node.endswith('_e'):
                extracellular += 1
            elif node.endswith('_c'):
                cytosol += 1
        for subnet in subnetworks:
            nodes = [n for n in graph.subgraph(subnet).nodes() if not n.startswith('R_') and not n.startswith('RD_')]
            reactions = [n for n in graph.subgraph(subnet).nodes() if n.startswith('R_') or n.startswith('RD_')]
            size_subnet_nodes = len(nodes)
            if size_subnet_nodes not in size_subnetworks_nodes: size_subnetworks_nodes[size_subnet_nodes] = 0
            size_subnetworks_nodes[size_subnet_nodes] += 1
            size_subnet_reactions = len(reactions)
            if size_subnet_reactions not in size_subnetworks_reactions: size_subnetworks_reactions[
                size_subnet_reactions] = 0
            size_subnetworks_reactions[size_subnet_reactions] += 1

        n_reactions = len([n for n in graph.nodes() if n.startswith('R_') or n.startswith('RD_')])
        percentage_reactions_largest_component = 100 * max(size_subnetworks_reactions.keys()) / n_reactions
        with open(self.output_report, 'a+') as file:
            outline = f'Percentage of reactions in {network_name} network in largest component: {percentage_reactions_largest_component}%\n'
            file.write(outline)

    ###### merge networks

    def read_unifuncnet_tsv(self, unifuncnet_tsv):
        res = {}
        with open(unifuncnet_tsv) as file:
            line = file.readline()
            while line:
                line = line.strip('\n').split('\t')
                internal_id = line.pop(0).split(':')[1]
                for i in line:
                    id_type = i.split(':')[0]
                    annot = i.replace(id_type + ':', '')
                    if internal_id not in res:
                        res[internal_id] = {}
                    if id_type not in res[internal_id] and id_type != 'reaction_compounds':
                        res[internal_id][
                        id_type] = set()
                    if id_type == 'reaction_compounds':
                        res[internal_id][id_type] = annot
                    else:
                        res[internal_id][id_type].add(annot)
                line = file.readline()
        return res

    def subset_unifuncnet_proteins(self, proteins_unifuncnet, mantis_annotations_tsv):
        res = {}
        mantis_annotations = self.extract_mantis_ids(mantis_annotations_tsv)
        for protein_id in proteins_unifuncnet:
            protein_info = proteins_unifuncnet[protein_id]
            passed = False
            if 'enzyme_ec' in protein_info and 'enzyme_ec' in mantis_annotations:
                if protein_info['enzyme_ec'].intersection(mantis_annotations['enzyme_ec']):
                    res[protein_id] = protein_info
            if 'metacyc' in protein_info and 'metacyc' in mantis_annotations:
                if protein_info['metacyc'].intersection(mantis_annotations['metacyc']):
                    res[protein_id] = protein_info
        return res

    def remove_proteins_without_reaction(self, proteins_dict):
        res = {}
        for i in proteins_dict:
            if 'reactions_connected' in proteins_dict[i]:
                res[i] = proteins_dict[i]
        return res

    def remove_proteins_unifuncnet_in_model(self, model_ids, proteins_dict):
        '''
        since mantis annotations sometimes only had KOs, we now check whether we already had the respective ECs in the model
        '''
        res = {}
        for i in proteins_dict:
            if 'enzyme_ec' in proteins_dict[i]:
                if not proteins_dict[i]['enzyme_ec'].intersection(model_ids['enzyme_ec']):
                    res[i] = proteins_dict[i]
        return res

    def remove_reactions_without_proteins_in_unifuncnet(self, proteins_dict, reactions_dict):
        '''
        some proteins were found by unifuncnet but were already in the model, so now we have to remove those and their respective reactions

        '''
        res = {}
        for i in proteins_dict:
            for r in proteins_dict[i]['reactions_connected']:
                res[r] = reactions_dict[r]
        return res

    def remove_reactions_unifuncnet_in_model_ids(self, model_ids, reactions_dict):
        '''
        some of the reactions from unifuncnet were already in the model, so now we remove them

        '''
        res = {}
        for i in reactions_dict:
            passed = False
            if 'kegg' in reactions_dict[i] and 'kegg_reaction' in model_ids:
                if reactions_dict[i]['kegg'].intersection(model_ids['kegg_reaction']):
                    passed = True
            if 'metacyc' in reactions_dict[i] and 'metacyc_reaction' in model_ids:
                if reactions_dict[i]['metacyc'].intersection(model_ids['metacyc_reaction']):
                    passed = True
            if not passed:
                if 'kegg' in reactions_dict[i] or 'metacyc' in reactions_dict[i]:
                    res[i] = reactions_dict[i]

        return res

    def match_compounds_unifuncnet_model(self, model_file, compounds_dict):
        model_compounds = self.extract_model_ids_compounds(model_file)
        for unifuncnet_internal_id in compounds_dict:
            for model_id in model_compounds:
                if model_id not in self.compounds_match['searched']:
                    for db in self.db_type_model_cpd:
                        if db in compounds_dict[unifuncnet_internal_id] and db in model_compounds[model_id]:
                            if compounds_dict[unifuncnet_internal_id][db].intersection(model_compounds[model_id][db]):
                                if unifuncnet_internal_id not in self.compounds_match['matched']:
                                    self.compounds_match['matched'][unifuncnet_internal_id] = set()
                                self.compounds_match['matched'][unifuncnet_internal_id].add(model_id)
                                self.compounds_match['searched'].add(model_id)

    def check_match_reactions(self, model_compounds, unifuncnet_compounds):
        res = 0
        for r1 in model_compounds:
            r1_set = {r1}
            for r2 in unifuncnet_compounds:
                if r1_set.intersection(r2):
                    res += 1
        return res

    def remove_reactions_unifuncnet_in_model_cpds(self, model_file, reactions_dict):
        compounds_match = self.compounds_match['matched']
        model_data = self.read_model(model_file)
        matched_reactions = {}
        res = {}
        matched = set()
        unmatched = set()
        for r_unifuncnet in reactions_dict:
            if 'reaction_compounds' in reactions_dict[r_unifuncnet]:
                reaction_cpds = reactions_dict[r_unifuncnet]['reaction_compounds']
                if ' => ' in reaction_cpds: reaction_cpds = reaction_cpds.replace(' => ', ' <=> ')
                if ' <= ' in reaction_cpds: reaction_cpds = reaction_cpds.replace(' <= ', ' <=> ')
                unifuncnet_reactants, unifuncnet_products = reaction_cpds.split('<=>')
                unifuncnet_reactants, unifuncnet_products = [j.strip() for j in unifuncnet_reactants.split('+')], [
                    j.strip() for j in unifuncnet_products.split('+')]
                matched_reactants, matched_products = [], []
                for cpd in unifuncnet_reactants:
                    if cpd in compounds_match:
                        matched_reactants.append(compounds_match[cpd])
                        matched.add(cpd)
                    else:
                        unmatched.add(cpd)
                for cpd in unifuncnet_products:
                    if cpd in compounds_match:
                        matched_products.append(compounds_match[cpd])
                        matched.add(cpd)
                    else:
                        unmatched.add(cpd)
                if len(matched_reactants) == len(unifuncnet_reactants) and len(matched_products) == len(
                        unifuncnet_products):
                    for r_model in model_data:
                        model_reactants, model_products = model_data[r_model]['reactants'], model_data[r_model][
                            'products']
                        match_reactants = self.check_match_reactions(model_reactants, matched_reactants)
                        match_products = self.check_match_reactions(model_products, matched_products)
                        if len(model_reactants) and len(model_products):
                            if match_reactants == len(model_reactants) and match_products == len(model_products):
                                matched_reactions[r_unifuncnet] = r_model

                            match_reactants = self.check_match_reactions(model_reactants, matched_products)
                            match_products = self.check_match_reactions(model_products, matched_reactants)
                            if match_reactants == len(model_reactants) and match_products == len(model_products):
                                matched_reactions[r_unifuncnet] = r_model
        # print('Matched compounds', len(matched))
        # print('Unmatched compounds', len(unmatched))
        for i in reactions_dict:
            if i not in matched_reactions:
                res[i] = reactions_dict[i]
        return res

    ###### generating report

    def read_model(self, model_file):
        reactions = {}
        with open(model_file) as file:
            line = file.readline()
            while line:
                line = line.strip('\n')
                if '<reaction metaid' in line:
                    reaction_id = line.split('metaid="')[1].split('"')[0]
                    reactions[reaction_id] = {'reactants': [], 'products': []}
                elif '<listOfReactants' in line:
                    while '</listOfReactants' not in line:
                        line = file.readline()
                        if '<speciesReference species=' in line:
                            reactant = line.split('<speciesReference species="')[1].split('stoichiometry=')[
                                0].strip().strip('"')
                            reactions[reaction_id]['reactants'].append(reactant)
                elif '<listOfProducts' in line:
                    while '</listOfProducts' not in line:
                        line = file.readline()
                        if '<speciesReference species=' in line:
                            reactant = line.split('<speciesReference species="')[1].split('stoichiometry=')[
                                0].strip().strip('"')
                            reactions[reaction_id]['products'].append(reactant)
                line = file.readline()

        return reactions

    ###### running tools

    def get_conda_prefix(self):
        current_prefix = str(subprocess.run('echo $CONDA_PREFIX', shell=True, stdout=subprocess.PIPE).stdout)
        current_prefix = current_prefix.split("'")[1].strip('\\n')
        try:
            base_prefix = current_prefix[:re.search('envs/', current_prefix).span()[0]]
        except:
            base_prefix = f'{current_prefix}{SPLITTER}'

        return base_prefix

    def create_mantis_config(self):
        with open(self.mantis_cfg, 'w+') as file:
            for db in self.unwanted_mantis_dbs:
                file.write(f'{db}_ref_folder=NA\n')

    def create_mantis_config_metacyc(self):
        with open(self.mantis_cfg, 'w+') as file:
            for db in self.unwanted_mantis_dbs:
                file.write(f'{db}_ref_folder=NA\n')
            file.write(f'custom_ref={self.metacyc_ref}\n')

    def run_carveme(self):
        print('Running CarveMe')
        activate_carveme_env = f'. {self.conda_prefix}/etc/profile.d/conda.sh && conda activate {self.carveme_env}'
        process = subprocess.run(activate_carveme_env, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if process.stderr:
            print('Could not find carveme environment!')
            raise Exception
        else:
            for model_fasta in os.listdir(self.input_folder):
                model_in_path = f'{self.input_folder}{model_fasta}'
                model_out = model_fasta.split('.')[0]
                model_out_path = f'{self.carveme_models}{model_out}.xml'
                model_tsv = f'{self.input_folder}{model_out}.tsv'
                if not os.path.exists(model_out_path):
                    print(f'Creating model {model_in_path}')
                    carveme_command = f'. {self.conda_prefix}/etc/profile.d/conda.sh && conda activate {self.carveme_env} && carve {model_in_path} -o {model_out_path}'
                    with open(self.workflow_console_out, 'a+') as file:
                        subprocess.run(carveme_command, shell=True, stdout=file)
                    os.remove(model_tsv)
                else:
                    print(f'Model already exists {model_out_path}')

    def run_mantis_setup(self):
        self.create_mantis_config()
        mantis_setup_command = f'mantis setup -mc {self.mantis_cfg}'
        subprocess.run(mantis_setup_command, shell=True)
        mantis_setup_command = f'mantis check -mc {self.mantis_cfg}'
        subprocess.run(mantis_setup_command, shell=True)
        self.create_mantis_config_metacyc()

    def create_mantis_input(self):
        res = 0
        with open(self.mantis_input, 'w+') as file:
            for model_fasta in os.listdir(self.input_folder):
                model_in_path = f'{self.input_folder}{model_fasta}'
                model_out = model_fasta.split('.')[0]
                mantis_output = f'{self.mantis_output}{model_out}{SPLITTER}consensus_annotation.tsv'
                if not os.path.exists(mantis_output):
                    line = [model_out, model_in_path]
                    line = '\t'.join(line) + '\n'
                    file.write(line)
                    res += 1
        return res

    def run_mantis(self):
        print('Running Mantis')
        n_input = self.create_mantis_input()
        if n_input:
            mantis_setup_command = f'mantis run -i {self.mantis_input} -o {self.mantis_output} -da heuristic -mc {self.mantis_cfg}'
            subprocess.run(mantis_setup_command, shell=True)
        else:
            print(f'Mantis already ran')

    def run_unifuncnet(self):
        print('Running UniFuncNet')
        self.create_mantis_input()
        if not os.listdir(self.unifuncnet_output):
            if self.database:
                run_unifuncnet_command = f'unifuncnet -i {self.unifuncnet_input} -o {self.unifuncnet_output} -db {self.database} -pt {self.politeness_timer} --ignore_cpd'
            else:
                run_unifuncnet_command = f'unifuncnet -i {self.unifuncnet_input} -o {self.unifuncnet_output} -pt {self.politeness_timer} --ignore_cpd'

            subprocess.run(run_unifuncnet_command, shell=True)
        else:
            print(f'We found files in {self.mantis_output}, so UniFuncNet will not run again')

    ###### main workflow

    def get_cpd_ids_to_ignore(self, compounds_unifuncnet):
        res = set()
        for cpd_id in compounds_unifuncnet:
            cpd_info = compounds_unifuncnet[cpd_id]
            for db in cpd_info:
                if db in self.cpds_to_ignore:
                    cpd_ids = set([i.lower() for i in cpd_info[db]])
                    if self.cpds_to_ignore[db].intersection(cpd_ids):
                        res.add(cpd_id)
        print(f'Compounds to ingore:\n{res}')
        return res

    def read_output_unifuncnet(self, model_file, unifuncnet_output, mantis_annotations_tsv):
        '''
        now we have information from unifuncnet where we used the extra mantis annotations as seed
        first we need to match the identifiers from unifuncnet/Proteins to the model ids, and exclude those, since it might be that we didn't find a match before due to the ids being different
        '''
        model_reactions = self.read_model(model_file)
        model_ids = self.extract_model_ids_proteins_and_reactions(model_file, verbose=False)
        compounds_unifuncnet = self.read_unifuncnet_tsv(unifuncnet_output + '/Compounds.tsv')
        cpd_ids_to_ignore = self.get_cpd_ids_to_ignore(compounds_unifuncnet)
        self.match_compounds_unifuncnet_model(model_file, compounds_unifuncnet)
        compounds_match = self.compounds_match['matched']
        with open(self.output_report, 'a+') as file:
            proteins_unifuncnet = self.read_unifuncnet_tsv(unifuncnet_output + '/Proteins.tsv')
            outline = f'All proteins found by UniFuncNet: {len(proteins_unifuncnet)}\n'
            file.write(outline)

            proteins_unifuncnet = self.subset_unifuncnet_proteins(proteins_unifuncnet, mantis_annotations_tsv)
            outline = f'All proteins in sample found by UniFuncNet: {len(proteins_unifuncnet)}\n'
            file.write(outline)

            # now we remove proteins without a reaction
            proteins_unifuncnet = self.remove_proteins_without_reaction(proteins_unifuncnet)
            outline = f'All proteins in sample connected to a reaction: {len(proteins_unifuncnet)}\n'
            file.write(outline)

            proteins_unifuncnet = self.remove_proteins_unifuncnet_in_model(model_ids, proteins_unifuncnet)
            outline = f'All proteins in sample absent in the model: {len(proteins_unifuncnet)}\n'
            file.write(outline)

            reactions_unifuncnet = self.read_unifuncnet_tsv(unifuncnet_output + '/Reactions.tsv')
            outline = f'All reactions found by UniFuncNet: {len(reactions_unifuncnet)}\n'
            file.write(outline)

            reactions_unifuncnet = self.remove_reactions_without_proteins_in_unifuncnet(proteins_unifuncnet,
                                                                                        reactions_unifuncnet)
            outline = f'All reactions connected to proteins in sample: {len(reactions_unifuncnet)}\n'
            file.write(outline)

            reactions_unifuncnet = self.remove_reactions_unifuncnet_in_model_ids(model_ids, reactions_unifuncnet)
            outline = f'All reactions in sample absent in the model (ID matching): {len(reactions_unifuncnet)}\n'
            file.write(outline)

            reactions_unifuncnet = self.remove_reactions_unifuncnet_in_model_cpds(model_file, reactions_unifuncnet)
            outline = f'All reactions in sample absent in the model (compound matching): {len(reactions_unifuncnet)}\n'
            file.write(outline)

            outline = f'UniFuncNet can potentially add {len(reactions_unifuncnet)} new reactions to the model\n'
            file.write(outline)

            model_network = self.create_network_model(model_reactions)
            self.evaluate_network(model_network, nx.weakly_connected_components, 'baseline')

            expanded_network, rejected_reactions = self.create_network_expanded(model_reactions, reactions_unifuncnet,
                                                                                cpd_ids_to_ignore)
            outline = f'We rejected {rejected_reactions} reactions so we only added {len(reactions_unifuncnet) - rejected_reactions} reactions\n'
            file.write(outline)
            self.evaluate_network(expanded_network, nx.weakly_connected_components, 'expanded')
        self.compare_dead_end_metabolites(model_network, expanded_network)

        return model_network, expanded_network

    def get_all_nodes(self, model_network, expanded_network):
        model_edges = set()
        expanded_edges = set()
        for edge in model_network.edges:
            n1, n2 = edge
            if n1.startswith('R_') or n1.startswith('RD_'):
                n1_type = 'r'
            else:
                n1_type = 'c'
            if n2.startswith('R_') or n2.startswith('RD_'):
                n2_type = 'r'
            else:
                n2_type = 'c'
            edge_str = f'{n1}\t{n1_type}{n2_type}\t{n2}'
            model_edges.add(edge_str)
        for edge in expanded_network.edges:
            n1, n2 = edge
            if n1.startswith('R_') or n1.startswith('RD_'):
                n1_type = 'r'
            else:
                n1_type = 'c'
            if n2.startswith('R_') or n2.startswith('RD_'):
                n2_type = 'r'
            else:
                n2_type = 'c'
            edge_str = f'{n1}\t{n1_type}{n2_type}\t{n2}'
            edge_str = edge_str.replace('RD_', '')
            if edge_str not in model_edges:
                expanded_edges.add(edge_str)
        return model_edges, expanded_edges

    def output_graph(self, output_sif, model_edges, expanded_edges):
        with open(output_sif, 'w+') as file:
            file.write('SOURCE\tINTERACTION\tTARGET\tEXPANSION\n')
            for model_e in model_edges:
                line = f'{model_e}\t0\n'
                file.write(line)
            for model_e in expanded_edges:
                line = f'{model_e}\t1\n'
                file.write(line)

    def output_results(self):
        print('Outputting results')
        for model in os.listdir(self.carveme_models):
            model_path = f'{self.carveme_models}{model}'
            model_name = model.replace('.xml', '')
            output_sif = f'{model_name}.sif'
            output_sif_path = f'{self.workflow_output}{model}'.replace('.xml', '.sif')
            mantis_annotations_tsv = f'{self.mantis_output}{model_name}{SPLITTER}consensus_annotation.tsv'
            if output_sif not in os.listdir(self.workflow_output):
                with open(self.output_report, 'a+') as file:
                    outline = f'####################################################################################################\nStarting analysis of {model}\n####################################################################################################\n'
                    file.write(outline)
                model_network, expanded_network = self.read_output_unifuncnet(model_path, self.unifuncnet_output,
                                                                              mantis_annotations_tsv)
                model_edges, expanded_edges = self.get_all_nodes(model_network, expanded_network)
                self.output_graph(output_sif_path, model_edges, expanded_edges)

    def workflow(self):
        self.run_mantis_setup()
        self.run_carveme()
        self.run_mantis()
        self.compile_input_unifuncnet()
        self.run_unifuncnet()
        self.output_results()
