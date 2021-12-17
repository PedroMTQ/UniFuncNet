import re
import sys
import os
import subprocess
from pathlib import Path
import argparse


if sys.platform.startswith('win'):
    SPLITTER = '\\'
else:
    SPLITTER = '/'

DRAX_FOLDER = os.path.abspath(os.path.dirname(__file__)).split(SPLITTER)[0:-2]
DRAX_FOLDER = SPLITTER.join(DRAX_FOLDER) + SPLITTER
RESOURCES_FOLDER=f'{DRAX_FOLDER}Resources{SPLITTER}'



class Compounds_to_Organisms_Mapping():

    def __init__(self,input_samples,metabolites,output_folder,database=None):
            if not input_samples.endswith('/'):input_samples+='/'
            self.input_folder=input_samples
            if not output_folder.endswith('/'):output_folder+='/'
            self.output_folder=output_folder
            self.metabolites=metabolites

            self.unwanted_mantis_dbs = ['nog', 'ncbi', 'tcdb']
            self.mantis_env='mantis_env'
            self.drax_env='drax_env'

            self.drax_input = f'{self.output_folder}drax_input.tsv'
            self.drax_output = f'{self.output_folder}drax_output{SPLITTER}'
            self.mantis_input = f'{self.output_folder}mantis_input.tsv'
            self.mantis_output = f'{self.output_folder}mantis_output{SPLITTER}'
            self.workflow_output = f'{self.output_folder}workflow_output{SPLITTER}'
            self.workflow_console_out = f'{self.output_folder}console.out'
            self.output_report = f'{self.workflow_output}Report.tsv'
            self.report={}

            self.conda_prefix = self.get_conda_prefix()
            self.database=database
            for p in [self.output_folder,self.drax_output,self.mantis_output,self.workflow_output]:
                Path(p).mkdir(parents=True, exist_ok=True)
            self.workflow()

    ###### running tools

    def get_conda_prefix(self):
        current_prefix = str(subprocess.run('echo $CONDA_PREFIX', shell=True, stdout=subprocess.PIPE).stdout)
        current_prefix = current_prefix.split("'")[1].strip('\\n')
        try:
            base_prefix = current_prefix[:re.search('envs/', current_prefix).span()[0]]
        except:
            base_prefix=f'{current_prefix}{SPLITTER}'
        return base_prefix

    def create_mantis_config(self,mantis_folder):
        with open(f'{mantis_folder}MANTIS.config','w+') as file:
            for db in self.unwanted_mantis_dbs:
                file.write(f'{db}_ref_folder=NA\n')

    def run_mantis_setup(self):
        mantis_folder=f'{RESOURCES_FOLDER}mantis{SPLITTER}'
        Path(mantis_folder).mkdir(parents=True, exist_ok=True)
        if not os.listdir(mantis_folder):
            mantis_url = 'https://github.com/PedroMTQ/mantis.git'
            #with open(self.workflow_console_out, 'a+') as file:
            #    subprocess.run(f'git clone {mantis_url} {mantis_folder}',shell=True, stdout=file,stderr=file)
            subprocess.run(f'git clone {mantis_url} {mantis_folder}',shell=True)
        activate_mantis_env = f'. {self.conda_prefix}/etc/profile.d/conda.sh && conda activate {self.mantis_env}'
        process = subprocess.run(activate_mantis_env, shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if process.stderr:
            print('Could not find mantis_env environment, creating environment')
            conda_create_env_command = f'conda env create -f {mantis_folder}mantis_env.yml'
            #with open(self.workflow_console_out, 'a+') as file:
            #    subprocess.run(conda_create_env_command, shell=True,stdout=file,stderr=file)
            subprocess.run(conda_create_env_command, shell=True)
        else:
            pass
        self.create_mantis_config(mantis_folder)
        mantis_setup_command = f'. {self.conda_prefix}/etc/profile.d/conda.sh && conda activate {self.mantis_env} && python {mantis_folder} setup_databases'
        #with open(self.workflow_console_out, 'a+') as file:
        #    subprocess.run(mantis_setup_command,shell=True,stdout=file,stderr=file)
        subprocess.run(mantis_setup_command,shell=True)
        mantis_setup_command = f'. {self.conda_prefix}/etc/profile.d/conda.sh && conda activate {self.mantis_env} && python {mantis_folder} check_installation'
        #with open(self.workflow_console_out, 'a+') as file:
        #    subprocess.run(mantis_setup_command,shell=True,stdout=file,stderr=file)
        subprocess.run(mantis_setup_command,shell=True)

    def create_mantis_input(self):
        res=0
        with open(self.mantis_input,'w+') as file:
            for model_fasta in os.listdir(self.input_folder):
                model_in_path=f'{self.input_folder}{model_fasta}'
                model_out=model_fasta.split('.')[0]
                mantis_output=f'{self.mantis_output}{model_out}{SPLITTER}consensus_annotation.tsv'
                if not os.path.exists(mantis_output):
                    line=[model_out,model_in_path]
                    line='\t'.join(line)+'\n'
                    file.write(line)
                    res+=1
        return res

    def run_mantis(self):
        mantis_folder=f'{RESOURCES_FOLDER}mantis{SPLITTER}'
        n_input=self.create_mantis_input()
        if n_input:
            mantis_setup_command = f'. {self.conda_prefix}/etc/profile.d/conda.sh && conda activate {self.mantis_env} && python {mantis_folder} run_mantis -t {self.mantis_input} -o {self.mantis_output} -da heuristic'
            subprocess.run(mantis_setup_command,shell=True)

    def run_drax(self):
        self.create_mantis_input()
        if not os.listdir(self.drax_output):
            mantis_setup_command = f'. {self.conda_prefix}/etc/profile.d/conda.sh && conda activate {self.drax_env} && python {DRAX_FOLDER} -i {self.drax_input} -o {self.drax_output}'
            subprocess.run(mantis_setup_command,shell=True)
        else:
            print(f'We found files in {self.mantis_output}, so DRAX will not run again')

    def compile_input_drax(self):
        with open(self.drax_input, 'w+') as out_file:
            with open(self.metabolites) as met_file:
                for line in met_file:
                    line=line.strip('\n').split('\t')
                    line_dict={}
                    for db_id in line:
                        try:
                            int(db_id)
                            db_type='chebi'
                        except:
                            db_type='synonyms'
                        outline=['compound','crp,prc',db_type,db_id,'\n']
                        outline='\t'.join(outline)
                        out_file.write(outline)


    ###### parsing functions

    def get_metabolites(self):
        res=[]
        with open(self.metabolites) as met_file:
            for line in met_file:
                line = line.strip('\n').split('\t')
                res.append(line)
        return res

    def read_drax_tsv(self,drax_tsv):
        res = {}
        with open(drax_tsv) as file:
            line = file.readline()
            while line:
                line = line.strip('\n').split('\t')
                internal_id = line.pop(0).split(':')[1]
                for i in line:
                    id_type = i.split(':')[0]
                    annot = i.replace(id_type + ':', '')
                    if internal_id not in res: res[internal_id] = {}
                    if id_type not in res[internal_id] and id_type != 'reaction_compounds': res[internal_id][
                        id_type] = set()
                    if id_type == 'reaction_compounds':
                        res[internal_id][id_type] = annot
                    else:
                        res[internal_id][id_type].add(annot)
                line = file.readline()
        return res

    def extract_mantis_ids(self,mantis_annotations):
        res = {}
        with open(mantis_annotations) as file:
            line = file.readline()
            while line:
                line = line.strip('\n')
                line = line.split('\t')
                protein_annotations = line[6:]
                for annot in protein_annotations:
                    if annot.startswith('enzyme_ec') or \
                            annot.startswith('kegg_reaction') or \
                            annot.startswith('kegg_ko'):
                        id_type, annotation = annot.split(':')
                        if id_type not in res: res[id_type] = set()
                        res[id_type].add(annotation)
                line = file.readline()
        return res

    ###### mapping functions

    def get_mapped_metabolites(self,compounds_info, metabolites):
        res = {}
        for cpd_id in compounds_info:
            cpd_info = set()
            if 'chebi' in compounds_info[cpd_id]:
                cpd_info.update(compounds_info[cpd_id]['chebi'])
            if 'synonyms' in compounds_info[cpd_id]:
                cpd_info.update(compounds_info[cpd_id]['synonyms'])
            for met in metabolites:
                if cpd_info.intersection(met):
                    res[cpd_id] = met
        mapped = res.values()
        unmapped = []
        for met in metabolites:
            if met not in mapped:
                unmapped.append(met)
        outline=f'DRAX managed to gather information on {len(mapped)} metabolites, but failed to do so for {len(unmapped)} metabolites.'
        print(outline)
        return res

    def get_mapped_reactions(self, reactions_info, mapped_internal):
        res = {}

        compounds = set()
        for reaction_id in reactions_info:
            reaction_cpds = reactions_info[reaction_id]['reaction_compounds']
            if ' => ' in reaction_cpds: reaction_cpds = reaction_cpds.replace(' => ', ' <=> ')
            if ' <= ' in reaction_cpds: reaction_cpds = reaction_cpds.replace(' <= ', ' <=> ')
            drax_reactants, drax_products = reaction_cpds.split('<=>')
            drax_reactants, drax_products = [j.strip() for j in drax_reactants.split('+')], [j.strip() for j in
                                                                                             drax_products.split('+')]
            reaction_cpds = set(drax_reactants + drax_products)
            metabolites_in_reactions = reaction_cpds.intersection(mapped_internal)
            if metabolites_in_reactions:
                res[reaction_id] = metabolites_in_reactions
            compounds.update(metabolites_in_reactions)
        outline=f'DRAX managed to link {len(compounds)} metabolites to a total of {len(res)} reactions.'
        print(outline)
        return res

    def get_mapped_proteins(self,linked_reactions, reactions_info):
        res = {}
        unconnected = set()
        total_proteins_connected=set()
        for reaction_id in linked_reactions:
            if 'proteins_connected' in reactions_info[reaction_id]:
                proteins_connected = set(reactions_info[reaction_id]['proteins_connected'])
                res[reaction_id]=proteins_connected
                total_proteins_connected.update(proteins_connected)
            else:
                unconnected.add(reaction_id)
        outline=f'DRAX managed to link {len(res)} reactions to {len(total_proteins_connected)} proteins, but failed to do so for {len(unconnected)} reactions.'
        print(outline)
        return res




    def get_mapped_annotations(self,linked_proteins, proteins_info):
        res = {}
        for reaction_id in linked_proteins:
            current_proteins=linked_proteins[reaction_id]
            for protein_id in current_proteins:
                res[protein_id] = set()
                if 'enzyme_ec' in proteins_info[protein_id]:
                    res[protein_id].update(proteins_info[protein_id]['enzyme_ec'])
                if 'kegg_ko' in proteins_info[protein_id]:
                    res[protein_id].update(proteins_info[protein_id]['kegg_ko'])
        return res

    def get_mapped_organisms(self,mapped_annotations, mantis_annotations):
        res = {}
        for ma in mantis_annotations:
            current_annotations = self.extract_mantis_ids(ma)
            kegg_ko = current_annotations['kegg_ko']
            enzyme_ec = current_annotations['enzyme_ec']
            organism = ma.split('/')[-2]
            res[organism] = {'annotations': set(), 'protein_ids': set()}
            for protein_id in mapped_annotations:
                protein_annotations = mapped_annotations[protein_id]
                ec_intersection = protein_annotations.intersection(enzyme_ec)
                ko_intersection = protein_annotations.intersection(kegg_ko)
                if ec_intersection:
                    res[organism]['annotations'].update(ec_intersection)
                    res[organism]['protein_ids'].add(protein_id)
                if ko_intersection:
                    res[organism]['annotations'].update(ko_intersection)
                    res[organism]['protein_ids'].add(protein_id)
        return res

    def get_unique_annotations(self,mapped_organisms):
        non_unique = set()
        unique = set()
        for organism in mapped_organisms:
            for annot in mapped_organisms[organism]['annotations']:
                if annot not in non_unique and annot not in unique:
                    unique.add(annot)
                elif annot in unique:
                    unique.remove(annot)
                    non_unique.add(annot)
        outline=f'After parsing the functional annotation from Mantis, we found that there were {len(unique)} unique annotations, and {len(non_unique)} non-unique annotations.'
        print(outline)

    def get_proteins_to_reactions(self,linked_reactions, organisms_protein_ids, reactions_info):
        res = set()
        for reaction_id in linked_reactions:
            if 'proteins_connected' in reactions_info[reaction_id]:
                proteins_connected = set(reactions_info[reaction_id]['proteins_connected'])
                if proteins_connected.intersection(organisms_protein_ids):
                    res.add(reaction_id)
        return res

    def get_organisms_to_metabolites(self,mapped_organisms, mapping, metabolites, linked_reactions, reactions_info):
        for organism in mapped_organisms:
            organisms_protein_ids = mapped_organisms[organism]['protein_ids']
            proteins_to_reactions = self.get_proteins_to_reactions(linked_reactions, organisms_protein_ids, reactions_info)
            for reaction_id in proteins_to_reactions:
                reaction_cpds = reactions_info[reaction_id]['reaction_compounds']
                if ' => ' in reaction_cpds: reaction_cpds = reaction_cpds.replace(' => ', ' <=> ')
                if ' <= ' in reaction_cpds: reaction_cpds = reaction_cpds.replace(' <= ', ' <=> ')
                drax_reactants, drax_products = reaction_cpds.split('<=>')
                drax_reactants, drax_products = [j.strip() for j in drax_reactants.split('+')], [j.strip() for j in drax_products.split('+')]
                reaction_cpds = set(drax_reactants + drax_products)
                for cpd_id in reaction_cpds:
                    if cpd_id in mapping:
                        metabolites_info = mapping[cpd_id]
                        metline=", ".join(metabolites_info)
                        if organism not in self.report: self.report[organism]={}
                        if metline not in self.report[organism]: self.report[organism][metline]=0
                        self.report[organism][metline]+=1
                        #print(organism, reaction_id,cpd_id,metabolites_info)


    def output_graph(self,linked_reactions,linked_proteins,mapped_organisms):
        output_sif = f'{self.workflow_output}Graph.sif'
        with open(output_sif,'w+') as file:
            file.write('SOURCE\tINTERACTION\tTARGET\n')

            for reaction_id in linked_reactions:
                metabolites=linked_reactions[reaction_id]
                for cpd_id in metabolites:
                    outline=f'{cpd_id}\tcr\t{reaction_id}\n'
                    file.write(outline)
            for reaction_id in linked_proteins:
                proteins=linked_proteins[reaction_id]
                for protein_id in proteins:
                    outline=f'{reaction_id}\trp\t{protein_id}\n'
                    file.write(outline)
            for organism in mapped_organisms:
                proteins = mapped_organisms[organism]['protein_ids']
                for protein_id in proteins:
                    outline = f'{protein_id}\tps\t{organism}\n'
                    file.write(outline)

    def output_report_tsv(self):
        headers=['Sample']
        for organism in self.report:
            for metline in self.report[organism]:
                if metline not in headers: headers.append(metline)
        headers='\t'.join(headers)+'\n'
        with open(self.output_report, 'w+') as file:
            file.write(headers)
            for organism in self.report:
                line=[organism]
                for metline in self.report[organism]:
                    line.append(str(self.report[organism][metline]))
                line='\t'.join(line)+'\n'
                file.write(line)

    def output_results(self):
        compounds = f'{self.drax_output}Compounds.tsv'
        genes = f'{self.drax_output}Genes.tsv'
        proteins = f'{self.drax_output}Proteins.tsv'
        reactions = f'{self.drax_output}Reactions.tsv'
        compounds_info = self.read_drax_tsv(compounds)
        reactions_info = self.read_drax_tsv(reactions)
        proteins_info = self.read_drax_tsv(proteins)
        metabolites = self.get_metabolites()
        organisms_annotations= [f'{self.mantis_output}{org}{SPLITTER}consensus_annotation.tsv' for org in os.listdir(self.mantis_output)]


        mapping = self.get_mapped_metabolites(compounds_info, metabolites)
        mapped_internal = set(mapping.keys())
        linked_reactions = self.get_mapped_reactions(reactions_info, mapped_internal)
        linked_proteins = self.get_mapped_proteins(linked_reactions, reactions_info)
        mapped_annotations = self.get_mapped_annotations(linked_proteins, proteins_info)
        mapped_organisms = self.get_mapped_organisms(mapped_annotations, organisms_annotations)

        self.get_unique_annotations(mapped_organisms)
        organisms_to_metabolites = self.get_organisms_to_metabolites(mapped_organisms, mapping, metabolites,linked_reactions, reactions_info)
        self.output_graph(linked_reactions,linked_proteins,mapped_organisms)
        self.output_report_tsv()


    ###### main workflow

    def workflow(self):
        self.run_mantis_setup()
        self.run_mantis()
        self.compile_input_drax()
        self.run_drax()
        self.output_results()

if __name__ == '__main__':
    if True:
        Compounds_to_Organisms_Mapping(input_samples='/home/pedroq/Desktop/test_mapping/samples/', metabolites='/home/pedroq/Desktop/test_mapping/metabolites.tsv',output_folder='/home/pedroq/Desktop/test_mapping/out',database=None)
    else:
        print('Executing command:\n', ' '.join(sys.argv))
        parser = argparse.ArgumentParser(description='This workflow suggests new connections for Carveme metabolic models\n',formatter_class=argparse.RawTextHelpFormatter)
        parser.add_argument('-i','--input_samples', help='[required]\tInput folder with protein sequences fastas')
        parser.add_argument('-m', '--metabolites', help='[required]\tMetabolites list (synonyms or ChEBI IDs), with each metabolite in a separate line')
        parser.add_argument('-o', '--output_folder', help='[required]\tOutput directory')
        parser.add_argument('-db','--database', help='[optional]\tDatabases to be used in DRAX')
        args = parser.parse_args()
        input_samples = args.input_samples
        metabolites = args.metabolites
        output_folder = args.output_folder
        database = args.database
        if input_folder and output_folder:
            Compounds_to_Organisms_Mapping(input_folder=input_folder,metabolites=metabolites,output_folder=output_folder,database=database)
        else:
            print('Missing input and output folders')