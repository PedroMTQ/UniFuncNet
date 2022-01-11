

from source.Utils.util import get_stoichiometry,standardize_reaction_str,SPLITTER,download_file_ftp,gunzip,RESOURCES_FOLDER,strip_tags
import re
import os
import sqlite3
from string import punctuation


class Metacyc_SQLITE_Connector():
    def __init__(self):
        self.insert_step=5000
        self.info_splitter='##'
        self.db_file = f'{RESOURCES_FOLDER}metacyc.db'
        self.metacyc_folder=f'{RESOURCES_FOLDER}metacyc{SPLITTER}'
        self.db_headers={}
        if os.path.exists(self.db_file):
            self.start_sqlite_cursor()
        else:
            check=self.check_all_resources()
            if check:
                self.create_sql_tables()
                self.start_sqlite_cursor()

            else:
                print(f'Missing metacyc files. Please download them and place them in {self.metacyc_folder}')

    def start_sqlite_cursor(self):
        self.sqlite_connection = sqlite3.connect(self.db_file)
        self.cursor = self.sqlite_connection.cursor()

    def get_headers_generator(self,generator):
        res=set()
        for i_dict in generator:
            res.update(i_dict.keys())
        return res

    def get_db_headers(self,target_table):
        res = set()
        if target_table in self.db_headers: return
        try:
            schema_command = f'PRAGMA table_info({target_table});'
            res_fetch = self.cursor.execute(schema_command).fetchall()
            res_fetch.pop(0)
            for line in res_fetch:
                link_type=line[1]
                res.add(link_type)
        except:
            extra_tables={'TABLECPDRXN':['METACYCRXN'],'TABLEINTRXNIDS':['METACYCPRT','METACYCRXN']}
            generator=None
            if target_table=='COMPOUNDS':
                generator=self.parse_compounds()
            elif target_table=='REACTIONS':
                generator=self.parse_reactions()
            elif target_table=='PROTEINS':
                generator=self.parse_proteins()
            elif target_table=='GENES':
                generator=self.parse_genes()
            if generator:
                res=self.get_headers_generator(generator)
            else:
                res=extra_tables[target_table]

            if 'metacyc' in res:
                res.remove('metacyc')

        self.db_headers[target_table]= sorted(list(res))

    def commit_and_close_sqlite_cursor(self):
        self.sqlite_connection.commit()
        self.sqlite_connection.close()

    def close_sql_connection(self):
        self.sqlite_connection.close()




    def generate_inserts(self, metadata_yielder):
        step=self.insert_step
        temp=[]
        for i in metadata_yielder:
            if len(temp)<step:
                temp.append(i)
            else:
                yield temp
                temp=[]
                temp.append(i)
        yield temp

    def generate_insert_command(self,table_name,main_id_str):
        self.get_db_headers(table_name)
        db_headers=self.db_headers[table_name]
        headers_str=', '.join(db_headers)
        headers_str=f'({main_id_str}, {headers_str})'
        n_values=['?' for i in range(len(db_headers)+1)]
        n_values_str=', '.join(n_values)
        n_values_str=f'({n_values_str})'
        insert_command = f'INSERT INTO {table_name} {headers_str} values {n_values_str}'
        return insert_command


    def store_main_data(self,target_table):
        exceptional_headers=[
            'reaction_stoichiometry',
                             ]

        if target_table == 'COMPOUNDS':
            generator = self.parse_compounds()
        elif target_table == 'REACTIONS':
            generator = self.parse_reactions()
        elif target_table == 'PROTEINS':
            generator = self.parse_proteins()
        elif target_table == 'GENES':
            generator = self.parse_genes()
        insert_command=self.generate_insert_command(target_table,'METACYC')
        metadata_yielder=self.yield_metadata(generator,target_table)
        if metadata_yielder:
            generator_insert = self.generate_inserts(metadata_yielder)
            for table_chunk in generator_insert:
                self.cursor.executemany(insert_command, table_chunk)
            self.sqlite_connection.commit()


    def convert_dict_to_sql(self,ref,row_dict,target_table):
        res=[ref]
        for db in self.db_headers[target_table]:
            if db in row_dict:
                if isinstance(row_dict[db],str): to_add=row_dict[db]
                elif db=='reaction_stoichiometry':
                    to_add=[]
                    stoichiometry=row_dict[db]
                    for i in stoichiometry:
                        if isinstance(i, list):
                            to_add.extend(i)
                        else:
                            to_add.append(i)
                    to_add=[str(i) for i in to_add]
                    to_add = self.info_splitter.join(to_add)
                else: to_add=self.info_splitter.join(row_dict[db])
                res.append(to_add)
            else:
                res.append(None)
        return res

    def yield_metadata(self,generator,target_table):
        for row_dict in generator:
            metacyc_ids=row_dict.pop('metacyc')
            for metacyc_id in metacyc_ids:
                yield self.convert_dict_to_sql(metacyc_id,row_dict,target_table)



    def yield_intermediate_reaction_ids(self):
        reactions_generator = self.parse_reactions()
        proteins_generator = self.parse_proteins()
        res={}
        for current_id,current_generator in [['METACYCPRT',proteins_generator],['METACYCRXN',reactions_generator]]:
            for current_dict in current_generator:
                metacyc_ids=current_dict['metacyc']
                if 'intermediate_reaction_ids' in current_dict:
                    intermediate_reaction_ids=current_dict['intermediate_reaction_ids']
                    for int_reaction_id in intermediate_reaction_ids:
                        if int_reaction_id not in res:
                            res[int_reaction_id]={}
                        if current_id not in res[int_reaction_id]:
                            res[int_reaction_id][current_id]=set()
                        res[int_reaction_id][current_id].update(metacyc_ids)
        #the keys need to be sorted here, so METACYCPRT comes before METACYCRXN
        for int_reaction_id in res:
            protein_ids_str=self.info_splitter.join(res[int_reaction_id]['METACYCPRT'])
            reaction_ids_str=self.info_splitter.join(res[int_reaction_id]['METACYCRXN'])
            yield [int_reaction_id,protein_ids_str,reaction_ids_str]


    def yield_cpd_to_rxn(self):
        generator = self.parse_reactions()
        res={}
        for reaction_dict in generator:
            reaction_stoichiometry=reaction_dict['reaction_stoichiometry']
            metacyc_ids=reaction_dict['metacyc']
            for i in reaction_stoichiometry:
                if isinstance(i,list):
                    cpd_id=i[1]
                    if cpd_id not in res: res[cpd_id]=set()
                    res[cpd_id].update(metacyc_ids)

        for cpd_id in res:
            reaction_ids_str=self.info_splitter.join(res[cpd_id])
            yield [cpd_id,reaction_ids_str]


    def create_main_tables(self):
        for current_table in ['COMPOUNDS','REACTIONS','PROTEINS','GENES']:
            #retrieving headers dynamicaly
            self.get_db_headers(current_table)
            current_db_headers=self.db_headers[current_table]
            #creating table
            create_table_command = f'CREATE TABLE {current_table} (METACYC TEXT, '
            for header in current_db_headers:
                create_table_command += f'{header} TEXT, '
            create_table_command = create_table_command.rstrip(', ')
            create_table_command += ')'
            self.cursor.execute(create_table_command)
            self.sqlite_connection.commit()
            #indexing table
            create_index_command = f'CREATE INDEX {current_table[0]}METACYC_IDX ON {current_table} (METACYC)'
            self.cursor.execute(create_index_command)
            self.sqlite_connection.commit()
            #storing data in table
            self.store_main_data(current_table)

    def create_intermediate_reaction_ids_table(self):
        target_table='TABLEINTRXNIDS'
        #creating new table to connect proteins and reactions
        create_table_command = f'CREATE TABLE {target_table} (' \
                                   f'INTRXNID TEXT,' \
                                   f'METACYCPRT TEXT,' \
                                   f'METACYCRXN  TEXT )'
        self.cursor.execute(create_table_command)
        # indexing table
        create_index_command = f'CREATE INDEX INTRXNID_IDX ON {target_table} (INTRXNID)'
        self.cursor.execute(create_index_command)
        self.sqlite_connection.commit()
        # storing data in table
        insert_command=self.generate_insert_command(target_table,'INTRXNID')
        metadata_yielder=self.yield_intermediate_reaction_ids()
        if metadata_yielder:
            generator_insert = self.generate_inserts(metadata_yielder)
            for table_chunk in generator_insert:
                self.cursor.executemany(insert_command, table_chunk)
            self.sqlite_connection.commit()


    def create_cpd_to_rxn_table(self):
        target_table='TABLECPDRXN'
        #creating new table to connect proteins and reactions
        create_table_command = f'CREATE TABLE {target_table} (' \
                                   f'METACYCCPD TEXT,' \
                                   f'METACYCRXN  TEXT )'
        self.cursor.execute(create_table_command)
        # indexing table
        create_index_command = f'CREATE INDEX METACYCCPD_IDX ON {target_table} (METACYCCPD)'
        self.cursor.execute(create_index_command)
        self.sqlite_connection.commit()
        # storing data in table
        insert_command=self.generate_insert_command(target_table,'METACYCCPD')
        metadata_yielder=self.yield_cpd_to_rxn()
        if metadata_yielder:
            generator_insert = self.generate_inserts(metadata_yielder)
            for table_chunk in generator_insert:
                self.cursor.executemany(insert_command, table_chunk)
            self.sqlite_connection.commit()


    def create_sql_tables(self):
        print('Creating SQL tables for Metacyc')
        if os.path.exists(self.db_file):
            os.remove(self.db_file)
        self.sqlite_connection = sqlite3.connect(self.db_file)
        self.cursor = self.sqlite_connection.cursor()

        self.create_cpd_to_rxn_table()
        self.create_intermediate_reaction_ids_table()
        self.create_main_tables()

        self.commit_and_close_sqlite_cursor()

    ######### fetching data

    def generate_fetch_command(self,wanted_id,table_name,main_id_str):
        self.get_db_headers(table_name)
        db_headers=self.db_headers[table_name]
        headers_str=', '.join(db_headers)
        headers_str=f'{main_id_str}, {headers_str}'.upper()
        fetch_command = f'SELECT {headers_str} FROM {table_name} WHERE {main_id_str}="{wanted_id}"'
        return fetch_command

    def convert_sql_to_dict(self,sql_result,table_name):
        sql_result=sql_result[1:]
        res={}
        for i in range(len(self.db_headers[table_name])):
            db=self.db_headers[table_name][i]
            db_res=sql_result[i]
            if db_res:
                db_res=db_res.split(self.info_splitter)
                if db=='reaction_stoichiometry':
                    temp=[]
                    for i in range(len(db_res)):
                        if db_res[i] in ['<-','<->','->']:
                            temp.append(db_res[i])
                        else:
                            try:
                                stoi=int(db_res[i])
                                cpd_id=db_res[i+1]
                                temp.append([stoi,cpd_id])
                            except: pass
                else:
                    if db not in res: res[db]=set()
                    res[db].update(db_res)
        return res

    def fetch_intermediate_reaction_ids(self,wanted_id):
        if not os.path.exists(self.db_file):
            return {}
        table_name='TABLEINTRXNIDS'
        main_id_str='INTRXNID'
        fetch_command=self.generate_fetch_command(wanted_id,table_name,main_id_str)
        try:
            res_fetch = self.cursor.execute(fetch_command).fetchone()
            temp=self.convert_sql_to_dict(res_fetch,table_name)
            res={}
            if 'METACYCPRT' in temp: res['protein_ids']=temp['METACYCPRT']
            if 'METACYCRXN' in temp: res['reaction_ids']=temp['METACYCPRT']
            return res
        except:
            print(f'Failed retrieving {wanted_id} in {self.db_file}.{table_name}')
            return {}

    def fetch_cpd_to_rxn(self,wanted_id):
        if not os.path.exists(self.db_file):
            return {}
        table_name='TABLECPDRXN'
        main_id_str='METACYCCPD'
        fetch_command=self.generate_fetch_command(wanted_id,table_name,main_id_str)
        try:
            res_fetch = self.cursor.execute(fetch_command).fetchone()
            res=self.convert_sql_to_dict(res_fetch,table_name)
            return res['METACYCRXN']
        except:
            print(f'Failed retrieving {wanted_id} in {self.db_file}.{table_name}')
            return {}

    def fetch_from_main_table(self,wanted_id,table_name):
        if not os.path.exists(self.db_file):
            return {}
        table_name=table_name.upper()+'S'
        main_id_str='METACYC'
        fetch_command=self.generate_fetch_command(wanted_id,table_name,main_id_str)
        try:
            res_fetch = self.cursor.execute(fetch_command).fetchone()
            res=self.convert_sql_to_dict(res_fetch,table_name)
            return res
        except:
            print(f'Failed retrieving {wanted_id} in {self.db_file}.{table_name}')
            return {}

    ######### the parsing happens below

    def check_all_resources(self):
        required_resources=[
            'compounds.dat',
            'proteins.dat',
            'reactions.dat',
            'gene-links.dat',
            'genes.dat',
        ]
        for file in os.listdir(self.metacyc_folder):
            if file not in required_resources:
                os.remove(f'{self.metacyc_folder}{file}')
        if len(required_resources)==len(os.listdir(self.metacyc_folder)):
            return True
        return False


    def process_chemical_formula(self,chemical_formula):
        res=[]
        for atom in chemical_formula:
            atom=atom.replace('(','')
            atom=atom.replace(')','')
            atom=atom.replace(' ','')
            res.append(atom)
        return ''.join(res)

    def parse_compounds(self):
        rejects=set()
        input_file=f'{self.metacyc_folder}compounds.dat'
        line_type = None
        temp={}
        with open(input_file,encoding="ISO-8859-1") as file:
            for line in file:
                line=line.strip('\n')
                if line.startswith('UNIQUE-ID'):
                    if 'chemical_formula' in temp:
                        temp['chemical_formula']=self.process_chemical_formula(temp['chemical_formula'])
                    if temp:
                        yield temp
                    line=line.replace('UNIQUE-ID - ','')
                    line=line.strip()
                    line=strip_tags(line)
                    temp={'metacyc':{line}}
                    line_type=None
                elif    line.startswith('COMMON-NAME') or\
                        line.startswith('ABBREV-NAME') or\
                        line.startswith('SYSTEMATIC-NAME') or\
                        line.startswith('SYNONYMS') or\
                        line.startswith('CHEMICAL-FORMULA') or\
                        line.startswith('INCHI-KEY') or\
                        line.startswith('SMILES') or\
                        line.startswith('DBLINKS'):
                    line_type=line.split()[0]
                else:
                    line_type=None

                if line_type:
                    line=line.replace(line_type,'').strip().strip('-').strip()
                    line=strip_tags(line)
                    if line_type in ['COMMON-NAME','ABBREV-NAME','SYSTEMATIC-NAME','SYNONYMS',]:
                        db_type='synonyms'
                    elif line_type in ['CHEMICAL-FORMULA']:
                        db_type = 'chemical_formula'
                    elif line_type in ['INCHI-KEY']:
                        db_type = 'inchi_key'
                    elif line_type in ['SMILES']:
                        db_type = 'smiles'
                    else:
                        line=line.split()
                        db_type=line[0][1:].strip()
                        line=line[1].strip('"').strip()
                        if db_type=='LIGAND-CPD':
                            db_type='kegg'
                        elif db_type=='KEGG-GLYCAN':
                            db_type='kegg'
                        elif db_type=='SEED':
                            db_type='seed'
                        elif db_type=='BRENDA-COMPOUND':
                            db_type='brenda'
                        elif db_type=='ECOCYC':
                            db_type='metacyc'
                        elif db_type=='BIGG':
                            db_type='bigg'
                        elif db_type=='PUBCHEM':
                            db_type='pubchem_cid'
                        elif db_type=='PUBCHEM-SID':
                            db_type='pubchem_sid'
                        elif db_type=='HMDB':
                            db_type='hmdb'
                        elif db_type=='CHEMSPIDER':
                            db_type='chemspider'
                        elif db_type=='CHEBI':
                            db_type='chebi'
                        elif db_type=='REACTOME-CPD':
                            db_type='reactome'
                        elif db_type=='CHEBI':
                            db_type='chebi'
                        elif db_type=='METANETX':
                            db_type='metanetx'
                        elif db_type=='DRUGBANK':
                            db_type='drugbank'
                        elif db_type=='CAS':
                            db_type='cas'
                        elif db_type=='METABOLIGHTS':
                            db_type='metabolights'
                        elif db_type=='REFMET':
                            db_type='synonyms'
                        elif db_type=='KNAPSACK':
                            db_type='knapsack'
                        elif db_type=='LIPID_MAPS':
                            db_type='lipid_maps'
                        else:
                            db_type=None
                    if db_type:
                        if db_type not in temp:
                            if db_type=='chemical_formula':
                                temp[db_type]=[]
                            else:
                                temp[db_type]=set()
                        if db_type == 'chemical_formula':
                            temp[db_type].append(line)
                        else:
                            temp[db_type].add(line)

        if 'chemical_formula' in temp:
            temp['chemical_formula'] = self.process_chemical_formula(temp['chemical_formula'])
        yield temp

    def get_genes_uniprot(self):
        res={}
        genes_to_uniprot_path=f'{self.metacyc_folder}gene-links.dat'
        with open(genes_to_uniprot_path) as file:
            for line in file:
                if not line.startswith('#'):
                    line = line.strip('\n')
                    line=line.split('\t')
                    gene_id,uniprot_id=line[0],line[1]
                    if uniprot_id:
                        res[gene_id]=uniprot_id
        return res

    def parse_genes(self):
        genes_to_uniprot=self.get_genes_uniprot()
        input_file=f'{self.metacyc_folder}genes.dat'
        line_type = None
        temp={}
        with open(input_file,encoding="ISO-8859-1") as file:
            for line in file:
                line=line.strip('\n')
                if line.startswith('UNIQUE-ID'):
                    if temp:
                        yield temp
                    line=line.replace('UNIQUE-ID - ','')
                    line=line.strip()
                    line=strip_tags(line)
                    temp={'metacyc':{line}}
                    line_type=None
                elif    line.startswith('ACCESSION-1') or \
                        line.startswith('ACCESSION-2') or \
                        line.startswith('COMMON-NAME') or \
                        line.startswith('SYNONYMS') or \
                        line.startswith('PRODUCT') or \
                        line.startswith('DBLINKS'):
                    line_type = line.split()[0]
                else:
                    line_type=None
                if line_type:
                    line=line.replace(line_type,'').strip().strip('-').strip()
                    line=strip_tags(line)
                    if line_type in ['COMMON-NAME','SYNONYMS',]:
                        db_type='synonyms'
                    elif line_type in ['ACCESSION-1','ACCESSION-2']:
                        db_type = 'metacyc'
                    elif line_type in ['PRODUCT']:
                        db_type = 'protein_ids'
                    else:
                        line = line.split()
                        db_type = line[0][1:].strip()
                        line = line[1].strip('"').strip()
                        if db_type=='REGULONDB':
                            db_type='regulondb'
                        elif db_type=='STRING':
                            db_type='string'
                        elif db_type=='NCBI-GENE':
                            db_type='ncbi_genbank_gene'
                        elif db_type=='REFSEQ':
                            db_type='refseq'
                        elif db_type=='KEGG':
                            db_type='kegg'
                        elif db_type=='UNIPROT':
                            db_type='uniprot'
                        elif db_type=='INTERPRO':
                            db_type='interpro'
                        elif db_type=='GENECARDS':
                            db_type='genecard'
                        elif db_type=='ENSEMBL':
                            db_type='ensembl'
                        elif db_type=='ENTREZ':
                            db_type='entrez'
                        else:
                            db_type=None
                    if db_type:
                        if db_type not in temp: temp[db_type]=set()
                        temp[db_type].add(line)

        yield temp

    def parse_proteins(self):
        input_file=f'{self.metacyc_folder}proteins.dat'
        line_type = None
        temp={}


        with open(input_file,encoding="ISO-8859-1") as file:
            for line in file:
                line=line.strip('\n')
                if line.startswith('UNIQUE-ID'):
                    if temp:
                        yield temp
                    line=line.replace('UNIQUE-ID - ','')
                    line=line.strip()
                    line=strip_tags(line)
                    temp={'metacyc':{line}}
                    line_type=None
                elif    line.startswith('ACCESSION-1') or \
                        line.startswith('COMMON-NAME') or\
                        line.startswith('ABBREV-NAME') or\
                        line.startswith('COMPONENT-OF') or\
                        line.startswith('COMPONENTS') or\
                        line.startswith('GENE') or\
                        line.startswith('CATALYZES') or\
                        line.startswith('TYPES') or\
                        line.startswith('DBLINKS'):
                    line_type = line.split()[0]
                else:
                    line_type=None
                if line_type:
                    line=line.replace(line_type,'').strip().strip('-').strip()
                    line=strip_tags(line)
                    if line_type in ['COMMON-NAME','ABBREV-NAME',]:
                        db_type='synonyms'
                    elif line_type in ['ACCESSION-1']:
                        db_type = 'metacyc'
                    elif line_type in ['COMPONENT-OF']:
                        db_type = 'complex_ids'
                    elif line_type in ['COMPONENTS']:
                        db_type = 'subunit_ids'
                    elif line_type in ['GENE']:
                        db_type = 'gene_ids'
                    elif line_type in ['CATALYZES']:
                        db_type = 'intermediate_reaction_ids'
                    elif line_type in ['TYPES']:
                        db_type = 'protein_class'

                    else:
                        line = line.split()
                        db_type = line[0][1:].strip()
                        line = line[1].strip('"').strip()
                        if db_type=='METANETX':
                            db_type='metanetx'
                        elif db_type=='PANTHER':
                            db_type='panther'
                        elif db_type=='ECOCYC':
                            db_type='metacyc'
                        elif db_type=='BIOCYC':
                            db_type='metacyc'
                        elif db_type=='TRANSPORTER_CLASSIFICATION_DATABASE':
                            db_type='tcdb'
                        elif db_type=='INTERPRO':
                            db_type='interpro'
                        elif db_type=='CAZY':
                            db_type='cazy'
                        elif db_type=='REFSEQ':
                            db_type='refseq'
                        elif db_type=='SEED':
                            db_type='seed'
                        elif db_type=='STRING':
                            db_type='string'
                        elif db_type=='PROSITE':
                            db_type='prosite'
                        elif db_type=='PFAM':
                            db_type='pfam'
                        elif db_type=='UNIPROT':
                            db_type='uniprot'
                        else:
                            db_type=None
                    if db_type:
                        if db_type not in temp: temp[db_type]=set()
                        temp[db_type].add(line)

        yield temp

    def process_reaction_stoichiometry(self,reaction_dict):
        direction=None
        sign_index=None
        reaction_stoichiometry=None
        if 'direction' in reaction_dict:
            direction=reaction_dict.pop('direction')
        if 'sign_index' in reaction_dict:
            sign_index=reaction_dict.pop('sign_index')
        if 'reaction_stoichiometry' in reaction_dict :
            reaction_stoichiometry=reaction_dict.pop('reaction_stoichiometry')
        if direction and sign_index and reaction_stoichiometry:
            reaction_stoichiometry.insert(sign_index, direction)
        if reaction_stoichiometry:
            reaction_dict['reaction_stoichiometry']=reaction_stoichiometry

    def parse_reactions(self):
        input_file=f'{self.metacyc_folder}reactions.dat'
        line_type = None
        temp={}
        with open(input_file,encoding="ISO-8859-1") as file:
            for line in file:
                line=line.strip('\n')
                if line.startswith('UNIQUE-ID'):
                    if 'reaction_stoichiometry' in temp:
                        self.process_reaction_stoichiometry(temp)
                    if temp:
                        yield temp
                    line=line.replace('UNIQUE-ID - ','')
                    line=line.strip()
                    line=strip_tags(line)
                    temp={'metacyc':{line}}
                    line_type=None
                elif line.startswith('TYPES') or\
                        line.startswith('EC-NUMBER') or\
                        line.startswith('ENZYMATIC-REACTION') or\
                        line.startswith('REACTION-DIRECTION') or\
                        line.startswith('LEFT') or\
                        line.startswith('LEFT') or\
                        line.startswith('RIGHT') or\
                        line.startswith('^COEFFICIENT') or\
                        line.startswith('DBLINKS'):
                    line_type = line.split()[0]
                else:
                    line_type=None
                if line_type:
                    line=line.replace(line_type,'').strip().strip('-').strip()
                    line=strip_tags(line)
                    if line_type in ['TYPES']:
                        db_type = 'reaction_class'
                    elif line_type in ['EC-NUMBER']:
                        db_type = 'enzyme_ec'
                    elif line_type in ['ENZYMATIC-REACTION']:
                        db_type = 'intermediate_reaction_ids'
                    elif line_type in ['REACTION-DIRECTION']:
                        db_type = 'direction'
                    elif line_type in ['LEFT','RIGHT']:
                        db_type = 'reaction_stoichiometry'
                    elif line_type in ['^COEFFICIENT']:
                        db_type = 'stoichiometry'

                    else:
                        line = line.split()
                        db_type = line[0][1:].strip()
                        line = line[1].strip('"').strip()
                        if db_type=='METANETX-RXN':
                            db_type='metanetx'
                        elif db_type=='LIGAND-RXN':
                            db_type='kegg'
                        elif db_type=='RHEA':
                            db_type='rhea'

                        else:
                            db_type=None
                    if db_type:
                        if db_type not in temp:
                            if db_type =='reaction_stoichiometry':
                                temp[db_type]=[]
                            elif db_type == 'direction':
                                temp[db_type]=None
                            elif db_type == 'stoichiometry':
                                pass
                            else:
                                temp[db_type]=set()
                        if db_type == 'reaction_stoichiometry':
                            if line_type == 'LEFT':
                                if 'sign_index' not in temp:
                                    temp['sign_index'] = 0
                                temp['sign_index'] += 1
                            temp[db_type].append([1,line])
                        elif db_type == 'stoichiometry':
                            last_cpd=temp['reaction_stoichiometry'][-1]
                            last_cpd[0]=line
                        elif db_type == 'direction':
                            if 'RIGHT-TO-LEFT' in line:
                                line='<-'
                            elif 'LEFT-TO-RIGHT' in line:
                                line='->'
                            elif 'REVERSIBLE' in line:
                                line='<->'
                            temp[db_type] = line
                        else:
                            temp[db_type].add(line)

        if 'reaction_stoichiometry' in temp:
            self.process_reaction_stoichiometry(temp)
        yield temp


if __name__ == '__main__':
    '''
    to parse:
            'reactions.dat',
    '''
    s=Metacyc_SQLITE_Connector()
    #s.create_sql_tables()
    #s.start_sqlite_cursor()
    print(s.fetch_cpd_to_rxn('GDP-MANNOSE'))
    print(s.fetch_from_main_table('GDP-MANNOSE','compound'))
    print(s.fetch_from_main_table('EG10368','gene'))
    print(s.fetch_from_main_table('CPLX-2401','protein'))
    print(s.fetch_from_main_table('MONOMER-2782','protein'))
    print(s.fetch_from_main_table('RXN-12763','reaction'))
    print(s.fetch_from_main_table('GDP-MANNOSE','compound'))
    print(s.fetch_intermediate_reaction_ids('ENZRXN-23736'))