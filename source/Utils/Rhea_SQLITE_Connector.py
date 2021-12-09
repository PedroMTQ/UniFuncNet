



from source.Utils.util import get_stoichiometry,standardize_reaction_str,SPLITTER,download_file_ftp,gunzip,RESOURCES_FOLDER
import re
import os
import sqlite3

class Rhea_SQLITE_Connector():
    def __init__(self):
        self.insert_step=5000
        self.db_file = f'{RESOURCES_FOLDER}rhea.db'
        if os.path.exists(self.db_file):
            self.start_sqlite_cursor()
        else:
            self.download_all_rhea_resources()
            self.create_sql_table()

    def start_sqlite_cursor(self):
        self.sqlite_connection = sqlite3.connect(self.db_file)
        self.cursor = self.sqlite_connection.cursor()

    def commit_and_close_sqlite_cursor(self):
        self.sqlite_connection.commit()
        self.sqlite_connection.close()

    def close_sql_connection(self):
        self.sqlite_connection.close()

    def check_table(self):
        self.cursor.execute("SELECT * FROM RHEAREACTIONS limit 10")
        res_fetch = self.cursor.fetchall()
        print(res_fetch)


    def download_all_rhea_resources(self):
        for url in ['https://ftp.expasy.org/databases/rhea/tsv/rhea2uniprot.tsv',
                  'https://ftp.expasy.org/databases/rhea/tsv/rhea2xrefs.tsv',
                  'https://ftp.expasy.org/databases/rhea/tsv/rhea-directions.tsv',
                  'https://ftp.expasy.org/databases/rhea/txt/rhea-reactions.txt.gz',
                  ]:
            uncompressed_file=RESOURCES_FOLDER+os.path.basename(url).replace('.gz','')
            compressed_file=RESOURCES_FOLDER+os.path.basename(url)
            if not os.path.exists(uncompressed_file) and not os.path.exists(compressed_file):
                download_file_ftp(url,compressed_file)

                if url.endswith('.gz'):
                    gunzip(compressed_file,remove_source=True)

    def parse_rhea2xrefs(self,rhea2xrefs_path):
        print('Parsing rhea2xrefs')
        rhea2ids={}
        with open(rhea2xrefs_path) as file:
            line=file.readline()
            line=file.readline()
            while line:
                line=line.strip('\n')
                if line:
                    rhea_id,direction,master_id,db_id,db_type=line.split('\t')
                    if master_id not in rhea2ids: rhea2ids[master_id]={}
                    if db_type=='EC': db_type='enzyme_ec'
                    elif db_type=='METACYC': db_type='biocyc'
                    elif db_type=='ECOCYC': db_type='biocyc'
                    elif db_type=='KEGG_REACTION': db_type='kegg'
                    elif db_type=='REACTOME': db_type=None
                    elif db_type=='MACIE': db_type=None
                    elif db_type=='GO': db_type=None
                    if db_type:
                        if db_type not in rhea2ids[master_id]: rhea2ids[master_id][db_type]=set()
                        rhea2ids[master_id][db_type].add(db_id)
                line=file.readline()
        return rhea2ids


    def parse_rhea_reactions(self,rhea_reactions_path):
        print('Parsing rheareactions')
        with open(rhea_reactions_path) as file:
            line=file.readline()
            rhea_id, reaction_str, chebi_equation = None, None, None
            while line:
                line=line.strip('\n')
                if line.startswith('///'): pass
                elif line.startswith('ENTRY'):
                    line=line.replace('ENTRY','')
                    line=line.strip()
                    rhea_id=line.replace('RHEA:','')
                elif line.startswith('DEFINITION'):
                    line=line.replace('DEFINITION','')
                    line=line.strip()
                    reaction_str=line
                elif line.startswith('EQUATION'):
                    line=line.replace('EQUATION','')
                    line=line.strip()
                    chebi_equation=line
                    if ',' in chebi_equation:
                        chebi_equation=chebi_equation.replace(',',' + ')
                if rhea_id and reaction_str and chebi_equation:
                    chebi_equation = standardize_reaction_str(chebi_equation)
                    chebi_equation=chebi_equation.replace('CHEBI:','')
                    yield rhea_id,reaction_str,chebi_equation
                    rhea_id,reaction_str,chebi_equation=None,None,None
                line=file.readline()

    def parse_rhea2uniprot(self,rhea2uniprot_path):
        print('Parsing rhea2uniprot')
        res={}
        with open(rhea2uniprot_path) as file:
            line=file.readline()
            line=file.readline()
            while line:
                line=line.strip('\n')
                if line:
                    rhea_id,direction,master_id,uniprot_id=line.split('\t')
                    if master_id not in res: res[master_id]=set()
                    res[master_id].add(uniprot_id)
                line=file.readline()
        return res

    def parse_rhea_directions(self,rhea_directions_path):
        res={}
        with open(rhea_directions_path) as file:
            line=file.readline()
            line=file.readline()
            while line:
                line=line.strip('\n')
                if line:
                    master_id,left_id,right_id,bidirectional_id=line.split('\t')
                    res[master_id]=left_id,right_id,bidirectional_id
                line=file.readline()
        return res

    def generate_alt_ids_yielder(self,rhea_masters):
        for master_id in rhea_masters:
            left_id, right_id, bidirectional_id= rhea_masters[master_id]
            yield left_id,master_id
            yield right_id,master_id
            yield bidirectional_id,master_id


    def generate_components_yielder(self,rhea2ids,rhea2uniprot,rhea_reactions,rhea_masters):
        '''
        rhea2ids -> protein (enzyme_ec and uniprot) and reaction(biocyc, kegg)
        rhea2uniprot -> protein
        '''
        for rhea_id,reaction_str,chebi_equation in rhea_reactions:
            if rhea_id in rhea_masters:
                alt_ids=rhea_masters[rhea_id]
                alt_ids=','.join(alt_ids)
                biocyc_ids=set()
                kegg_ids=set()
                enzyme_ec_ids=set()
                uniprot_ids=set()
                if rhea_id in rhea2ids:
                    if 'enzyme_ec' in rhea2ids[rhea_id]: enzyme_ec_ids.update(rhea2ids[rhea_id]['enzyme_ec'])
                    if 'biocyc' in rhea2ids[rhea_id]: biocyc_ids.update(rhea2ids[rhea_id]['biocyc'])
                    if 'kegg' in rhea2ids[rhea_id]: kegg_ids.update(rhea2ids[rhea_id]['kegg'])
                if rhea_id in rhea2uniprot:
                    if 'uniprot' in rhea2uniprot[rhea_id]: uniprot_ids.update(rhea2uniprot[rhea_id]['uniprot'])
                biocyc_ids=','.join(biocyc_ids)
                kegg_ids=','.join(kegg_ids)
                enzyme_ec_ids=','.join(enzyme_ec_ids)
                uniprot_ids=','.join(uniprot_ids)
                yield rhea_id,alt_ids,biocyc_ids,kegg_ids,enzyme_ec_ids,uniprot_ids,reaction_str,chebi_equation


    def scrape_rhea(self):
        rhea_reactions_path=RESOURCES_FOLDER+'rhea-reactions.txt'
        rhea2uniprot_path=RESOURCES_FOLDER+'rhea2uniprot.tsv'
        rhea2xrefs_path=RESOURCES_FOLDER+'rhea2xrefs.tsv'
        rhea_directions_path=RESOURCES_FOLDER+'rhea-directions.tsv'


        rhea2ids = self.parse_rhea2xrefs(rhea2xrefs_path)
        rhea_reactions=self.parse_rhea_reactions(rhea_reactions_path)
        rhea2uniprot=self.parse_rhea2uniprot(rhea2uniprot_path)
        rhea_masters=self.parse_rhea_directions(rhea_directions_path)


        rhea_yielder= self.generate_components_yielder(rhea2ids,rhea2uniprot,rhea_reactions,rhea_masters)
        self.store_rheadb(rhea_yielder)

        rhea_yielder= self.generate_alt_ids_yielder(rhea_masters)
        self.store_alt_ids(rhea_yielder)
        for i in [rhea_reactions_path,rhea2uniprot_path,rhea2xrefs_path,rhea_directions_path]:
            os.remove(i)


    def create_sql_table(self):
        '''
        reaction table
        rhea id| ids| str | cpd chebi ids

        alt_ids table
        rhea alt id | master id
        '''
        if os.path.exists(self.db_file):
            os.remove(self.db_file)
        self.sqlite_connection = sqlite3.connect(self.db_file)
        self.cursor = self.sqlite_connection.cursor()
        create_reaction_table_command = f'CREATE TABLE RHEAREACTIONS (' \
                            f'RHEA INTEGER,' \
                            f'ALTIDS TEXT,' \
                            f'BIOCYC TEXT,' \
                            f'KEGG TEXT,' \
                            f'ENZYMEEC TEXT,' \
                            f'UNIPROT TEXT,' \
                            f'EQUATIONSTR TEXT,' \
                            f'EQUATIONCHEBI  TEXT )'
        self.cursor.execute(create_reaction_table_command)
        create_index_command = f'CREATE INDEX RHEA_IDX ON RHEAREACTIONS (RHEA)'
        self.cursor.execute(create_index_command)

        create_alt_table_command = f'CREATE TABLE RHEAALTIDS (' \
                            f'ALTID INTEGER,' \
                            f'MASTERID  INTEGER )'
        self.cursor.execute(create_alt_table_command)

        create_index_command = f'CREATE INDEX ALTID_IDX ON RHEAALTIDS (ALTID)'
        self.cursor.execute(create_index_command)


        self.sqlite_connection.commit()
        self.scrape_rhea()


    def generate_inserts(self, input_generator):
        step=self.insert_step
        temp=[]
        for i in input_generator:
            if len(temp)<step:
                temp.append(i)
            elif len(temp)==step:
                yield temp
                temp=[]
        yield temp

    def store_alt_ids(self,rhea_yielder):
        generator_insert = self.generate_inserts(rhea_yielder)
        for table_chunk in generator_insert:
            insert_command = f'INSERT INTO RHEAALTIDS (ALTID, MASTERID) values (?,?)'
            self.cursor.executemany(insert_command, table_chunk)
        self.sqlite_connection.commit()

    def store_rheadb(self,rhea_yielder):
        generator_insert = self.generate_inserts(rhea_yielder)
        for table_chunk in generator_insert:
            insert_command = f'INSERT INTO RHEAREACTIONS (RHEA,ALTIDS, BIOCYC, KEGG, ENZYMEEC, UNIPROT, EQUATIONSTR, EQUATIONCHEBI) values (?,?,?,?,?,?,?,?)'
            self.cursor.executemany(insert_command, table_chunk)
        self.sqlite_connection.commit()

    def fetch_rhea_id_info(self,rhea_id):
        res={}
        try:    rhea_id=int(rhea_id)
        except: return res
        fetch_master_id=f'SELECT ALTID, MASTERID FROM RHEAALTIDS WHERE ALTID="{rhea_id}"'
        res_fetch=self.cursor.execute(fetch_master_id).fetchone()
        if res_fetch:
            alt_id,master_id=res_fetch
        else:
            master_id=rhea_id
        fetch_command = f'SELECT RHEA, ALTIDS, BIOCYC, KEGG, ENZYMEEC, UNIPROT, EQUATIONSTR, EQUATIONCHEBI FROM RHEAREACTIONS WHERE RHEA = "{master_id}"'
        res_fetch=self.cursor.execute(fetch_command).fetchone()
        if not res_fetch: return res
        master_id,alt_ids,biocyc_ids,kegg_ids,enzyme_ec_ids,uniprot_ids,reaction_str,chebi_equation = res_fetch
        res['alt_ids']=[i for i in alt_ids.split(',') if i]
        res['biocyc']=[i for i in biocyc_ids.split(',') if i]
        res['kegg']=[i for i in kegg_ids.split(',') if i]
        res['enzyme_ec']=[i for i in enzyme_ec_ids.split(',') if i]
        res['uniprot']=[i for i in uniprot_ids.split(',') if i]
        res['reaction_str']=reaction_str
        res['chebi_equation']=chebi_equation

        return res

    def fetch_rhea_from_id(self,id_type,input_id):
        res = []
        id_type_sql=None
        if id_type=='enzyme_ec': id_type_sql='ENZYMEEC'
        elif id_type in ['biocyc','kegg','uniprot']: id_type_sql=id_type.upper()
        if not id_type_sql or not input_id: return res
        fetch_command = f"SELECT RHEA,ALTIDS FROM RHEAREACTIONS WHERE {id_type_sql} = '{input_id}'"
        res_fetch=self.cursor.execute(fetch_command).fetchall()
        if not res_fetch: return res
        for rhea_id,alt_ids in res_fetch:
            alt_ids = alt_ids.split(',')
            alt_ids.insert(0,rhea_id)
            res.append(alt_ids)
        return rhea_id

    def find_reactions_chebi(self,chebi_id):
        res=set()
        try:    chebi_id=int(chebi_id)
        except: return res
        fetch_command = f'SELECT RHEA, ALTIDS, BIOCYC, KEGG, ENZYMEEC, UNIPROT, EQUATIONSTR, EQUATIONCHEBI FROM RHEAREACTIONS WHERE EQUATIONCHEBI LIKE "%{chebi_id}%"'
        res_fetch=self.cursor.execute(fetch_command).fetchall()
        for master_id,alt_ids,biocyc_ids,kegg_ids,enzyme_ec_ids,uniprot_ids,reaction_str,chebi_equation in res_fetch:
            all_chebi=re.findall('\d+',chebi_equation)
            all_chebi=[int(i) for i in all_chebi]
            if chebi_id in all_chebi:
               res.add(master_id)
        return res



if __name__ == '__main__':
    s=Rhea_SQLITE_Connector()

    r=s.find_reactions_chebi('16459')
    print(f'found {len(r)} for this id')
    r=s.fetch_rhea_id_info('10000')
    print(r)
    #this is an alternative id for 10000
    #r=s.fetch_rhea_id_info('10001')
    #print(r)
    #r=s.fetch_rhea_id_info('11210000')
    #print(r)
    #r=s.fetch_rhea_from_id('enzyme_ec','3.5.1.50')
    #print(r)

