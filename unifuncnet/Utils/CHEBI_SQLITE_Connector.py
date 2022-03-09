import os
import sqlite3
from unifuncnet.Utils.util import RESOURCES_FOLDER,download_file_ftp,gunzip



class CHEBI_SQLITE_Connector():
    '''
    this just creates an sql database from two chebi files. first we download them and create the db
    then anytime we need to fetch info we just open the db, fetch the info, and close the connection
    '''
    def __init__(self):
        self.insert_step=10000
        self.chebi_db = f'{RESOURCES_FOLDER}chebi.db'
        self.download_chebi()

        if os.path.exists(self.chebi_db):
            self.chebi_start_sqlite_cursor()
        else:
            self.chebi_create_db()

    def chebi_start_sqlite_cursor(self):
        self.chebi_sqlite_connection = sqlite3.connect(self.chebi_db)
        self.chebi_cursor = self.chebi_sqlite_connection.cursor()

    def chebi_commit(self):
        self.chebi_sqlite_connection.commit()

    def chebi_execute(self, command):
        return self.chebi_cursor.execute(command)

    def chebi_executemany(self, command, chunk):
        return self.chebi_cursor.executemany(command, chunk)

    def chebi_commit_and_close_sqlite_cursor(self):
        self.chebi_commit()
        self.chebi_close_sql_connection()

    def chebi_close_sql_connection(self):
        self.chebi_sqlite_connection.close()

    def trim_chebi_obo(self,infile_path,outfile_path):
        a=set()
        with open(infile_path) as infile:
            with open(outfile_path,'a+') as outfile:
                line=infile.readline()
                while line:
                    line=line.strip('\n')
                    if line.startswith('id: CHEBI:'):
                        main_id=line.split('CHEBI:')[1].strip()
                    elif line.startswith('alt_id:'):
                        current_info=line.split('CHEBI:')[1].strip()
                        outline=f'{main_id}\tchebi\t{current_info}'
                        outfile.write(f'{outline}\n')
                    elif line.startswith('property_value:'):
                        line=line.split()
                        if len(line)==4:
                            link_type,current_info=line[1],line[2]
                            current_info=current_info.strip('\"')
                            if link_type.endswith('formula'): link_type='chemical_formula'
                            elif link_type.endswith('smiles'): link_type='smiles'
                            elif link_type.endswith('inchikey'): link_type='inchi_key'
                            else: link_type=None
                            if link_type:
                                outline=f'{main_id}\t{link_type}\t{current_info}'
                                outfile.write(f'{outline}\n')
                                a.add(link_type)
                    line=infile.readline()
        print(a)

    def trim_chebi_accession(self,infile_path,outfile_path):
        res=set()
        with open(infile_path) as infile:
            with open(outfile_path,'w+') as outfile:
                line=infile.readline()
                while line:
                    line=line.strip('\n')
                    #ID	COMPOUND_ID	SOURCE	TYPE	ACCESSION_NUMBER
                    _,chebi_id,_,id_type,secondary_id= line.split('\t')
                    outline=None
                    if id_type=='KEGG COMPOUND':
                        outline=f'{chebi_id}\tkegg\t{secondary_id}'
                    elif id_type=='KEGG DRUG':
                        outline=f'{chebi_id}\tkegg\t{secondary_id}'
                    elif id_type=='KEGG DRUG accession':
                        outline=f'{chebi_id}\tkegg\t{secondary_id}'
                    elif id_type=='KEGG COMPOUND accession':
                        outline=f'{chebi_id}\tkegg\t{secondary_id}'
                    elif id_type=='MetaCyc accession':
                        outline=f'{chebi_id}\tmetacyc\t{secondary_id}'
                    elif id_type=='HMDB accession':
                        outline=f'{chebi_id}\thmdb\t{secondary_id}'
                    elif id_type=='Chemspider accession':
                        outline=f'{chebi_id}\tchemspider\t{secondary_id}'
                    else:
                        res.add(id_type)
                    if outline:
                        outfile.write(f'{outline}\n')

                    line=infile.readline()
        os.remove(infile_path)

    def download_chebi_obo(self):
        url='https://ftp.ebi.ac.uk/pub/databases/chebi/ontology/chebi.obo.gz'
        infile_path=f'{RESOURCES_FOLDER}chebi.obo.gz'
        outfile_path=f'{RESOURCES_FOLDER}chebi2others.tsv'
        download_file_ftp(url,infile_path)
        gunzip(infile_path)
        self.trim_chebi_obo(infile_path.replace('.gz',''),outfile_path)
        os.remove(infile_path)
        os.remove(infile_path.replace('.gz',''))

    def download_chebi_to_others(self):
        url='https://ftp.ebi.ac.uk/pub/databases/chebi/Flat_file_tab_delimited/database_accession.tsv'
        infile_path=f'{RESOURCES_FOLDER}database_accession.tsv'
        outfile_path=f'{RESOURCES_FOLDER}chebi2others.tsv'
        download_file_ftp(url,infile_path)
        self.trim_chebi_accession(infile_path,outfile_path)

    def download_chebi(self):
        if not os.path.exists(self.chebi_db):
            self.download_chebi_to_others()
            self.download_chebi_obo()

    def get_chebi_to_others(self):
        input_path = f'{RESOURCES_FOLDER}chebi2others.tsv'
        with open(input_path) as file:
            line = file.readline()
            while line:
                line = line.strip('\n')
                if line:
                    current_chebi_id, db, db_id = line.split('\t')
                    yield current_chebi_id, db, db_id
                line = file.readline()
        os.remove(input_path)

    def chebi_create_db(self):
        #this will probably need to be changed to an output_folder provided by the user
        outfile_path=f'{RESOURCES_FOLDER}chebi2others.tsv'
        if os.path.exists(self.chebi_db):
            os.remove(self.chebi_db)
        self.chebi_start_sqlite_cursor()
        create_table_command = f'CREATE TABLE CHEBI2OTHERS (' \
                            f'CHEBI INTEGER,' \
                            f'DATABASE TEXT,' \
                            f'ALTID  TEXT )'
        self.chebi_execute(create_table_command)

        create_index_command = f'CREATE INDEX CHEBI_IDX ON CHEBI2OTHERS (CHEBI)'
        self.chebi_execute(create_index_command)

        self.chebi_commit()
        self.store_chebi2others()

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

    def store_chebi2others(self):
        chebi2others=self.get_chebi_to_others()
        generator_insert = self.generate_inserts(chebi2others)
        for table_chunk in generator_insert:
            insert_command = f'INSERT INTO CHEBI2OTHERS (CHEBI, DATABASE, ALTID) values (?,?,?)'
            self.chebi_executemany(insert_command, table_chunk)
        self.chebi_commit()

    def fetch_chebi_alt_id_info(self,chebi_id):
        res={}
        try:    chebi_id=int(chebi_id)
        except: return res
        fetch_command = f"SELECT CHEBI,DATABASE, ALTID FROM CHEBI2OTHERS WHERE ALTID = {chebi_id} AND DATABASE='chebi'"
        res_fetch = self.chebi_execute(fetch_command).fetchall()
        main_chebi_id=set()
        for i in res_fetch:
            chebi_id,db,alt_id=i
            main_chebi_id.add(chebi_id)
        if len(main_chebi_id)>1: raise Exception
        elif not main_chebi_id: return res
        else:
            return main_chebi_id.pop()

    def fetch_chebi_id_info(self,chebi_id):
        res={}
        try:    chebi_id=int(chebi_id)
        except: return res
        fetch_command = f"SELECT CHEBI,DATABASE, ALTID FROM CHEBI2OTHERS WHERE CHEBI = {chebi_id}"
        res_fetch=self.chebi_execute(fetch_command).fetchall()
        if res_fetch:
            main_chebi_id=str(chebi_id)
        for i in res_fetch:
            chebi_id,db,alt_id=i
            if db not in res: res[db]=set()
            res[db].add(alt_id)
        if not res:
            main_chebi_id=self.fetch_chebi_alt_id_info(chebi_id)
            main_chebi_id=str(main_chebi_id)
            return self.fetch_chebi_id_info(main_chebi_id)
        return main_chebi_id,res


if __name__ == '__main__':
    sql=CHEBI_SQLITE_Connector()
    alt_ids=sql.fetch_chebi_id_info('15377')
    print(alt_ids)