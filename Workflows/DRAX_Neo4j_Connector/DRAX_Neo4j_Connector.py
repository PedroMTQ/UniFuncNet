import re
import sys
import os
import subprocess
from pathlib import Path
import argparse
from neo4j import GraphDatabase


if sys.platform.startswith('win'):
    SPLITTER = '\\'
else:
    SPLITTER = '/'

DRAX_FOLDER = os.path.abspath(os.path.dirname(__file__)).split(SPLITTER)[0:-2]
DRAX_FOLDER = SPLITTER.join(DRAX_FOLDER) + SPLITTER
RESOURCES_FOLDER=f'{DRAX_FOLDER}Resources{SPLITTER}'


class DRAX_Neo4j_Connector():
    def __init__(self,username=None,password=None,try_limit=10,db_name=None,uri=None,):
        # Database Credentials
        port = '7687'
        if not uri: uri = 'bolt://localhost:'+port
        self.try_limit=try_limit
        self.uri=uri
        self.username=username
        self.password=password
        self.first_connection=True
        self.bolt_connections=0
        if not db_name: self.db_name='neo4j'
        else: self.db_name=db_name
        if not hasattr(self,'neo4j_driver'): self.neo4j_driver=None
        self.connect_to_neo4j()
        self.insert_step=2500

    def create_node_labels(self,drax_output_folder):
        drax_output_files=os.listdir(drax_output_folder)
        #this is for the drax neo4j db, a general db wont have these
        exceptional_labels=['_connected','reaction_compounds','has_subunits','in_complex','stoichiometry']
        self.node_labels={}
        self.main_labels=set()
        all_indexes=set()
        for drax_output_file in drax_output_files:
            node_type=drax_output_file.split('.')[0]
            #node_type for Proteins.tsv becomes Protein
            node_type=f'{node_type[0].upper()}{node_type[1:]}'
            node_type=node_type.strip('s')
            drax_output_file_path = f'{drax_output_folder}{drax_output_file}'
            generator_drax_output = self.yield_tsv(drax_output_file_path)
            for internal_id,node_info in generator_drax_output:
                for db_type in node_info:
                    added=False
                    for exc_label in exceptional_labels:
                        if db_type.endswith(exc_label):
                            if exc_label=='_connected':
                                entity_type=db_type.replace('s_connected','')
                                entity_type = f'{entity_type[0].upper()}{entity_type[1:]}'
                                self.node_labels[db_type] = entity_type
                            added=True
                            break
                    if not added:
                        self.node_labels[db_type]=db_type
                self.node_labels[node_type]=node_type
                self.main_labels.add(node_type)
        for i in self.node_labels:
            all_indexes.add(self.node_labels[i])
        #self.add_constraints(all_indexes)
        self.reset_db()
        self.add_indexes(all_indexes)

    def reset_db(self):
        reset_db_command = 'MATCH (n) DETACH DELETE n'
        self.run_command_neo4j(reset_db_command)
        #self.drop_constraints()
        self.drop_indexes()
        print('DATABASE RESET')




    def create_indexes(self,all_indexes):
        res=[]
        for i in all_indexes:
            if i not in self.main_labels:
                command=f' INDEX ON :{i}(node_info)'
                res.append(command)
        for i in self.main_labels:
            command=f' INDEX ON :{i}(drax_id)'
            res.append(command)
        return res

    def add_indexes(self,all_indexes):
        active_indexes = [i['labelsOrTypes'] for i in self.run_command_neo4j('CALL db.indexes() YIELD labelsOrTypes RETURN labelsOrTypes')]
        list_commands=[]
        for c in self.create_indexes(all_indexes):
            processed_c=c.split(':')[1].split('(')[0]
            if [processed_c] not in active_indexes:
                list_commands.append(f'CREATE {c}')
        if list_commands:

            for command in list_commands:
                self.run_command_neo4j(command)

    def drop_indexes(self):
        for i in self.run_command_neo4j("CALL db.indexes",return_output=True):
            try:
                label=i['labelsOrTypes'][0]
                properties=i['properties'][0]
                #print(f'Dropping index :{label}({properties})')
                self.run_command_neo4j(f"DROP INDEX ON :{label}({properties})")
            except:
                pass

    def get_available_indexes(self):
        res=set()
        for i in self.run_command_neo4j("CALL db.indexes",return_output=True):
            try:
                label=i['labelsOrTypes'][0]
                res.add(label)
            except:
                pass
        return res



    def create_constraints(self, all_indexes):
        res = []
        for i in all_indexes:
            if i not in self.main_labels:
                command = f' CONSTRAINT {i} FOR (c:{i}) REQUIRE c.node_info IS UNIQUE '
                res.append(command)
        for i in self.main_labels:
            command = f' CONSTRAINT {i} FOR (c:{i}) REQUIRE c.drax_id IS UNIQUE '
            res.append(command)

        return res

    def add_constraints(self,all_indexes):
        active_constraints = [i['description'] for i in self.run_command_neo4j('CALL db.constraints() YIELD description RETURN description')]
        active_constraints=[i.lower().split('assert')[1].split('is unique')[0] for i in active_constraints]
        active_constraints=[i.strip().strip('(').strip(')').strip() for i in active_constraints]
        list_commands=[]
        for c in self.create_constraints(all_indexes):
            processed_c=c.lower().split('require')[1].split('is unique')[0]
            processed_c=processed_c.strip().strip('(').strip(')').strip()
            processed_c=processed_c.split('.')[1]
            processed_c+='.'+processed_c
            if processed_c not in active_constraints:
                list_commands.append(f'CREATE {c}')
        if list_commands:
            for command in list_commands:
                self.run_command_neo4j(command)


    def drop_constraints(self):
        for constraint in self.run_command_neo4j("CALL db.constraints",return_output=True):
            print('Dropping constraint',constraint['name'])
            self.run_command_neo4j("DROP " + constraint['description'])





    def connect_to_databases(self,connect_ES=True):
        if connect_ES: self.connect_to_ElasticSearch()
        self.connect_to_neo4j()

    def connect_to_ElasticSearch(self):
        self.ElasticSearch_driver=ElasticSearch_Connector()
        self.ElasticSearch_driver.connect_to_ElasticSearch()

    def connect_to_neo4j(self):
        #classes are inheritable and work by themselves as well. so to avoid reconnecting all the time, we limit the initial connection
        if self.first_connection:
            try:
                print(f'Connecting to Neo4j database {repr(self.db_name)}')
                self.bolt_connections = 0
                # Connect to the neo4j database server
                self.neo4j_driver = GraphDatabase.driver(self.uri, auth=(self.username, self.password),encrypted=False)
                self.first_connection=False
            except:
                print('Couldn\'t connect Neo4j! Crashing here')
        else:
            self.bolt_connections = 0
            self.neo4j_driver = GraphDatabase.driver(self.uri,auth=(self.username, self.password),encrypted=False)

    def close_connection_neo4j(self):
        self.neo4j_driver.close()

    def reconnect_neo4j(self):
        c=0
        while c<=self.try_limit:
            try:
                self.close_connection_neo4j()
                self.connect_to_neo4j()
                return
            except:
                c+=1
        print('Couldn"t reconnect to neo4j! Check this urgently')



    def is_broken_connection(self,connection_error):
        if 'Failed to read from defunct connection Address' in connection_error or\
                'NoneType" object has no attribute "acquire"' in connection_error or\
            '[SSL] internal error' in connection_error:
            return True
        return False



    def run_command_neo4j2(self,command_to_run,return_output=False,batch=None):
        self.bolt_connections += 1
        c = 0
        # the bolt driver only allows around 100 connections before it needs to be restarted. so we preemptively restart the connection so it doesnt lag
        if self.bolt_connections > 98: self.reconnect_neo4j()
        while c <= self.try_limit:
            try:
                with self.neo4j_driver.session(database=self.db_name) as session:
                    res=session.run(command_to_run,batch=batch)
                    if 'return' in command_to_run.lower() or return_output:
                        return res.data()
                    else:
                        res.consume()
                        return
                    session.commit()
            except Exception as e:
                if self.is_broken_connection(str(e)):
                    self.reconnect_neo4j()
                else:
                    c+=1
                    print('Error,  ' + str(e) + '   ,while executing neo4j command\n', command_to_run)
        print('neo4j didnt run command so it will crash here!\n',command_to_run)
        raise Exception


    def run_command_neo4j(self,command_to_run,return_output=False,batch=None):

                with self.neo4j_driver.session(database=self.db_name) as session:
                    res=session.run(command_to_run,batch=batch)
                    if 'return' in command_to_run.lower() or return_output:
                        return res.data()
                    else:
                        res.consume()
                        return
                    session.commit()



    def yield_tsv(self, drax_tsv):
        try:
            if os.path.exists(drax_tsv):
                with open(drax_tsv) as file:
                    for line in file:
                        line = line.strip('\n').split('\t')
                        internal_id = line.pop(0).split(':')[1]
                        temp={}
                        for i in line:
                            id_type = i.split(':')[0]
                            annot = i.replace(id_type + ':', '')
                            if id_type not in temp: temp[id_type] = set()
                            temp[id_type].add(annot)

                        yield internal_id,temp
        except:
            return {}

    def replace_detail(self,current_detail):
        res=str(current_detail)
        res = res.replace('\\', '\\\\')
        res = res.replace('"', '\\"')
        res = res.replace('\'', '\\\'')
        return res


    def process_node_info(self,node_info):
        res={}
        c=0
        for detail in node_info:
            detail_neo4j=self.node_labels.get(detail)
            if detail_neo4j:
                if detail_neo4j not in self.main_labels:
                    for current_detail in node_info[detail]:
                        temp_command=[]
                        if isinstance(current_detail,str):  current_detail_str=self.replace_detail(current_detail)
                        else: current_detail_str=current_detail
                        if detail_neo4j not in res: res[detail_neo4j]=[]
                        temp = {'db_id': current_detail_str}
                        res[detail_neo4j].append(current_detail_str)
        return res

    def fill_out_nodes(self,chunk):
        node_types=set()
        for node in chunk:
            for sub_node in node['node_data']:
                node_types.add(sub_node)
        for nt in node_types:
            for node in chunk:
                if nt not in node['node_data']:
                    node['node_data'][nt]=[]
        return node_types


    def generate_command_subnodes(self,node_types):
        res=[]
        for nt in node_types:
            others_command = f'FOREACH (db_id in node_data.{nt} | MERGE (sn:{nt} {{node_info: db_id}}) CREATE (n)-[:HAS_{nt}]-> (sn)  ) '
            res.append(others_command)
        return ' '.join(res)

    def generate_chunks_drax_output(self, drax_wielder,node_type):
        step=self.insert_step
        temp=[]
        for node_id,node_info in drax_wielder:
            sub_nodes=self.process_node_info(node_info)
            node_data = {'drax_id': node_id, 'node_type': node_type,'node_data':sub_nodes}
            if node_type=='Reaction':
                sign=self.get_sign(node_info['reaction_compounds'])
                if sign=='<=>': reversible=True
                else: reversible=False
                node_data['reversible']=reversible
            if len(temp)<step:
                temp.append(node_data)
            else:
                yield temp
                temp = []
                temp.append(node_data)
        yield temp

    def create_nodes_drax(self,drax_output_folder):
        drax_output_files=os.listdir(drax_output_folder)
        for drax_output_file in drax_output_files:
            node_type=drax_output_file.split('.')[0]
            node_type=f'{node_type[0].upper()}{node_type[1:]}'
            node_type=node_type.strip('s')
            drax_output_file_path=f'{drax_output_folder}{drax_output_file}'
            generator_drax_output= self.yield_tsv(drax_output_file_path)
            generator_insert = self.generate_chunks_drax_output(generator_drax_output,node_type)
            for chunk in generator_insert:
                if chunk:
                    if node_type == 'Reaction':
                        command_to_run=f'WITH $batch as chunk UNWIND chunk as chunk_data ' \
                                       f'CREATE (n:{node_type}  {{node_type: chunk_data.node_type ,drax_id: chunk_data.drax_id, reversible: chunk_data.reversible }}) ' \
                                       f'WITH n, chunk_data.node_data as node_data '
                    else:
                        command_to_run=f'WITH $batch as chunk UNWIND chunk as chunk_data ' \
                                       f'CREATE (n:{node_type}  {{node_type: chunk_data.node_type ,drax_id: chunk_data.drax_id}} ) ' \
                                       f'WITH n, chunk_data.node_data as node_data '
                    node_types=self.fill_out_nodes(chunk)
                    subnode_commands=self.generate_command_subnodes(node_types)
                    command_to_run=f'{command_to_run}  {subnode_commands}'
                    self.run_command_neo4j(command_to_run=command_to_run,batch=chunk)


    def yield_list(self,list_to_yield):
        step=self.insert_step
        temp=[]
        for item in list_to_yield:
            if len(temp)<step:
                temp.append(item)
            else:
                yield temp
                temp = []
                temp.append(item)
        yield temp

    def get_sign(self,reaction_str):
        if '<=>' in reaction_str: return '<=>'
        elif '=>' in reaction_str: return '=>'
        elif '<=' in reaction_str: return '<='

    def get_reaction_compounds_neo4j(self,drax_id,reaction_str,stoichiometry):
        res=[]
        stoichiometry = stoichiometry.pop()
        stoichiometry = stoichiometry.split(',')
        reaction_sign = self.get_sign(reaction_str)
        left_side,right_side = reaction_str.split(reaction_sign)
        left_side=left_side.strip()
        left_side=left_side.split('+')
        left_side=[i.strip() for i in left_side]
        right_side=right_side.strip()
        right_side=right_side.split('+')
        right_side=[i.strip() for i in right_side]
        res={'IS_SUBSTRATE':[],'IS_PRODUCT':[]}
        for cpd_id in left_side:
            current_stoichiometry=stoichiometry.pop(0)
            item={'start_drax_id':drax_id,'end_drax_id':cpd_id,'stoichiometry':current_stoichiometry}
            res['IS_SUBSTRATE'].append(item)
        for cpd_id in right_side:
            current_stoichiometry=stoichiometry.pop(0)
            item={'start_drax_id':drax_id,'end_drax_id':cpd_id,'stoichiometry':current_stoichiometry}
            res['IS_PRODUCT'].append(item)
        return res

    def get_reaction_compounds_tsv(self,nodes_dict):
        for node in nodes_dict:
            reaction_compounds=[]
            stoichiometry=[]
            node_info=nodes_dict[node]
            substrates= node_info.pop('IS_SUBSTRATE')
            products= node_info.pop('IS_PRODUCT')
            reaction_str= node_info.pop('reaction_str').pop()
            reaction_sign = self.get_sign(reaction_str)
            for stoi,cpd in substrates:
                stoichiometry.append(stoi)
                reaction_compounds.append(cpd)
                reaction_compounds.append('+')
            reaction_compounds.pop(-1)
            reaction_compounds.append(reaction_sign)
            for stoi,cpd in products:
                stoichiometry.append(stoi)
                reaction_compounds.append(cpd)
                reaction_compounds.append('+')
            reaction_compounds.pop(-1)
            reaction_compounds=' '.join(reaction_compounds)
            stoichiometry=','.join(stoichiometry)
            nodes_dict[node]['reaction_compounds']=reaction_compounds
            nodes_dict[node]['stoichiometry']=stoichiometry




    def get_connections(self, generator_drax_output):
        res={}
        exceptional_connections=['reaction_compounds','has_subunits','in_complex','stoichiometry']
        for i in self.main_labels:
            exceptional_connections.append(f'{i.lower()}s_connected')
        for drax_id,node_info in generator_drax_output:
            for connection_type in node_info:
                if connection_type in exceptional_connections:
                    for connection_id in node_info[connection_type]:
                        connection_label = None
                        if connection_type =='reaction_compounds':
                            if 'stoichiometry' in node_info:
                                stoichiometry=node_info['stoichiometry']
                                reaction_compounds=self.get_reaction_compounds_neo4j(drax_id,connection_id,stoichiometry)
                                for connection_label in reaction_compounds:
                                    if connection_label not in res: res[connection_label] = []
                                    for item in reaction_compounds[connection_label]:
                                        res[connection_label].append(item)
                                connection_label=None
                        elif connection_type=='stoichiometry':
                            pass
                        elif connection_type=='has_subunits':
                            item = {'start_drax_id': drax_id,'end_drax_id': connection_id}
                            connection_label='HAS_SUBUNIT'
                        elif connection_type=='in_complex':
                            item = {'start_drax_id': drax_id,'end_drax_id': connection_id}
                            connection_label='IN_COMPLEX'
                        else:
                            item = {'start_drax_id': drax_id,'end_drax_id': connection_id}
                            connection_label='CONNECTED_TO'
                        #change this check
                        if connection_label:
                            if connection_label not in res: res[connection_label]=[]
                            res[connection_label].append(item)
        return res


    def connect_nodes_drax(self,drax_output_folder):
        drax_output_files = os.listdir(drax_output_folder)
        for drax_output_file in drax_output_files:
            node_type = drax_output_file.split('.')[0]
            node_type = f'{node_type[0].upper()}{node_type[1:]}'
            node_type = node_type.strip('s')
            drax_output_file_path = f'{drax_output_folder}{drax_output_file}'
            generator_drax_output= self.yield_tsv(drax_output_file_path)
            all_connections = self.get_connections(generator_drax_output)
            for connection_label in all_connections:
                for chunk in self.yield_list(all_connections[connection_label]):
                    if chunk:
                        if connection_label in ['IS_SUBSTRATE','IS_PRODUCT']:
                            command_to_run = f'WITH $batch as chunk UNWIND chunk as chunk_data ' \
                                             f'MATCH (n1:{node_type}  {{drax_id: chunk_data.start_drax_id}}) MATCH (n2  {{drax_id: chunk_data.end_drax_id}}) ' \
                                             f'CREATE (n1)-[r:{connection_label}]->(n2) SET r.stoichiometry=chunk_data.stoichiometry'
                        else:
                            command_to_run = f'WITH $batch as chunk UNWIND chunk as chunk_data ' \
                                             f'MATCH (n1:{node_type}  {{drax_id: chunk_data.start_drax_id}}) MATCH (n2  {{drax_id: chunk_data.end_drax_id}})' \
                                             f'CREATE (n1)-[:{connection_label}]->(n2)'
                        self.run_command_neo4j(command_to_run=command_to_run, batch=chunk)

    def remove_nas(self):
        command_to_run='MATCH (n) WHERE n.node_info="NA" DETACH DELETE n'
        self.run_command_neo4j(command_to_run=command_to_run)

    def export_drax_to_neo4j(self,drax_output_folder):
        self.start_time=time.time()
        if not drax_output_folder.endswith(SPLITTER):
            drax_output_folder+=SPLITTER
        self.create_node_labels(drax_output_folder)
        self.create_nodes_drax(drax_output_folder)
        self.connect_nodes_drax(drax_output_folder)

    def parse_consensus_tsv(self,input_file, wanted_ids):
        res = {i:set() for i in wanted_ids}
        with open(input_file) as file:
            file.readline()
            for line in file:
                line = line.strip('\n')
                line = line.split('\t')
                separator = line.index('|')
                annotations = line[separator + 1:]
                for db_annot in annotations:
                    db = db_annot.split(':')[0]
                    # to avoid bad splitting when dealing with descriptions
                    annot = db_annot[len(db) + 1:]
                    if db in res:
                        res[db].add(annot)
        res={i:res[i] for i in res if res[i]}
        unwind_res={i:[] for i in res}
        for db in res:
            for db_id in res[db]:
                item={'node_info':db_id}
                unwind_res[db].append(item)
        return unwind_res

    def parse_fetch_annotations_results(self,db,fetch_results,res):
        for i in fetch_results:
            #if db not in res: res[db]={}
            db_id=i['db_id']
            drax_id=i['drax_id']
            #if db_id not in res[db]: res[db][db_id]=set()
            #res[db][db_id].add(drax_id)
            res.add(drax_id)

    def parse_fetch_node_results(self,fetch_results,res):
        for i in fetch_results:
            relationship=i['r'][1]
            stoichiometry=i['r.stoichiometry']
            n1_node_type=i['n1']['node_type']
            n1_drax_id=i['n1']['drax_id']
            n2_node_type = i['n2'].get('node_type')
            n2_drax_id = i['n2'].get('drax_id')
            if not n2_drax_id:
                n2_drax_id = i['n2'].get('node_info')

            if relationship =='CONNECTED_TO':
                relationship=f'{n2_node_type.lower()}s_connected'
            elif relationship=='HAS_SUBUNIT':
                relationship='has_subunits'
            elif relationship=='IN_COMPLEX':
                relationship='in_complex'
            elif relationship.startswith('HAS'):
                relationship=relationship.replace('HAS_','')
            if n1_drax_id not in res:
                res[n1_drax_id]={}

            if relationship not in res[n1_drax_id]:
                if relationship in ['IS_SUBSTRATE','IS_PRODUCT']:
                    if n1_node_type!='Compound':
                        res[n1_drax_id][relationship]=[]
                else:
                    res[n1_drax_id][relationship]=set()
            if relationship in ['IS_SUBSTRATE', 'IS_PRODUCT']:
                if n1_node_type != 'Compound':
                    res[n1_drax_id][relationship].append([stoichiometry,n2_drax_id])
            else:
                res[n1_drax_id][relationship].add(n2_drax_id)


    def get_nodes_info(self,node_type,nodes_list):
        unwind_items=[]
        res={}
        for drax_id in nodes_list:
            unwind_items.append({'drax_id':drax_id})
        for chunk in self.yield_list(unwind_items):
            command_to_run = f'WITH $batch as chunk UNWIND chunk as chunk_data ' \
                             f'MATCH (n1:{node_type} {{drax_id: chunk_data.drax_id}})-[r]-(n2) RETURN n1,n2,r,r.stoichiometry'
            fetch_results = self.run_command_neo4j(command_to_run=command_to_run, batch=chunk)
            self.parse_fetch_node_results(fetch_results, res)
        if node_type=='Reaction':
            self.get_reaction_compounds_tsv(res)
        return res

    def export_neo4j_to_tsv(self,output_tsv_folder,genes_info,proteins_info,reactions_info,compounds_info):
        if not output_tsv_folder.endswith(SPLITTER):
            output_tsv_folder=f'{output_tsv_folder}{SPLITTER}'
        if not os.path.exists(output_tsv_folder):
            os.mkdir(output_tsv_folder)

        output_tsvs=[['g',genes_info],['p',proteins_info],['r',reactions_info],['c',compounds_info]]
        top_line='SOURCE\tINTERACTION\tTARGET\n'

        with open(f'{output_tsv_folder}Graph.sif', 'w+') as file:

            for source_type,info_dict in output_tsvs:
                for drax_id in info_dict:
                    line_info=info_dict[drax_id]
                    for db in line_info:
                        if db=='reaction_compounds':
                            for target_id in line_info[db].split():
                                try:
                                    int(target_id)
                                    target_type='c'
                                    line=f'{drax_id}\t{source_type}{target_type}\t{target_id}\n'
                                    file.write(line)
                                except:
                                    pass
                        elif db.endswith('_connected'):
                            target_type=db.lower()[0]
                            for target_id in line_info[db]:
                                line = f'{drax_id}\t{source_type}{target_type}\t{target_id}\n'
                                file.write(line)

                        elif db in ['in_complex','has_subunits']:
                            if db=='in_complex': target_type='cplx'
                            elif db=='has_subunits': target_type='subu'
                            for target_id in line_info[db]:
                                line = f'{drax_id}\t{source_type}{target_type}\t{target_id}\n'
                                file.write(line)

    def export_neo4j_to_sif(self,output_tsv_folder,genes_info,proteins_info,reactions_info,compounds_info):
        if not output_tsv_folder.endswith(SPLITTER):
            output_tsv_folder=f'{output_tsv_folder}{SPLITTER}'
        if not os.path.exists(output_tsv_folder):
            os.mkdir(output_tsv_folder)

        output_tsvs=[['Genes',genes_info],['Proteins',proteins_info],['Reactions',reactions_info],['Compounds',compounds_info]]
        for tsv_name,info_dict in output_tsvs:
            with open(f'{output_tsv_folder}{tsv_name}.tsv','w+') as file:
                for drax_id in info_dict:
                    line=[f'internal_id:{drax_id}']
                    line_info=info_dict[drax_id]
                    for db in line_info:
                        if db in ['reaction_compounds','stoichiometry']:
                            line.append(f'{db}:{line_info[db]}')
                        else:
                            for db_id in line_info[db]:
                                line.append(f'{db}:{db_id}')
                    line.append('\n')
                    line='\t'.join(line)
                    file.write(line)


    def get_all_nodes_info(self,protein_ids):
        proteins_info = self.get_nodes_info('Protein', protein_ids)
        reaction_ids = set()
        compound_ids = set()
        gene_ids = set()
        for node in proteins_info:
            if 'reactions_connected' in proteins_info[node]:
                reaction_ids.update(proteins_info[node].get('reactions_connected'))
        for node in proteins_info:
            if 'genes_connected' in proteins_info[node]:
                gene_ids.update(proteins_info[node]['genes_connected'])
        reactions_info = self.get_nodes_info('Reaction', reaction_ids)
        for node in reactions_info:
            if 'reaction_compounds' in reactions_info[node]:
                reaction_compounds=reactions_info[node].get('reaction_compounds')
                reaction_compounds=reaction_compounds.split()
                for i in reaction_compounds:
                    try:
                        int(i)
                        compound_ids.add(i)
                    except: pass
        compounds_info=self.get_nodes_info('Compound', compound_ids)
        genes_info=self.get_nodes_info('Gene', gene_ids)
        return genes_info,proteins_info,reactions_info,compounds_info

    def get_mantis_network(self,mantis_annotations,output_tsv_folder):
        self.start_time=time.time()
        available_indexes=self.get_available_indexes()
        mantis_annotations=self.parse_consensus_tsv(mantis_annotations,available_indexes)
        #res={}
        protein_ids=set()
        for db in mantis_annotations:
            for chunk in self.yield_list(mantis_annotations[db]):
                command_to_run = f'WITH $batch as chunk UNWIND chunk as chunk_data ' \
                                 f'MATCH (n:{db} {{node_info: chunk_data.node_info}})<--(n2) RETURN chunk_data.node_info as db_id,n2.drax_id as drax_id'
                fetch_results=self.run_command_neo4j(command_to_run=command_to_run, batch=chunk)
                self.parse_fetch_annotations_results(db,fetch_results,protein_ids)
        genes_info,proteins_info,reactions_info,compounds_info=self.get_all_nodes_info(protein_ids)
        self.export_neo4j_to_tsv(output_tsv_folder,genes_info,proteins_info,reactions_info,compounds_info)
        self.export_neo4j_to_sif(output_tsv_folder,genes_info,proteins_info,reactions_info,compounds_info)



    ####ELASTIC SEARCH####
    def save_detail_ElasticSearch(self,detail_type,detail):
        return self.ElasticSearch_driver.save_detail(detail_type,detail)

    def get_match_detail_ElasticSearch(self, detail_type, detail):
        return self.ElasticSearch_driver.get_match_detail(detail_type,detail)

    def get_detail_ElasticSearch(self, detail_type, detail_ES_ID):
        return self.ElasticSearch_driver.get_detail(detail_type,detail_ES_ID)

    def get_all_details_by_type_ElasticSearch(self, detail_type):
        return self.ElasticSearch_driver.get_all_details_by_type(detail_type)

    def reset_db_ElasticSearch(self):
        return self.ElasticSearch_driver.reset_db()






if __name__ == '__main__':
    import time
    # Database Credentials
    uri = 'bolt://localhost:7687'
    username = 'neo4j'
    password = 'drax_neo4j'
    neo4j_driver = DRAX_Neo4j_Connector(username, password,db_name='neo4j',uri=uri)
    drax_output_folder='/home/pedroq/Desktop/test_drax/out/'
    mantis_input_tsv='/home/pedroq/Desktop/test_drax/test_mantis1.tsv'
    output_tsv_folder='/home/pedroq/Desktop/test_drax/network_mantis/'
    #neo4j_driver.reset_db()

    #neo4j_driver.export_drax_to_neo4j(drax_output_folder)
    neo4j_driver.get_mantis_network(mantis_input_tsv,output_tsv_folder)
    #nodes_info = neo4j_driver.get_nodes_info('Reaction', ['24'])
    #print(nodes_info)

    #print('time',time.time()-neo4j_driver.start_time)
