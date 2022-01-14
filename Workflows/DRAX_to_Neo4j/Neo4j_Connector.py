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


class Neo4j_Connector():
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
        self.reset_db()
        self.create_node_labels()
        self.start_time=time.time()
        self.drax_connections={}
        self.insert_step=10

    def create_node_labels(self):
        self.node_labels={
             None:'identifier',
            'synonyms':'Synonym',
            'pathway':'Pathway',
            'reaction_lass':'Class',
            'protein_class':'Class',
            'location':'Location',
            'chemical_formula':'Chemical_formula',
            'reaction_str':'Reaction_str',
            'reaction_compounds':'compound',
            'in_complex':'protein',
            'has_subunits':'protein',
            'genes_connected':'gene',
            'proteins_connected':'protein',
            'reactions_connected':'reaction',

        }

    def reset_db(self):
        reset_db_command = 'MATCH (n) DETACH DELETE n'
        self.run_command_neo4j(reset_db_command)
        print('DATABASE RESET')

    def create_constraints(self):
        return [ ' CONSTRAINT ON (Synonym:Synonym) ASSERT Synonym.Synonym IS UNIQUE',
                 ' CONSTRAINT ON (Pathway:Pathway) ASSERT Pathway.Pathway IS UNIQUE',
                 ' CONSTRAINT ON (Class:Class) ASSERT Class.Class IS UNIQUE',
                 ' CONSTRAINT ON (Location:Location) ASSERT Location.Location IS UNIQUE',
                 ' CONSTRAINT ON (Chemical_formula:Chemical_formula) ASSERT Chemical_formula.Chemical_formula IS UNIQUE',
                 ' CONSTRAINT ON (Reaction_str:Reaction_str) ASSERT Reaction_str.Reaction_str IS UNIQUE',
        ]

    def create_indexes(self):
        return [
                #main nodes
                ' INDEX ON :Protein(node_type)',
                ' INDEX ON :Reaction(node_type)',
                ' INDEX ON :Gene(node_type)',
                ' INDEX ON :Compound(node_type)',
                #sub nodes
                ' INDEX ON :Identifier(Identifier)',
                ' INDEX ON :Synonym(Synonym)',
                ' INDEX ON :Pathway(Pathway)',
                ' INDEX ON :Class(Class)',
                ' INDEX ON :Location(Location)',
                ' INDEX ON :Chemical_formula(Chemical_formula)',
                ' INDEX ON :Reaction_str(Reaction_str)',
                ]

    def drop_constraints(self):
        for constraint in self.run_command_neo4j("CALL db.constraints",return_output=True):
            print('Dropping constraint',constraint)
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
                self.add_constraints()
                self.add_indexes()
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

    def add_indexes(self):
        active_indexes = [i['labelsOrTypes'] for i in self.run_command_neo4j('CALL db.indexes() YIELD labelsOrTypes RETURN labelsOrTypes')]
        list_commands=[]
        for c in self.create_indexes():
            processed_c=c.split(':')[1].split('(')[0]
            if [processed_c] not in active_indexes:
                list_commands.append(f'CREATE {c}')
        if list_commands:

            list_commands=' '.join(list_commands)
            self.run_command_neo4j(list_commands)


    def add_constraints(self):
        active_constraints = [i['description'] for i in self.run_command_neo4j('CALL db.constraints() YIELD description RETURN description')]
        active_constraints=[i.lower().split('assert')[1].split('is unique')[0] for i in active_constraints]
        active_constraints=[i.strip().strip('(').strip(')').strip() for i in active_constraints]
        list_commands=[]
        for c in self.create_constraints():
            processed_c=c.lower().split('assert')[1].split('is unique')[0]
            processed_c=processed_c.strip().strip('(').strip(')').strip()
            processed_c=processed_c.split('.')[1]
            processed_c+='.'+processed_c
            if processed_c not in active_constraints:
                list_commands.append(f'CREATE {c}')
        if list_commands:
            list_commands = ' '.join(list_commands)
            self.run_command_neo4j(list_commands)

    def is_broken_connection(self,connection_error):
        if 'Failed to read from defunct connection Address' in connection_error or\
                'NoneType" object has no attribute "acquire"' in connection_error or\
            '[SSL] internal error' in connection_error:
            return True
        return False

    def run_command_neo4j_legacy(self,command_to_run,return_output=False):
        if command_to_run:
            if isinstance(command_to_run,list):
                to_run=''.join(command_to_run)
            else: to_run=command_to_run
            self.bolt_connections+=1
            c=0
            #the bolt driver only allows around 100 connections before it needs to be restarted. so we preemptively restart the connection so it doesnt lag
            if self.bolt_connections>98: self.reconnect_neo4j()
            while c<=self.try_limit:
                try:
                    #NEED TO MAKE A PERMANENT SESSION
                    with self.neo4j_driver.session(database=self.db_name) as session:
                        res=session.run(to_run)
                        if c>0: print('neo4j managed to run command!')
                        if 'return' in to_run.lower() or return_output:
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
                        print('Error,  ' + str(e) + '   ,while executing neo4j command\n', to_run)
            print('neo4j didnt run command so it will crash here!\n',to_run)
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


    def yield_drax_tsv(self, drax_tsv):
        if os.path.exists(drax_tsv):
            with open(drax_tsv) as file:
                line = file.readline()
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

    def process_node_info(self,node_info):
        res={}
        c=0
        for detail in node_info:
            detail_neo4j=self.node_labels.get(detail)
            if not detail_neo4j: detail_neo4j=self.node_labels[detail_neo4j]
            if detail_neo4j not in ['gene','protein','reaction','compound']:
                for current_detail in node_info[detail]:
                    temp_command=[]
                    if isinstance(current_detail,str):  current_detail_str=self.replace_detail(current_detail)
                    else: current_detail_str=current_detail
                    if detail_neo4j=='Identifier':
                        temp={'database':detail,'id':current_detail_str}
                        res[detail_neo4j].append(temp)

                    elif detail_neo4j=='Synonym':
                        temp = {detail_neo4j: current_detail_str}
                        if detail_neo4j not in res: res[detail_neo4j]=[]
                        res[detail_neo4j].append(temp)

                    elif detail_neo4j=='Pathway':
                        temp = {detail_neo4j: current_detail_str}
                        if detail_neo4j not in res: res[detail_neo4j]=[]
                        res[detail_neo4j].append(temp)

                    elif detail_neo4j=='Class':
                        temp = {detail_neo4j: current_detail_str}
                        if detail_neo4j not in res: res[detail_neo4j]=[]
                        res[detail_neo4j].append(temp)

                    elif detail_neo4j=='Location':
                        temp = {detail_neo4j: current_detail_str}
                        if detail_neo4j not in res: res[detail_neo4j]=[]
                        res[detail_neo4j].append(temp)

                    elif detail_neo4j=='Chemical_formula':
                        temp = {detail_neo4j: current_detail_str}
                        if detail_neo4j not in res: res[detail_neo4j]=[]
                        res[detail_neo4j].append(temp)

                    elif detail_neo4j=='Reaction_str':
                        temp = {detail_neo4j: current_detail_str}
                        if detail_neo4j not in res: res[detail_neo4j]=[]
                        res[detail_neo4j].append(temp)

        return res

    def fill_out_nodes(self,chunk):
        for node in chunk:
            for sub_node in node['node_data']:
                print(sub_node)
        raise Exception


    def generate_chunks_drax_output(self, drax_wielder,node_type):
        step=self.insert_step
        temp=[]
        for node_id,node_info in drax_wielder:
            sub_nodes=self.process_node_info(node_info)
            node_data = {'drax_id': node_id, 'node_type': node_type,'node_data':sub_nodes}
            if len(temp)<step:
                temp.append(node_data)
            else:
                yield temp
                temp = []
                temp.append(node_data)

        yield temp


    def export_drax_to_neo4j(self,drax_output_folder):
        if not drax_output_folder.endswith(SPLITTER):
            drax_output_folder+=SPLITTER
        drax_output_files=['gene']
        #drax_output_files=['gene','protein','reaction','compound']
        for node_type in drax_output_files:
            output_str=f'{node_type[0].upper()}{node_type[1:]}s'
            node_type=f'{node_type[0].upper()}{node_type[1:]}'
            drax_tsv_file=f'{drax_output_folder}{output_str}.tsv'
            generator_drax_output= self.yield_drax_tsv(drax_tsv_file)
            generator_insert = self.generate_chunks_drax_output(generator_drax_output,node_type)
            for chunk in generator_insert:
                available_data=[]
                print(chunk)
                self.fill_out_nodes(chunk)

                command_to_run=f'WITH $batch as chunk UNWIND chunk as chunk_data ' \
                               f'CREATE (n:{node_type}  {{node_type: chunk_data.node_type ,drax_id: chunk_data.drax_id}}) ' \
                               f'WITH n, chunk_data.node_data as node_data_chunk '

                identifiers_command= f'UNWIND node_data_chunk.Identifier as ids ' \
                               f'MERGE (n)-[:HAS_ID]->(sn:Identifier {{id: ids.id, database: ids.database}}) '
                others_command=f'WITH n, node_data_chunk ' \
                               f'UNWIND node_data_chunk.Synonym as synonyms ' \
                               f'MERGE (n)-[:HAS_SYNONYM]->(sn:Synonym {{node_type: synonyms.Synonym}}) ' \


                print(command_to_run+others_command)
                res=self.run_command_neo4j(command_to_run=command_to_run,batch=chunk)
                print(res)


                #self.generate_node_neo4j(node_type,node_id,node_info)
                #self.create_dict_neo4j(node_type,node_id,node_info)

    def create_dict_neo4j(self,node_type,node_id,node_info):
        res={}
        for detail in node_info:
            detail_neo4j=self.node_labels.get(detail)
            if not detail_neo4j: detail_neo4j=self.node_labels[detail_neo4j]
            if detail_neo4j not in ['gene','protein','reaction','compound']:
                for current_detail in node_info[detail]:
                    if isinstance(current_detail,str):  current_detail_str=self.replace_detail(current_detail)
                    else: current_detail_str=current_detail
                    if detail_neo4j=='Identifier':
                        temp_command.append(f' ')

        create_command=f'CREATE (n0:{node_type} {{node_type: "{node_type}", drax_id: "{node_id}"}}) WITH n0 '
        unwind_command = f'UNWIND $SUBNODES as subnodes MERGE (n{c}:{detail_neo4j} {{ Database:"{detail}", {detail_neo4j}: "{current_detail_str}" }}) '


    def generate_node_neo4j(self,node_type,node_id,node_info):
        list_commands=[]
        creation_main_node =f'CREATE (n0:{node_type} {{node_type: "{node_type}", drax_id: "{node_id}"}}) WITH n0 '
        list_commands.append(creation_main_node)
        c=1
        for detail in node_info:
            detail_neo4j=self.node_labels.get(detail)
            if not detail_neo4j: detail_neo4j=self.node_labels[detail_neo4j]
            if detail_neo4j not in ['gene','protein','reaction','compound']:
                for current_detail in node_info[detail]:
                    temp_command=[]
                    if isinstance(current_detail,str):  current_detail_str=self.replace_detail(current_detail)
                    else: current_detail_str=current_detail
                    if detail_neo4j=='Identifier':
                        temp_command.append(f' MERGE (n{c}:{detail_neo4j} {{ Database:"{detail}", {detail_neo4j}: "{current_detail_str}" }}) ')
                    else:
                        temp_command.append(f' MERGE (n{c}:{detail_neo4j} {{ {detail_neo4j}: "{current_detail_str}" }}) ')
                    temp_command.append(f' WITH n0,n{c} MERGE (n0)-[r:HAS_{detail_neo4j.upper()}]->(n{c}) ')
                    if temp_command:
                        temp_command=''.join(temp_command)
                        list_commands.append(temp_command)
                        c+=1
            else:
                if node_id not in self.drax_connections:
                    self.drax_connections[node_id]={}
                if detail_neo4j not in self.drax_connections[node_id]:
                    self.drax_connections[node_id][detail_neo4j]=set()
                self.drax_connections[node_id][detail_neo4j].update(node_info[detail])
        return_id=' RETURN ID(n0) AS neo4j_id '
        list_commands.append(return_id)
        merged_command=' '.join(list_commands)
        print(merged_command)
        self.run_command_neo4j(merged_command)






    #legacy?

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





    def add_set_convergence(self,node_instance,detail,current_detail):
        res=[]
        for c in node_instance.get_convergence_db(detail, current_detail):
            converged_db, converged_time = c
            res.append(' SET rc.Converged_in_' + converged_db + '="' + converged_time + '" ')
        if res :
            res.insert(0,'WITH rc ')
            return ''.join(res)
        return ''

    def load_convergence(self,node_instance,bio_instance,relationship_properties,instances_type):
        convergence = list(relationship_properties.keys())
        if convergence:
            convergence_dbs = [i.split('_')[-1] for i in convergence]
        else:
            convergence_dbs = []
        for db in convergence_dbs:
            node_instance.set_detail(instances_type, bio_instance, converged_in=db)

    def replace_detail(self,current_detail):
        res=str(current_detail)
        res = res.replace('\\', '\\\\')
        res = res.replace('"', '\\"')
        res = res.replace('\'', '\\\'')
        return res

    def get_instance_node_id_neo4j(self,node_id,node_type):
        if not hasattr(self,'neo4j_connector_memory' ):
            self.neo4j_connector_memory={}
        if node_type not in self.neo4j_connector_memory:
            self.neo4j_connector_memory[node_type]=[]
        for n in self.neo4j_connector_memory[node_type]:
            if n.get_detail('internal_id')== node_id:
                return n


    def add_instance(self,instance_to_add,node_type):
        if not hasattr(self,'neo4j_connector_memory' ):
            self.neo4j_connector_memory={}
        if node_type not in self.neo4j_connector_memory:
            self.neo4j_connector_memory[node_type]=[]
        self.neo4j_connector_memory[node_type].append(instance_to_add)




####CREATE NODE FROM INSTANCE####



    def generate_neo4j_cpd(self,node_instance):
        node_instance.saved()
        node_id=node_instance.get_detail('internal_id')
        match_main_node ='MATCH (n1:Compound) WHERE ID(n1)='+str(node_id)+' WITH n1 '
        for detail in node_instance.get_details_list():
            all_current_details= node_instance.get_detail(detail,all_possible=True)
            for current_detail in copy(all_current_details):
                res=[]
                if isinstance(current_detail,str):  current_detail_str=self.replace_detail(current_detail)
                else: current_detail_str=current_detail

                if detail=='synonyms':
                    c=str(all_current_details[current_detail])
                    res.append(' MERGE (n2:Synonym { Synonym: "' + current_detail_str + '" }) ')
                    res.append(' WITH n1,n2 MERGE (n1)-[r:HAS_SYNONYM]->(n2) ' \
                                    'SET r.Counter = "'+c+'" ')
                else:
                    c=str(all_current_details[current_detail])
                    res.append(' MERGE (n2:Identifier {Database: "'+detail+'", ID: "'+current_detail_str+'"}) ')
                    res.append(' WITH n1,n2 MERGE (n1)-[r:HAS_IDENTIFIER]->(n2) ' \
                                   'SET r.Counter = '+c+' ')
                if res:
                    res.insert(0,match_main_node)
                    self.run_command_neo4j(res)

    def create_node_neo4j_cpd(self,node_instance):
        creation_main_node ='CREATE (n1:Compound {node_type: "Compound"}) RETURN ID(n1) AS id_instance'
        command_run = self.run_command_neo4j(creation_main_node)
        id_instance = command_run[0]['id_instance']
        if not node_instance.get_detail('internal_id'):
            node_instance.set_detail('internal_id',id_instance)
        self.generate_neo4j_cpd(node_instance)


    def update_node_neo4j_cpd(self,node_instance):
        self.generate_neo4j_cpd(node_instance)




    def generate_neo4j_gene(self,node_instance):
        node_instance.saved()
        node_id=node_instance.get_detail('internal_id')
        match_main_node ='MATCH (n1:Gene) WHERE ID(n1)='+str(node_id)+' WITH n1 '
        for detail in node_instance.get_details_list():
            all_current_details = node_instance.get_detail(detail, all_possible=True)
            for current_detail in copy(all_current_details):
                res = []
                if isinstance(current_detail,str):  current_detail_str=self.replace_detail(current_detail)
                else: current_detail_str=current_detail
                if detail == 'synonyms':
                    c=str(all_current_details[current_detail])
                    res.append('  MERGE (n2:Synonym { Synonym: "' + current_detail_str + '" }) ')
                    res.append(' WITH n1,n2 MERGE (n1)-[r:HAS_SYNONYM]->(n2) ' \
                           'SET  r.Counter="' + c + '" ')

                elif detail=='aa_seq':
                    ES_ID=self.save_detail_ElasticSearch(detail,current_detail_str)
                    res.append(' MERGE (n2:Aminoacid_sequence { ES_ID: "'+ES_ID+'" }) ')
                    res.append(' WITH n1,n2 MERGE (n1)-[:HAS_AA_SEQ]->(n2)')

                elif detail=='nt_seq':
                    ES_ID=self.save_detail_ElasticSearch(detail,current_detail_str)
                    res.append('MERGE (n2:Nucleotide_sequence { ES_ID: "'+ES_ID+'" }) ')
                    res.append('WITH n1,n2 MERGE (n1)-[:HAS_NT_SEQ]->(n2)')

                elif detail.endswith('_instances'):

                    instance_to_add=self.save_instance_neo4j(current_detail)
                    if instance_to_add:

                        instance_id=str(instance_to_add.get_detail('internal_id'))
                        convergence_str=self.add_set_convergence(node_instance,detail,instance_to_add)
                        instance_type=get_instance_type(instance_to_add)
                        if convergence_str:
                            res.append(f' MATCH (n2:{instance_type}) WHERE ID(n2)={instance_id} ' \
                                      f' WITH n1,n2 MERGE (n1)-[rc:HAS_{instance_type.upper()}]->(n2) ')
                            res.append(convergence_str)

                else:
                    c=str(all_current_details[current_detail])
                    res.append('MERGE (n2:Identifier {Database: "'+detail+'", ID: "'+current_detail_str+'"}) ')
                    res.append(' WITH n1,n2  MERGE (n1)-[r:HAS_IDENTIFIER]->(n2) ' \
                                   'SET r.Counter = '+c+' ')
                if res:
                    res.insert(0,match_main_node)
                    self.run_command_neo4j(res)


    def create_node_neo4j_gene(self, node_instance):
        creation_main_node = 'CREATE (n1:Gene {node_type: "Gene"})  RETURN ID(n1) AS id_instance'
        command_run = self.run_command_neo4j(creation_main_node)
        id_instance = command_run[0]['id_instance']
        if not node_instance.get_detail('internal_id'):
            node_instance.set_detail('internal_id',id_instance)
        self.generate_neo4j_gene(node_instance)


    def update_node_neo4j_gene(self, node_instance):
        self.generate_neo4j_gene(node_instance)


    def generate_neo4j_protein(self, node_instance):
        node_instance.saved()
        node_id=node_instance.get_detail('internal_id')
        match_main_node ='MATCH (n1:Protein) WHERE ID(n1)='+str(node_id)+' WITH n1'
        for detail in node_instance.get_details_list():
            all_current_details = node_instance.get_detail(detail, all_possible=True)
            for current_detail in copy(all_current_details):
                if current_detail in node_instance.get_detail(detail, all_possible=True):

                    res=[]
                    if isinstance(current_detail,str):  current_detail_str=self.replace_detail(current_detail)
                    else: current_detail_str=current_detail
                    if detail == 'synonyms':
                        c=str(all_current_details[current_detail])
                        res.append(' MERGE (n2:Synonym { Synonym: "'+current_detail_str+'" }) ')
                        res.append(' WITH n1,n2 MERGE (n1)-[r:HAS_SYNONYM]->(n2) ' \
                                         'SET  r.Counter="' + c + '" ')

                    elif detail=='pathways':
                        res.append(' MERGE (n2:Pathway { Pathway: "'+current_detail_str+'" }) ')
                        res.append(' WITH n1,n2 MERGE (n1)-[:IN_PATHWAY]->(n2)')

                    elif detail=='aa_seq':
                        ES_ID=self.save_detail_ElasticSearch(detail,current_detail_str)
                        converged_in= node_instance.get_convergence_db('aa_seq', current_detail_str)
                        res.append(' MERGE (n2:Aminoacid_sequence { ES_ID: "'+ES_ID+'" }) ')
                        res.append(' WITH n1,n2 MERGE (n1)-[:HAS_AA_SEQ]->(n2)')


                    elif detail.endswith('_instances'):

                        instance_to_add=self.save_instance_neo4j(current_detail)

                        if instance_to_add:
                            instance_id=str(instance_to_add.get_detail('internal_id'))
                            convergence_str=self.add_set_convergence(node_instance,detail,instance_to_add)
                            instance_type=get_instance_type(instance_to_add)
                            if convergence_str:
                                res.append(f' MATCH (n2:{instance_type}) WHERE ID(n2)={instance_id} ' \
                                          f' WITH n1,n2 MERGE (n1)-[rc:HAS_{instance_type.upper()}]->(n2) ')
                                res.append(convergence_str)

                    else:
                        c=str(all_current_details[current_detail])
                        res.append(' MERGE (n2:Identifier {Database: "'+detail+'", ID: "'+current_detail_str+'"}) ')
                        res.append(' WITH n1,n2 MERGE (n1)-[r:HAS_IDENTIFIER]->(n2) ' \
                               ' SET r.Counter="' + c+ '" ')
                    if res:
                        res.insert(0,match_main_node)
                        self.run_command_neo4j(res)



    def create_node_neo4j_protein(self, node_instance):
        creation_main_node ='CREATE (n1:Protein {node_type: "Protein"}) RETURN ID(n1) AS id_instance'
        command_run = self.run_command_neo4j(creation_main_node)
        id_instance = command_run[0]['id_instance']
        if not node_instance.get_detail('internal_id'):
            node_instance.set_detail('internal_id',id_instance)
        self.generate_neo4j_protein(node_instance)


    def update_node_neo4j_protein(self,node_instance):
        self.generate_neo4j_protein(node_instance)



    def generate_neo4j_reaction(self,node_instance):
        node_instance.saved()
        node_id = node_instance.get_detail('internal_id')
        match_main_node = 'MATCH (n1:Reaction) WHERE ID(n1)=' + str(node_id) + ' WITH n1 '
        for detail in node_instance.get_details_list():
            if detail == 'reaction_with_instances':
                # current_detail will be left or right
                reaction_with_instances=node_instance.get_detail(detail)
                left=reaction_with_instances['left']
                right=reaction_with_instances['right']
                left_instances=[]
                right_instances=[]
                for stoi_n, cpd_inst in left:
                    cpd_inst = self.save_instance_neo4j(cpd_inst)
                    if not cpd_inst:
                        print('ERROR HERE', node_id, current_detail)
                        raise Exception
                    cpd_id=cpd_inst.get_detail('internal_id')
                    left_instances.append(cpd_id)
                for stoi_n, cpd_inst in right:
                    cpd_inst = self.save_instance_neo4j(cpd_inst)
                    if not cpd_inst:
                        print('ERROR HERE', node_id, current_detail)
                        raise Exception
                    cpd_id=cpd_inst.get_detail('internal_id')
                    right_instances.append(cpd_id)
                    side_dict={'left':left_instances,'right':right_instances}
                    for side in side_dict:
                        for cpd_id in side_dict[side]:
                            res = [match_main_node]
                            res.append(f' MATCH (cpd:Compound) WHERE ID(cpd)={cpd_id} ')
                            res.append(' WITH n1,cpd MERGE (n1)-[:HAS_CPD {N_molecules: ' + str(stoi_n) + ', side: "' + side + '" }]->(cpd) ')
                            self.run_command_neo4j(res)
                            res = [match_main_node]

            else:
                all_current_details = node_instance.get_detail(detail, all_possible=True)
                for current_detail in copy(all_current_details):
                    if current_detail in node_instance.get_detail(detail, all_possible=True):
                        res=[]
                        if isinstance(current_detail,str):  current_detail_str=self.replace_detail(current_detail)
                        else: current_detail_str=current_detail
                        if detail=='reaction_str':
                            c=str(all_current_details[current_detail])
                            sign=find_sign(current_detail_str)
                            if sign:
                                sign=sign.group()
                                sign=uniform_sign(sign)
                                if is_reversible(sign):
                                    set_reversibility_command='MATCH (n1:Reaction) WHERE ID(n1)=' + str(node_id) + ' AND (NOT EXISTS(n1.Reversible) OR n1.Reversible="False") set n1.Reversible="True"'
                                else:
                                    set_reversibility_command='MATCH (n1:Reaction) WHERE ID(n1)=' + str(node_id) + ' AND NOT EXISTS(n1.Reversible) SET n1.Reversible="False"'
                                self.run_command_neo4j(set_reversibility_command)
                            res.append('MERGE (n2:Reaction_str {Reaction_str: "'+current_detail_str+'"}) ')

                            res.append(' WITH n1,n2 MERGE (n1)-[r:HAS_REACTION_STR]->(n2) ' \
                                   ' ON CREATE SET ' \
                                   ' n2.Reaction_str="' + current_detail_str + '", ' \
                                   ' r.Counter="'+c+'" ' \
                                   ' ON MATCH SET ' \
                                   ' r.Counter="'+c+'" ')
                        elif detail=='pathways':
                            res.append('MERGE (n2:Pathway {Pathway: "'+current_detail_str+'"}) ')
                            res.append(' WITH n1,n2 MERGE (n1)-[:IN_PATHWAY]->(n2) ')
                        elif detail=='rn_with_ids':
                            res.append('MERGE (n2:rn_with_ids {rn_with_ids_str: "'+self.replace_detail(current_detail[0])+'", rn_with_ids_ids: "'+self.replace_detail(str(current_detail[1]))+'", rn_with_ids_db: "'+current_detail[2]+'"}) ')
                            res.append(' WITH n1,n2 MERGE (n1)-[:HAS_RN_WITH_IDS]->(n2) ')

                        elif detail.endswith('_instances') and detail!='reaction_with_instances':

                            instance_to_add=self.save_instance_neo4j(current_detail)
                            if instance_to_add:

                                instance_id=str(instance_to_add.get_detail('internal_id'))
                                convergence_str=self.add_set_convergence(node_instance,detail,instance_to_add)
                                instance_type=get_instance_type(instance_to_add)
                                if convergence_str:
                                    res.append(f' MATCH (n2:{instance_type}) WHERE ID(n2)={instance_id} ' \
                                              f' WITH n1,n2 MERGE (n1)-[rc:HAS_{instance_type.upper()}]->(n2) ')
                                    res.append(convergence_str)


                        else:
                            c=str(all_current_details[current_detail])
                            res.append('MERGE (n2:Identifier {Database: "'+detail+'", ID: "'+current_detail_str+'"}) ')
                            res.append(' WITH n1,n2 MERGE (n1)-[r:HAS_IDENTIFIER]->(n2) ' \
                                   'SET  r.Counter="' + c+ '" ')
                        if res:
                            res.insert(0,match_main_node)
                            self.run_command_neo4j(res)

    def create_node_neo4j_reaction(self, node_instance):
        creation_main_node ='CREATE (n1:Reaction {node_type: "Reaction"}) RETURN ID(n1) AS id_instance'
        command_run = self.run_command_neo4j(creation_main_node)
        id_instance = command_run[0]['id_instance']
        if not node_instance.get_detail('internal_id'):
            node_instance.set_detail('internal_id',id_instance)
        self.generate_neo4j_reaction(node_instance)


    def update_node_neo4j_reaction(self, node_instance):
        self.generate_neo4j_reaction(node_instance)

    def generate_neo4j_organism(self, node_instance):
        node_instance.saved()
        node_id = node_instance.get_detail('internal_id')
        match_main_node = 'MATCH (n1:Organism) WHERE ID(n1)=' + str(node_id) + ' WITH n1 '

        for detail in node_instance.get_details_list():
            all_current_details = node_instance.get_detail(detail, all_possible=True)
            for current_detail in copy(all_current_details):
                if isinstance(current_detail, str):
                    current_detail_str = self.replace_detail(current_detail)
                else:
                    current_detail_str = current_detail
                res = []
                if detail == 'synonyms':
                    c = str(all_current_details[current_detail])
                    res.append(' MERGE (n2:Synonym {Synonym: "' + current_detail_str + '"}) ')
                    res.append(' WITH n1,n2 MERGE (n1)-[r:HAS_SYNONYM]->(n2) ' \
                               'SET r.Counter="' + c + '" ')
                elif detail.endswith('_instances'):
                    instance_to_add=self.save_instance_neo4j(current_detail)
                    if instance_to_add:
                        instance_id=str(instance_to_add.get_detail('internal_id'))
                        convergence_str=self.add_set_convergence(node_instance,detail,instance_to_add)
                        instance_type=get_instance_type(instance_to_add)
                        if convergence_str:
                            res.append(f' MATCH (n2:{instance_type}) WHERE ID(n2)={instance_id} ' \
                                      f' WITH n1,n2 MERGE (n1)-[rc:HAS_{instance_type.upper()}]->(n2) ')
                            res.append(convergence_str)
                else:
                    c = str(all_current_details[current_detail])
                    res.append(' MERGE (n2:Identifier {Database: "' + detail + '", ID: "' + current_detail_str + '"}) ')
                    res.append(' WITH n1,n2 MERGE (n1)-[r:HAS_IDENTIFIER]->(n2) ' \
                               'SET  r.Counter="' + c + '" ')
                if res:
                    res.insert(0, match_main_node)
                    self.run_command_neo4j(res)

    def create_node_neo4j_organism(self, node_instance):
        creation_main_node = 'CREATE (n1:Organism {node_type: "Organism"}) RETURN ID(n1) AS id_instance'
        command_run = self.run_command_neo4j(creation_main_node)
        id_instance = command_run[0]['id_instance']
        if not node_instance.get_detail('internal_id'):
            node_instance.set_detail('internal_id', id_instance)
        self.generate_neo4j_organism(node_instance)


    def update_node_neo4j_organism(self, node_instance):
        self.generate_neo4j_organism(node_instance)

####FIND QUERY MATCHES####

    def find_matches_neo4j_str(self,bio_query,bio_db,node_type,best_matches=True,load_instances=True,get_all_instances_connections=False):
        """
        :best_matches will benefit the nodes which have a higher counter in the respective bio_query
        """
        if not bio_query: return []
        # sometimes we try to match strings which may have regex special characters
        bio_query = regex_escape(bio_query)
        if bio_db=='synonyms':
            total_command=['MATCH (n1:'+node_type+')-[r:HAS_SYNONYM]->(n2:Synonym) WHERE n2.Synonym=~"(?i)'+bio_query+'" WITH n1,r,n2 ']
        elif bio_db=='reaction_str':
            total_command=['MATCH (n1:'+node_type+')-[r:HAS_REACTION_STR]->(n2:Reaction_str) WHERE n2.Reaction_str=~"(?i)'+bio_query+'" WITH n1,r,n2 ']
        else:
            total_command=['MATCH (n1:'+node_type+')-[r:HAS_IDENTIFIER]->(n2) WHERE n2.Database=~"(?i)'+bio_db+'" AND n2.ID=~"(?i)'+bio_query+'"  WITH n1,r,n2 ']
        total_command.append('MATCH (n1:'+node_type+')-[all_r]->() RETURN id(n1) as id,r.Counter as Counter,COUNT(r) as N_Relationships ORDER BY COUNT(r) DESC')
        matches=self.run_command_neo4j(total_command)
        res=[]
        best_res={}
        for r in matches:
            res.append(r['id'])
            if best_matches:
                c=int(r['Counter'])
                #we benefit nodes with lots of relationships
                n= ceil(int(r['N_Relationships'])/len(matches))
                if n>1: c+=n
                if c in best_res: best_res[c].append(r['id'])
                else: best_res[c]=[r['id']]
        if best_matches and res:
            max_counter = max(best_res, key=int)
            res=best_res[max_counter]
        if load_instances:res= list(self.load_instance_neo4j_list(nodes_id=res,node_type=node_type,get_all_instances_connections=get_all_instances_connections))
        return res




####FIND INSTANCE MATCHES####

    def check_match_reaction_with_instances(self,reaction_with_instances):
        res=[]
        var_counter=0
        for side in reaction_with_instances:
            for n_stoi,cpd in reaction_with_instances[side]:
                db_instance= self.find_matches_neo4j_instance(cpd)
                if db_instance: cpd=db_instance
                else:           cpd=self.save_instance_neo4j(cpd)
                cpd_id=cpd.get_detail('internal_id')
                if var_counter==0:  res.append('MATCH (rn:Reaction)-')
                else:               res.append('MATCH (rn)-')
                res.append('[:HAS_CPD {N_molecules:'+str(n_stoi)+'}]->(cpd'+str(var_counter)+':Compound)  WHERE ID(cpd'+str(var_counter)+')='+str(cpd_id)+' ')
                res.append(' WITH rn, cpd'+str(var_counter)+' ')
                var_counter+=1
        res.append(' MATCH (rn)-[r:HAS_CPD]->() with rn,COUNT(r) as r_count WHERE r_count='+str(var_counter))
        res.append(' RETURN ID(rn) as id_instance')
        data=self.run_command_neo4j(res)
        if not data: return None
        id_instance = data[0]['id_instance']
        return id_instance



    def clean_up_matches(self,instance_matches,node_type):
        # if there's a match in the graph we take the first match -> since this is the one with the most relationships
        # could be any other but usually the number of matches should be low - if everything is working properly
        main_node = instance_matches[0]
        for n in instance_matches[1:]:
            id_to_remove=n.get_detail('internal_id')
            main_node.unite_instances(n)
            #so objects pointers change accordingly
            n=main_node
            # we unite the instances and remove the extra nodes
            if id_to_remove:
                self.remove_node_neo4j(id_to_remove, node_type=node_type)
        #main_node.require_save()
        self.remove_unconnected_nodes()
        return main_node

    def find_matches_neo4j_instance(self,instance_to_match,get_all_instances_connections=False):
        if instance_to_match.get_detail('internal_id'): return instance_to_match
        node_type=get_instance_type(instance_to_match)
        possible_matches=set()
        for unique_detail in get_unique_details(bio_instance_type=node_type):
            all_current_details = instance_to_match.get_detail(unique_detail, all_possible=True)
            if all_current_details:
                if unique_detail == 'reaction_with_instances':
                    match_cpd_instances=self.check_match_reaction_with_instances(all_current_details)
                    if match_cpd_instances: possible_matches.add(match_cpd_instances)
                else:
                    for current_detail in copy(all_current_details):
                        if current_detail:
                            matches=self.find_matches_neo4j_str(current_detail,unique_detail,node_type=node_type,load_instances=False)
                            for m in matches:
                                possible_matches.add(m)
        to_merge=[]
        for n in self.load_instance_neo4j_list(nodes_id=possible_matches,node_type=node_type,get_all_instances_connections=get_all_instances_connections):
            if n.is_match_instances(instance_to_match):
                to_merge.append(n)
        if not to_merge: return None
        #to_merge.append(instance_to_match)
        main_node=self.clean_up_matches(to_merge,node_type)
        return main_node


####LOAD INSTANCES FROM NODE IDS####

    def load_instance_neo4j_list(self,nodes_id,node_type,get_all_instances_connections=False):
        for n in nodes_id:
            node_instance= self.load_instance_neo4j(node_id=n,node_type=node_type,get_all_instances_connections=get_all_instances_connections)
            if node_instance: yield node_instance

    def load_instance_neo4j(self,node_id,node_type='',get_all_instances_connections=False):
        node_instance=self.get_instance_node_id_neo4j(node_id,node_type)
        if node_instance:
            return node_instance
        if node_type:c_node_type=':'+node_type
        else: c_node_type=''
        if get_all_instances_connections:
            total_command='MATCH (n1'+c_node_type+')-[r]-(n2) WHERE ID(n1)='+str(node_id)+' RETURN properties(n1) as n1_properties, ' \
                                                                                'type(r) AS relationship_type, ' \
                                                                                'properties(r) AS r_properties, ' \
                                                                                'properties(n2) AS n2_properties, ' \
                                                                                'ID(n2) as n2_id'
        else:
            #returning central node (node id) and connecting nodes and the relationships that connect them
            total_command='MATCH (n1'+c_node_type+')-[r]->(n2) WHERE ID(n1)='+str(node_id)+' RETURN properties(n1) as n1_properties, ' \
                                                                                'type(r) AS relationship_type, ' \
                                                                                'properties(r) AS r_properties, ' \
                                                                                'properties(n2) AS n2_properties, ' \
                                                                                'ID(n2) as n2_id'

        node_relationships=self.run_command_neo4j(total_command)
        node_properties=node_relationships[0]['n1_properties']
        if not node_relationships: return
        if node_type: n1_type=node_type
        else: n1_type = node_properties['node_type']
        if n1_type== 'Compound':
            node_instance= self.load_instance_neo4j_cpd(node_id,node_relationships)
        elif n1_type== 'Gene':
            node_instance= self.load_instance_neo4j_gene(node_id,node_relationships)
        elif n1_type== 'Protein':
            node_instance= self.load_instance_neo4j_protein(node_id,node_relationships)
        elif n1_type== 'Reaction':
            node_instance= self.load_instance_neo4j_reaction(node_id,node_relationships,node_properties)
        elif n1_type == 'Organism':
            node_instance = self.load_instance_neo4j_organism(node_id, node_relationships)
        node_instance.saved()
        return node_instance


    def load_instance_neo4j_cpd(self,node_id,node_relationships):
        node_instance=Compound({'internal_id':node_id})
        #we add it here to avoid recursivity
        self.add_instance(node_instance)
        for nr in node_relationships:
            relationship_type = nr['relationship_type']
            relationship_properties = nr['r_properties']
            n2_properties = nr['n2_properties']

            if relationship_type == 'HAS_SYNONYM':
                syn = n2_properties['Synonym']
                c = int(relationship_properties['Counter'])
                node_instance.set_detail('synonyms',{syn:c})

            elif relationship_type == 'HAS_IDENTIFIER':
                db = n2_properties['Database']
                db_id = n2_properties['ID']
                c = int(relationship_properties['Counter'])
                node_instance.set_detail(db,{db_id:c})

        return node_instance


    def load_instance_neo4j_gene(self,node_id,node_relationships):
        node_instance=Gene({'internal_id':node_id})
        #we add it here to avoid recursivity
        self.add_instance(node_instance)
        node_instance_type=get_instance_type(node_instance)
        for nr in node_relationships:
            relationship_type = nr['relationship_type']
            relationship_properties = nr['r_properties']
            n2_properties = nr['n2_properties']
            n2_id = nr['n2_id']

            if relationship_type == 'HAS_SYNONYM':
                syn = n2_properties['Synonym']
                c = int(relationship_properties['Counter'])
                node_instance.set_detail('synonyms',{syn:c})

            elif relationship_type == 'HAS_AA_SEQ':
                list_detail = n2_properties['ES_ID']
                detail_str=self.get_detail_ElasticSearch('aa_seq',list_detail)
                node_instance.set_detail('aa_seq',detail_str)
            elif relationship_type == 'HAS_NT_SEQ':
                list_detail = n2_properties['ES_ID']
                detail_str=self.get_detail_ElasticSearch('nt_seq',list_detail)
                node_instance.set_detail('nt_seq',detail_str)


            elif relationship_type in ['HAS_PROTEIN','HAS_REACTION','HAS_GENE','HAS_ORGANISM']:
                node_type = nr['n2_properties']['node_type']
                instances_type=f'{node_type.lower()}_instances'
                bio_instance=self.load_instance_neo4j(n2_id,node_type=node_type)
                # here we dont load convergence since this only appears if we are not looking for directed graph connections
                if relationship_type==f'HAS_{node_instance_type.upper()}':
                    node_instance.set_detail(instances_type, bio_instance)
                else:
                    self.load_convergence(node_instance=node_instance, bio_instance=bio_instance,
                                          relationship_properties=relationship_properties,
                                          instances_type=instances_type)


            elif relationship_type == 'HAS_IDENTIFIER':
                db = n2_properties['Database']
                db_id = n2_properties['ID']
                c = int(relationship_properties['Counter'])
                node_instance.set_detail(db,{db_id:c})

        return node_instance

    def load_instance_neo4j_protein(self,node_id,node_relationships):
        node_instance=Protein({'internal_id':node_id})
        #we add it here to avoid recursivity
        self.add_instance(node_instance)
        eocs_to_add=[]
        node_instance_type=get_instance_type(node_instance)
        for nr in node_relationships:
            relationship_type = nr['relationship_type']
            relationship_properties = nr['r_properties']
            n2_properties = nr['n2_properties']
            n2_id = nr['n2_id']

            if relationship_type == 'HAS_SYNONYM':
                syn = n2_properties['Synonym']
                c = int(relationship_properties['Counter'])
                node_instance.set_detail('synonyms',{syn:c})

            elif relationship_type == 'HAS_AA_SEQ':
                list_detail = n2_properties['ES_ID']
                detail_str=self.get_detail_ElasticSearch('aa_seq',list_detail)
                node_instance.set_detail('aa_seq',detail_str)


            elif relationship_type == 'IN_PATHWAY':
                list_detail = n2_properties['Pathway']
                node_instance.set_detail('pathways',list_detail)


            elif relationship_type in ['HAS_PROTEIN','HAS_REACTION','HAS_GENE','HAS_ORGANISM']:
                node_type = nr['n2_properties']['node_type']
                instances_type=f'{node_type.lower()}_instances'
                bio_instance=self.load_instance_neo4j(n2_id,node_type=node_type)
                # here we dont load convergence since this only appears if we are not looking for directed graph connections
                if relationship_type==f'HAS_{node_instance_type.upper()}':
                    node_instance.set_detail(instances_type, bio_instance)
                else:
                    self.load_convergence(node_instance=node_instance, bio_instance=bio_instance,
                                          relationship_properties=relationship_properties,
                                          instances_type=instances_type)

            elif relationship_type == 'HAS_IDENTIFIER':
                db = n2_properties['Database']
                db_id = n2_properties['ID']
                c = int(relationship_properties['Counter'])
                node_instance.set_detail(db,{db_id:c})


        return node_instance

    def load_instance_neo4j_reaction(self,node_id,node_relationships,node_properties):
        node_instance=Reaction({'internal_id':node_id})
        reversibility=bool(node_properties['Reversible'])
        node_instance.is_reversible=reversibility
        #we add it here to avoid recursivity
        self.add_instance(node_instance)
        reaction_with_instances={'left':[],'right':[]}
        node_instance_type=get_instance_type(node_instance)

        for nr in node_relationships:
            relationship_type = nr['relationship_type']
            relationship_properties = nr['r_properties']
            n2_properties = nr['n2_properties']
            n2_id = nr['n2_id']

            if relationship_type == 'HAS_REACTION_STR':
                rn_str = n2_properties['Reaction_str']
                c = int(relationship_properties['Counter'])
                node_instance.set_detail('reaction_str',{rn_str:c})

            elif relationship_type == 'IN_PATHWAY':
                list_detail = n2_properties['Pathway']
                node_instance.set_detail('pathways',list_detail)

            elif relationship_type == 'HAS_RN_WITH_IDS':
                rn_with_ids_str = n2_properties['rn_with_ids_str']
                rn_with_ids_ids = n2_properties['rn_with_ids_ids']
                rn_with_ids_db = n2_properties['rn_with_ids_db']
                node_instance.set_detail('rn_with_ids',[rn_with_ids_str,rn_with_ids_ids,rn_with_ids_db])


            elif relationship_type in ['HAS_PROTEIN','HAS_REACTION','HAS_GENE','HAS_ORGANISM']:
                node_type = nr['n2_properties']['node_type']
                instances_type=f'{node_type.lower()}_instances'
                bio_instance=self.load_instance_neo4j(n2_id,node_type=node_type)
                # here we dont load convergence since this only appears if we are not looking for directed graph connections
                if relationship_type==f'HAS_{node_instance_type.upper()}':
                    node_instance.set_detail(instances_type, bio_instance)
                else:
                    self.load_convergence(node_instance=node_instance, bio_instance=bio_instance,
                                          relationship_properties=relationship_properties,
                                          instances_type=instances_type)


            elif relationship_type == 'HAS_CPD':
                side=relationship_properties['side']
                n_molecules=relationship_properties['N_molecules']
                bio_instance=self.load_instance_neo4j(n2_id,node_type='Compound')
                reaction_with_instances[side].append([n_molecules,bio_instance])
            elif relationship_type == 'HAS_IDENTIFIER':
                db = n2_properties['Database']
                db_id = n2_properties['ID']
                c = int(relationship_properties['Counter'])
                node_instance.set_detail(db,{db_id:c})

        node_instance.set_detail('reaction_with_instances', reaction_with_instances)


        return node_instance

    def load_instance_neo4j_organism(self,node_id,node_relationships):
        node_instance=Organism({'internal_id':node_id})
        #we add it here to avoid recursivity
        self.add_instance(node_instance)
        for nr in node_relationships:
            relationship_type = nr['relationship_type']
            relationship_properties = nr['r_properties']
            n2_properties = nr['n2_properties']
            if relationship_type == 'HAS_SYNONYM':
                syn = n2_properties['Synonym']
                c = int(relationship_properties['Counter'])
                node_instance.set_detail('synonyms',{syn:c})
            elif relationship_type == 'HAS_IDENTIFIER':
                db = n2_properties['Database']
                db_id = n2_properties['ID']
                c = int(relationship_properties['Counter'])
                node_instance.set_detail(db,{db_id:c})
        return node_instance


def test(neo4_driver):
    c1={'chebi': {'38665'}, 'chemical_formula': {'C8H19O2PS2', 'C8H19O2P1S2'}, 'inchi': {'1S/C8H19O2PS2/c1-4-7-12-11(9,10-6-3)13-8-5-2/h4-8H2,1-3H3'}, 'inchi_key': {'VJYFKVYYMZPMAB-UHFFFAOYSA-N'}, 'kegg': {'C18687'}, 'metacyc': {'CPD-22347'}, 'pubchem_cid': {'3289'}, 'smiles': {'CCCSP(=O)(OCC)SCCC', 'CCCSP(SCCC)(=O)OCC'}, 'synonyms': {'o-ethyl s,s-dipropyl phosphorodithioate', 'ethoprophos', 'ethoprop', 'o-ethyl s,s-dipropyl dithiophosphate'}}
    c2={'chebi': {'8082'}, 'chemical_formula': {'C8H8O2'}, 'inchi': {'1S/C8H8O2/c1-7(9)10-8-5-3-2-4-6-8/h2-6H,1H3'}, 'inchi_key': {'IPBVNPXQWQGGJP-UHFFFAOYSA-N'}, 'kegg': {'C15583', 'C00548'}, 'metacyc': {'CPD-22363'}, 'pubchem_cid': {'31229'}, 'smiles': {'CC(=O)Oc1ccccc1', 'CC(OC1(\\C=C/C=C\\C=1))=O'}, 'synonyms': {'phenol acetate', 'phenyl acetate', 'acetic acid,phenyl ester', 'acetylphenol'}}
    list_instances={'compound':[c1,c2]}





if __name__ == '__main__':
    import time
    # Database Credentials
    uri = 'bolt://localhost:7687'
    username = 'neo4j'
    password = 'drax_neo4j'
    neo4j_driver = Neo4j_Connector(username, password,db_name='neo4j',uri=uri)
    drax_output_folder='/home/pedroq/Desktop/test_drax/out/'
    neo4j_driver.export_drax_to_neo4j(drax_output_folder)

