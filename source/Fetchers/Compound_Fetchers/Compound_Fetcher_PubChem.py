
from source.Fetchers.Compound_Fetchers.Compound_Fetcher import *

class Compound_Fetcher_PubChem(Compound_Fetcher):
    def __init__(self,compound_id,db,memory_storage=None):
        Compound_Fetcher.__init__( self, compound_id,memory_storage)
        self.db=db
        self.compound=self.get_compound_pubchem()
        self.add_compound()

    def convert_inchi_to_pubchem_id(self):
        json_result=None
        if self.db=='inchi_key':
            webpage=self.get_with_fetcher('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/cids/JSON/',
                                          data={'inchikey': f'{self.compound_id}'},original_response=True)
            if webpage:
                json_result=webpage.json()
        elif self.db=='inchi':
            webpage=self.get_with_fetcher('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchi/cids/JSON/',
                                          data={'inchi': f'InChI={self.compound_id}'},original_response=True)
            if webpage:
                json_result=webpage.json()
        if json_result:
            return 'cid',json_result['IdentifierList']['CID'][0]
        else:
            if self.db == 'inchi_key':
                webpage = self.get_with_fetcher('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/sids/JSON/',
                                                data={'inchikey': f'{self.compound_id}'}, original_response=True)
                if webpage:
                    json_result = webpage.json()
            elif self.db == 'inchi':
                webpage = self.get_with_fetcher('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchi/sids/JSON/',
                                                data={'inchi': f'InChI={self.compound_id}'}, original_response=True)
                if webpage:
                    json_result = webpage.json()
        if json_result:
            return 'sid',json_result['IdentifierList']['SID'][0]
        else:
            return None,None

    def get_synonyms_pubchem(self,cpd_dict,cpd_type,id_type,pubchem_id):
        synonyms_url=f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/{cpd_type}/{id_type}/{pubchem_id}/synonyms/json'
        #https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/10920476/synonyms/json
        webpage = self.get_with_fetcher(synonyms_url, original_response=True)
        if webpage:
            json_result=webpage.json()
            cpd_dict['synonyms']=set()
            cpd_dict['synonyms'].update(json_result['InformationList']['Information'][0]['Synonym'])

    def get_ids_pubchem(self,cpd_dict,cpd_type,id_type,pubchem_id):
        db_ids_url=f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/{cpd_type}/{id_type}/{pubchem_id}/xrefs/SBURL/JSON'
        webpage = self.get_with_fetcher(db_ids_url, original_response=True)
        if webpage:
            json_result=webpage.json()
            all_db_links=json_result['InformationList']['Information'][0]['SBURL']
            wanted_dbs={'www.chemspider.com/Chemical-Structure.':'chemspider',
                        'chebiID=CHEBI:':'chebi',
                        'genome.jp/dbget-bin/www_bget?cpd:':'kegg',
                        'genome.jp/dbget-bin/www_bget?dr:':'kegg',
                        'www.hmdb.ca/metabolites/':'hmdb',
                        'biocyc.org/compound?orgid=META&id=':'metacyc',
                        'drugbank.ca/drugs/':'drugbank',
                        'ebi.ac.uk/chembl/compound_report_card/':'chembl',
                        'inchikey=':'inchi_key',
                        'InChI%3D':'inchi',
                        'InChI=':'inchi',
                        }
            for db_link in all_db_links:
                for w_db in wanted_dbs:
                    if w_db in db_link:
                        db_type=wanted_dbs[w_db]
                        if db_type not in cpd_dict: cpd_dict[db_type]=set()
                        db_link=db_link.split(w_db)[1]
                        db_link=db_link.replace('.html','')
                        #html encoding for "line feed"
                        db_link=db_link.replace('%0A','')
                        db_link=strip_tags(db_link)
                        cpd_dict[db_type].add(db_link)
                        break

    def get_compound_pubchem(self):
        cpd_dict={self.db:{self.compound_id}}
        if self.db in ['inchi_key','inchi']:
            id_type,pubchem_id=self.convert_inchi_to_pubchem_id()
        else:
            id_type = self.db.split('_')[1]
            pubchem_id = self.compound_id
        if not id_type or not pubchem_id: return None
        cpd_type='substance' if id_type=='sid' else 'compound'
        #has too much junk in the synonyms
        #self.get_synonyms_pubchem(cpd_dict,cpd_type,id_type,pubchem_id)
        self.get_ids_pubchem(cpd_dict,cpd_type,id_type,pubchem_id)
        compound_instance = Compound(cpd_dict)
        return compound_instance

    def converge_compound_global(self):
        pass


if __name__ == '__main__':
    search= Compound_Fetcher_PubChem('HONVBDHDLJUKLX-UHFFFAOYSA-N','inchi_key')
    search= Compound_Fetcher_PubChem('962','pubchem_cid')
    search= Compound_Fetcher_PubChem('1S/C5H7NO2/c6-5-3(7)1-2-4(5)8/h5H,1-2,6H2','inchi')
    search.get_compound().get_all_info()
    print('###')
    search.get_compound().get_top_n_synonyms()
