
from source.Fetchers.Compound_Fetchers.Compound_Fetcher import *

class Compound_Fetcher_Chemspider(Compound_Fetcher):
    def __init__(self,compound_id,memory_storage=None,search_by_inchi_key=False):
        Compound_Fetcher.__init__( self, compound_id,memory_storage)
        if search_by_inchi_key:
            self.db='inchi_key'
        else:
            self.db='chemspider'
        self.compound=self.get_compound_chemspider()
        self.add_compound()

    def get_compound_chemspider(self):
        xpath_list = ['//*[@id="btnTabListModalTrigger1"]',
                      '//*[@id="ctl00_ctl00_ContentSection_ContentPlaceHolder1_RecordViewTabDetailsControl_data_sources_more_menu"]']
        dbs_to_look_for = {'hmdb.ca/metabolites/': 'hmdb',
                           #'foodb.ca/compounds/': 'food_db',
                           'ebi.ac.uk/chebi/searchId': 'chebi',
                           'ebi.ac.uk/chembl/compound': 'chembl',
                           'drugbank.ca/': 'drugbank',
                           'genome.jp.*cpd:': 'kegg',
                           'biocyc.org/compound': 'biocyc',
                           'rcsb.org/ligand/': 'pdb',
                           'pubchem.ncbi.nlm.nih.gov/summary': 'pubchem_cid'}
        if self.db=='inchi_key':
            url = f'http://www.chemspider.com/inchikey={self.compound_id}'
        else:
            url = f'http://www.chemspider.com/Chemical-Structure.{self.compound_id}.html'
        webpage = self.get_with_fetcher(url,selenium=True,xpath=xpath_list,timer=5)
        if not webpage: return None
        soup = BeautifulSoup(webpage, 'lxml')
        title_details = ['Molecular Formula', 'ChemSpider ID']
        res = {}
        res['synonyms']= set()
        for detail in title_details:
            found_detail = soup.find_all('span', class_='prop_title', text=detail)
            if found_detail:
                found_detail = strip_tags(str(found_detail[0].nextSibling))
                if detail == 'Molecular Formula':     res['chemical_formula'] = found_detail
                if detail == 'ChemSpider ID':         res['chemspider'] = found_detail
        main_name_tag = 'ctl00_ctl00_ContentSection_ContentPlaceHolder1_RecordViewDetails_rptDetailsView_ctl00_WrapTitle'
        name = soup.find_all('span', {'id': main_name_tag})
        if name:
            name = name[0].text
            res['synonyms'] = {name}
        more_details = ['Systematic name', 'SMILES']
        for detail in more_details:
            found_detail = soup.find_all('span', class_='prop_title', text=detail)
            if found_detail:
                found_detail = strip_tags(str(found_detail[0].findNext('p')))
                found_detail = found_detail.strip('\n').split('\n')[0]
                if detail == 'Systematic name':
                    res['synonyms'].add(found_detail)
                elif detail == 'SMILES':
                    res['smiles'] = found_detail
        # Now to find the DBs IDs
        for db in dbs_to_look_for:
            db_ids = soup.find_all('a', {'href': re.compile(db)})
            if db_ids:
                db_id = db_ids[0].text
                if dbs_to_look_for[db] == 'chebi':
                    chebi_search = re.search('CHEBI:', db_id)
                    if chebi_search: db_id = db_id[chebi_search.span()[1]:]
                    else: db_id=None
                elif dbs_to_look_for[db] == 'drugbank':
                    if db_id.lower()=='drugbank':
                        db_id = db_ids[1].text
                if db_id:
                    res[dbs_to_look_for[db]] = db_id
        res=self.remove_unwanted_info(res)
        if number_of_nones_dict(res)==len(res): return None
        compound_instance = Compound(res)
        return compound_instance

    def converge_compound_global(self):
        pass


if __name__ == '__main__':
    search= Compound_Fetcher_Chemspider('58575')
    search.get_compound().get_all_info()
