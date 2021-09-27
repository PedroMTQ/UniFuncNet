
from source.Fetchers.Compound_Fetchers.Compound_Fetcher import *
import xml.etree.cElementTree as celementTree


class Compound_Fetcher_Biocyc(Compound_Fetcher):
    def __init__(self,compound_id,memory_storage=None):
        Compound_Fetcher.__init__( self, compound_id,memory_storage)
        self.db='biocyc'
        self.compound=self.get_compound_biocyc()
        self.add_compound()

    def parse_biocyc_xml(self,xml):
        syns = []
        res ={}
        res['biocyc'] = self.compound_id
        if 'selenium' in str(type(xml)):content=xml.page_source
        else: content=xml.content
        context = celementTree.fromstring(content)
        if  self.compound_id == 'E-':
            res['synonyms'] = ['electron', 'e-']
            return res

        context = context[1]  # compound part
        for i in context:
            if i.tag == 'cml':
                for first_nest in i:
                    for second_nest in first_nest:
                        if 'title' in second_nest.attrib:
                            if second_nest.attrib['title'] == 'smiles':
                                res['SMILES'] = second_nest.text
                        if second_nest.tag == 'formula':
                            res['Chemical_formula'] = second_nest.attrib['concise'].replace(' ', '')
        
            if i.tag == 'synonym':            syns.append(strip_tags(i.text).lower())
            if i.tag == 'common-name':        syns.insert(0, strip_tags(i.text).lower())
            if i.tag == 'dblink':  # get the other DBs IDs
                for j in i:
                    if j.tag == 'dblink-db':
                        temp_db = j.text
                    if j.tag == 'dblink-oid':
                        if temp_db == 'HMDB':
                            temp_db = 'hmdb'
                        elif temp_db == 'CHEMSPIDER':
                            temp_db = 'chemspider'
                        elif temp_db == 'LIGAND-CPD':
                            temp_db = 'kegg'
                        elif temp_db == 'CHEBI':
                            temp_db = 'chebi'
                        elif temp_db == 'PUBCHEM':
                            temp_db = 'pubchem_cid'
                        elif temp_db == 'CAS':
                            temp_db = 'cas'
                        elif temp_db == 'KNAPSACK':
                            temp_db = 'knapsack'
                        elif temp_db == 'BIGG':
                            temp_db = 'bigg'
                        res[temp_db] = j.text
        res['synonyms'] = syns
        res=self.remove_unwanted_info(res)
        return res

    def get_compound_biocyc(self):
        url = 'https://websvc.biocyc.org/getxml?id=META:' + self.compound_id
        xml = self.get_with_fetcher(url,original_response=True)
        if not xml: return None
        parsed_info= self.parse_biocyc_xml(xml)
        if 'selenium.webdriver' in repr(xml): xml.quit()
        if number_of_nones_dict(parsed_info)==len(parsed_info): return None
        return Compound(parsed_info)


if __name__ == '__main__':
    search= Compound_Fetcher_Biocyc('ATP')
