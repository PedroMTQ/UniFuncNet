
from source.Fetchers.Compound_Fetchers.Compound_Fetcher import *
import xml.etree.cElementTree as celementTree


class Compound_Fetcher_Biocyc(Compound_Fetcher):
    def __init__(self,compound_id,memory_storage=None):
        Compound_Fetcher.__init__( self, compound_id,memory_storage)
        self.db='biocyc'
        self.set_convergence_args()

        self.compound=self.get_compound_biocyc()
        self.add_compound()


    def get_reactions(self):
        res=set()
        link = f'https://biocyc.org/compound?orgid=META&id={self.compound_id}#RXNS=&tab=RXNS'
        source = self.get_with_fetcher(link,selenium=True)
        if not source: return None
        soup = BeautifulSoup(source, 'lxml')
        reactions = soup.find_all('td', class_='reactionEquation')
        for r in reactions:
            r_link = r.findChild()['href']
            rn_id = r_link[re.search('object=', r_link).span()[1]:]
            res.add(rn_id)
        return res

    def set_convergence_args(self):
        #args for convergence
        self.convergence_args['reactions'] = set()


    def parse_biocyc_xml(self,xml):
        syns = []
        res ={}
        res['biocyc'] = self.compound_id
        if 'selenium' in str(type(xml)):content=xml.page_source
        else: content=xml.content
        try:
            context = celementTree.fromstring(content)
        except: return None
        if  self.compound_id == 'E-':
            res['synonyms'] = ['electron', 'e-']
            return res
        try:
            context = context[1]  # compound part
        except:
            return None
        for i in context:
            if i.tag == 'cml':
                for first_nest in i:
                    for second_nest in first_nest:
                        if 'title' in second_nest.attrib:
                            if second_nest.attrib['title'] == 'smiles':
                                res['smiles'] = second_nest.text
                        if second_nest.tag == 'formula':
                            res['chemical_formula'] = second_nest.attrib['concise'].replace(' ', '')
        
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
        url = f'https://websvc.biocyc.org/getxml?id=META:{self.compound_id}'
        xml = self.get_with_fetcher(url,original_response=True)
        if not xml: return None
        parsed_info= self.parse_biocyc_xml(xml)
        if not parsed_info: return None
        if 'selenium.webdriver' in repr(xml): xml.quit()
        if number_of_nones_dict(parsed_info)==len(parsed_info): return None
        return Compound(parsed_info)

    def converge_compound_global(self):
        self.converge_compound_to_reaction()


    #RP with enzyme ec
    def converge_compound_to_reaction(self):
        self.convergence_args['reactions']=self.get_reactions()

        if self.convergence_args['reactions']:
            for reaction_id in self.convergence_args['reactions']:
                print(f'Linking from compound {self.compound_id} in {self.db} to reaction {reaction_id}')
                self.find_reaction(query_id=reaction_id)



if __name__ == '__main__':
    search= Compound_Fetcher_Biocyc('CPD-15291')
    search.get_compound().get_all_info()
