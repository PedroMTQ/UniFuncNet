

# External modules
import ast
import types
import difflib
import io
import os
import re
import sys
from html.parser import HTMLParser
from io import StringIO
from urllib.parse import quote_plus

from types import GeneratorType as generator
from sys import platform

# to use pubchem API you can use requests.post
# inchi='InChI=1S/C15H12N2O2/c16-15(18)17-11-7-3-1-5-9(11)13-14(19-13)10-6-2-4-8-12(10)17/h1-8,13-14H,(H2,16,18)'
# r=requests.post('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchi/cids/JSON/',data={'inchi':inchi})
# print(r.text)

SCRAPPABLE_DBS=['biocyc','kegg','hmdb']
SCRAPPABLE_DBS.extend(['uniprot', 'chemspider', 'inchi_key'])
SCRAPPABLE_DBS.extend(['uniprot', 'inchi_key'])


if platform.startswith('win'):
    SPLITTER = '\\'
else:
    SPLITTER = '/'

DRAX_FOLDER = os.path.abspath(os.path.dirname(__file__)).split(SPLITTER)[0:-2]
DRAX_FOLDER = SPLITTER.join(DRAX_FOLDER) + SPLITTER

def set_scrappable_dbs(user_databases):
    if user_databases:
        while SCRAPPABLE_DBS: SCRAPPABLE_DBS.pop()
        for i in user_databases: SCRAPPABLE_DBS.append(i)
        SCRAPPABLE_DBS.extend(['uniprot', 'chemspider', 'inchi_key'])
        SCRAPPABLE_DBS.extend(['uniprot', 'inchi_key'])
        return SCRAPPABLE_DBS


class MLStripper(HTMLParser):
    def __init__(self):
        super().__init__()
        self.reset()
        self.strict = False
        self.convert_charrefs = True
        self.fed = []

    def handle_data(self, d):
        self.fed.append(d)

    def get_data(self):
        return ''.join(self.fed)

def strip_tags(html):
    s = MLStripper()
    s.feed(html)
    return s.get_data()

def remove_inchi_key_equal(inchi_key):
    search_inchi=re.search('=',inchi_key)
    if search_inchi:
        inchi_key=inchi_key[search_inchi.span()[1]:]
    return inchi_key

def xstr(s):
    if not s:
        return ''
    if isinstance(s,generator): s=list(s)
    if isinstance(s,list):
        return ' , '.join(s)
    return str(s)

def is_ec(enz_id,required_level=3):
    if enz_id:
        ec_pattern = re.compile('^\d+\.\d+\.\d+(\.(-|\d+|([a-zA-Z]\d+)|[a-zA-Z]))?')
        if re.search(ec_pattern, enz_id):
            enz_id_copy=str(enz_id).replace('.-','')
            if len(enz_id_copy.split('.'))>=required_level:
                return True
    return False

def clean_ec(original_ec ):
    temp_ec = str(original_ec)
    if '-' in temp_ec[-1]:
        temp_ec = temp_ec[:-2]
    return temp_ec

def find_ecs(string_to_search, required_level=3):
    res = set()
    # greedy match of confounders
    ec_pattern = re.compile('\d(\.(-|\d{1,3}|([a-zA-Z]\d{1,3}))){2,3}')
    search = re.finditer(ec_pattern, string_to_search)
    for i in search:
        ec = i.group()
        passed = False
        start = i.span()[0]
        end = i.span()[1]
        if len(string_to_search) > end:
            if string_to_search[start - 1] != '.' and\
                    string_to_search[end] != '.' \
                    and not re.match('\.|[a-zA-Z]|\d{1,3}', string_to_search[end]) and not re.match('-', string_to_search[end]):
                passed = True
        else:
            if string_to_search[start - 1] != '.':
                passed = True
        if passed:
            if ec.count('.') >= required_level - 1:
                if ec.count('.') + 1 - ec.count('-') >= required_level:
                    res.add(ec)
    return res


def join_or_not(res,sep=' ; ',nones='N_A',include_dict_numbers=True):
    if isinstance(res,generator): res=list(res)
    if res and isinstance(res, list):
        temp = []  # remove equal values - happens in organisms
        for i in res:
            if i not in temp: temp.append(i)
        return sep.join([i for i in temp if i != 'Not Available'])
    elif res and isinstance(res, dict):
        if include_dict_numbers:
            return sep.join([i + ': ' + str(res[i]) for i in res])
        else:
            return sep.join([i for i in res])
    elif not res:
        return nones
    else:
        return res




def sub_prod_to_reaction(substrates, products, arrow=' → '):
    res = ''
    for i in substrates:
        if not i: return None
        res += i + ' + '
    res = res[0:-3]
    res += arrow
    for i in products:
        if not i: return None
        res += i + ' + '
    res = res[0:-3]
    return res




def find_path(to_search_name,to_search='file',end_dir='DRAX',parent_directory=None):
    past_path=str(os.getcwd())
    current_dir=os.getcwd()
    current_dir_limit=3
    c=0
    while end_dir not in str(os.getcwd()).split(SPLITTER)[-1] and c<= current_dir_limit:
        if current_dir==os.getcwd():
            c+=1
        else:
            current_dir = os.getcwd()
            c = 0
        os.chdir("..")
    print('changed dir',to_search_name)
    for root, dirs, files in os.walk(os.getcwd()):
        parent_dir=root.split(SPLITTER)[-1]
        if to_search=='directory': target=dirs
        else: target=files
        if to_search_name in target:
            if parent_directory:
                if parent_dir==parent_directory:
                    os.chdir(past_path)
                    return os.path.join(root, to_search_name)
            else:
                os.chdir(past_path)
                return os.path.join(root, to_search_name)

def split_by_plus(str_to_split,split_str):
    return [i.strip() for i in str_to_split.split(split_str)]

# this allows the inclusion of stoichiometry in the reaction with ids and instances
def get_stoichiometry_reaction(reaction,ignore_numbers=False,split_plus_str=' + '):
    find = find_sign(reaction)
    if not find: return None,None
    sub = reaction[0:find.span()[0]]
    prod = reaction[find.span()[1]:]
    l_sub = split_by_plus(sub,split_plus_str)
    l_prod = split_by_plus(prod,split_plus_str)
    complete_l = l_sub + l_prod
    for i in range(len(complete_l)):
        #we use [b-z] since some reactions have an undefined amount of compounds, they usually use n but
        #we dont use 'a' because 'a compound' refers to a unit
        numbers = re.match('(\d+n)|(\d+ )|(a )|(an )|([b-z]\s)|(\(\d?n[\+\-]?\d?\))', complete_l[i].lower().strip())
        if numbers and not ignore_numbers:
            to_add= complete_l[i][numbers.span()[1]:].strip()
            n=numbers.group().strip()
            if re.match('(a\s?)|(an\s?)',n.lower()):
                n=1
            #undefined (e.g. n molecules of something) numbers will be -1
            elif re.match('(\(\d?n[\+\-]?\d?\))|([b-z]\s?)|(\d+n)',n.lower()):
                #print('POSSIBLE UNDEFINED NUMBER OF COMPOUNDS, CHECK THIS TO BE SURE NO ERRORS ARE OCCURRING',reaction)
                n=-1
            else: n=int(n)
            complete_l[i] = [int(n), to_add]
        else:
            complete_l[i] = [1, complete_l[i]]
    return complete_l,len(l_sub)

# this allows the inclusion of stoichiometry in the reaction with ids and instances

def get_stoichiometry(reaction,reaction_with_ids,ignore_numbers=False):

    if not reaction and not reaction_with_ids: return None
    if not reaction: reaction=reaction_with_ids
    if not reaction_with_ids: reaction_with_ids=reaction
    complete_l,len_sub=get_stoichiometry_reaction(reaction,ignore_numbers)
    #when the reaction is incomplete (either no products or no substrates, we won't get stoichiometry
    if not complete_l: return None
    already_split=False
    if not isinstance(reaction_with_ids,list):
        complete_l_id,len_sub=get_stoichiometry_reaction(reaction_with_ids,ignore_numbers)
    else:
        already_split=True
        complete_l_id=list(reaction_with_ids)
        if len(complete_l)!=len(complete_l_id):
            temp_complete_l, temp_len_sub = get_stoichiometry_reaction(reaction, ignore_numbers,split_plus_str=' +')
            if len(temp_complete_l)!=len(complete_l_id):
                temp_complete_l, temp_len_sub = get_stoichiometry_reaction(reaction, ignore_numbers,split_plus_str='+ ')
            if len(temp_complete_l) == len(complete_l_id):
                complete_l=temp_complete_l
                len_sub=temp_len_sub
    res=[]
    c=0
    unite=True
    #depending on the reaction html code, it may have a different length between its html IDS list and the reaction string
    if len(complete_l)==len(complete_l_id): subtractor=1
    else: subtractor=0
    if already_split:
        for combo in complete_l:
            stoi=combo[0]
            #biocyc sometimes has the same ID twice (instead of merging it into one and adding a multiplier)
            while c<len(complete_l_id)-1:
                # we need to take this into account when defining whether we should unite or not - this avoids uniting when the last substrate is the same as the first product
                if complete_l_id[c+1]==complete_l_id[c] and unite:
                    if c == len_sub-subtractor:  unite = False
                    if unite:                    c+=1
                else:           break
            res.append([stoi,complete_l_id[c]])
            c+=1
            unite=True
    else: res=complete_l_id
    for i in range(len(res)):
        res[i][0]=complete_l[i][0]
    return res,complete_l,len_sub



def standardize_reaction_str(reaction_str):
    sign=find_sign(reaction_str).group()
    uniformed_sign=uniform_sign(sign)
    res = reaction_str.replace(sign,uniformed_sign)
    return res

def find_sign(string_to_search):
    #first search should be for one of these signs, which means directionality is set
    compiled_sign = re.compile('( = | <=> | <= | => | → | ← | ↔ | ⇄ | <-> | <--> | <- | -> | [=-]?&gt; | &lt;=&gt; | &lt;[=-]? )')
    sign = re.search(compiled_sign, string_to_search)
    if sign:
        return sign
    else:
        #the second search should be for unknown directionality
        #this way it avoids recognizing unknown compounds "?" as the direction sign of a reaciton
        interrogation_sign=re.compile('\?')
        sign = re.search(interrogation_sign, string_to_search)
        if sign: return sign
    return None

def uniform_sign(sign):
    both = re.compile('( = | <=> | ↔ | ⇄ | <-> | <--> | &lt;=&gt; )')
    right = re.compile('( => | → | -> | [=-]?&gt; )')
    left = re.compile('( <= | ← | <- | &lt;[=-]? )')
    if re.search(both,sign):  return ' <=> '
    if re.search(right,sign):  return ' => '
    if re.search(left,sign):  return ' <= '

def fix_html_sign(string_to_fix):
    res=str(string_to_fix)
    if re.search('[=-]?&gt;',res):
        res=res.replace(re.search('[=-]?&gt;',res).group(),'=>')
    if re.search('&lt;[=-]?',res):
        res=res.replace(re.search('&lt;[=-]?',res).group(),'=>')
    if re.search('&lt;[=-]?',res):
        res=string_to_fix.replace(re.search('&lt;=&gt;',res).group(),'<=>')
    if re.search('%2b',res):
        res=string_to_fix.replace(re.search('%2b',res).group(),'%2b')
    return res
    
    
def get_id_or_str(rn_cpd):
    if not isinstance(rn_cpd,str): return rn_cpd.get_most_common_synonym()
    else: return rn_cpd





# transforms list of compound ids into a string with the stoichiometry and the compound ids`
def l_rn_ids_to_str(rn, rn_ids,without_stoichiometry=False):
    res = ''
    if not rn_ids: return 'N_A'
    rn_ids = any_nones_to_na(rn_ids)
    sign = find_sign(rn).group()
    temp = rn.split(sign)
    len_sub, len_prod = len(temp[0].split(' + ')), len(temp[1].split(' + '))
    for i in range(len_sub):
        if without_stoichiometry:
            res += get_id_or_str(rn_ids[i][1])
        else:
            res += str(rn_ids[i][0]) + ' ' + get_id_or_str(rn_ids[i][1])
        if i < len_sub - 1: res += ' + '
    sign=uniform_sign(sign).strip()
    res += ' ' + sign + ' '
    for j in range(len_sub, len_sub + len_prod):
        if without_stoichiometry:
            res += get_id_or_str(rn_ids[j][1])
        else:
            res += str(rn_ids[j][0]) + ' ' + get_id_or_str(rn_ids[j][1])
        if j < len_prod + len_sub - 1: res += ' + '
    return res



def any_nones_to_na(list_to_convert):
    for i in range(len(list_to_convert)):
        if not list_to_convert[i][0]: list_to_convert[i][0] = '1'
        if not list_to_convert[i][1]: list_to_convert[i][1] = 'N_A'
    return list_to_convert


#some reactions have an undefined extra substrates
def match_plus(separated_reaction,sign):
    expression=re.compile('\s\+\s+{}'.format(sign))
    search=re.search(expression,separated_reaction)
    if search:
        splitter0=search.span()[0]
        sub=separated_reaction[:splitter0]
        splitter1=search.span()[1]
        prod=separated_reaction[splitter1-2:]
        return sub+prod
    return separated_reaction


def test_match_possible_ids(possible_ids_1,possible_ids_2):
    if not possible_ids_1 or not possible_ids_2: return False
    if isinstance(possible_ids_1,str) or isinstance(possible_ids_1,int):
        temp_possible_ids_1= {possible_ids_1:1}
    elif isinstance(possible_ids_1,dict):
        temp_possible_ids_1=set(possible_ids_1.keys())
    else: temp_possible_ids_1=possible_ids_1


    if isinstance(possible_ids_2,str) or isinstance(possible_ids_2,int):
        temp_possible_ids_2= {possible_ids_2:1}
    elif isinstance(possible_ids_2,dict):
        temp_possible_ids_2=set(possible_ids_2.keys())
    else: temp_possible_ids_2=possible_ids_2

    for i in temp_possible_ids_1:
        for j in temp_possible_ids_2:
            if i and j:
                if i == j:
                    return True
    return False


def score_match_possible_ids(possible_ids_1,possible_ids_2,match_c=1,mismatch_c=-0.25):
    """
    just a scoring system. if one instance doesnt have an id we dont penalize
    if the ids match we benefit the match, if they dont we add a minor penalty.
    The penalty is lower since we assume databases already have some curation, thus the amount of errors shouldnt be too high
    """

    if not possible_ids_1 or not possible_ids_2: return 0
    if isinstance(possible_ids_1,set) or isinstance(possible_ids_1,list):
        temp_possible_ids_1= {i:1 for i in possible_ids_1}
    elif not isinstance(possible_ids_1,dict):
        temp_possible_ids_1= {possible_ids_1:1}
    else:
        temp_possible_ids_1=dict(possible_ids_1)


    if isinstance(possible_ids_2,set) or isinstance(possible_ids_2,list):
        temp_possible_ids_2= {i:1 for i in possible_ids_2}
    elif not isinstance(possible_ids_2, dict):
        temp_possible_ids_2 = {possible_ids_2: 1}
    else:
        temp_possible_ids_2=dict(possible_ids_2)
    for i in temp_possible_ids_1:
        for j in temp_possible_ids_2:
            if i and j:
                if i == j:
                    if not is_ec(i):
                        return match_c
    return mismatch_c

def unite_possible_ids(self_instance,instance_2,detail_type):
    if instance_2.get_detail(detail_type):
        possible_ids = instance_2.get_detail(detail_type,all_possible=True)
        if possible_ids:
            for possible_id in possible_ids:
                possible_id_count = possible_ids[possible_id]
                self_instance.set_detail(detail_type, {possible_id: possible_id_count})


def number_of_nones_dict(test_dict):
    res=0
    for i in test_dict:
        if not test_dict[i]: res+=1
    return res


def check_if_any_none(list_to_check,pos=None):
    if not list_to_check: return True
    for i in list_to_check:
        if pos:
            if not i[pos]: return True
        else:
            if not i : return True
    return False

#only used when we invoke a specific detail type
def list_has_common_items(iterable1,iterable2):
    if iterable1 and iterable2:
        for i in iterable1:
            for j in iterable2:
                if i.lower()==j.lower(): return True
    return False

def unite_instance_list(list_of_instances):
    for p1 in list_of_instances:
        for p2 in list_of_instances:
            if p1 is not p2: p1.unite_instances(p2,always_unite=True)
    if list_of_instances:
        return list_of_instances[0]
    else: return list_of_instances

def get_instance_type(instance_to_test):
    res=str(type(instance_to_test))
    res=res.split('.')[-1].replace('\'>','')
    return res

def check_match_in_list(list_to_check,instance_to_check):
    for i in list_to_check:
        if instance_to_check.is_match_instances(i): return True
    return False






def regex_escape(regex_string):
    regex_characters=[
        "\\",
        ".",
        "+",
        "*",
        "?",
        "[",
        "^",
        "]",
        "$",
        "(",
        ")",
        "{",
        "}",
        "=",
        "!",
        "<",
        ">",
        "|",
        "\'",
        "\"",
        ":",
        "-"]
    l_regex_string=[i for i in regex_string]
    for i in range(len(l_regex_string)):
        if l_regex_string[i] in regex_characters:
            if l_regex_string[i]=='\\': replacer='\\\\'
            else: replacer='\\'
            l_regex_string[i]=replacer+l_regex_string[i]
    return ''.join(l_regex_string)



if __name__ == '__main__':
    rn='Acetyl-CoA + n Malonyl-CoA + 2n NADPH + 2n H+ <=> Long-chain fatty acid + n CO2 + 2n NADP+ + (n+1) CoA + n H2O'
    rn_ids='C00024 + n C00083 + 2n C00005 + 2n C00080 <=> C00638 + n C00011 + 2n C00006 + (n+1) C00010 + n C00001'
    get_stoichiometry(rn,rn_ids)
