from source.Utils.util import get_instance_type,SCRAPPABLE_DBS

###############GENERAL FUNCTION###############
def generate_details_dict(bio_instance):
    return {i: None for i in get_details_list(bio_instance)}


def get_details_list(bio_instances,remove_from_list=[],append_to_list=[]):
    details=set()
    if not isinstance(bio_instances,list):
        bio_instances=[bio_instances]
    for bio_inst in bio_instances:
        vars_instance= vars(bio_inst)
        for v in vars_instance:
            if v=='identifiers':
                for i in vars_instance[v]:details.add(i.replace('_id',''))
            elif v=='instances':
                for i in vars_instance[v]:details.add(i)
            else:
                if v not in [
                            'verbose',
                            'need_to_save',
                             ]:
                    details.add(v)
    for v in remove_from_list:
        if v in details:
            details.remove(v)
    for v in append_to_list:
        details.add(v)

    return details


def get_unique_details(bio_instance=None,bio_instance_type=None,remove_from_list=[],append_to_list=[]):
    res=[]
    if type(remove_from_list)==str:remove_from_list=[remove_from_list]
    if type(append_to_list)==str:append_to_list=[append_to_list]
    if bio_instance:
        bio_instance_type=get_instance_type(bio_instance)
    if bio_instance_type=='Compound':
        res= get_unique_details_compound()
    elif bio_instance_type=='Gene':
        res= get_unique_details_gene()
    elif bio_instance_type=='Protein':
        res= get_unique_details_protein()
    elif bio_instance_type=='Reaction':
        res= get_unique_details_reaction()
    res.extend(SCRAPPABLE_DBS)
    for detail in append_to_list:
        if detail not in res:res.append(detail)
    for detail in remove_from_list:
        if detail in res:
            res.remove(detail)
    return res



def get_unique_details_compound():
    res=[
        'chebi',
        'chemspider',
        'pubchem_cid',
        'pubchem_sid',
        'inchi_key',
        'inchi',
                       ]
    return res

def get_unique_details_gene():
    res=[
        'uniprot',
        'uniprot_name',
               ]
    return res


def get_unique_details_protein():
    res=[
         ]
    return res



def get_unique_details_reaction():
    res = [
            'reaction_str',
            'reaction_with_instances',
           ]
    return res

