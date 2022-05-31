'''
Obtaining the most common compounds in all reactions from a unifuncnet output

'''
import os
import re


def get_cpd_counts(unifuncnet_tsv):
    res = {}
    digit_pattern = re.compile('\d+')
    try:
        if os.path.exists(unifuncnet_tsv):
            with open(unifuncnet_tsv) as file:
                for line in file:
                    line = line.strip('\n').split('\t')
                    wanted_info = [i for i in line if i.startswith('reaction_compounds')]
                    if wanted_info:
                        wanted_info = wanted_info[0]
                        wanted_info = wanted_info.replace('reaction_compounds:', '')
                        cpd_ids = digit_pattern.findall(wanted_info)
                        for c_id in cpd_ids:
                            if c_id not in res: res[c_id] = 0
                            res[c_id] += 1
        return res

    except:
        return {}


def get_description_cpds(unifuncnet_tsv, cpd_ids):
    if os.path.exists(unifuncnet_tsv):
        with open(unifuncnet_tsv) as file:
            for line in file:
                line = line.strip('\n').split('\t')
                internal_id = line[0].split(':')[1]
                wanted_info = [i for i in line if i.startswith('synonyms')]
                wanted_info = [i.replace('synonyms:', '') for i in wanted_info]
                if wanted_info:
                    if internal_id in cpd_ids:
                        # print(internal_id,cpd_ids[internal_id],wanted_info)
                        print('\t'.join(wanted_info))


def get_top_n_count(cpd_counts, unifuncnet_tsv, top_n=100):
    res = sorted(cpd_counts.items(), key=lambda item: item[1], reverse=True)
    res = res[:top_n]
    cpd_ids = {i[0]: i[1] for i in res}
    print(cpd_ids)
    get_description_cpds(unifuncnet_tsv, cpd_ids)


reactions_tsv = '/home/pedroq/PycharmProjects/unifuncnet_resources/universal_input/Reactions.tsv'
compounds_tsv = '/home/pedroq/PycharmProjects/unifuncnet_resources/universal_input/Compounds.tsv'
cpd_counts = get_cpd_counts(reactions_tsv)
top_n = get_top_n_count(cpd_counts, compounds_tsv)
