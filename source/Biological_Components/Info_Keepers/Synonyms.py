
import re
from types import GeneratorType as generator

#This class saves synonyms of a specific compound, useful to know which synonyms are the most common.
# This way, matches start by the most common synonyms, instead of less common synonyms.

class Synonyms():
    def __init__(self, synonyms,bio_type):
        self.bio_type=bio_type
        if synonyms:
            if isinstance(synonyms, list) or isinstance(synonyms,generator) or isinstance(synonyms,set):
                self.possible_synonyms = {self.fix_synonym(syn):1 for syn in synonyms}
            elif isinstance(synonyms, dict):
                self.possible_synonyms = {self.fix_synonym(syn):int(synonyms[syn]) for syn in synonyms}
            elif isinstance(synonyms,str):
                self.possible_synonyms={self.fix_synonym(synonyms):1}
        else: self.possible_synonyms={}

    def __str__(self):
        return    str(self.possible_synonyms)

    def get_most_common_synonym(self):
        try:    return next(self.get_synonyms())
        except: return ''

    def remove_synonym(self,syn):
        if syn in self.possible_synonyms:
            self.possible_synonyms.pop(syn)

    def get_synonyms(self):
        res= list(sorted(self.get_possible_synonyms(), key=self.get_possible_synonyms().__getitem__, reverse=True))
        for r in res: yield r

    def get_possible_synonyms(self):
        return dict(self.possible_synonyms)

    def append(self,synonyms,count=1):
        if synonyms:
            if isinstance(synonyms, list) or isinstance(synonyms,generator) or isinstance(synonyms,set):
                [self.inner_append(i,count) for i in synonyms]
            elif isinstance(synonyms, dict):
                [self.inner_append(i,synonyms[i]) for i in synonyms]
            elif isinstance(synonyms,str):
                self.inner_append(synonyms,count)

    def inner_append(self, new_syn,count=1):
        if not isinstance(count,int): count=int(count)
        new_syn=self.fix_synonym(new_syn)
        if new_syn:
            if new_syn in self.possible_synonyms:
                self.possible_synonyms[new_syn] += count
            else:
                self.possible_synonyms[new_syn] = count
        if None in self.possible_synonyms.keys() and len(self.possible_synonyms)>1: self.possible_synonyms.pop(None)

    def extend(self,list_synonyms):
        for syn in list_synonyms:
            self.append(self.fix_synonym(syn))

    def fix_synonym(self,synonym):
        if not synonym: return synonym
        res=synonym
        if   self.bio_type=='Compound':
            res=self.fix_synonym_compound(synonym)
        elif self.bio_type=='Gene':
            res=self.fix_synonym_gene(synonym)
        elif self.bio_type=='Protein':
            res=self.fix_synonym_protein(synonym)
        return res


    def fix_synonym_generic(self,synonym):
        syn=synonym.lower()
        syn=syn.replace('\n','')
        syn=syn.strip()
        for s in ['(non-preferred name)',
                  '(brand name)'
                  '(ambiguous)'
                  ]:
            if s in syn:
                syn=syn.replace('(ambiguous)','')
                syn=syn.strip()

        return syn

    def fix_synonym_compound(self,synonym):
        syn=self.fix_synonym_generic(synonym)
        match = re.match('a[n]?\s',synonym)
        if match:
            syn=syn[match.span()[1]:]
        return syn

    def fix_synonym_gene(self,synonym):
        syn=self.fix_synonym_generic(synonym)
        return syn

    def fix_synonym_protein(self,synonym):
        syn=self.fix_synonym_generic(synonym)
        return syn


if __name__ == '__main__':
    pass