
import datetime
from types import GeneratorType as generator


class Item_Set():
    def __init__(self,item_or_list=None,converged_in=None):
        self.set_of_items=set()
        self.convergence={}
        self.set_item_set(item_or_list,converged_in)

    def __str__(self):
        return str(i for i in self.get_item_set())

    def set_item_set(self,item_or_list,converged_in=None):
        if converged_in:
            if isinstance(converged_in,str):
                converged_in=[converged_in]
        if item_or_list:
            if not (isinstance(item_or_list, list) and not isinstance(item_or_list,generator)
                    and not isinstance(item_or_list,dict)) and not isinstance(item_or_list,set):
                item_or_list = [item_or_list]
            for i in item_or_list:
                if i and i not in self.set_of_items:  self.set_of_items.add(i)
                if i and converged_in:
                    for db in converged_in:
                        self.add_convergence_properties(i,db)
        if None in self.set_of_items: self.set_of_items.remove(None)

    def remove_item(self,item_to_remove):
        if item_to_remove in self.set_of_items:
            self.set_of_items.remove(item_to_remove)

    def clean_item_set(self):
        self.set_of_items=set()

    def get_item_set(self):
        return self.set_of_items
        #if not self.set_of_items: return []
        #for i in self.set_of_items

    def add_convergence_properties(self,i,converged_in):
        if converged_in:
            if converged_in not in self.convergence: self.convergence[converged_in]={}
            self.convergence[converged_in][i]=str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

    def get_convergence_db(self,i):
        for db in self.convergence:
            if i in self.convergence[db]:
                yield [db,self.convergence[db][i]]

