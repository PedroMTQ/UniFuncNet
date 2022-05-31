from types import GeneratorType as generator


class Counter:
    def __init__(self, string):
        if isinstance(string, list) or isinstance(string, generator) or isinstance(string, set):
            self.possible_strings = {db_string.strip(): 1 for db_string in string}
        elif isinstance(string, dict):
            self.possible_strings = {db_string.strip(): int(string[db_string]) for db_string in string}
        else:
            self.possible_strings = {string.strip(): 1}

    def __str__(self):
        return str(self.possible_strings)

    def get_most_common_string(self):
        if self.get_strings():
            try:
                return next(self.get_strings())
            except:
                return None
        else:
            return None

    def remove_string(self, string_to_remove):
        if string_to_remove in self.possible_strings:
            self.possible_strings.pop(string_to_remove)

    def get_possible_strings(self):
        if None in self.possible_strings and len(self.possible_strings) == 1:
            return {}
        return dict(self.possible_strings)

    def get_strings(self):
        res = list(sorted(self.get_possible_strings(), key=self.get_possible_strings().__getitem__, reverse=True))
        for r in res:
            yield r

    def append(self, new_string, count=1):
        if isinstance(count, str): count = int(count)
        if new_string:
            if isinstance(new_string, str): new_string = [new_string]
            for nstring in new_string:
                if nstring:
                    nstring = nstring.strip()
                    if nstring in self.possible_strings:
                        self.possible_strings[nstring] += count
                    else:
                        self.possible_strings[nstring] = count
            if None in self.possible_strings.keys() and len(self.possible_strings) > 1: self.possible_strings.pop(None)
