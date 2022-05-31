from unifuncnet.fetchers.protein_fetchers.protein_fetcher import *


class ProteinFetcherMetacyc(ProteinFetcher):
    def __init__(self, protein_id=None, extra_args={}, memory_storage=None):
        ProteinFetcher.__init__(self, protein_id=protein_id, memory_storage=memory_storage)
        self.db = 'metacyc'
        self.set_convergence_args()
        self.protein = self.get_protein_metacyc()
        self.add_protein()

    def set_convergence_args(self):
        self.convergence_args['genes'] = set()
        self.convergence_args['enzyme_ecs'] = set()
        self.convergence_args['proteins'] = set()
        self.convergence_args['complex'] = set()
        self.convergence_args['subunits'] = set()
        self.convergence_args['reactions'] = set()

    def converge_ec_to_protein(self):
        # first we get reactions
        reactions = self.fetch_metacyc_rxn_from_ec(self.protein_id)
        for reaction_id in reactions:
            reaction_info = self.fetch_metacyc_id_info(reaction_id, 'reaction')
            # intermediate ids to connect to proteins
            if 'intermediate_reaction_ids' in reaction_info:
                intermediate_reaction_ids = reaction_info['intermediate_reaction_ids']
                for int_reaction_id in intermediate_reaction_ids:
                    # each intermediate id will connect to reactions and proteins
                    fetched_info = self.fetch_metacyc_intermediate_rxn_ids(int_reaction_id)
                    # now we get the proteins and check what they are
                    if 'protein_ids' in fetched_info:
                        for protein_id in fetched_info['protein_ids']:
                            protein_info = self.fetch_metacyc_id_info(protein_id, 'protein')
                            # if they are subunits they will have complex_ids
                            if 'complex_ids' in protein_info:
                                entry_type = 'subunits'
                            # if they are complexes they will have subunits
                            elif 'subunit_ids' in protein_info:
                                entry_type = 'complex'
                            else:
                                # but they can be just normal proteins
                                entry_type = 'proteins'
                            self.convergence_args[entry_type].add(protein_id)
                    if 'reaction_ids' in fetched_info:
                        self.convergence_args['reactions'].update(fetched_info['reaction_ids'])

    def converge_non_ec_to_protein(self):
        fetched_info = self.fetch_metacyc_id_info(self.protein_id, 'protein')
        if 'complex_ids' in fetched_info:
            self.convergence_args['complex'].update(fetched_info.pop('complex_ids'))
        if 'subunit_ids' in fetched_info:
            self.convergence_args['subunits'].update(fetched_info.pop('subunit_ids'))
        if 'gene_ids' in fetched_info:
            self.convergence_args['genes'].update(fetched_info.pop('gene_ids'))

        if 'intermediate_reaction_ids' in fetched_info:
            intermediate_reaction_ids = fetched_info.pop('intermediate_reaction_ids')
            for int_reaction_id in intermediate_reaction_ids:
                int_fetched_info = self.fetch_metacyc_intermediate_rxn_ids(int_reaction_id)
                if 'reaction_ids' in int_fetched_info:
                    reactions = int_fetched_info['reaction_ids']
                    self.convergence_args['reactions'].update(reactions)
                    for reaction_id in reactions:
                        reaction_info = self.fetch_metacyc_id_info(reaction_id, 'reaction')
                        # intermediate ids to connect to proteins
                        if 'enzyme_ec' in reaction_info:
                            self.convergence_args['enzyme_ecs'].update(reaction_info['enzyme_ec'])

        return fetched_info

    def get_protein_metacyc(self):
        if not self.protein_id: return None
        if is_ec(self.protein_id):
            self.converge_ec_to_protein()
            fetched_info = {'enzyme_ec': {self.protein_id}}
        else:
            fetched_info = self.converge_non_ec_to_protein()
        fetched_info['metacyc'] = self.protein_id
        return Protein(fetched_info)

    def converge_protein_to_protein(self):
        convergence_dict = {'enzyme_ecs': 'protein_instances',
                            'proteins': 'protein_instances',
                            'complex': 'complex_instances',
                            'subunits': 'subunit_instances',
                            }

        for entry_type in convergence_dict:
            for protein_id in self.convergence_args[entry_type]:
                print(f'Linking from protein {self.protein_id} in {self.db} to {entry_type} {protein_id}')
                if entry_type in ['complex', 'subunits']:
                    convergence_search = True
                else:
                    convergence_search = False
                protein_instance = self.find_protein(query_id=protein_id, convergence_search=convergence_search)
                if protein_instance:
                    self.get_protein().set_detail(convergence_dict[entry_type], protein_instance, converged_in=self.db)

    def converge_protein_global(self):
        self.converge_protein_to_gene()
        self.converge_protein_to_reaction()

    def converge_protein_gpr(self):
        self.converge_protein_to_reaction()

    def converge_protein_rpg(self):
        # PG part of the CRPG pipeline
        self.converge_protein_to_gene()

    def converge_protein_to_reaction(self):
        for reaction_id in self.convergence_args['reactions']:
            print(f'Linking from protein {self.protein_id} in {self.db} to reaction {reaction_id}')
            reaction_instance = self.find_reaction(query_id=reaction_id)
            if reaction_instance:
                self.get_protein().set_detail('reaction_instances', reaction_instance, converged_in=self.db)

    def converge_protein_to_gene(self):
        for gene_id in self.convergence_args['genes']:
            print(f'Linking from protein {self.protein_id} in {self.db} to gene {gene_id}')
            gene_instance = self.find_gene(query_id=gene_id)
            if gene_instance:
                self.get_protein().set_detail('gene_instances', gene_instance, converged_in=self.db)


if __name__ == '__main__':
    from unifuncnet.biological_components.Protein import Protein
    from unifuncnet.fetchers.Reaction_Fetchers.Reaction_Fetcher import *

    protein = ProteinFetcherMetacyc('CYT-D-UBIOX-CPLX')
    # protein=ProteinFetcherMetacyc('CPLX-2401')
    # protein=ProteinFetcherMetacyc('MONOMER-2782')
    # protein=ProteinFetcherMetacyc('5.3.1.26')
