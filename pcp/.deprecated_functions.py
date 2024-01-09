    def create_initial_separations(self):
        ''' this is to be run before you start any simulations
        as we superimpose the molecules on each other in the cell this can give us strange intermolecular distances so we estimate
        '''
        elems = []
        for x in self.molecular_units:
            elems.append(
                list(
                    dict.fromkeys(x.get_chemical_symbols())
                    )
                    )
        self.elems = elems
        # for the elements within a molecule we can easily define this
        dict_of_separations = {}
        for i,elem_list in enumerate(elems):
            combinations = list(it.combinations_with_replacement(elem_list,2))
            analysis = Analysis(self.molecular_units[i])
            for combination in combinations:
                a1,a2 = combination
                try:
                    dict_of_separations['{}-{}'.format(a1,a2)] = np.min(
                        analysis.get_values(
                            analysis.get_bonds(a1,a2,unique=True)
                            )
                            )

                except:
                    if a1 == a2:
                        dict_of_separations['{}-{}'.format(a1,a2)]  = self.min_sep*2
                    else:
                        dict_of_separations['{}-{}'.format(a1,a2)]  = self.min_sep

        missing = it.chain.from_iterable(elems)
        combinations = list(it.combinations_with_replacement(missing,2))
        existing = list(dict_of_separations)
        for combination in combinations:
            a1,a2 = combination
            if not '{}-{}'.format(a1,a2) in existing:
                dict_of_separations['{}-{}'.format(a1,a2)] = self.init_sep_val
        
        #self.dict_of_separations = {k:float("{:.2f}".format(v)) for k,v in dict_of_separations.items() if not v == None}
        self.dict_of_separations = {k:self.min_sep for k,v in dict_of_separations.items() if not v == None} # temp.


    def create_initial_separations_from_seed_new(self,seed):
        all_distances = seed.get_all_distances()
        
        self.elems = list(
                    dict.fromkeys(seed.get_chemical_symbols())
                    )
        
        element_indices = {}
        for e in self.elems:
            element_indices[e] = [i for i,x in enumerate(seed) if x.symbol == e] 
        
        combinations = list(it.combinations_with_replacement(self.elems,2))
        
        dict_of_separations = {}
        for combination in combinations:
            a1,a2 = combination
            a1_ind = element_indices[a1] 
            a2_ind = element_indices[a2]
            products = it.product(a1_ind,a2_ind)
            dict_of_separations['{}-{}'.format(a1,a2)] = float("{:.2f}".format(np.min([all_distances[x[0]][x[1]] 
                                                                                       for x in products 
                                                                                       if not all_distances[x[0]][x[1]] == 0])))
        self.dict_of_separations = dict_of_separations