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