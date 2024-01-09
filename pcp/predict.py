from smilesbox.smilesbox import SMILESbox
from airsspy import SeedAtoms, Buildcell
from chgnet.model.model import CHGNet
from chgnet.model import StructOptimizer

# Import ase assets
from ase.optimize import BFGSLineSearch, sciopt, FIRE, BFGS, precon
from ase.calculators.lj import LennardJones
from ase.constraints import UnitCellFilter, ExpCellFilter
from ase.io import read,write

import copy 
import itertools as it 
from monty.serialization import loadfn,dumpfn

from ase.geometry.analysis import Analysis
import numpy as np 
import re, os, sys,io 
from tqdm import tqdm
from datetime import  datetime as dt
from pymatgen.io.ase import AseAtomsAdaptor
import pandas as pd 
from pymatgen.analysis.structure_prediction.volume_predictor import DLSVolumePredictor

import warnings 
warnings.filterwarnings("ignore", module="pymatgen")
warnings.filterwarnings("ignore", module="ase")

import pathos.multiprocessing as multiprocessing 
from pathos.helpers import mp as pmp 
from pathos.pools import ProcessPool 
from pathos.pp import ParallelPool 
from pathos.serial import SerialPool
import math


class PredictStructure:
    def __init__(self,
                 boxsize=10,
                 init_sep_val=0.5,
                 num_units=1,
                 units=True,
                 min_sep = 0.5,
                 max_failures = 10,
                 nprocs=2,
                 use_device='cpu'):
        
        self.boxsize = boxsize
        self.init_sep_val = init_sep_val
        self.num_units = num_units
        self.min_sep = float(min_sep)
        self.units = units
        self.data = {}
        self.max_failures = max_failures
        self.energies = []
        self.max_forces = []
        self.convergence = []
        self.nprocs = nprocs
        self.use_device = use_device # cpu, cuda, mps
        

    def show(self,seed):
        from pyxtal import pyxtal
        cryst = pyxtal()
        AseAtomsAdaptor.get_structure(seed)
        cryst._from_pymatgen(AseAtomsAdaptor.get_structure(seed))
        return(cryst)
    
    def create_initial_seed(self,smileses=None,dls=True,**kws):
        '''
        need a better way of doing the box size as it makes the airss configurations more sensitive
        '''

        molecular_units = []
        for mol in smileses:
            sb = SMILESbox()
            sb.smiles_to_atoms(mol)
            if not self.boxsize == None:
                sb.add_box([self.boxsize,self.boxsize,self.boxsize])
            molecule = sb.atoms
            molecular_units.append(molecule)
        self.molecular_units = molecular_units
        initial_seed = copy.deepcopy(molecular_units)
        for i,j in enumerate(initial_seed):
            if not i == 0:
                initial_seed[0]+=j
        if dls == True and not self.boxsize==None:
            self.initial_seed = AseAtomsAdaptor().get_atoms(
                DLSVolumePredictor(
                    cutoff=self.boxsize+5,**kws).get_predicted_structure(
                        AseAtomsAdaptor().get_structure(initial_seed[0])
            ))
            self.seed = copy.deepcopy(self.initial_seed)


        else:
            self.initial_seed = initial_seed[0]
            self.seed = copy.deepcopy(self.initial_seed)


    def create_initial_separations_from_seed(self,seed): # this needs changing
        '''seed must be Atoms object'''
        
        # lowest_distance = np.min(seed.get_cell_lengths_and_angles()[0:3])

        self.elems = list(
            dict.fromkeys(seed.get_chemical_symbols())
            )
        dict_of_separations = {}
        combinations = list(it.combinations_with_replacement(self.elems,2))
        analysis = Analysis(seed)
        for combination in combinations:
            a1,a2 = combination
            try:
                dict_of_separations['{}-{}'.format(a1,a2)] = np.min(
                    analysis.get_values(
                        analysis.get_bonds(a1,a2,unique=True)
                        )
                        )

            except:
                dict_of_separations['{}-{}'.format(a1,a2)]  = None #lowest_distance

        self.dict_of_separations = {k:float("{:.2f}".format(v)) for k,v in dict_of_separations.items() if not v == None}
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
            dict_of_separations['{}-{}'.format(a1,a2)] = [float("{:.2f}".format(np.min(all_distances[x[0]][x[1]]))) 
                                                             for x in products 
                                                             if not all_distances[x[0]][x[1]] == 0]
        self.dict_of_separations = dict_of_separations
    def generate_airss_input(self,
                             targvol = None,
                             system = None,
                             ): # add some more
        '''
        this generates a SeedAtoms atoms object which can then be used to generate random configurations
        '''

        self.seed = SeedAtoms(self.seed)
        try:
            self.seed.gentags.minsep = [self.min_sep,self.dict_of_separations]
        except:
            self.seed.gentags.minsep = self.min_sep
        self.seed.gentags.targvol = targvol
        self.seed.gentags.system = system
        # need to add more gentags

        # for now this is highly specialised towards my own system:
        for i,mol in enumerate(self.molecular_units):
            elems = list(
                dict.fromkeys(mol.get_chemical_symbols())
                )
            for x in self.seed:
                elem = re.split('(\d+)', x.tagname)
                if any(y for y in elem if y in elems):
                    if self.units == True:
                        x.tagname = '{}-mol'.format(i+1)
                    x.num = self.num_units

        self.airrs_input_file = '\n'.join(self.seed.get_cell_inp_lines()) # incase you need an input file

    def generate_random_cells(self,num_cells,**kws):

        try:
            random_cells = [self.seed.build_random_atoms(**kws) for x in tqdm(range(num_cells),
                                                                              desc='Building Randomised Cells',
                                                                              leave=True,
                                                                              bar_format='{l_bar}{bar:10}{r_bar}')]
        except:
            random_cells = None
            print(' -> unable to build random atoms from seed.')
        return(random_cells)
    
    def gen_and_relax_mp(self,i_chunk,chunk,relaxer,directory,steps,out_q=None):
        _data = {}
        for i,c in enumerate(chunk):
            #i+=1
            structure_number = i_chunk+i+len(chunk)*i_chunk
            destination = os.path.join('.',directory,'structure-{}'.format(structure_number))
            os.makedirs(destination,exist_ok=True)
            try:
                write('{}/input.vasp'.format(destination),c,vasp5=True,sort=True)
                #logfile = open(sys.stdout.fileno(),'wb',0)
                sys.stdout = io.TextIOWrapper(open(os.path.join(destination,'out.log'),'wb',0),write_through=True) #logging the output 
                print(self.airrs_input_file)
                print('\n---\n')
                result = relaxer.relax(c,steps=steps,verbose=True)
                fmax = np.max(
                    [np.linalg.norm(x) for x in result['trajectory'].forces[-1]]
                    )
                final_energy = result['trajectory'].energies[-1]
                final_structure = result['final_structure']
                final_structure.sort()
                final_structure.to(filename='{}/output.vasp'.format(destination),fmt='poscar')
                _data[structure_number] = {
                    'final_structure':final_structure,
                    'final_energy':final_energy,
                    'max_force':fmax}
            except:
                pass
        out_q.put(_data)

    def _mp_function(self,run_num,random_cells,relaxer,steps,dls=False,):
        '''this is the multiprocessing function'''

        data_chunks = [random_cells[chunksize*i:chunksize*(i+1)]
                            for i in range(self.nprocs)
                            for chunksize in [int(math.ceil(len(random_cells)/float(self.nprocs)))]]
        
        manager = pmp.Manager()
        out_queue = manager.Queue()
        #run = 0        
        data = {}
        jobs = []
        directory = os.path.join('runs','run_{}'.format(run_num))
        for i,chunk in enumerate(data_chunks):
            process = pmp.Process(target=self.gen_and_relax_mp,
                                  args=(i,chunk,relaxer,directory,steps,out_queue))
            jobs.append(process)
            process.start()

        for proc in jobs:
            data.update(out_queue.get())
        
        for proc in jobs:
            proc.terminate()

        for proc in jobs:
            proc.join()
        
        return(data)
        
    
    @staticmethod
    def check_convergence(points,n_points,energy_convergence):
        '''checks if the gradient of the line is under convergence criterion'''

        points_considered = list(reversed(points))[0:n_points]
        xrange = range(len(points_considered))
        a, b = np.polyfit(xrange, points_considered, 1)

        if abs(a) <= energy_convergence:
            if points_considered[-1] == np.min(points):
                return(1,abs(a))
            else:
                return(0,abs(a))
        else:
            return(0,abs(a))
        
        
    def run_seeds(self,
                  num_seeds=100,
                  optimizer_class='FIRE',
                  energy_convergence=0.1,
                  num_points_to_converge=3,
                  dir='.',
                  dls=False,
                  steps=100,
                  **kws):
        
        chgnetrelaxer = StructOptimizer(optimizer_class=optimizer_class,use_device=self.use_device)
        #initially

        run = 0 # should rename to "generation"

        self.energy_convergence = energy_convergence
        self.generate_airss_input() # needs more options
        print('\nGeneration {}:'.format(run))
        random_atoms = self.generate_random_cells(num_cells=num_seeds) # add kws
        
        start = dt.now()
        data = self._mp_function(run,random_atoms,chgnetrelaxer,steps=steps,dls=dls,)
        total = dt.now() - start
        

        df = pd.DataFrame(data).T.sort_values(by='final_energy',ascending=True)
        df.reset_index(inplace=True)

        self.data[run] = copy.deepcopy(df)
        self.energies.append(df.T[0]['final_energy'])
        self.max_forces.append(df.T[0]['max_force'])
        self.seed = SeedAtoms(AseAtomsAdaptor().get_atoms(df.T[0]['final_structure']))
        self.seed.write('{}/run_{}.vasp'.format(dir,run))

        print('''   -> victor:  structure-{}, 
              energy: {:.2F},
              fmax: {:.2F},
              time: {}s'''.format(
                  self.data[run].T[0]['index'],
                  self.energies[-1],
                  self.max_forces[-1],
                  total.total_seconds()
                  )
                  )

        #now to converge the energies by looping airss
        convergence = 0
        while convergence == 0: # perhaps this should be a function
            run+=1
            print('\nGeneration {}:'.format(run))

            data = {}
            #self.create_initial_separations_from_seed(self.seed)
            self.create_initial_separations_from_seed_new(self.seed)
            self.num_units = 1 # to avoid exponentially increasing the structure
            self.generate_airss_input() # need to have some kws
            random_atoms = self.generate_random_cells(num_cells=num_seeds) # add kws

            start = dt.now()
            data = self._mp_function(run,random_atoms,chgnetrelaxer,steps=steps,dls=dls)
            total = dt.now() - start

            df = pd.DataFrame(data).T.sort_values(by='final_energy',ascending=True)
            df.reset_index(inplace=True)

            self.data[run] = copy.deepcopy(df)

            if df.T[0]['final_energy'] <= np.min(self.energies):
                self.seed = AseAtomsAdaptor().get_atoms(df.T[0]['final_structure']) 
                self.seed.write('{}/run_{}.vasp'.format(dir,run))
            else:
                self.seed = AseAtomsAdaptor().get_atoms(df.T[0]['final_structure']) 
                self.seed.write('{}/run_{}.vasp'.format(dir,run))
                self.seed = AseAtomsAdaptor().get_atoms(self.data[np.argmin(self.energies)].T[0]['final_structure'])
                print('warning: lower energy seed was found previously - restarting from that')
                
            self.energies.append(df.T[0]['final_energy'])
            self.max_forces.append(df.T[0]['max_force'])
            
            if run >= num_points_to_converge:

                convergence,absey = self.check_convergence(points=self.energies,
                                                           n_points=num_points_to_converge,
                                                           energy_convergence=energy_convergence)
                
                self.convergence.append(absey)

                print('''   -> victor:  structure-{}, 
                      energy: {:.2F},
                      convergence:{:.2F} 
                      fmax: {:.2F},
                      time: {:.2F}s'''.format(
                          self.data[run].T[0]['index'],
                          self.energies[-1],
                          absey,
                          self.max_forces[-1],
                          total.total_seconds()
                          )
                          )
                
            else:
                self.convergence.append(None)
                print('''   -> victor:  structure-{}, 
                      energy: {:.2F},
                      fmax: {:.2F},
                      time: {:.2F}s'''.format(
                          self.data[run].T[0]['index'],
                          self.energies[-1],
                          self.max_forces[-1],
                          total.total_seconds()
                          )
                          )

    def plot_trajectory(self):
        import matplotlib.pyplot as plt 
        fig,(ax1,ax2,ax3) = plt.subplots(ncols=3,figsize=(18,6),dpi=100)
        pd.Series(self.max_forces).plot.line(ax=ax1,color='tab:blue')
        pd.Series(self.energies).plot.line(ax=ax2,color='tab:green')
        pd.Series(self.convergence).plot.line(ax=ax3,color='tab:red')
        ax3.axhline(self.energy_convergence,linestyle='--',color='black')
        ax3.set_xlim(0)
        ax1.set_title('max force')
        ax2.set_title('total energy')
        ax3.set_title('convergence')

        ax1.set_ylabel('max force (eV $\\AA^{-1}$)')
        ax2.set_ylabel('energy (eV)')
        ax3.set_ylabel('convergence (arb.units)')
        [ax.set_xlabel('run') for ax in [ax1,ax2,ax3]]
        plt.tight_layout()
        plt.savefig('test.pdf')    
