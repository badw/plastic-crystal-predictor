from pcp.predict import PredictStructure
from datetime import datetime as dt 
import pickle

def main():
    ps = PredictStructure(
            num_units=1,
            units=True,
            boxsize=10,
            min_sep=2,
            nprocs=32)
    
    ps.create_initial_seed(smileses=['CCC[N+](CCC)(CCC)CCC','Cl[Fe-](Cl)(Cl)Cl'],
                           dls = True,
                           **{'min_scaling':0.1,'max_scaling':0.4})
    #ps.initial_seed.write('runs/init.vasp')
    start = dt.now()
    ps.run_seeds(num_seeds=96,
                 convergence=0.1,
                 optimizer_class='FIRE',
                 dir='./runs',
                 num_points_to_converge=3,
                 steps=200)
    
    end= dt.now()
    total = end - start
    print('time took: {} seconds'.format(total.total_seconds()))
    import pickle
    ps.plot_trajectory() 
    pickle.dump(ps,open('data.p','wb'))

if __name__ == "__main__":
    main()
