from pcp.predict import PredictStructure
import datetime as dt 
import pickle

def main():
    ps = PredictStructure(
            num_units=1,
            units=True,
            boxsize=10,
            nprocs=2)
    
    ps.create_initial_seed(smileses=['CCC[N+](CCC)(CCC)CCC','Cl[Fe-](Cl)(Cl)Cl'],
                           dls = True,
                           **{'min_scaling':0.1,'max_scaling':0.4})
    ps.initial_seed.write('testing/init.vasp')
    start = dt.now()
    ps.run_seeds(num_seeds=100,
                 convergence=0.1,
                 optimizer_class='FIRE',
                 dir='./testing',
                 num_points_to_converge=3,
                 **{'verbose':True,'steps':200})
    
    end= dt.now()
    total = end - start
    print('time took: {} seconds'.format(total.total_seconds()))
    import pickle
    
    pickle.dump(ps,open('data.p','wb'))

if __name__ == "__main__":
    main()
