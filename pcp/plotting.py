import pandas as pd 
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
        plt.show()    