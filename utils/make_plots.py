import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import os
from matplotlib.ticker import AutoMinorLocator

plot_dir = "plots/"


sns.set_context("paper", 
                rc={"font.size":20,
                    "axes.titlesize":20,
                    "axes.labelsize":20, 
                    "axes.linewidth":1.25,
                    "lines.linewidth":2.5,
                    "lines.markersize":7.5},
                font_scale=1.5)
sns.set_style({"xtick.direction": "in","ytick.direction": "in"})

def make_energy_plots(s_range, predict_data, true_data, exec_args):

    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)

    plot_fname = plot_dir+"E_flow_{}_{}_{}_{:0.2f}_{:0.2f}_{:0.2f}.pdf".format(exec_args["nobs"], 
                                                                               int(exec_args["exact"]), 
                                                                               exec_args["trunc"], 
                                                                               exec_args["t0"], 
                                                                               exec_args["t1"], 
                                                                               exec_args["dt"])
    
    fig,axes = plt.subplots(figsize=[8,4])
    for a in np.ravel(axes):
        a.tick_params(which="both", bottom=True)
        a.xaxis.set_minor_locator(AutoMinorLocator())
        a.yaxis.set_minor_locator(AutoMinorLocator())
        

    sns.lineplot(x=s_range, y=true_data, ax=axes)        
    sns.lineplot(x=s_range, y=predict_data, ax=axes)
    axes.axvline(s_range[exec_args["nobs"]], color='blue', linestyle='--')
    axes.axhline(true_data[0], color='red', linestyle='--')
    axes.axhline(true_data[-1], color='red', linestyle='--')

    axes.set_ylabel("E(s)")
    axes.set_xlabel("s")
    
    axes.legend(["IMSRG E flow", "DMD E flow", r"Max observed $s$"])
    
    fig.tight_layout()
    fig.savefig(plot_fname)
    
def make_correlation_plots(s_range, predict_data, true_data, exec_args):
    if not os.path.exists(plot_dir):
        os.mkdir(plot_dir)

    plot_fname = plot_dir+"Corr_flow_{}_{}_{}_{:0.2f}_{:0.2f}_{:0.2f}.pdf".format(exec_args["nobs"], 
                                                                                  int(exec_args["exact"]), 
                                                                                  exec_args["trunc"], 
                                                                                  exec_args["t0"], 
                                                                                  exec_args["t1"], 
                                                                                  exec_args["dt"])
    fig,axes = plt.subplots(figsize=[8,4])
    for a in np.ravel(axes):
        a.tick_params(which="both", bottom=True)
        a.xaxis.set_minor_locator(AutoMinorLocator())
        a.yaxis.set_minor_locator(AutoMinorLocator())
        

    sns.lineplot(x=s_range, y=[np.linalg.norm(predict_data[:,i] - true_data[:,i])/np.linalg.norm(true_data[:,i]) for i in range(true_data.shape[1])], ax=axes)        

    axes.set_ylabel(r"$\frac{||H_{DMD}(s) - H_{IMSRG}(s)||_F}{||H_{IMSRG}(s)||_F}$")
    axes.set_xlabel("s")
        
    fig.tight_layout()
    fig.savefig(plot_fname)
    
