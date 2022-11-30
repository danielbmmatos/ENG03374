# -*- coding: utf-8 -*-
"""
                             Plot config
                     
                     
                     
Author: Daniel B. M. Matos and Augusto Bopsin Borges
Criation: 05/09/2022
Modified: 05/09/2022
"""
def plot_config():
    import matplotlib
    from matplotlib import rc
    rc('font',**{'family':'Times New Roman','size' :12})
    matplotlib.rcParams['axes.unicode_minus'] = False
    matplotlib.rc('text', usetex=True)
    matplotlib.rc('axes', labelsize=12)
    matplotlib.rc('xtick', labelsize=12) 
    matplotlib.rc('ytick', labelsize=12)
    matplotlib.rc('lines', lw=1.0,color='k')
    matplotlib.rc('axes',lw=0.75)
    matplotlib.rc('legend', fontsize=12)