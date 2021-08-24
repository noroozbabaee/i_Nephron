import numpy as np
import pandas as pd
import altair as alt
import scipy
import scipy.optimize
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.pyplot import figure
import math
import pickle

Figure_4 = 1
Figure_5 = 1
Figure_7 = 1
Figure_9_10 = 1
# To produce Figures 6A, 6B, 8A and 8B, we use Altair which is an open-source python library.
# We can create just one figure using Altair.
Figure_6A = 1
Figure_6B = 0
Figure_8A = 0
Figure_8B = 0

if Figure_4:
    f = open('Data_Figure_4a.py', 'rb')
    read_file = f.read()
    loaded_list_4a = pickle.loads(read_file)
    f.close()
    fig4 = plt.figure(constrained_layout=False,  figsize=(10, 10))
    spec2 = gridspec.GridSpec(ncols=1, nrows=3, figure=fig4)
    ax1 = fig4.add_subplot(spec2 [ 0, 0 ])
    ax1.set_title('Figure_4', color='black')
    # ax1.set_ylim([ -20, 45])
    x1 = ax1.plot(loaded_list_4a['mylist_hco3'][52:-2], loaded_list_4a['flux_cl'][52:-2], 'black')
    x2 = ax1.plot(loaded_list_4a['mylist_hco3'][52:-2], loaded_list_4a['fikm_cl'][52:-2], 'g')
    x3 = ax1.plot(loaded_list_4a['mylist_hco3'][52:-2], loaded_list_4a['fekm_cl'][52:-2], 'blue')
    ax1.set_xlabel('Peritubular & Luminal HCO3 Concentration [M]', color='black')
    plt.setp(plt.gca(), xticklabels=[ ])
    f = open('Data_Figure_4b.py', 'rb')
    read_file = f.read()
    loaded_list_4b = pickle.loads(read_file)
    print('loaded_list_4b',loaded_list_4b)
    f.close()
    ax2 = fig4.add_subplot(spec2 [1, 0 ])
    ax2.set_ylabel('flux_cl [mmol/sec]', color='black')
    # ax2.set_ylim([ -20, 45 ])
    x1 = ax2.plot(loaded_list_4b [ 'mylist_hco2' ] [ 52:-2 ], loaded_list_4b [ 'flux_cl' ] [ 52:-2 ], 'black')
    x2 = ax2.plot(loaded_list_4b [ 'mylist_hco2' ] [ 52:-2 ], loaded_list_4b [ 'fikm_cl' ] [ 52:-2 ], 'g')
    x3 = ax2.plot(loaded_list_4b [ 'mylist_hco2' ] [ 52:-2 ], loaded_list_4b [ 'fekm_cl' ] [ 52:-2 ], 'blue')
    ax2.legend([ x1 [ 0 ], x2 [ 0 ], x3 [ 0 ] ], [ 'Total Epithelial', 'Cellular', 'Junctional' ], frameon=False)
    plt.setp(plt.gca(), xticklabels=[ ])
    f = open('Data_Figure_4c.py', 'rb')
    read_file = f.read()
    loaded_list_4c = pickle.loads(read_file)
    f.close()
    ax3= fig4.add_subplot(spec2 [2, 0 ])
    ax3.set_xlabel('Peritubular & Luminal HCO2 Concentration [M]', color='black')
    # ax3.set_ylim([ -20, 45])
    x1 = ax3.plot(loaded_list_4c['mylist_hco2'][52:-2], loaded_list_4c['flux_cl'][52:-2], 'black')
    x2 = ax3.plot(loaded_list_4c['mylist_hco2'][52:-2], loaded_list_4c['fikm_cl'][52:-2], 'g')
    x3 = ax3.plot(loaded_list_4c['mylist_hco2'][52:-2], loaded_list_4c['fekm_cl'][52:-2], 'blue')
    plt.show()
if Figure_5:
    f = open('Data_Figure_5.py', 'rb')
    read_file = f.read()
    loaded_list_5 = pickle.loads(read_file)
    print(loaded_list_5)
    f.close()
    fig3 = plt.figure(constrained_layout=True, figsize=(10, 10))
    spec2 = gridspec.GridSpec(ncols=1, nrows=5, figure=fig3)
    ax1 = fig3.add_subplot(spec2 [ 0, 0 ])
    ax1.plot(loaded_list_5 ['t'] [ 0:-2 ], loaded_list_5['mylist_na'] [ 0:-2 ], color='#007282', label="Na ")
    ax1.plot(loaded_list_5 ['t'] [ 0:-2 ], loaded_list_5['mylist_cl'] [ 0:-2 ], color='#b66dff', label="Cl ")
    ax1.set_ylabel('Cm, Cs', color='black')
    ax1.set_title('Figure_5', color='black')
    ax1.set_xlim(-5, loaded_list_5 ['t'] [ -2 ])
    ax1.legend(frameon=False)
    ax0 = fig3.add_subplot(spec2 [ 1, 0 ])
    ax0.set_ylabel('Vm ', color='black')
    ax0.plot(loaded_list_5 ['t'] [ 0:-2 ], loaded_list_5['vm'] [ 0:-2 ], 'black')
    ax0.set_xlim(-5,loaded_list_5 ['t'] [ -2 ])
    ax0.set_ylim(-0.0, 0.2)
    ax2 = fig3.add_subplot(spec2 [ 2, 0 ])
    ax2.set_ylabel('ci_na ', color='black')
    ax2.plot(loaded_list_5 ['t'] [ 0:-2 ], loaded_list_5['ci_na'] [ 0:-2 ], 'black')
    ax2.set_xlim(-5,loaded_list_5 ['t'] [ -2 ])

    ax3 = fig3.add_subplot(spec2 [ 3, 0 ])
    ax3.plot(loaded_list_5 ['t'] [ 0:-2 ], loaded_list_5['ci_k'] [ 0:-2 ], 'black')
    ax3.set_ylabel('ci_k ', color='black')
    ax3.set_xlim(-5,loaded_list_5 ['t'] [ -2 ])

    ax4 = fig3.add_subplot(spec2 [ 4, 0 ])
    ax4.plot(loaded_list_5 ['t'] [ 0:-2 ], loaded_list_5['ci_cl'] [ 0:-2 ], 'black')
    ax4.set_ylabel('ci_cl ', color='black')
    ax4.set_xlim(-5,loaded_list_5 ['t'] [ -2 ])
    plt.show()
if Figure_7:
    f = open('Data_Figure_7.py', 'rb')
    read_file = f.read()
    my_loaded_list0 = pickle.loads(read_file)
    f.close()
    print('my_loaded_list0', my_loaded_list0['epithelial_flx_variation_default'])
    fig = plt.figure(num=None, figsize=(12, 12))
    x1=['ES', 'IE', 'IS', 'ME', 'MI','EPTL']
    ax1 = fig.add_subplot(3, 3, 1)
    ax1.bar(x1, np.array(my_loaded_list0['flx_na_mem_default']), color=[ 'black', '#b66dff', '#24ff24', '#007282', '#2B9F78', '#920000'])
    # ax1.set_ylim(-0.2,7.5)
    plt.setp(plt.gca(), xticklabels=[])

    # Defines values for the second plot
    x2 = ['CSF', 'GLD', 'ECHMCL']
    print(type(my_loaded_list0['epithelial_flx_variation_default']),my_loaded_list0['epithelial_flx_variation_default'])
    ax2 = fig.add_subplot(3, 3, 2)

    Epithelial = my_loaded_list0['epithelial_flx_variation_default']
    ax2.bar(x2, Epithelial, color=['black', '#b66dff','#2B9F78'])
    # ax2.set_ylim(-0.2,8.4)
    plt.setp(plt.gca(), xticklabels=[])
    x3 = ['NHE3', 'SGLT','NaH2PO4']
    ax3 = fig.add_subplot(3, 3, 3)
    ax3.bar(x3, my_loaded_list0['epithelial_flx_elct_chmcl_default'], color=['#db6d00', '#b66dff', '#2B9F78'])
    # ax3.set_ylim(-0.2,3.8)
    plt.setp(plt.gca(), xticklabels=[])
    ax4 = fig.add_subplot(3, 3, 4)
    labels = ['Interspace Basement', 'Cell Lateral', 'Basal Cell', 'Tight Junction', 'Cell Apical', 'Epithelial']
    colors = ['black', '#b66dff', '#24ff24', '#007282', '#2B9F78', '#920000']
    for tmpX, tmpY, color, label in zip(x1, my_loaded_list0['flx_na_mem_nak'], colors, labels):
        ax4.bar(tmpX, tmpY, color=color, label=label)
    # ax4.set_ylim(-0.2, 7.5)
    plt.setp(plt.gca(), xticklabels=[])
    ax4.legend(shadow=False, fancybox=False, frameon=False)
    ax5 = fig.add_subplot(3, 3, 5)
    labels = ['Covective ', 'Passive', 'Electochemical']
    colors = ['black', '#b66dff', '#2B9F78']
    for tmpX, tmpY, color, label in zip(x2, my_loaded_list0['epithelial_flx_variation_nak'], colors, labels):
        ax5.bar(tmpX, tmpY, color=color, label=label)
    plt.setp(plt.gca(), xticklabels=[])
    # ax5.set_ylim(-0.2, 8.4)
    ax5.legend(shadow=False, fancybox=False, frameon=False)
    ax6=fig.add_subplot(3, 3, 6)
    colors=['#db6d00', '#b66dff','#2B9F78']
    labels=['NHE3','SGLT','NaH2PO4']
    for tmpX, tmpY, color, label in zip(x3, my_loaded_list0['epithelial_flx_elct_chmcl_nak'], colors, labels):
        ax6.bar(tmpX, tmpY, color=color, label=label)
    plt.setp(plt.gca(), xticklabels=[])
    # ax6.set_ylim(-0.1,3.8)
    ax6.legend(shadow=False, fancybox=False, frameon=False)
    ax7=fig.add_subplot(3, 3, 7)
    ax7.bar(x1, my_loaded_list0['flx_na_mem_nhe3'], color=['black', '#b66dff','#24ff24', '#007282','#2B9F78','#920000'])
    ax7.set_xlabel('Membrane Flux_Na [mmol/s.cm2]', color='black')
    # ax7.set_ylim(-0.2,7.5)

    plt.setp(plt.gca(), xticklabels=[])
    # ax7.legend(shadow=False, fancybox=False, frameon=False)
    ax8=fig.add_subplot(3, 3, 8)
    ax8.bar(x2,my_loaded_list0['epithelial_flx_variation_nhe3'], color=['black', '#b66dff','#2B9F78'])
    ax8.set_xlabel('Epithelial Flux_Na [mmol/s.cm2]', color='black')
    plt.setp(plt.gca(), xticklabels=[])
    # ax8.set_ylim(-0.2,8.4)
    # ax8.legend(shadow=False, fancybox=False, frameon=False)
    # label=['ES', 'IE', 'IS', 'ME', 'MI','EPTL']
    ax9 = fig.add_subplot(3, 3, 9)
    ax9.bar(x3, my_loaded_list0['epithelial_flx_elct_chmcl_nhe3'], color=['#db6d00', '#b66dff', '#2B9F78'])
    plt.setp(plt.gca(), xticklabels=[])
    ax9.set_xlabel('Electrochemical Flux_Na [mmol/s.cm2]', color='black')
    # ax9.set_ylim(-0.1,3.8)
    # ax9.legend(shadow=False, fancybox=False, frameon=False)

    for ax, color in zip([ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9], ['black', '#920000', '#2B9F78', 'black', '#920000', '#2B9F78', 'black', '#920000','#2B9F78']):
        for ticks in ax.xaxis.get_ticklines() + ax.yaxis.get_ticklines():
            ticks.set_color(color)
        for pos in ['top', 'bottom', 'right', 'left']:
            ax.spines[pos].set_edgecolor(color)

        plt.show()
if Figure_9_10:
    f = open('Data_Figure_9_10.py', 'rb')
    read_file = f.read()
    loaded_list = pickle.loads(read_file)
    f.close()
    print(loaded_list)
    for key in loaded_list:
        print([key])
    fig8 = plt.figure(constrained_layout=True)
    spec2 = gridspec.GridSpec(ncols=1, nrows=4, figure=fig8)
    ax1 = fig8.add_subplot(spec2 [ 0, 0 ])
    ax1.plot(loaded_list['t'], loaded_list['L_NHE3'], 'r-')
    ax1.set_ylabel('L_NHE3', color='blue')

    ax4 = fig8.add_subplot(spec2 [ 1, 0 ])
    ax4.plot(loaded_list['t'], loaded_list['flux_na_nah'], 'b-')
    ax4.set_ylabel('flux_na_nah', color='blue')

    ax2 = fig8.add_subplot(spec2 [ 2, 0 ])
    ax2.plot(loaded_list['t'], loaded_list['flux_hco3_clcho3'], 'b-')
    ax2.set_ylabel('flux_hco3_clcho3 ', color='blue')

    ax3 = fig8.add_subplot(spec2 [ 3, 0 ])
    ax3.plot(loaded_list['t'], loaded_list['flux_Na_nak'], 'b-')
    ax3.set_ylabel('flux_na_nak ', color='blue')

    fig3 = plt.figure(constrained_layout=True)
    spec2 = gridspec.GridSpec(ncols=1, nrows =  6, figure=fig3)
    ax1 = fig3.add_subplot(spec2 [ 0, 0 ])
    ax1.plot(loaded_list['t'], loaded_list['L_NHE3'], 'r-')
    ax1.set_ylabel('L_NHE3', color='blue')

    ax2 = fig3.add_subplot(spec2 [ 1, 0 ])
    ax2.plot(loaded_list['t'], loaded_list['ci_na'], 'b-')
    ax2.set_ylabel(' ci_na ', color='blue')

    ax3 = fig3.add_subplot(spec2 [ 2, 0 ])
    ax3.plot(loaded_list['t'], loaded_list['ci_cl'], 'b-')
    ax3.set_ylabel('ci_cl ', color='blue')

    ax4 = fig3.add_subplot(spec2 [ 3, 0 ])
    ax4.plot(loaded_list['t'], loaded_list['ci_hco3'], 'b-')
    ax4.set_ylabel('ci_hco3 ', color='blue')
    ax5 = fig3.add_subplot(spec2 [ 4, 0 ])
    ax5.plot(loaded_list['t'], loaded_list['ci_hco2'], 'b-')
    ax5.set_ylabel('ci_hco2 ', color='blue')

    ax6 = fig3.add_subplot(spec2 [ 5, 0 ])
    ax6.plot(loaded_list['t'], loaded_list['clvl_imp'], 'b-')
    ax6.set_ylabel('clvl ', color='blue')
    ax6.set_xlabel('Time (s) ', color='blue')
    plt.show()
if Figure_6A:
    f = open('Data_Figure_6A.py', 'rb')
    read_file = f.read()
    my_loaded_list0 = pickle.loads(read_file)
    f.close()
    scale_factor = 1
    flux_me = np.array(my_loaded_list0 [ 'fekm' ])
    flux_mi = np.array(my_loaded_list0 [ 'fikm' ])
    flux_is = np.array(my_loaded_list0 [ 'fiks' ])
    flux_ie = np.array(my_loaded_list0 [ 'fike' ])
    flux_es = np.array(my_loaded_list0 [ 'feks' ])
    df1 = pd.DataFrame(scale_factor * flux_me, index=[ "Original Setup", "NaK = 0", "KCl = 0", "NaHCO3 = 0" ],
                       columns=[ "Na", "K", "Cl", "Gluc" ])
    df2 = pd.DataFrame(scale_factor * flux_is, index=[ "Original Setup", "NaK = 0", "KCl = 0", "NaHCO3 = 0" ],
                       columns=[ "Na", "K", "Cl", "Gluc" ])
    df3 = pd.DataFrame(scale_factor * flux_ie, index=[ "Original Setup", "NaK = 0", "KCl = 0", "NaHCO3 = 0" ],
                       columns=[ "Na", "K", "Cl", "Gluc" ])
    df4 = pd.DataFrame(scale_factor * flux_es, index=[ "Original Setup", "NaK = 0", "KCl = 0", "NaHCO3 = 0" ],
                       columns=[ "Na", "K", "Cl", "Gluc" ])
    df5 = pd.DataFrame(scale_factor * flux_mi, index=[ "Original Setup", "NaK = 0", "KCl = 0", "NaHCO3 = 0" ],
                       columns=[ "Na", "K", "Cl", "Gluc" ])
    def prep_df(df, name):
        df = df.stack().reset_index()
        df.columns = [ 'c1', 'c2', 'values' ]
        df [ 'Fluxes' ] = name
        return df
    df1 = prep_df(df1, 'ME')
    df2 = prep_df(df2, 'IS')
    df3 = prep_df(df3, 'IE')
    df4 = prep_df(df4, 'ES')
    df5 = prep_df(df5, 'MI')
    df = pd.concat([ df1, df2, df3, df4, df5 ])
    print(df)

    alt.renderers.enable('altair_viewer')
    chart = alt.Chart(df).mark_bar().encode(
        alt.X('c2', title=None, sort=[ "Na", "K", "Cl", "Gluc" ]),
        alt.Y('sum(values)',
              axis=alt.Axis(
                  grid=False,
                  title=None)),

        alt.Column('c1', title=None, sort=[ "Original Setup", "NaK = 0", "KCl = 0", "NaHCO3 = 0" ]),
        alt.Color('Fluxes',
                  scale=alt.Scale(
                      range=[ 'black', '#D55E00', '#F0E442', '#007282', '#2B9F78' ]
                  ),
                  )) \
        .configure_view(
    ).configure_axis(
        grid=False
    )
    chart.show()

if Figure_6B:
    f = open('Data_Figure_6B.py', 'rb')
    read_file = f.read()
    my_loaded_list1 = pickle.loads(read_file)
    f.close()
    Na = np.array(my_loaded_list1 [ 'Na' ])
    K = np.array(my_loaded_list1 [ 'K' ])
    Cl = np.array(my_loaded_list1 [ 'Cl' ])
    Gluc = np.array(my_loaded_list1 [ 'Gluc' ])
    df = pd.DataFrame({
        'index': [ "A Original Setup", "B NaK = 0", "C KCl = 0", "D NaHCO3 = 0" ],
        'Na': Na,
        'K': K,
        'Cl': Cl,
        'Gluc': Gluc
    })
    chart = alt.Chart(df.melt('index')).mark_bar().encode(
        alt.X('variable:N', axis=alt.Axis(title=''), sort=alt.Sort(None)),
        alt.Y('value:Q', axis=alt.Axis(title='', grid=False)),
        column='index:N'
    ).configure_view(
    ).encode(
        color=alt.Color('variable:N', scale=alt.Scale(range=[ 'black', '#D55E00', '#007282', '#2B9F78' ]))
    )
    chart.show()

if Figure_8A:
    f = open('Data_Figure_8A.py', 'rb')
    read_file = f.read()
    my_loaded_list0 = pickle.loads(read_file)
    f.close()
    scale_factor = 1
    flux_me = np.array(my_loaded_list0 [ 'fekm' ])
    flux_mi = np.array(my_loaded_list0 [ 'fikm' ])
    flux_is = np.array(my_loaded_list0 [ 'fiks' ])
    flux_ie = np.array(my_loaded_list0 [ 'fike' ])
    flux_es = np.array(my_loaded_list0 [ 'feks' ])
    df1 = pd.DataFrame(flux_me,
                       index=[ "Original Setup", "NHE3 = 0", "SGLT = 0", "NaH2PO4 = 0" ],
                       columns=[ "Na", "K", "Cl", "Gluc" ])
    df2 = pd.DataFrame(flux_is,
                       index=[ "Original Setup", "NHE3 = 0", "SGLT = 0", "NaH2PO4 = 0" ],
                       columns=[ "Na", "K", "Cl", "Gluc" ])
    df3 = pd.DataFrame(flux_ie,
                       index=[ "Original Setup", "NHE3 = 0", "SGLT = 0", "NaH2PO4 = 0" ],
                       columns=[ "Na", "K", "Cl", "Gluc" ])
    df4 = pd.DataFrame(flux_es,
                       index=[ "Original Setup", "NHE3 = 0", "SGLT = 0", "NaH2PO4 = 0" ],
                       columns=[ "Na", "K", "Cl", "Gluc" ])
    df5 = pd.DataFrame(flux_mi,
                       index=[ "Original Setup", "NHE3 = 0", "SGLT = 0", "NaH2PO4 = 0" ],
                       columns=[ "Na", "K", "Cl", "Gluc" ])
    def prep_df(df, name):
        df = df.stack().reset_index()
        df.columns = [ 'c1', 'c2', 'values' ]
        df [ 'Fluxes' ] = name
        return df
    df1 = prep_df(df1, 'ME')
    df2 = prep_df(df2, 'IS')
    df3 = prep_df(df3, 'IE')
    df4 = prep_df(df4, 'ES')
    df5 = prep_df(df5, 'MI')
    df = pd.concat([ df1, df2, df3, df4, df5 ])
    print(df)

    alt.renderers.enable('altair_viewer')
    chart = alt.Chart(df).mark_bar().encode(
        alt.X('c2', title=None, sort=[ "Na", "K", "Cl", "Gluc" ]),
        alt.Y('sum(values)',
              axis=alt.Axis(
                  grid=False,
                  title=None)),
        alt.Column('c1', title=None, sort=[ "Original Setup", "NHE3 = 0", "SGLT = 0", "NaH2PO4 = 0" ]),
        alt.Color('Fluxes',
                  scale=alt.Scale(
                       range=[ 'black', '#D55E00', '#F0E442', '#007282', '#2B9F78' ]
                    ),
                  )) \
        .configure_view(
    ).configure_axis(
        grid=False
    )
    chart.show()

if Figure_8B:
    f = open('Data_Figure_8B.py', 'rb')
    read_file = f.read()
    my_loaded_list1 = pickle.loads(read_file)

    f.close()
    print('my_loaded_list1', my_loaded_list1)
    Na = np.array(my_loaded_list1 [ 'Na' ])
    K = np.array(my_loaded_list1 [ 'K' ])
    Cl = np.array(my_loaded_list1 [ 'Cl' ])
    Gluc = np.array(my_loaded_list1 [ 'Gluc' ])

    df = pd.DataFrame({
        'index': [ "A Original Setup", "B NHE3 = 0", "C SGLT = 0", "D NaH2PO4 = 0" ],
        'Na': Na,
        'K': K,
        'Cl': Cl,
        'Gluc': Gluc
    })

    chart = alt.Chart(df.melt('index')).mark_bar().encode(
        alt.X('variable:N', axis=alt.Axis(title=''), sort=alt.Sort(None)),
        alt.Y('value:Q', axis=alt.Axis(title='', grid=False)),
        column='index:N'
    ).configure_view(
    ).encode(
        color=alt.Color('variable:N', scale=alt.Scale(range=[ 'black', '#D55E00', '#007282', '#2B9F78' ]))
    )
    chart.show()
    
