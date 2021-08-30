import numpy as np
import scipy
import scipy.optimize
from scipy import optimize
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import math
from W_PCT_E_Glob import *
import pandas as pd
import matplotlib.cm as cm
import pickle

def plot_clustered_stacked(dfall, labels=None, title="multiple stacked bar plot", H="/",**kwargs):
    """Given a list of dataframes, with identical columns and index, create a clustered stacked bar plot.
labels is a list of the names of the dataframe, used for the legend
title is a string for the title of the plot
H is the hatch used for identification of the different dataframe"""

    n_df = len(dfall)
    n_col = len(dfall[0].columns)
    n_ind = len(dfall[0].index)

    axe = plt.subplot(111)

    for df in dfall : # for each data frame
        axe = df.plot(kind="bar",
                      linewidth=0,
                      stacked=True,
                      colormap="Oranges_r",
                      ax=axe,
                      legend = False,
                      grid = False,
                      **kwargs)  # make bar plots

    h,l = axe.get_legend_handles_labels() # get the handles we want to modify
    for i in range(0, n_df * n_col, n_col): # len(h) = n_col * n_df
        for j, pa in enumerate(h[i:i+n_col]):
            for rect in pa.patches: # for each index
                rect.set_x(rect.get_x() + 1 / float(n_df + 1) * i / float(n_col))
                rect.set_hatch(H * int(i / n_col)) #edited part
                rect.set_width(1 / float(n_df + 1))

    axe.set_xticks((np.arange(0, 2 * n_ind, 2) + 1 / float(n_df + 1)) / 2.)
    axe.set_xticklabels(df.index, rotation = 0)
    axe.set_title(title)

    # Add invisible data to add another legend
    n= []
    for i in range(n_df):
        n.append(axe.bar(0, 0, color="gray", hatch=H * i))

    l1 = axe.legend(h[:n_col], l[:n_col], loc= [1.01, 0.5])
    if labels is not None:
        l2 = plt.legend(n, labels, loc= [1.01, 0.1])
    axe.add_artist(l1)
    return axe


def Extract(lst):
    return[item[-1] for item in lst]


def csf(ca, cb, flux, s, param_csf):
    # Convective Solute Fluxes, see Eqs: (37, 38)
    if param_csf == 0:
        return 0

    def lmmsc(ca, cb):
        if ca > 0 and cb > 0 and ca - cb != 0 and cb != 0 and np.log10(ca / cb) != 0:
            return (ca - cb) / (np.log10(ca / cb))
        else:
            return cb
    return flux * (1.00 - s) * lmmsc(ca, cb)


def goldman(hab, a, z, va, vb, ca, cb, param_goldman):
    # see  Eqs: (39), (41)
    zab = (z * f * (va - vb) * 1.e-6) / rte
    if param_goldman == 0:
        return[0]
    elif zab == 0 or va == vb and ca > 0 and cb > 0:
        return hab * a * (ca - cb)
    elif zab > 0:
        return hab * a * zab * (ca - cb * np.exp(-zab)) / (1 - np.exp(-zab))
    else:
        return hab * a * zab * (ca * np.exp(zab) - cb) / (np.exp(zab) - 1)


def k_cl(ca_k, cb_k, ca_cl, cb_cl, z_k, z_cl, va, vb, a, l_kcl, param_kcl):
    # see Eq: (42)
    if param_kcl == 0:
        return[0, 0]
    else:
        xa_k = f_eps(ca_k, z_k, va)
        xa_cl = f_eps(ca_cl, z_cl, va)
        xb_k = f_eps(cb_k, z_k, vb)
        xb_cl = f_eps(cb_cl, z_cl, vb)
        k_kcl = l_kcl * a * (xa_k - xb_k + xa_cl - xb_cl)
        cl_kcl = l_kcl * a * (xa_k - xb_k + xa_cl - xb_cl)
    return[k_kcl, cl_kcl]


def sglt_mi(cm_na, ci_na, cm_gluc, ci_gluc, z_na, z_gluc, vm, vi, ami, lmi_nagluc, param_sglt_mi):
    #  Eq: (43)
    if param_sglt_mi == 0:
        return[0, 0]
    else:
        xm_na = f_eps(cm_na, z_na, vm)
        xm_gluc = f_eps(cm_gluc, z_gluc, vm)
        xi_na = f_eps(ci_na, z_na, vi)
        xi_gluc = f_eps(ci_gluc, z_gluc, vi)
        na_mi_nagluc = lmi_nagluc * ami * (xm_na - xi_na + xm_gluc - xi_gluc)
        gluc_mi_nagluc = lmi_nagluc * ami * (xm_na - xi_na + xm_gluc - xi_gluc)
    return[na_mi_nagluc, gluc_mi_nagluc]


def nah2po4_mi(cm_na, ci_na, cm_h2po4, ci_h2po4, z_na, z_h2po4, vm, vi, ami, lmi_nah2po4, param_nah2po4_mi):
    #  Eq: (44)
    if param_nah2po4_mi == 0:
        return[0, 0]
    else:
        xm_na = f_eps(cm_na, z_na, vm)
        xm_h2po4 = f_eps(cm_h2po4, z_h2po4, vm)
        xi_na = f_eps(ci_na, z_na, vi)
        xi_h2po4 = f_eps(ci_h2po4, z_h2po4, vi)
        na_mi_nah2po4 = lmi_nah2po4 * ami * (xm_na - xi_na + xm_h2po4 - xi_h2po4)
        h2po4_mi_nah2po4 = lmi_nah2po4 * ami * (xm_na - xi_na + xm_h2po4 - xi_h2po4)
    return[na_mi_nah2po4, h2po4_mi_nah2po4]


def nah(cm_na, ci_na, cm_h, ci_h, z_na, z_h, vm, vi, ami, lmi_nah, param_nah_mi):
    # Eq: (45)
    if param_nah_mi == 0:
        return[0, 0]
    else:
        xm_na = f_eps(cm_na, z_na, vm)
        xm_h = f_eps(cm_h, z_h, vm)
        xi_na = f_eps(ci_na, z_na, vi)
        xi_h = f_eps(ci_h, z_h, vi)
        na_mi_nah = lmi_nah * ami * (xm_na - xi_na - xm_h + xi_h)
        h_mi_nah = - lmi_nah * ami * (xm_na - xi_na - xm_h + xi_h)
    return[na_mi_nah, h_mi_nah]


def nanh4(cm_na, ci_na, cm_nh4, ci_nh4 ,z_na, z_nh4, vm, vi, ami, lmi_nanh4, param_nanh4_mi):
    # Eq: (46)
    if param_nanh4_mi == 0:
        return[0, 0]
    else:
        xm_na = f_eps(cm_na, z_na, vm)
        xm_nh4= f_eps(cm_nh4, z_nh4, vm)
        xi_na = f_eps(ci_na, z_na, vi)
        xi_nh4 = f_eps(ci_nh4, z_nh4, vi)
        na_mi_nanh4 = lmi_nanh4 * ami * (xm_na - xi_na - xm_nh4 + xi_nh4)
        h_mi_nanh4 = - lmi_nanh4 * ami * (xm_na - xi_na - xm_nh4 + xi_nh4)
    return[na_mi_nanh4, h_mi_nanh4]


def clhco2_mi(cm_cl, ci_cl, cm_hco2, ci_hco2, z_cl, z_hco2, vm, vi, ami, lmi_clhco2, param_clhco2_mi):
    # Eq: (47)
    if param_clhco2_mi == 0:
        return[0, 0]
    else:
        xm_cl = f_eps(cm_cl, z_cl, vm)
        xm_hco2 = f_eps(cm_hco2, z_hco2, vm)
        xi_cl = f_eps(ci_cl, z_cl, vi)
        xi_hco2 = f_eps(ci_hco2, z_hco2, vi)
        cl_mi_clhco2 = lmi_clhco2 * ami * (xm_cl - xi_cl - xm_hco2 + xi_hco2)
        hco2_mi_clhco2 = - lmi_clhco2 * ami * (xm_cl - xi_cl - xm_hco2 + xi_hco2)
        return[cl_mi_clhco2, hco2_mi_clhco2]


def clhco3_mi(cm_cl, ci_cl, cm_hco3, ci_hco3, z_cl, z_hco3, vm, vi, ami, lmi_clhco3, param_clhco3_mi):
    # See Eq: (48)
    if param_clhco3_mi == 0:
        return[0, 0]
    else:
        xm_cl = f_eps(cm_cl, z_cl, vm)
        xm_hco3 = f_eps(cm_hco3, z_hco3, vm)
        xi_cl = f_eps(ci_cl, z_cl, vi)
        xi_hco3 = f_eps(ci_hco3, z_hco3, vi)
        cl_mi_clhco3 = lmi_clhco3 * ami * (xm_cl - xi_cl - xm_hco3 + xi_hco3)
        hco3_mi_clhco3 = - lmi_clhco3 * ami * (xm_cl - xi_cl - xm_hco3 + xi_hco3)
    return[cl_mi_clhco3, hco3_mi_clhco3]


def na_hco3(ci_na, cs_na, ci_hco3, cs_hco3, z_na, z_hco3, vi, vs, ais, lis_nahco3, param_na_hco3):
    # See Eq: (49)
    if param_na_hco3 == 0:
        return [ 0, 0]
    else:
        xi_na = f_eps(ci_na, z_na, vi)
        xi_hco3 = f_eps(ci_hco3, z_hco3, vi)
        xs_na = f_eps(cs_na, z_na, vs)
        xs_hco3 = f_eps(cs_hco3, z_hco3, vs)
        na_is_nahco3 = lis_nahco3 * ais * (xi_na - xs_na + 3 * (xi_hco3 - xs_hco3))
        hco3_is_nahco3 =  3 * lis_nahco3 * ais * (xi_na - xs_na + 3 * (xi_hco3 - xs_hco3))
        return[na_is_nahco3, hco3_is_nahco3]


def na_cl_hco3(ci_na, cs_na, ci_cl, cs_cl, ci_hco3, cs_hco3, z_na, z_cl, z_hco3, vi, vs, ais, lis_na_clhco3,
                 param_na_cl_hco3):
    # See Eq: (50)
    if param_na_cl_hco3 == 0:
        return[0, 0, 0]
    else:
        xi_na = f_eps(ci_na, z_na, vi)
        xi_cl = f_eps(ci_cl, z_cl, vi)
        xi_hco3 = f_eps(ci_hco3, z_hco3, vi)
        xs_na = f_eps(cs_na, z_na, vs)
        xs_cl = f_eps(cs_cl, z_cl, vs)
        xs_hco3 = f_eps(cs_hco3, z_hco3, vs)
        na_na_clhco3 = + ais * lis_na_clhco3 * (xi_na - xs_na - xi_cl + xs_cl + 2 * (xi_hco3 - xs_hco3))
        cl_na_clhco3 = - ais * lis_na_clhco3 * (xi_na - xs_na - xi_cl + xs_cl + 2 * (xi_hco3 - xs_hco3))
        hco3_na_clhco3 = 2 * ais * lis_na_clhco3 * (xi_na - xs_na - xi_cl + xs_cl + 2 * (xi_hco3 - xs_hco3))
        return[na_na_clhco3, cl_na_clhco3, hco3_na_clhco3]


def h_atp_mi(cm_h, ci_h, vm, vi, z_h, ami, parama_mi_h):
    # See Eq: (53)
    xm_h = f_eps(cm_h, z_h, vm)
    xi_h = f_eps(ci_h, z_h, vi)
    gamma = -(xihp * (xm_h - xi_h - xhp))
    if parama_mi_h == 0:
        return 0
    elif gamma < 0:
        return -lhp * ami * (1 - 1 / (1 + np.exp(gamma)))
    return -lhp * ami * (1 / (1 + np.exp(-gamma)))


def nak_atp(ci_na, ce_na, ci_k,  ce_k, ce_nh4, param_nak_atp):
    # See Eqs: (54, 55, 56, 57, 58)
    if param_nak_atp == 0:
        return[0, 0, 0]
    else:
        knpn = 0.0002 * (1.0 + ci_k / .00833)
        knpk = 0.0001 * (1.0 + ce_na / .0185)
        a = (ci_na / (knpn + ci_na)) ** 3
        b = ((ce_k + ce_nh4) / (knpk + ce_k + ce_nh4)) ** 2
        atp_na = n_p * a * b
        atp_k = -atp_na * (2/3) * ce_k / (ce_k + ce_nh4 / knh4)
        atp_nh4 = -atp_na * (2/3) * ce_nh4 / (knh4 * ce_k + ce_nh4)
        return[atp_na, atp_k, atp_nh4]


def lch(ca, cb):
    if ca > 0 and cb > 0 and (ca / cb) != 0 and cb != 0:
        return np.log10(ca / cb)
    else:
        return np.log10(abs(ca / cb))


def f_eps(c, z, v):
    # see Eq: (52)
    if c > 0 and c != 0:
        return rte * np.log(c) + z * f * v * 1.e-6
    else:
        print('uy')
        return rte * np.log(abs(c)) + z * f * v * 1.e-6


def ebuf(lch, pk, ca, cb, param_ebuf):
    #  see Eqs: (15), (16)
    if param_ebuf == 0:
        return 0
    if ca > 0 and cb > 0 and (ca / cb) != 0 and cb != 0:
        return lch - pk - np.log10(ca / cb)

    else:
        return lch - pk - np.log10(abs(ca / cb))


def buff_activation_co2_formate_phosphate_ammonia(q_h, q_nh4, q_hco3, q_h2co2, q_h2po4, q_co2, q_h2co3, q_cb, c_co2, c_h2co3, volume, scale,flow_dependent, Co2_Progressive_Activation_Param):
    # See Eqs: (28), (29)
    if Co2_Progressive_Activation_Param == 1:
        # print('q_hco3,  q_h2co3, q_co2', str(q_hco3), str(q_h2co3), str(q_co2))

        if flow_dependent == 1:
            q_hco3 = -1 * (q_h + q_nh4 - q_hco3 + q_h2co2 + q_h2po4 - q_cb)
            q_h2co3 = -1 * (- q_h2co3 - q_h2co3 + scale * volume * (khy * c_co2 - kdhy * c_h2co3))
            q_co2  = q_hco3 + q_h2co3 + q_co2
        else:
            q_hco3 = -1 * ( q_h + q_nh4 - q_hco3 + q_h2co2 + q_h2po4)
            q_h2co3 = -1 * (- q_h2co3 - q_h2co3 + scale * volume * (khy * c_co2 - kdhy * c_h2co3))
            q_co2 = q_hco3 + q_h2co3 + q_co2

    return q_hco3,  q_h2co3, q_co2


def buff_activation(q_h_vary, q_h2_vary, lc_h, pk, c_h_vary, c_h2_vary, qi_nh3_cell, Buffer_Activation_Param):
    # See Eqs:(20), (21), (25)
    if Buffer_Activation_Param == 1:
        if qi_nh3_cell == 1:
           q_h_vary = q_h_vary + q_h2_vary - 1.e6 * qiamm
           q_h2_vary = lc_h - pk - lch(c_h_vary ,c_h2_vary)
        else:
           q_h_vary = q_h_vary + q_h2_vary
           q_h2_vary = lc_h - pk - lch(c_h_vary, c_h2_vary)
    else:
        print('No Buffer Effect for phosphate')
        pass
    return q_h_vary, q_h2_vary


def zab(z, va, vb):
    return z * f * (va - vb) * 1.e-6 / rte


def phi_scale(phi, scale):
    return phi * scale


def nhe3(ci_h, ci_na, ci_nh4, cm_h, cm_na, cm_nh4, param_nah):
     if param_nah == 0:
        return[0, 0, 0]
     else:
        cxt = 1.00
        knah_na = 0.3000e-01
        knah_h = 0.7200e-07
        knah_nh4 = 0.2700e-01
        knah_i = 0.1000e-05
        pnah_na = 0.792000e+04
        pnah_h = 0.23800e+04
        pnah_nh4 = 0.792000e+04
        # paramters from original work NHE_1995
        # pnah_na = 0.160000e-02
        # pnah_h = 0.4800e-03
        # pnah_nh4 = 0.160000e-2
        pnah_mm = 0.2000e+01
        pnah_m = 0.0000
        psnah_na = pnah_na * (pnah_mm * ci_h + pnah_m * knah_i) / (ci_h + knah_i)
        psnah_h = pnah_h * (pnah_mm * ci_h + pnah_m * knah_i) / (ci_h + knah_i)
        psnah_nh4 = pnah_nh4 * (pnah_mm * ci_h + pnah_m * knah_i) / (ci_h + knah_i)
        cisnah_na = ci_na / knah_na
        cmsnah_na = cm_na / knah_na
        cisnah_h = ci_h / knah_h
        cmsnah_h = cm_h / knah_h
        cisnah_nh4 = ci_nh4 / knah_nh4
        cmsnah_nh4 = cm_nh4 / knah_nh4
        delta_i = 1.0
        delta_m = 1.0
        delta_i = delta_i + cisnah_na + cisnah_h + cisnah_nh4
        delta_m = delta_m + cmsnah_na + cmsnah_h + cmsnah_nh4
        epsilon_i = psnah_na * cisnah_na + psnah_h * cisnah_h + psnah_nh4 * cisnah_nh4
        epsilon_m = psnah_na * cmsnah_na + psnah_h * cmsnah_h + psnah_nh4 * cmsnah_nh4
        sigma_nhe3 = delta_i * epsilon_m + delta_m * epsilon_i
        jnah_na = 0.0
        jnah_h = 0.0
        jnah_nh4 = 0.0
        jnah_na = jnah_na + (psnah_na * psnah_h * cxt / sigma_nhe3) * (cmsnah_na * cisnah_h - cisnah_na * cmsnah_h) + (
                    psnah_na * psnah_nh4 * cxt / sigma_nhe3) * (cmsnah_na * cisnah_nh4 - cisnah_na * cmsnah_nh4)

        jnah_h = jnah_h + (psnah_h * psnah_na * cxt / sigma_nhe3) * (cmsnah_h * cisnah_na - cisnah_h * cmsnah_na) + (
                    psnah_h * psnah_nh4 * cxt / sigma_nhe3) * (cmsnah_h * cisnah_nh4 - cisnah_h * cmsnah_nh4)
        jnah_nh4 = jnah_nh4 + (psnah_nh4 * psnah_na * cxt / sigma_nhe3) * (
                    cmsnah_nh4 * cisnah_na - cisnah_nh4 * cmsnah_na) + (psnah_nh4 * psnah_h * cxt / sigma_nhe3) * (
                               cmsnah_nh4 * cisnah_h - cmsnah_h * cisnah_nh4)
        # jnah_na_max = cxt * psnah_na * psnah_h / (psnah_na + psnah_h)
        return[jnah_na, jnah_h, jnah_nh4]


def matrix(init, h):
    return[init for x in range(h)]


def eQs(guess, solver):
    
    for i in range(len(guess)):
        vars[i][t] = guess[i]
    ae[t] = ae0 * (1 + mua * (pe[t] - pm))
    ae[t] = ae0 if ae[t] < ae0 else ae[t]

    chvl[t] = chvl0 * (1.0 + muv * (pe[t] - pm))
    chvl[t] = chvl0 if chvl[t] < chvl0 else chvl[t]
    l[t] = chvl[t]

    lche = pkc + lch(ce_hco3[t], ce_h2co3[t])
    ce_h[t] = 10 ** (-lche)

    clvl[t] = clvl0 * imp0 / imp[t]
    clvl[t] = clvl0 if imp[t] == 0 else clvl[t]
    clvl_imp = clvl0 * imp0 / imp[t]
    
    l[t] = l[t] + clvl[t]
    p_i = pm

    lchi = pkc + lch(ci_hco3[t], ci_h2co3[t])
    ci_h[t] = 10 ** (-lchi)

    # The Water Fluxes, See Eq: (36)
    fevm = ame * lpme * (pm - pe[t] - rt * impm) / rt
    fevs = ae[t] * lpes * (rt * imps + pe[t] - ps) / rt
    fivm = ami * lpmi * (rt * imp[t] - rt * impm + pm - p_i) / rt
    fivs = ais * lpis * (p_i - ps + rt * imps - rt * imp[t]) / rt
    five = aie * lpie * (p_i - pe[t] - rt * imp[t]) / rt
    
    fevm = fevm + ame * lpme * (
            sme_na * (ce_na[t] - cm_na) + sme_k * (
                    ce_k[t] - cm_k) + sme_cl * (
                    ce_cl[t] - cm_cl) + sme_hco3 * (
                    ce_hco3[t] - cm_hco3) + sme_h2co3 * (
                    ce_h2co3[t] - cm_h2co3) + sme_co2 * (
                    ce_co2[t] - cm_co2) + sme_hpo4 * (
                    ce_hpo4[t] - cm_hpo4) + sme_h2po4 * (
                    ce_h2po4[t] - cm_h2po4) + sme_urea * (
                    ce_urea[t] - cm_urea) + sme_nh3 * (
                    ce_nh3[t] - cm_nh3) + sme_nh4 * (
                    ce_nh4[t] - cm_nh4) + sme_h * (
                    ce_h[t]- cm_h) + sme_hco2 * (
                    ce_hco2[t] - cm_hco2) + sme_h2co2 * (
                    ce_h2co2[t] - cm_h2co2) + sme_gluc * (
                    ce_gluc[t] - cm_gluc))

    fevs = fevs + ae[t] * lpes * (
            ses_na * (cs_na - ce_na[t]) + ses_k * (
                    cs_k - ce_k[t]) + ses_cl * (
                    cs_cl - ce_cl[t]) + ses_hco3 * (
                    cs_hco3 - ce_hco3[t]) + ses_h2co3 * (
                    cs_h2co3 - ce_h2co3[t]) + ses_co2 * (
                    cs_co2 - ce_co2[t]) + ses_hpo4 * (
                    cs_hpo4 - ce_hpo4[t]) + ses_h2po4 * (
                    cs_h2po4 - ce_h2po4[t]) + ses_urea * (
                    cs_urea - ce_urea[t]) + ses_nh3 * (
                    cs_nh3 - ce_nh3[t]) + ses_nh4 * (
                    cs_nh4 - ce_nh4[t]) + ses_h * (
                    cs_h - ce_h[t]) + ses_hco2 * (
                    cs_hco2 - ce_hco2[t]) + ses_h2co2 * (
                    cs_h2co2 - ce_h2co2[t]) + ses_gluc * (
                    cs_gluc - ce_gluc[t]))

    fivm = fivm + ami * lpmi * (
            smi_na * (ci_na[t] - cm_na) + smi_k * (
                    ci_k[t] - cm_k) + smi_cl * (
                    ci_cl[t] - cm_cl) + smi_hco3 * (
                    ci_hco3[t] - cm_hco3) + smi_h2co3 * (
                    ci_h2co3[t] - cm_h2co3) + smi_co2 * (
                    ci_co2[t] - cm_co2) + smi_hpo4 * (
                    ci_hpo4[t] - cm_hpo4) + smi_h2po4 * (
                    ci_h2po4[t] - cm_h2po4) + smi_urea * (
                    ci_urea[t] - cm_urea) + smi_nh3 * (
                    ci_nh3[t] - cm_nh3) + smi_nh4 * (
                    ci_nh4[t] - cm_nh4) + smi_h * (
                    ci_h[t] - cm_h) + smi_hco2 * (
                    ci_hco2[t] - cm_hco2) + smi_h2co2 * (
                    ci_h2co2[t] - cm_h2co2) + smi_gluc * (
                    ci_gluc[t] - cm_gluc))

    fivs = fivs + ais * lpis * (
            sis_na * (cs_na - ci_na[t]) + sis_k * (
                    cs_k - ci_k[t]) + sis_cl * (
                    cs_cl - ci_cl[t]) + sis_hco3 * (
                    cs_hco3 - ci_hco3[t]) + sis_h2co3 * (
                    cs_h2co3 - ci_h2co3[t]) + sis_co2 * (
                    cs_co2 - ci_co2[t]) + sis_hpo4 * (
                    cs_hpo4 - ci_hpo4[t]) + sis_h2po4 * (
                    cs_h2po4 - ci_h2po4[t]) + sis_urea * (
                    cs_urea - ci_urea[t]) + sis_nh3 * (
                    cs_nh3 - ci_nh3[t]) + sis_nh4 * (
                    cs_nh4 - ci_nh4[t]) + sis_h * (
                    cs_h - ci_h[t]) + sis_hco2 * (
                    cs_hco2 - ci_hco2[t]) + sis_h2co2 * (
                    cs_h2co2 - ci_h2co2[t]) + sis_gluc * (
                    cs_gluc - ci_gluc[t]))

    five = five + lpie * aie * (
            sis_na * (ce_na[t] - ci_na[t]) + sis_k * (
                    ce_k[t] - ci_k[t]) + sis_cl * (
                    ce_cl[t] - ci_cl[t]) + sis_hco3 * (
                    ce_hco3[t] - ci_hco3[t]) + sis_h2co3 * (
                    ce_h2co3[t] - ci_h2co3[t]) + sis_co2 * (
                    ce_co2[t] - ci_co2[t]) + sis_hpo4 * (
                    ce_hpo4[t] - ci_hpo4[t]) + sis_h2po4 * (
                    ce_h2po4[t] - ci_h2po4[t]) + sis_urea * (
                    ce_urea[t] - ci_urea[t]) + sis_nh3 * (
                    ce_nh3[t] - ci_nh3[t]) + sis_nh4 * (
                    ce_nh4[t] - ci_nh4[t]) + sis_h * (
                    ce_h[t] - ci_h[t]) + sis_hco2 * (
                    ce_hco2[t] - ci_hco2[t]) + sis_h2co2 * (
                    ce_h2co2[t] - ci_h2co2[t]) + sis_gluc * (
                    ce_gluc[t] - ci_gluc[t]))

    Cellular_water_fluxes = fivm - fivs - five

    # Convective Fluxes, See Eqs: (37), (38)
    param_csf = 1
    fekm_na = csf(ce_na[t], cm_na, fevm, sme_na, param_csf)
    fekm_na_csf = csf(ce_na[t], cm_na, fevm, sme_na, param_csf)
    fekm_k = csf(ce_k[t], cm_k, fevm, sme_k, param_csf)
    fekm_cl = csf(ce_cl[t], cm_cl, fevm, sme_cl, param_csf)
    fekm_hco3 = csf(ce_hco3[t], cm_hco3, fevm, sme_hco3, param_csf)
    fekm_h2co3 = csf(ce_h2co3[t], cm_h2co3, fevm, sme_h2co3, param_csf)
    fekm_co2 = csf(ce_co2[t], cm_co2, fevm, sme_co2, param_csf)
    fekm_hpo4 = csf(ce_hpo4[t], cm_hpo4, fevm, sme_hpo4, param_csf)
    fekm_h2po4 = csf(ce_h2po4[t], cm_h2po4, fevm, sme_h2po4, param_csf)
    fekm_urea = csf(ce_urea[t], cm_urea, fevm, sme_urea, param_csf)
    fekm_nh3 = csf(ce_nh3[t], cm_nh3, fevm, sme_nh3, param_csf)
    fekm_nh4 = csf(ce_nh4[t], cm_nh4, fevm, sme_nh4, param_csf)
    fekm_h = csf(ce_h[t], cm_h, fevm, sme_h, param_csf)
    fekm_hco2 = csf(ce_hco2[t], cm_hco2, fevm, sme_hco2, param_csf)
    fekm_h2co2 = csf(ce_h2co2[t], cm_h2co2, fevm, sme_h2co2, param_csf)
    fekm_gluc = csf(ce_gluc[t], cm_gluc, fevm, sme_gluc, param_csf)

    feks_na = csf(ce_na[t], cs_na, fevs, ses_na, param_csf)
    feks_na_csf = csf(ce_na[t], cs_na, fevs, ses_na, param_csf)
    feks_k = csf(ce_k[t], cs_k, fevs, ses_k, param_csf)
    feks_cl = csf(ce_cl[t], cs_cl, fevs, ses_cl, param_csf)
    feks_hco3 = csf(ce_hco3[t], cs_hco3, fevs, ses_hco3, param_csf)
    feks_h2co3 = csf(ce_h2co3[t], cs_h2co3, fevs, ses_h2co3, param_csf)
    feks_co2 = csf(ce_co2[t], cs_co2, fevs, ses_co2, param_csf)
    feks_hpo4 = csf(ce_hpo4[t], cs_hpo4, fevs, ses_hpo4, param_csf)
    feks_h2po4 = csf(ce_h2po4[t], cs_h2po4, fevs, ses_h2po4, param_csf)
    feks_urea = csf(ce_urea[t], cs_urea, fevs, ses_urea, param_csf)
    feks_nh3 = csf(ce_nh3[t], cs_nh3, fevs, ses_nh3, param_csf)
    feks_nh4 = csf(ce_nh4[t], cs_nh4, fevs, ses_nh4, param_csf)
    feks_h = csf(ce_h[t], cs_h, fevs, ses_h, param_csf)
    feks_hco2 = csf(ce_hco2[t], cs_hco2, fevs, ses_hco2, param_csf)
    feks_h2co2 = csf(ce_h2co2[t], cs_h2co2, fevs, ses_h2co2, param_csf)
    feks_gluc = csf(ce_gluc[t], cs_gluc, fevs, ses_gluc, param_csf)

    fikm_na = csf(ci_na[t], cm_na, fivm, smi_na, param_csf)
    fikm_na_csf = csf(ci_na[t], cm_na, fivm, smi_na, param_csf)
    fikm_k = csf(ci_k[t], cm_k, fivm, smi_k, param_csf)
    fikm_cl = csf(ci_cl[t], cm_cl, fivm, smi_cl, param_csf)
    fikm_hco3 = csf(ci_hco3[t], cm_hco3, fivm, smi_hco3, param_csf)
    fikm_h2co3 = csf(ci_h2co3[t], cm_h2co3, fivm, smi_h2co3, param_csf)
    fikm_co2 = csf(ci_co2[t], cm_co2, fivm, smi_co2, param_csf)
    fikm_hpo4 = csf(ci_hpo4[t], cm_hpo4, fivm, smi_hpo4, param_csf)
    fikm_h2po4 = csf(ci_h2po4[t], cm_h2po4, fivm, smi_h2po4, param_csf)
    fikm_urea = csf(ci_urea[t], cm_urea, fivm, smi_urea, param_csf)
    fikm_nh3 = csf(ci_nh3[t], cm_nh3, fivm, smi_nh3, param_csf)
    fikm_nh4 = csf(ci_nh4[t], cm_nh4, fivm, smi_nh4, param_csf)
    fikm_h = csf(ci_h[t], cm_h, fivm, smi_h, param_csf)
    fikm_hco2 = csf(ci_hco2[t], cm_hco2, fivm, smi_hco2, param_csf)
    fikm_h2co2 = csf(ci_h2co2[t], cm_h2co2, fivm, smi_h2co2, param_csf)
    fikm_gluc = csf(ci_gluc[t], cm_gluc, fivm, smi_gluc, param_csf)

    fiks_na = csf(ci_na[t], cs_na, fivs, sis_na, param_csf)
    fiks_na_csf = csf(ci_na[t], cs_na, fivs, sis_na, param_csf)
    fiks_k = csf(ci_k[t], cs_k, fivs, sis_k, param_csf)
    fiks_cl = csf(ci_cl[t], cs_cl, fivs, sis_cl, param_csf)
    fiks_hco3 = csf(ci_hco3[t], cs_hco3, fivs, sis_hco3, param_csf)
    fiks_h2co3 = csf(ci_h2co3[t], cs_h2co3, fivs, sis_h2co3, param_csf)
    fiks_co2 = csf(ci_co2[t], cs_co2, fivs, sis_co2, param_csf)
    fiks_hpo4 = csf(ci_hpo4[t], cs_hpo4, fivs, sis_hpo4, param_csf)
    fiks_h2po4 = csf(ci_h2po4[t], cs_h2po4, fivs, sis_h2po4, param_csf)
    fiks_urea = csf(ci_urea[t], cs_urea, fivs, sis_urea, param_csf)
    fiks_nh3 = csf(ci_nh3[t], cs_nh3, fivs, sis_nh3, param_csf)
    fiks_nh4 = csf(ci_nh4[t], cs_nh4, fivs, sis_nh4, param_csf)
    fiks_h = csf(ci_h[t], cs_h, fivs, sis_h, param_csf)
    fiks_hco2 = csf(ci_hco2[t], cs_hco2, fivs, sis_hco2, param_csf)
    fiks_h2co2 = csf(ci_h2co2[t], cs_h2co2, fivs, sis_h2co2, param_csf)
    fiks_gluc = csf(ci_gluc[t], cs_gluc, fivs, sis_gluc, param_csf)

    fike_na = csf(ci_na[t], ce_na[t], five, sis_na, param_csf)
    fike_na_csf = csf(ci_na[t], ce_na[t], five, sis_na, param_csf)
    fike_k = csf(ci_k[t], ce_k[t], five, sis_k, param_csf)
    fike_cl = csf(ci_cl[t], ce_cl[t], five, sis_cl, param_csf)
    fike_hco3 = csf(ci_hco3[t], ce_hco3[t], five, sis_hco3, param_csf)
    fike_h2co3 = csf(ci_h2co3[t], ce_h2co3[t], five, sis_h2co3, param_csf)
    fike_co2 = csf(ci_co2[t], ce_co2[t], five, sis_co2, param_csf)
    fike_hpo4 = csf(ci_hpo4[t], ce_hpo4[t], five, sis_hpo4, param_csf)
    fike_h2po4 = csf(ci_h2po4[t], ce_h2po4[t], five, sis_h2po4, param_csf)
    fike_urea = csf(ci_urea[t], ce_urea[t], five, sis_urea, param_csf)
    fike_nh3 = csf(ci_nh3[t], ce_nh3[t], five, sis_nh3, param_csf)
    fike_nh4 = csf(ci_nh4[t], ce_nh4[t], five, sis_nh4, param_csf)
    fike_h = csf(ci_h[t], ce_h[t], five, sis_h, param_csf)
    fike_hco2 = csf(ci_hco2[t], ce_hco2[t], five, sis_hco2, param_csf)
    fike_h2co2 = csf(ci_h2co2[t], ce_h2co2[t], five, sis_h2co2, param_csf)
    fike_gluc = csf(ci_gluc[t], ce_gluc[t], five, sis_gluc, param_csf)

    # Goldman Fluxes, See Eqs: (39), (40), (41)
    param_goldman = 1
    fekm_na = fekm_na + goldman(hme_na, ame, z_na, vm[t], ve[t], cm_na,
                                ce_na[t], param_goldman)
    fekm_na_goldman = goldman(hme_na, ame, z_na, vm[t], ve[t], cm_na,
                      ce_na[t], param_goldman)
    fekm_k = fekm_k + goldman(hme_k, ame, z_k, vm[t], ve[t], cm_k, ce_k[t], param_goldman)
    fekm_cl = fekm_cl + goldman(hme_cl, ame, z_cl, vm[t], ve[t], cm_cl, ce_cl[t], param_goldman)
    fekm_hco3 = fekm_hco3 + goldman(hme_hco3, ame, z_hco3, vm[t], ve[t], cm_hco3, ce_hco3[t], param_goldman)
    fekm_h2co3 = fekm_h2co3 + goldman(hme_h2co3, ame, z_h2co3, vm[t], ve[t], cm_h2co3,
                                      ce_h2co3[t], param_goldman)
    fekm_co2 = fekm_co2 + goldman(hme_co2, ame, z_co2, vm[t], ve[t], cm_co2, ce_co2[t], param_goldman)
    fekm_hpo4 = fekm_hpo4 + goldman(hme_hpo4, ame, z_hpo4, vm[t], ve[t], cm_hpo4,
                                    ce_hpo4[t], param_goldman)
    fekm_h2po4 = fekm_h2po4 + goldman(hme_h2po4, ame, z_h2po4, vm[t], ve[t], cm_h2po4,
                                      ce_h2po4[t], param_goldman)
    fekm_urea = fekm_urea + goldman(hme_urea, ame, z_urea, vm[t], ve[t], cm_urea, ce_urea[t], param_goldman)
    fekm_nh3 = fekm_nh3 + goldman(hme_nh3, ame, z_nh3, vm[t], ve[t], cm_nh3, ce_nh3[t], param_goldman)
    fekm_nh4 = fekm_nh4 + goldman(hme_nh4, ame, z_nh4, vm[t], ve[t], cm_nh4, ce_nh4[t], param_goldman)
    fekm_h = fekm_h + goldman(hme_h, ame, z_h, vm[t], ve[t], cm_h, ce_h[t], param_goldman)
    fekm_hco2 = fekm_hco2 + goldman(hme_hco2, ame, z_hco2, vm[t], ve[t], cm_hco2, ce_hco2[t], param_goldman)
    fekm_h2co2 = fekm_h2co2 + goldman(hme_h2co2, ame, z_h2co2, vm[t], ve[t], cm_h2co2, ce_h2co2[t],
                                      param_goldman)
    fekm_gluc = fekm_gluc + goldman(hme_gluc, ame, z_gluc, vm[t], ve[t], cm_gluc, ce_gluc[t], param_goldman)

    feks_na = feks_na + goldman(hes_na, ae[t], z_na, ve[t], vs, ce_na[t],
                                cs_na, param_goldman)
    feks_na_goldman = goldman(hes_na, ae[t], z_na, ve[t], vs, ce_na[t],
                                cs_na, param_goldman)
    feks_k = feks_k + goldman(hes_k, ae[t], z_k, ve[t], vs, ce_k[t], cs_k, param_goldman)
    feks_cl = feks_cl + goldman(hes_cl, ae[t], z_cl, ve[t], vs, ce_cl[t], cs_cl, param_goldman)
    feks_hco3 = feks_hco3 + goldman(hes_hco3, ae[t], z_hco3, ve[t], vs, ce_hco3[t], cs_hco3, param_goldman)
    feks_h2co3 = feks_h2co3 + goldman(hes_h2co3, ae[t], z_h2co3, ve[t], vs, ce_h2co3[t],
                                      cs_h2co3, param_goldman)
    feks_co2 = feks_co2 + goldman(hes_co2, ae[t], z_co2, ve[t], vs, ce_co2[t], cs_co2, param_goldman)
    feks_hpo4 = feks_hpo4 + goldman(hes_hpo4, ae[t], z_hpo4, ve[t], vs, ce_hpo4[t],
                                    cs_hpo4, param_goldman)
    feks_h2po4 = feks_h2po4 + goldman(hes_h2po4, ae[t], z_h2po4, ve[t], vs, ce_h2po4[t],
                                      cs_h2po4, param_goldman)
    feks_urea = feks_urea + goldman(hes_urea, ae[t], z_urea, ve[t], vs, ce_urea[t], cs_urea, param_goldman)
    feks_nh3 = feks_nh3 + goldman(hes_nh3, ae[t], z_nh3, ve[t], vs, ce_nh3[t], cs_nh3, param_goldman)
    feks_nh4 = feks_nh4 + goldman(hes_nh4, ae[t], z_nh4, ve[t], vs, ce_nh4[t], cs_nh4, param_goldman)
    feks_h = feks_h + goldman(hes_h, ae[t], z_h, ve[t], vs, ce_h[t], cs_h, param_goldman)
    feks_hco2 = feks_hco2 + goldman(hes_hco2, ae[t], z_hco2, ve[t], vs, ce_hco2[t], cs_hco2, param_goldman)

    feks_h2co2 = feks_h2co2 + goldman(hme_h2co2, ame, z_h2co2, ve[t], vs, ce_h2co2[t], cs_h2co2, param_goldman)
    feks_gluc = feks_gluc + goldman(hes_gluc, ae[t], z_gluc, ve[t], vs, ce_gluc[t], cs_gluc, param_goldman)

    fikm_na = fikm_na + goldman(hmi_na, ami, z_na, vm[t], vi[t], cm_na,
                                ci_na[t], param_goldman)
    fikm_na_goldman = goldman(hmi_na, ami, z_na, vm[t], vi[t], cm_na,
                      ci_na[t], param_goldman)
    fikm_k = fikm_k + goldman(hmi_k, ami, z_k, vm[t], vi[t], cm_k, ci_k[t], param_goldman)
    fikm_cl = fikm_cl + goldman(hmi_cl, ami, z_cl, vm[t], vi[t], cm_cl, ci_cl[t], param_goldman)
    fikm_hco3 = fikm_hco3 + goldman(hmi_hco3, ami, z_hco3, vm[t], vi[t], cm_hco3, ci_hco3[t], param_goldman)
    fikm_h2co3 = fikm_h2co3 + goldman(hmi_h2co3, ami, z_h2co3, vm[t], vi[t], cm_h2co3,
                                      ci_h2co3[t], param_goldman)
    fikm_co2 = fikm_co2 + goldman(hmi_co2, ami, z_co2, vm[t], vi[t], cm_co2, ci_co2[t], param_goldman)
    fikm_hpo4 = fikm_hpo4 + goldman(hmi_hpo4, ami, z_hpo4, vm[t], vi[t], cm_hpo4,
                                    ci_hpo4[t], param_goldman)
    fikm_h2po4 = fikm_h2po4 + goldman(hmi_h2po4, ami, z_h2po4, vm[t], vi[t], cm_h2po4,
                                      ci_h2po4[t], param_goldman)
    fikm_urea = fikm_urea + goldman(hmi_urea, ami, z_urea, vm[t], vi[t], cm_urea, ci_urea[t], param_goldman)
    fikm_nh3 = fikm_nh3 + goldman(hmi_nh3, ami, z_nh3, vm[t], vi[t], cm_nh3, ci_nh3[t], param_goldman)
    fikm_nh4 = fikm_nh4 + goldman(hmi_nh4, ami, z_nh4, vm[t], vi[t], cm_nh4, ci_nh4[t], param_goldman)
    fikm_h = fikm_h + goldman(hmi_h, ami, z_h, vm[t], vi[t], cm_h, ci_h[t], param_goldman)
    fikm_hco2 = fikm_hco2 + goldman(hmi_hco2, ami, z_hco2, vm[t], vi[t], cm_hco2, ci_hco2[t], param_goldman)
    fikm_h2co2 = fikm_h2co2 + goldman(hmi_h2co2, ami, z_h2co2, vm[t], vi[t], cm_h2co2,
                                      ci_h2co2[t], param_goldman)
    fikm_gluc = fikm_gluc + goldman(hmi_gluc, ami, z_gluc, vm[t], vi[t], cm_gluc, ci_gluc[t], param_goldman)

    fike_na = fike_na + goldman(his_na, aie, z_na, vi[t], ve[t], ci_na[t], ce_na[t], param_goldman)
    fike_na_goldman = goldman(his_na, aie, z_na, vi[t], ve[t], ci_na[t], ce_na[t], param_goldman)
    fike_k = fike_k + goldman(his_k, aie, z_k, vi[t], ve[t], ci_k[t], ce_k[t], param_goldman)
    fike_cl = fike_cl + goldman(his_cl, aie, z_cl, vi[t], ve[t], ci_cl[t], ce_cl[t], param_goldman)
    fike_hco3 = fike_hco3 + goldman(his_hco3, aie, z_hco3, vi[t], ve[t], ci_hco3[t], ce_hco3[t], param_goldman)
    fike_h2co3 = fike_h2co3 + goldman(his_h2co3, aie, z_h2co3, vi[t], ve[t], ci_h2co3[t],
                                  ce_h2co3[t], param_goldman)
    fike_co2 = fike_co2 + goldman(his_co2, aie, z_co2, vi[t], ve[t], ci_co2[t], ce_co2[t], param_goldman)
    fike_hpo4 = fike_hpo4 + goldman(his_hpo4, aie, z_hpo4, vi[t], ve[t], ci_hpo4[t],
                                ce_hpo4[t], param_goldman)
    fike_h2po4 = fike_h2po4 + goldman(his_h2po4, aie, z_h2po4, vi[t], ve[t], ci_h2po4[t],
                                  ce_h2po4[t], param_goldman)
    fike_urea = fike_urea + goldman(his_urea, aie, z_urea, vi[t], ve[t], ci_urea[t], ce_urea[t], param_goldman)
    fike_nh3 = fike_nh3 + goldman(his_nh3, aie, z_nh3, vi[t], ve[t], ci_nh3[t], ce_nh3[t], param_goldman)
    fike_nh4 = fike_nh4 + goldman(his_nh4, aie, z_nh4, vi[t], ve[t], ci_nh4[t], ce_nh4[t], param_goldman)
    fike_h = fike_h + goldman(his_h, aie, z_h, vi[t], ve[t], ci_h[t], ce_h[t], param_goldman)
    fike_hco2 = fike_hco2 + goldman(his_hco2, aie, z_hco2, vi[t], ve[t], ci_hco2[t], ce_hco2[t], param_goldman)
    fike_h2co2 = fike_h2co2 + goldman(his_h2co2, aie, z_h2co2, vi[t], ve[t], ci_h2co2[t], ce_h2co2[t],
                                  param_goldman)
    fike_gluc = fike_gluc + goldman(his_gluc, aie, z_gluc, vi[t], ve[t], ci_gluc[t], ce_gluc[t], param_goldman)

    fiks_na = fiks_na + goldman(his_na, ais, z_na, vi[t], vs, ci_na[t], cs_na, param_goldman)
    fiks_na_goldman = goldman(his_na, ais, z_na, vi[t], vs, ci_na[t], cs_na, param_goldman)
    fiks_k = fiks_k + goldman(his_k, ais, z_k, vi[t], vs, ci_k[t], cs_k, param_goldman)
    fiks_cl = fiks_cl + goldman(his_cl, ais, z_cl, vi[t], vs, ci_cl[t], cs_cl, param_goldman)
    fiks_hco3 = fiks_hco3 + goldman(his_hco3, ais, z_hco3, vi[t], vs, ci_hco3[t], cs_hco3, param_goldman)
    fiks_h2co3 = fiks_h2co3 + goldman(his_h2co3, ais, z_h2co3, vi[t], vs, ci_h2co3[t],
                                      cs_h2co3, param_goldman)
    fiks_co2 = fiks_co2 + goldman(his_co2, ais, z_co2, vi[t], vs, ci_co2[t], cs_co2, param_goldman)
    fiks_hpo4 = fiks_hpo4 + goldman(his_hpo4, ais, z_hpo4, vi[t], vs, ci_hpo4[t],
                                    cs_hpo4, param_goldman)
    fiks_h2po4 = fiks_h2po4 + goldman(his_h2po4, ais, z_h2po4, vi[t], vs, ci_h2po4[t],
                                      cs_h2po4, param_goldman)
    fiks_urea = fiks_urea + goldman(his_urea, ais, z_urea, vi[t], vs, ci_urea[t], cs_urea, param_goldman)
    fiks_nh3 = fiks_nh3 + goldman(his_nh3, ais, z_nh3, vi[t], vs, ci_nh3[t], cs_nh3, param_goldman)
    fiks_nh4 = fiks_nh4 + goldman(his_nh4, ais, z_nh4, vi[t], vs, ci_nh4[t], cs_nh4, param_goldman)
    fiks_h = fiks_h + goldman(his_h, ais, z_h, vi[t], vs, ci_h[t], cs_h, param_goldman)
    fiks_hco2 = fiks_hco2 + goldman(his_hco2, ais, z_hco2, vi[t], vs, ci_hco2[t], cs_hco2, param_goldman)
    fiks_h2co2 = fiks_h2co2 + goldman(his_h2co2, ais, z_h2co2, vi[t], vs, ci_h2co2[t],
                                      cs_h2co2, param_goldman)
    fiks_gluc = fiks_gluc + goldman(his_gluc, ais, z_gluc, vi[t], vs, ci_gluc[t], cs_gluc, param_goldman)

    # net transporters on mi bourder
    sglt = sglt_mi(cm_na, ci_na[t], cm_gluc, ci_gluc[t], z_na, z_gluc, vm[t], vi[t], ami, lmi_nagluc,
                   param_sglt)
    na_mi_nagluc = sglt[0]
    gluc_mi_nagluc = sglt[1]

    nah2po4 = nah2po4_mi(cm_na, ci_na[t], cm_h2po4, ci_h2po4[t], z_na, z_h2po4, vm[t], vi[t], ami,
                         lmi_nah2po4, param_h2po4)
    na_mi_nah2po4 = nah2po4[0]
    h2po4_mi_nah2po4 = nah2po4[1]

    clhco3 = clhco3_mi(cm_cl, ci_cl[t], cm_hco3, ci_hco3[t], z_cl, z_hco3, vm[t], vi[t], ami,
                       lmi_clhco3, param_clhco3)
    cl_mi_clhco3 = clhco3[0]
    hco3_mi_clhco3 = clhco3[1]

    clhco2 = clhco2_mi(cm_cl, ci_cl[t], cm_hco2, ci_hco2[t], z_cl, z_hco2, vm[t], vi[t], ami,
                       lmi_clhco2, param_clhco2)
    cl_mi_clhco2 = clhco2[0]
    hco2_mi_clhco2 = clhco2[1]

    nah_mi = nah(cm_na, ci_na[t], cm_h, ci_h[t], z_na, z_h, vm[t], vi[t], ami, lmi_nah, param_nah_mi)
    na_mi_nah = nah_mi[0]
    h_mi_nah = nah_mi[1]

    nanh4_mi = nanh4(cm_na, ci_na[t], cm_nh4, ci_nh4[t], z_na, z_nh4, vm[t], vi[t], ami, lmi_nanh4, param_nanh4_mi)
    na_mi_nanh4 = nanh4_mi[0]
    nh4_mi_nanh4 = nanh4_mi[1]

    mynah = nhe3(ci_h[t], ci_na[t], ci_nh4[t], cm_h, cm_na, cm_nh4, param_nhe3)
    jnah_na = mynah[0]
    jnah_h = mynah[1]
    jnah_nh4 = mynah[2]
    jnhe3_na = nnhe3 * ami * jnah_na
    jnhe3_h = nnhe3 * ami * jnah_h
    jnhe3_nh4 = nnhe3 * ami * jnah_nh4

    fikm_na = fikm_na + na_mi_nagluc + na_mi_nah2po4 + jnhe3_na + na_mi_nanh4 + na_mi_nah
    fikm_cl = fikm_cl + cl_mi_clhco2 + cl_mi_clhco3
    fikm_hco3 = fikm_hco3 + hco3_mi_clhco3
    fikm_h2po4 = fikm_h2po4 + h2po4_mi_nah2po4
    fikm_hco2 = fikm_hco2 + hco2_mi_clhco2
    fikm_gluc = fikm_gluc + gluc_mi_nagluc
    fikm_h = fikm_h + jnhe3_h + h_mi_nah
    fikm_nh4 = fikm_nh4 + jnhe3_nh4 + nh4_mi_nanh4

    nahco3 = na_hco3(ci_na[t], cs_na, ci_hco3[t], cs_hco3, z_na, z_hco3, vi[t], vs, ais,
                     lis_nahco3, param_nahco3)
    na_is_nahco3 = nahco3[0]
    # print('na_is_nahco3', na_is_nahco3)
    hco3_is_nahco3 = nahco3[1]
    # print('hco3_is_nahco3', hco3_is_nahco3)
    kcl = k_cl(ci_k[t], cs_k, ci_cl[t], cs_cl, z_k, z_cl, vi[t], vs, ais, lis_kcl, param_kcl)
    k_is_kcl = kcl[0]
    cl_is_kcl = kcl[1]

    na_clhco3 = na_cl_hco3(ci_na[t], cs_na, ci_cl[t], cs_cl, ci_hco3[t], cs_hco3, z_na, z_cl, z_hco3,
                           vi[t], vs, ais, lis_na_clhco3, param_na_cl_hco3)
    na_is_na_clhco3 = na_clhco3[0]
    cl_is_na_clhco3 = na_clhco3[1]
    hco3_is_na_clhco3 = na_clhco3[2]
    #fiks_na = fiks_na  + na_is_na_clhco3
    fiks_na = fiks_na + na_is_nahco3 + na_is_na_clhco3
    fiks_k = fiks_k + k_is_kcl
    fiks_cl = fiks_cl + cl_is_kcl + cl_is_na_clhco3
    fiks_hco3 = fiks_hco3 + hco3_is_na_clhco3 + hco3_is_nahco3

    nahco3_aie = na_hco3(ci_na[t], ce_na[t], ci_hco3[t], ce_hco3[t], z_na, z_hco3, vi[t], ve[t],
                           aie,
                           lis_nahco3, param_nahco3)
    na_aie_nahco3 = nahco3_aie[0]
    hco3_aie_nahco3 = nahco3_aie[1]

    kcl_aie = k_cl(ci_k[t], ce_k[t], ci_cl[t], ce_cl[t], z_k, z_cl, vi[t], ve[t], aie, lis_kcl, param_kcl)
    k_aie_kcl = kcl_aie[0]
    cl_aie_kcl = kcl_aie[1]

    na_clhco3_aie = na_cl_hco3(ci_na[t], ce_na[t], ci_cl[t], ce_cl[t], ci_hco3[t], ce_hco3[t], z_na,
                                 z_cl, z_hco3,
                                 vi[t], ve[t], aie, lis_na_clhco3, param_na_cl_hco3)
    na_aie_na_clhco3 = na_clhco3_aie[0]
    cl_aie_na_clhco3 = na_clhco3_aie[1]
    hco3_aie_na_clhco3 = na_clhco3_aie[2]
    fike_na = fike_na + na_aie_nahco3 + na_aie_na_clhco3
    #fike_na = fike_na + na_aie_na_clhco3
    fike_k = fike_k + k_aie_kcl
    fike_cl = fike_cl + cl_aie_kcl + cl_aie_na_clhco3
    fike_hco3 = fike_hco3 + hco3_aie_na_clhco3 + hco3_aie_nahco3
    # fike_hco3 = fike_hco3 + hco3_aie_na_clhco3
    nak = nak_atp(ci_na[t], cs_na, ci_k[t], cs_k,  cs_nh4, param_sodium_pumps)
    atis_na = nak[0]
    atis_k = nak[1]
    atis_nh4 = nak[2]
    fiks_na = fiks_na + ais * atis_na
    fiks_k = fiks_k + ais * atis_k
    fiks_nh4 = fiks_nh4 + ais * atis_nh4

    nak = nak_atp(ci_na[t], ce_na[t], ci_k[t], ce_k[t],  ce_nh4[t], param_sodium_pumps)
    atie_na = nak[0]
    atie_k = nak[1]
    atie_nh4 = nak[2]
    fike_na = fike_na + aie * atie_na
    fike_k = fike_k + aie * atie_k
    fike_nh4 = fike_nh4 + aie * atie_nh4

    atmi_h = h_atp_mi(cm_h, ci_h[t], vm[t], vi[t], z_h, ami, param_H_pumps)
    fikm_h = fikm_h + atmi_h
    # See Eq:(30), (31)
    phie_en = 0
    phie_en = phie_en + z_na * ce_na[t] + z_k * ce_k[t] + z_cl * ce_cl[
        t] + z_hco3 * ce_hco3[t] + z_h2co3 * ce_h2co3[t] + z_co2 * ce_co2[t] + z_hpo4 * ce_hpo4[
                  t] + z_h2po4 * ce_h2po4[t] + z_urea * ce_urea[t] + z_nh3 * ce_nh3[t] + z_nh4 * ce_nh4[
                  t] + z_h * ce_h[t] + z_hco2 * ce_hco2[t] + z_h2co2 * ce_h2co2[t] + z_gluc * ce_gluc[t]

    phii_en = imp[t] * zimp
    phii_en = phii_en - cbuf[t] + z_na * ci_na[t] + z_k * ci_k[t] \
              + z_cl * ci_cl[t] + z_hco3 * ci_hco3[t] + z_h2co3 * ci_h2co3[t] + z_co2 * ci_co2[t] + z_hpo4 * \
              ci_hpo4[t] + z_h2po4 * ci_h2po4[t] + z_urea * ci_urea[t] + z_nh3 * ci_nh3[t] + z_nh4 * \
              ci_nh4[t] + z_h * ci_h[t] + z_hco2 * ci_hco2[t] + z_h2co2 * ci_h2co2[t] + z_gluc * ci_gluc[
                  t]

    # Mass conservation, See Eqs: (6), (7), (8), (9)
    phie_vlm = fevs - fevm - five + rtau * (chvl[t] - chvl[t - 1])
    qe_na = feks_na - fekm_na - fike_na + rtau * (ce_na[t] * chvl[t] - ce_na[t - 1] * chvl[t - 1])
    qe_k = feks_k - fekm_k - fike_k + rtau * (ce_k[t] * chvl[t] - ce_k[t - 1] * chvl[t - 1])
    qe_cl = feks_cl - fekm_cl - fike_cl + rtau * (ce_cl[t] * chvl[t] - ce_cl[t - 1] * chvl[t - 1])
    qe_hco3 = feks_hco3 - fekm_hco3 - fike_hco3 + rtau * (ce_hco3[t] * chvl[t] - ce_hco3[t - 1] * chvl[t - 1])
    qe_h2co3 = feks_h2co3 - fekm_h2co3 - fike_h2co3 + rtau * (
            ce_h2co3[t] * chvl[t] - ce_h2co3[t - 1] * chvl[t - 1])
    qe_co2 = feks_co2 - fekm_co2 - fike_co2 + rtau * (ce_co2[t] * chvl[t] - ce_co2[t - 1] * chvl[t - 1])
    qe_hpo4 = feks_hpo4 - fekm_hpo4 - fike_hpo4 + rtau * (ce_hpo4[t] * chvl[t] - ce_hpo4[t - 1] * chvl[t - 1])
    qe_h2po4 = feks_h2po4 - fekm_h2po4 - fike_h2po4 + rtau * (
            ce_h2po4[t] * chvl[t] - ce_h2po4[t - 1] * chvl[t - 1])
    qe_urea = feks_urea - fekm_urea - fike_urea + rtau * (
            ce_urea[t] * chvl[t] - ce_urea[t - 1] * chvl[t - 1])
    qe_nh3 = feks_nh3 - fekm_nh3 - fike_nh3 + rtau * (ce_nh3[t] * chvl[t] - ce_nh3[t - 1] * chvl[t - 1])
    qe_nh4 = feks_nh4 - fekm_nh4 - fike_nh4 + rtau * (ce_nh4[t] * chvl[t] - ce_nh4[t - 1] * chvl[t - 1])
    qe_h = feks_h - fekm_h - fike_h + rtau * (ce_h[t] * chvl[t] - ce_h[t - 1] * chvl[t - 1])
    qe_hco2 = feks_hco2 - fekm_hco2 - fike_hco2 + rtau * (ce_hco2[t] * chvl[t] - ce_hco2[t - 1] * chvl[t - 1])
    qe_h2co2 = feks_h2co2 - fekm_h2co2 - fike_h2co2 + rtau * (
            ce_h2co2[t] * chvl[t] - ce_h2co2[t - 1] * chvl[t - 1])

    qe_gluc = feks_gluc - fekm_gluc - fike_gluc + rtau * (
            ce_gluc[t] * chvl[t] - ce_gluc[t - 1] * chvl[t - 1])

    phii_vlm = fivs - fivm + five + rtau * (clvl[t] - clvl[t - 1])

    qi_na = fiks_na - fikm_na + fike_na + rtau * (
            ci_na[t] * clvl[t] - ci_na[t - 1] * clvl[t - 1])

    qi_k = fiks_k - fikm_k + fike_k + rtau * (
            ci_k[t] * clvl[t] - ci_k[t - 1] * clvl[t - 1])
    qi_cl = fiks_cl - fikm_cl + fike_cl + rtau * (
            ci_cl[t] * clvl[t] - ci_cl[t - 1] * clvl[t - 1])

    qi_hco3 = fiks_hco3 - fikm_hco3 + fike_hco3 + rtau * (
            ci_hco3[t] * clvl[t] - ci_hco3[t - 1] * clvl[t - 1])
    qi_h2co3 = fiks_h2co3 - fikm_h2co3 + fike_h2co3 + rtau * (
            ci_h2co3[t] * clvl[t] - ci_h2co3[t - 1] * clvl[t - 1])
    qi_co2 = fiks_co2 - fikm_co2 + fike_co2 + rtau * (
            ci_co2[t] * clvl[t] - ci_co2[t - 1] * clvl[t - 1])
    qi_hpo4 = fiks_hpo4 - fikm_hpo4 + fike_hpo4 + rtau * (
            ci_hpo4[t] * clvl[t] - ci_hpo4[t - 1] * clvl[t - 1])
    qi_h2po4 = fiks_h2po4 - fikm_h2po4 + fike_h2po4 + rtau * (
            ci_h2po4[t] * clvl[t] - ci_h2po4[t - 1] * clvl[t - 1])

    qi_urea = fiks_urea - fikm_urea + fike_urea + rtau * (
            ci_urea[t] * clvl[t] - ci_urea[t - 1] * clvl[t - 1])
    qi_nh3 = fiks_nh3 - fikm_nh3 + fike_nh3 + rtau * (
            ci_nh3[t] * clvl[t] - ci_nh3[t - 1] * clvl[t - 1])
    qi_nh4 = fiks_nh4 - fikm_nh4 + fike_nh4 + rtau * (
            ci_nh4[t] * clvl[t] - ci_nh4[t - 1] * clvl[t - 1])
    qi_h = fiks_h - fikm_h + fike_h + rtau * (
            ci_h[t] * clvl[t] - ci_h[t - 1] * clvl[t - 1])
    qi_hco2 = fiks_hco2 - fikm_hco2 + fike_hco2 + rtau * (
            ci_hco2[t] * clvl[t] - ci_hco2[t - 1] * clvl[t - 1])
    qi_h2co2 = fiks_h2co2 - fikm_h2co2 + fike_h2co2 + rtau * (
            ci_h2co2[t] * clvl[t] - ci_h2co2[t - 1] * clvl[t - 1])
    qi_gluc = fiks_gluc - fikm_gluc + fike_gluc + rtau * (
            ci_gluc[t] * clvl[t] - ci_gluc[t - 1] * clvl[t - 1])

    qi_h = qi_h - rtau * (cbuf[t] * clvl[t] - cbuf[t - 1] * clvl[t - 1])
    scale = 1e6
    phie_vlm = phi_scale(phie_vlm, scale)
    qe_na = phi_scale(qe_na, scale)
    qe_k = phi_scale(qe_k, scale)
    qe_cl = phi_scale(qe_cl, scale)
    qe_hco3 = phi_scale(qe_hco3, scale)
    qe_h2co3 = phi_scale(qe_h2co3, scale)
    qe_co2 = phi_scale(qe_co2, scale)
    qe_hpo4 = phi_scale(qe_hpo4, scale)
    qe_h2po4 = phi_scale(qe_h2po4, scale)
    qe_urea = phi_scale(qe_urea, scale)
    qe_nh3 = phi_scale(qe_nh3, scale)
    qe_nh4 = phi_scale(qe_nh4, scale)
    qe_h = phi_scale(qe_h, scale)
    qe_hco2 = phi_scale(qe_hco2, scale)
    qe_h2co2 = phi_scale(qe_h2co2, scale)
    qe_gluc = phi_scale(qe_gluc, scale)

    phii_vlm = phi_scale(phii_vlm, scale)
    qi_na = phi_scale(qi_na, scale)
    qi_k = phi_scale(qi_k, scale)
    qi_cl = phi_scale(qi_cl, scale)
    qi_hco3 = phi_scale(qi_hco3, scale)
    qi_h2co3 = phi_scale(qi_h2co3, scale)
    qi_co2 = phi_scale(qi_co2, scale)
    qi_hpo4 = phi_scale(qi_hpo4, scale)
    qi_h2po4 = phi_scale(qi_h2po4, scale)
    qi_urea = phi_scale(qi_urea, scale)
    qi_nh3 = phi_scale(qi_nh3, scale)
    qi_nh4 = phi_scale(qi_nh4, scale)
    qi_h = phi_scale(qi_h, scale)
    qi_hco2 = phi_scale(qi_hco2, scale)
    qi_h2co2 = phi_scale(qi_h2co2, scale)
    qi_gluc = phi_scale(qi_gluc, scale)

    # solver = 1
    if solver == 1:
        phi[0] = phie_en
        phi[1] = phie_vlm
        phie_na = qe_na
        phi[2] = phie_na
        phie_k = qe_k
        phi[3] = phie_k
        phie_cl = qe_cl
        phi[4] = phie_cl
        phi[10] = qe_urea
        phie_gluc = qe_gluc
        phi[15] = phie_gluc
        phi[16] = phii_en
        phi[17] = phii_vlm
        phii_na = qi_na
        phi[18] = phii_na
        phii_k = qi_k
        phi[19] = phii_k
        phii_cl = qi_cl
        phi[20] = phii_cl
        phii_urea = qi_urea
        phi[26] = phii_urea
        phii_gluc = qi_gluc
        phi[31] = phii_gluc

        Buffer_Co2_LIS_Param = 1
        Buffer_Co2_Cell_Param = 1
        Buffer_HPO4_LIS_Param = 1
        Buffer_HPO4_CELL_Param = 1
        Buffer_NH3_LIS_Param = 1
        Buffer_NH3_CELL_Param = 1
        Buffer_HCO2_LIS_Param = 1
        Buffer_HCO2_CELL_Param = 1

        phi[8], phi[9] = buff_activation(qe_hpo4, qe_h2po4, lche, pkp, ce_hpo4[t], ce_h2po4[t],
                                               0, Buffer_HPO4_LIS_Param)
        phi[24], phi[25] = buff_activation(qi_hpo4, qi_h2po4, lchi, pkp, ci_hpo4[t], ci_h2po4[t],
                                                 0, Buffer_HPO4_CELL_Param)
        phi[11], phi[12] = buff_activation(qe_nh3, qe_nh4, lche, pkn, ce_nh3[t], ce_nh4[t],
                                                 0, Buffer_NH3_LIS_Param)
        phi[27], phi[28] = buff_activation(qi_nh3, qi_nh4, lchi, pkn, ci_nh3[t], ci_nh4[t],
                                                 1, Buffer_NH3_CELL_Param)
        phi[13], phi[14] = buff_activation(qe_hco2, qe_h2co2, lche, pkf, ce_hco2[t], ce_h2co2[t],
                                                 0, Buffer_HCO2_LIS_Param)
        phi[29], phi[30] = buff_activation(qi_hco2, qi_h2co2, lchi, pkf, ci_hco2[t], ci_h2co2[t],
                                                 0, Buffer_HCO2_CELL_Param)

        qi_cb = cbuf[t] + hcbuf[t] - (tbuf * clvl0 / clvl[t])
        phi[32] = qi_cb
        c_eQ_h2po4 = np.log10(abs(
            (ci_hpo4[t] * hcbuf[t]) / (ci_h2po4[t] * cbuf[t]))) if ci_h2po4[t] * cbuf[t] == 0 or (
                (ci_hpo4[t] * hcbuf[t]) / (ci_h2po4[t] * cbuf[t]) <= 0) else np.log10(
            (ci_hpo4[t] * hcbuf[t]) / (ci_h2po4[t] * cbuf[t]))
        phi_ph_eQ = pkb - pkp - c_eQ_h2po4
        phi[33] = phi_ph_eQ
        phi[5], phi[6], phi[7] = buff_activation_co2_formate_phosphate_ammonia(qe_h, qe_nh4, qe_hco3, qe_h2co2,
                                                                                        qe_h2po4, qe_co2, qe_h2co3, 0,
                                                                                        ce_co2[t],
                                                                                        ce_h2co3[t], chvl[t],
                                                                                        scale, 0,
                                                                                        Buffer_Co2_LIS_Param)
        phi[21], phi[22], phi[23] = buff_activation_co2_formate_phosphate_ammonia(qi_h, qi_nh4, qi_hco3,
                                                                                           qi_h2co2,
                                                                                           qi_h2po4, qi_co2, qi_h2co3,
                                                                                           qi_cb,
                                                                                           ci_co2[t],
                                                                                           ci_h2co3[t], clvl[t],
                                                                                           scale,
                                                                                           1,
                                                                                           Buffer_Co2_Cell_Param)
        # See Eqs (33), (34)
        cure = f * (
                z_na * fekm_na + z_k * fekm_k + z_cl * fekm_cl + z_hco3 * fekm_hco3 + z_h2co3 * fekm_h2co3 + z_co2
                * fekm_co2 + z_hpo4 * fekm_hpo4 + z_h2po4 * fekm_h2po4 + z_urea * fekm_urea + z_nh3 * fekm_nh3 + z_nh4
                * fekm_nh4 + z_h * fekm_h + z_hco2 * fekm_hco2 + z_h2co2 * fekm_h2co2 + z_gluc *
                fekm_gluc)

        curi = f * (
                z_na * fikm_na + z_k * fikm_k + z_cl * fikm_cl + z_hco3 * fikm_hco3 + z_h2co3 * fikm_h2co3 + z_co2
                * fikm_co2 + z_hpo4 * fikm_hpo4 + z_h2po4 * fikm_h2po4 + z_urea * fikm_urea + z_nh3 * fikm_nh3 + z_nh4 *
                fikm_nh4 + z_h * fikm_h + z_hco2 * fikm_hco2 + z_h2co2 * fikm_h2co2 + z_gluc *
                fikm_gluc)
        phie_cur = cure + curi
        phi[34] = phie_cur
        return phi
    else:
        flux = 1
        if flux:
            flux_na = fikm_na + fekm_na
            na_nak =  aie * atie_na
            flux_na_goldman = fekm_na_goldman + fikm_na_goldman
            hco3_clhco3 = clhco3[1]
    return[fekm_na, fekm_cl, fekm_k, fekm_gluc,
           fikm_na, fikm_cl, fikm_k, fikm_gluc,
           fiks_na, fiks_cl, fiks_k, fiks_gluc,
           fike_na, fike_cl, fike_k, fike_gluc,
           feks_na, feks_cl, feks_k, feks_gluc, fekm_na_csf, fekm_na_goldman, fikm_na, fikm_na_csf,
           fikm_na_goldman, jnhe3_na, na_mi_nagluc, na_mi_nah2po4, feks_na_csf, feks_na_goldman, flux_na,
           fekm_na_csf + fikm_na_csf, flux_na_goldman,
           jnhe3_na + na_mi_nagluc + na_mi_nah2po4, gluc_mi_nagluc, na_nak, hco3_clhco3, Cellular_water_fluxes,
           clvl_imp]


t0 = 0
tf = 2000
T = 20000
dt = (tf - t0) / (T - 1)
rtau = 1 / dt

ve = matrix(-0.894329e-02, T)
pe = matrix(-0.230937e+02, T)
ce_na = matrix(0.1404e+00, T)
ce_k = matrix(0.46537e-02, T)
ce_cl = matrix(0.1120e+00, T)
ce_hco3 = matrix(0.2560787e-01, T)
ce_h2co3 = matrix(0.4375464e-05, T)
ce_co2 = matrix(0.14987447089e-02, T)
ce_hpo4 = matrix(0.29852338e-02, T)
ce_h2po4 = matrix(0.86622197e-03, T)
ce_urea = matrix(0.490469598e-02, T)
ce_nh3 = matrix(0.269149195e-05, T)
ce_nh4 = matrix(0.17484103e-03, T)
ce_hco2 = matrix(0.77584319e-03, T)
ce_h2co2 = matrix(0.20531685e-06, T)
ce_gluc = matrix(0.772369e-02, T)
vi = matrix(-0.549457859e+02, T)
p_i = matrix(0.0e-01, T)
ci_na = matrix(0.204436152649e-01, T)
ci_k = matrix(0.1374233529e+00, T)
ci_cl = matrix(0.16901004e-01, T)
ci_hco3 = matrix(0.250905e-01, T)
ci_h2co3 = matrix(0.4375052e-05, T)
ci_co2 = matrix(0.1498841868e-02, T)
ci_hpo4 = matrix(0.94254445e-02, T)
ci_h2po4 = matrix(0.279109e-02, T)
ci_urea = matrix(0.4954689e-02, T)
ci_nh3 = matrix(0.352918888e-05, T)
ci_nh4 = matrix(0.23395e-03, T)
ci_hco2 = matrix(0.53552e-03, T)
ci_h2co2 = matrix(0.9337925e-07, T)
ci_gluc = matrix(0.1487617424e-01, T)
hcbuf = matrix(0.269599e-01, T)
cbuf = matrix(0.400120786e-01, T)
vm = matrix(-0.185736e+00, T)
ce_h = matrix(4.59e-8, T)
ci_h = matrix(4.69e-8, T)
ae = matrix(0.02000, T)


chvl = matrix(0.7000e-04, T)
clvl = matrix(0.1000e-02, T)

l = matrix(0.1000e-02, T)
imp = matrix(0.6000e-01, T)
rm = matrix(0.1250e-02, T)
am = matrix(0, T)
phi = matrix(0, 35)

lchm = pkp + np.log10(cm_hpo4 / cm_h2po4)
lchs = pkp + np.log10(cs_hpo4 / cs_h2po4)
cm_h = 10. ** (-lchm)
cs_h = 10. ** (-lchs)

clvl_imp = matrix(0, T)
flux_na = matrix(0, T)
# flux_k = matrix(0, T)
flux_cl = matrix(0, T)
cl_clhco3 = matrix(0, T)
cl_kcl = matrix(0, T)
cl_clhco2 = matrix(0, T)
cl_na_clhco3 = matrix(0, T)
na_nak = matrix(0, T)
k_nak = matrix(0, T)
na_nah = matrix(0, T)

fikm_na_csf = matrix(0, T)
fikm_na_goldman = matrix(0, T)
fekm_na_goldman = matrix(0, T)
fekm_na_csf=matrix(0, T)
flux_na_goldman = matrix(0, T)
feks_na_csf=matrix(0, T)
feks_na_goldman=matrix(0, T)
Cellular_water_fluxes = matrix(0, T)


fikm_na = matrix(0, T)
fikm_k = matrix(0, T)
fikm_cl = matrix(0, T)
fikm_gluc = matrix(0, T)

fiks_na = matrix(0, T)
fiks_k = matrix(0, T)
fiks_cl = matrix(0, T)
fiks_gluc = matrix(0, T)

fekm_na = matrix(0, T)
fekm_k = matrix(0, T)
fekm_cl = matrix(0, T)
fekm_gluc = matrix(0, T)

fike_na = matrix(0, T)
fike_k = matrix(0, T)
fike_cl = matrix(0, T)
fike_gluc = matrix(0, T)

feks_na = matrix(0, T)
feks_k = matrix(0, T)
feks_cl = matrix(0, T)
feks_gluc = matrix(0, T)
hco3_clhco3 = matrix(0, T)

na_mi_nagluc = matrix(0, T)
gluc_mi_nagluc = matrix(0, T)
na_mi_nah2po4 = matrix(0, T)
jnhe3_na = matrix(0, T)
cl_mi_clhco2 = matrix(0, T)
cl_mi_clhco3 = matrix(0, T)
vars = [ve, pe, ce_na, ce_k,
         ce_cl, ce_hco3, ce_h2co3, ce_co2, ce_hpo4, ce_h2po4, ce_urea, ce_nh3, ce_nh4, ce_hco2, ce_h2co2, ce_gluc,
         vi, imp, ci_na, ci_k, ci_cl, ci_hco3, ci_h2co3, ci_co2, ci_hpo4, ci_h2po4, ci_urea, ci_nh3, ci_nh4, ci_hco2,
         ci_h2co2, ci_gluc, cbuf, hcbuf, vm]

Figure_4a = 0
Figure_4b = 0
Figure_4c = 0
Figure_5 = 0
Figure_9_10 = 0
Figure_6_7_8 = 1

flx = [fekm_na, fekm_cl, fekm_k, fekm_gluc, fikm_na, fikm_cl, fikm_k, fikm_gluc,
      fiks_na, fiks_cl, fiks_k, fiks_gluc, fike_na, fike_cl, fike_k, fike_gluc,
      feks_na, feks_cl, feks_k, feks_gluc, fekm_na_csf, fekm_na_goldman, fikm_na, fikm_na_csf,
      fikm_na_goldman, jnhe3_na, na_mi_nagluc, na_mi_nah2po4, feks_na_csf, feks_na_goldman, flux_na,
      fekm_na_csf + fikm_na_csf, flux_na_goldman,
      jnhe3_na + na_mi_nagluc + na_mi_nah2po4, gluc_mi_nagluc, na_nak, hco3_clhco3, Cellular_water_fluxes,
      clvl_imp]
if Figure_9_10:
    from PCT_GLOB_new import *

    t = 0
    guess = [var[t] for var in vars]
    from scipy import optimize

    mylist = [0.0 for i in range(T)]
    pm = 15.000
    cm_cl = 0.11322499
    cm_na = 0.14
    cm_gluc = 0.005
    stepwise = 1
    if stepwise:
        cmna = [0.100, 0.110, 0.120, 0.13, 0.140]
        n_nhe3 = [0.27500e-08, 0.2000e-08, 0.1500e-08, 0.1000e-08, 0.01000e-08]
        while True:
            t += 1
            if t == T:
                break
            elif 0 < t <= T / len(cmna):
                nnhe3 = n_nhe3[0]
                mylist[t] = n_nhe3[0]
                result = scipy.optimize.root(eQs, np.array(guess), 1)
                guess = [var[t] for var in vars]
                result_flx = eQs(guess, 0)
                for i in range(len(result_flx)):
                    flx[i][t] = result_flx[i]
            elif 1 * T / len(cmna) < t <= 2 * T / len(cmna):
                nnhe3 = n_nhe3[1]
                mylist[t] = n_nhe3[1]
                result = scipy.optimize.root(eQs, np.array(guess), 1)
                guess = [var[t] for var in vars]
                result_flx = eQs(guess, 0)
                for i in range(len(result_flx)):
                    flx[i][t] = result_flx[i]
            elif 2 * T / len(cmna) < t <= 3 * T / len(cmna):
                nnhe3 = n_nhe3[2]
                mylist[t] = n_nhe3[2]
                result = scipy.optimize.root(eQs, np.array(guess), 1)
                guess = [var[t] for var in vars]
                result_flx = eQs(guess, 0)
                for i in range(len(result_flx)):
                    flx[i][t] = result_flx[i]
            elif 3 * T / len(cmna) < t <= 4 * T / len(cmna):
                nnhe3 = n_nhe3[3]
                mylist[t] = n_nhe3[3]
                result = scipy.optimize.root(eQs, np.array(guess), 1)
                guess = [var[t] for var in vars]
                result_flx = eQs(guess, 0)
                for i in range(len(result_flx)):
                    flx[i][t] = result_flx[i]
            elif 4 * T / len(cmna) < t <= 5 * T / len(cmna):
                # lmi_nagluc  = lmi_na_gluc [4]
                # mylist[t] =lmi_na_gluc [4]
                # lmi_nah = lmi_na_h[4]
                # mylist[t] = lmi_na_h[4]
                nnhe3 = n_nhe3[4]
                mylist[t] = n_nhe3[4]

                # n_p = l_na_k[4]
                # mylist[t] = l_na_k[4]

                result = scipy.optimize.root(eQs, np.array(guess), 1)
                guess = [var[t] for var in vars]
                result_flx = eQs(guess, 0)
                for i in range(len(result_flx)):
                    flx[i][t] = result_flx[i]
    else:
        print('Error')

    t = np.linspace(t0, tf, T)
    t_t = np.transpose(t)
    mylist = np.asarray(mylist)
    fig8 = plt.figure(constrained_layout=True)
    spec2 = gridspec.GridSpec(ncols=1, nrows=4, figure=fig8)
    ax1 = fig8.add_subplot(spec2[0, 0])
    ax1.plot(t[0:-2], mylist[0:-2], 'r-')
    ax1.set_ylabel('L_NHE3', color='blue')
    ax1.set_title('Figure_10', color='black')
    ax4 = fig8.add_subplot(spec2[1,0])
    ax4.plot(t, jnhe3_na, 'b-')
    ax4.set_ylabel('flux_na_nhe3', color='blue')

    ax2 = fig8.add_subplot(spec2[2, 0])
    ax2.plot(t, hco3_clhco3, 'b-')
    ax2.set_ylabel(' flux_hco3_clcho3 ', color='blue')

    ax3 = fig8.add_subplot(spec2[3, 0])
    ax3.plot(t, na_nak, 'b-')
    ax3.set_ylabel('flux_Na_nak ', color='blue')
    ax3.set_xlabel('Time[s] ', color='black')
    fig3 = plt.figure(constrained_layout=True)
    spec2 = gridspec.GridSpec(ncols=1, nrows=6, figure=fig3)
    ax1 = fig3.add_subplot(spec2[0, 0])
    ax1.plot(t[0:-2], mylist[0:-2], 'r-')
    ax1.set_ylabel('L_NHE3', color='blue')
    ax1.set_title('Figure_9', color='black')
    ax2 = fig3.add_subplot(spec2[1, 0])
    ax2.plot(t, ci_na, 'b-')
    ax2.set_ylabel(' ci_na ', color='blue')

    ax3 = fig3.add_subplot(spec2[2, 0])
    ax3.plot(t, ci_cl, 'b-')
    ax3.set_ylabel('ci_cl ', color='blue')

    ax4 = fig3.add_subplot(spec2[3, 0])
    ax4.plot(t, ci_hco3, 'b-')
    ax4.set_ylabel('ci_hco3 ', color='blue')
    ax5 = fig3.add_subplot(spec2[4, 0])
    ax5.plot(t, ci_hco2, 'b-')
    ax5.set_ylabel('ci_hco2 ', color='blue')

    ax6 = fig3.add_subplot(spec2[5, 0])
    ax6.plot(t, clvl_imp, 'b-')
    ax6.set_ylabel('clvl ', color='blue')
    ax6.set_xlabel('Time[s] ', color='blue')
    plt.show()

    my_list_Figure_9_10 = {'t': t, 'L_NHE3': mylist, 'flux_na_nah': jnhe3_na, 'flux_hco3_clcho3': hco3_clhco3,
                           'flux_Na_nak': na_nak, 'ci_na': ci_na, 'ci_cl': ci_cl, 'ci_hco2': ci_hco2,
                           'ci_hco3': ci_hco3, 'clvl_imp': clvl_imp}
    print(my_list_Figure_9_10)
    pickled_list = pickle.dumps(my_list_Figure_9_10)
    of = open('Data_Figure_9_10.py', 'wb')
    of.write(pickled_list)
    of.close()
if Figure_5:
    from PCT_GLOB_new import *
    t = 0
    guess = [var[t] for var in vars]
    NaCl = 1
    scale_per = 0.12
    cmna = [0.14 + NaCl * j * 0.14 * scale_per for j in range(4)]
    cmna.append(0.14)
    cmcl = [0.1132 + NaCl * j * 0.1132 * scale_per for j in range(4)]
    cmcl.append(0.1132)
    mylist_na = [0.1 for i in range(T)]
    mylist_cl = [0.1 for i in range(T)]
    while True:
        t += 1
        if t == T:
            break
        elif 0 < t <= T / len(cmna):
            cm_na = cmna[0]
            cm_cl = cmcl[0]
            cs_na = cmna[0]
            cs_cl = cmcl[0]
            mylist_na[t] = cmna[0]
            mylist_cl[t] = cmcl[0]
            result = scipy.optimize.root(eQs, np.array(guess), 1)
            guess = [var[t] for var in vars]
            result_flx=eQs(guess, 0)
            for i in range(len(result_flx)):
                flx[i][t]=result_flx[i]
        elif 1 * T / len(cmna) < t <= 2 * T / len(cmna):
            mylist_na[t]=cmna[1]
            mylist_cl[t]=cmcl[1]
            cm_na=cmna[1]
            cm_cl=cmcl[1]
            cs_na=cmna[1]
            cs_cl=cmcl[1]
            result=scipy.optimize.root(eQs, np.array(guess), 1)
            guess = [var[t] for var in vars]
            result_flx=eQs(guess, 0)
            for i in range(len(result_flx)):
                flx[i][t]=result_flx[i]
        elif 2 * T / len(cmna) < t <= 3 * T / len(cmna):
            cm_na=cmna[2]
            cm_cl=cmcl[2]
            cs_na=cmna[2]
            cs_cl=cmcl[2]
            mylist_na[t]=cmna[2]
            mylist_cl[t]=cmcl[2]
            result=scipy.optimize.root(eQs, np.array(guess), 1)
            guess = [var[t] for var in vars]
            result_flx=eQs(guess, 0)
            for i in range(len(result_flx)):
                flx[i][t]=result_flx[i]
        elif 3 * T / len(cmna) < t <= 4 * T / len(cmna):
            cm_na=cmna[3]
            cm_cl=cmcl[3]
            cs_na=cmna[3]
            cs_cl=cmcl[3]
            mylist_na[t]=cmna[3]
            mylist_cl[t]=cmcl[3]
            result=scipy.optimize.root(eQs, np.array(guess), 1)
            guess = [var[t] for var in vars]
            result_flx=eQs(guess, 0)
            for i in range(len(result_flx)):
                flx[i][t]=result_flx[i]
        elif 4 * T / len(cmna) < t <= 5 * T / len(cmna):
            cm_na=cmna[4]
            cm_cl=cmcl[4]
            cs_na=cmna[4]
            cs_cl=cmcl[4]
            mylist_na[t]=cmna[4]
            mylist_cl[t]=cmcl[4]
            result=scipy.optimize.root(eQs, np.array(guess), 1)
            guess = [var[t] for var in vars]
            result_flx=eQs(guess, 0)
            for i in range(len(result_flx)):
                flx[i][t] = result_flx[i]
        else:
            print('t is not in a defined range.')
    t=np.linspace(t0, tf, T)
    my_list_Figure_5 = {'t': t, 'mylist_na': mylist_na, 'mylist_cl': mylist_cl, 'vm': vm,
                            'ci_na': ci_na, 'ci_cl': ci_cl, 'ci_k': ci_k }
    pickled_list = pickle.dumps(my_list_Figure_5)
    of = open('Data_Figure_5_v1.py', 'wb')
    of.write(pickled_list)
    of.close()
    fig3 = plt.figure(constrained_layout=True)
    spec2 = gridspec.GridSpec(ncols=1, nrows=5, figure=fig3)
    ax1 = fig3.add_subplot(spec2[0, 0])
    ax1.plot(t[0:-2], mylist_na[0:-2], color='#007282', label="Na ")
    ax1.plot(t[0:-2], mylist_cl[0:-2], color='#b66dff', label="Cl ")
    ax1.set_ylabel('Cm, Cs', color='black')
    ax1.set_xlim(-5, tf)
    ax1.legend(frameon=False)
    ax0 = fig3.add_subplot(spec2[1, 0])
    ax0.set_ylabel('Vm', color='black')
    ax0.plot(t[0:-2], vm[0:-2], 'black')
    ax0.set_xlim(-5, tf)
    ax0.set_ylim(0.0, 0.15)
    ax2 = fig3.add_subplot(spec2[2, 0])
    ax2.set_ylabel('ci_na', color='black')
    ax2.plot(t[0:-2], ci_na[0:-2], 'black')
    ax2.set_xlim(-5, tf)
    ax3 = fig3.add_subplot(spec2[3, 0])
    ax3.plot(t[0:-2], ci_k[0:-2], 'black')
    ax3.set_ylabel('ci_k', color='black')
    ax3.set_xlim(-5, tf)
    ax4 = fig3.add_subplot(spec2[4, 0])
    ax4.plot(t[0:-2], ci_cl[0:-2], 'black')
    ax4.set_ylabel('ci_cl', color='black')
    ax4.set_xlim(-5, tf)
    plt.show()
if Figure_4a:
    from PCT_GLOB_new import *
    mylist_na = [0.1 for i in range(T)]
    mylist_cl = [0.1 for i in range(T)]
    mylist_hco3 = [0.1 for i in range(T)]
    cm_hco3_0 = 0.004
    cm_hco3_f = 0.05
    slope = (cm_hco3_f - cm_hco3_0) / (tf - t0)
    mylist_hco3 = [cm_hco3_0 + 0.1 * slope * i for i in range(T)]
    t = 0
    guess = [var[t] for var in vars]

    while True:
        t += 1
        if t == T:
            break
        else:
            lmi_clhco3 = 0.8000e-08
            lmi_clhco2 = 2.0000e-08
            hmi_h2co2 = 0.05000e00
            his_h2co2 = 0.0600
            hie_h2co2 = 0.0600
            cm_hco3 = mylist_hco3[t]
            result = scipy.optimize.root(eQs, np.array(guess), 1)
            guess = [var[t] for var in vars]
            result_flx = eQs(guess, 0)
            for i in range(len(result_flx)):
                flx[i][t] = result_flx[i]

    my_list = {'fekm_cl': fekm_cl, 'fikm_cl': fikm_cl, 'flux_cl': flux_cl, 'mylist_hco3': mylist_hco3}
    pickled_list = pickle.dumps(my_list)
    of = open('Data_Figure_4a_v1.py', 'wb')
    of.write(pickled_list)
    of.close()
    fig4 = plt.figure(constrained_layout=True)
    spec2 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig4)
    ax3 = fig4.add_subplot(spec2[0, 0])
    ax3.set_ylabel('flux_cl[mmol/sec]', color='black')
    ax3.set_xlabel('cm_hco3[M]', color='black')
    ax3.set_title('Figure_4a', color='black')
    x1 = ax3.plot(mylist_hco3[30:-2], flux_cl[30:-2], 'black')
    x2 = ax3.plot(mylist_hco3[30:-2], fikm_cl[30:-2], 'g')
    x3 = ax3.plot(mylist_hco3[30:-2], fekm_cl[30:-2], 'blue')
    ax3.legend([x1[0], x2[0], x3[0]],['Total Epithelial', 'Cellular', 'Junctional'])
if Figure_4b:
    from PCT_GLOB_new import *
    mylist_hco2 = [0.1 for i in range(T)]
    cm_hco2_0 = 0.004
    cm_hco2_f = 0.05
    slope = (cm_hco2_f - cm_hco2_0) / (tf - t0)
    mylist_hco2 = [cm_hco2_0 + 0.1 * slope * i for i in range(T)]
    t = 0
    guess = [var[t] for var in vars]

    while True:
        t += 1
        if t == T:
            break
        else:
            lmi_clhco3 = 0.5000e-08
            lmi_clhco2 = 2.5000e-08
            hmi_h2co2 = 0.05000e00
            hie_h2co2 = 0.0600
            his_h2co2 = 0.0600
            cm_hco2 = mylist_hco2[t]
            cs_hco2 = mylist_hco2[t]
            cm_hco3 = 0.004
            result = scipy.optimize.root(eQs, np.array(guess), 1)
            guess = [var[t] for var in vars]
            result_flx = eQs(guess, 0)
            for i in range(len(result_flx)):
                flx[i][t] = result_flx[i]
    fig4b = plt.figure(constrained_layout=True)
    spec2 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig4b)
    ax1 = fig4b.add_subplot(spec2[0, 0])
    ax1.set_title('Figure_4b', color='black')
    ax1.set_ylabel('flux_cl[mmol/sec]', color='black')
    ax1.set_xlabel('cm_hco2[M] & cs_hco2[M] ', color='black')
    x1 = ax1.plot(mylist_hco2[30:-2], flux_cl[30:-2], 'black')
    x2 = ax1.plot(mylist_hco2[30:-2], fikm_cl[30:-2], 'g')
    x3 = ax1.plot(mylist_hco2[30:-2], fekm_cl[30:-2], 'blue')
    ax1.legend([x1[0], x2[0], x3[0]],['Total Epithelial', 'Cellular', 'Junctional'])
    plt.show()
    my_list = {'fekm_cl': fekm_cl, 'fikm_cl': fikm_cl, 'flux_cl': flux_cl, 'mylist_hco2': mylist_hco2}
    pickled_list = pickle.dumps(my_list)
    of = open('Data_Figure_4b_v1.py', 'wb')
    of.write(pickled_list)
    of.close()
if Figure_4c:
    from PCT_GLOB_new import *
    t = 0
    guess = [var[t] for var in vars]
    cm_hco3 = 0.004
    cm_hco2_0 = 0.004
    cm_hco2_f = 0.05
    slope = (cm_hco2_f - cm_hco2_0) / (tf - t0)
    mylist_hco2 = [cm_hco2_0 + 0.1 * slope * i for i in range(T)]
    while True:
        t += 1
        if t == T:
            break
        else:
            lmi_clhco3 = 0.0
            lmi_clhco2 = 5.0000e-08
            hmi_h2co2 = 0.1 * 0.05000e00
            hie_h2co2 = 0.1 * 0.0600
            his_h2co2 = 0.1 * 0.0600
            cm_hco2 = mylist_hco2[t]
            cs_hco2 = mylist_hco2[t]
            result = scipy.optimize.root(eQs, np.array(guess), 1)
            guess = [var[t] for var in vars]
            result_flx = eQs(guess, 0)
            for i in range(len(result_flx)):
                flx[i][t] = result_flx[i]
    my_list = {'fekm_cl': fekm_cl, 'fikm_cl': fikm_cl, 'flux_cl': flux_cl, 'mylist_hco2': mylist_hco2}
    pickled_list = pickle.dumps(my_list)
    of = open('Data_Figure_4c_v1.py', 'wb')
    of.write(pickled_list)
    of.close()
    fig4 = plt.figure(constrained_layout=True)
    spec2 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig4)
    ax1 = fig4.add_subplot(spec2[0, 0])
    ax1.set_title('Figure_4c', color='black')
    ax1.set_ylabel('flux_cl[mmol/sec]', color='black')
    ax1.set_xlabel('cm_hco2[M] & cs_hco2[M] ', color='black')
    x1 = ax1.plot(mylist_hco2[30:-2], flux_cl[30:-2], 'black')
    x2 = ax1.plot(mylist_hco2[30:-2], fikm_cl[30:-2], 'g')
    x3 = ax1.plot(mylist_hco2[30:-2], fekm_cl[30:-2], 'blue')
    ax1.legend([x1[0], x2[0], x3[0]],['Total Epithelial', 'Cellular', 'Junctional'])
if Figure_6_7_8:

    from PCT_GLOB_new import *
    t = 0
    guess = [var[t] for var in vars]
    while True:
        t += 1
        if t == T:
            break
        else:
            result = scipy.optimize.root(eQs, np.array(guess), 1)
            guess = [var[t] for var in vars]
            result_flx = eQs(guess, 0)
            for i in range(len(result_flx)):
                flx[i][t] = result_flx[i]
    my_list_default = {'na': [fekm_na[-1], fikm_na[-1], fiks_na[-1], fike_na[-1], feks_na[-1]],
                       'k': [fekm_k[-1], fikm_k[-1], fiks_k[-1], fike_k[-1], feks_k[-1]],
                       'cl': [fekm_cl[-1], fikm_cl[-1], fiks_cl[-1], fike_cl[-1], feks_cl[-1]],
                       'gluc': [fekm_gluc[-1], fikm_gluc[-1], fiks_gluc[-1],fike_gluc[-1], feks_gluc[-1]],
                       'fekm': [fekm_na[-1], fekm_k[-1], fekm_cl[-1], fekm_gluc[-1]],
                       'fikm': [fikm_na[-1], fikm_k[-1], fikm_cl[-1],  fikm_gluc[-1]],
                       'fiks': [fiks_na[-1], fiks_k[-1], fiks_cl[-1], fiks_gluc[-1]],
                       'fike': [fike_na[-1], fike_k[-1], fike_cl[-1],  fike_gluc[-1]],
                       'feks': [feks_na[-1], feks_k[-1], feks_cl[-1], feks_gluc[-1]]}

    slt_con_default = {'Na': ci_na[-1], 'K': ci_k[-1], 'Cl': ci_cl[-1],  'Gluc': ci_gluc[-1]}
    flx_na_mem_default = [feks_na[-1], fike_na[-1], fiks_na[-1], fekm_na[-1], fikm_na[-1], flux_na[-1]]

    epithelial_flx_variation_default = [fekm_na_csf[-1] + fikm_na_csf[-1],
                                        fekm_na_goldman[-1] + fikm_na_goldman[-1],
                                        jnhe3_na[-1] + na_mi_nagluc[-1] + na_mi_nah2po4[-1]]

    epithelial_flx_elct_chmcl_default = [jnhe3_na[-1], na_mi_nagluc[-1], na_mi_nah2po4[-1]]

    t = 0
    guess = [var[t] for var in vars]
    while True:
        t += 1
        if t == T:
            break
        else:
            param_sodium_pumps = 0
            result = scipy.optimize.root(eQs, np.array(guess), 1)
            guess = [var[t] for var in vars]
            result_flx = eQs(guess, 0)
            for i in range(len(result_flx)):
                flx[i][t] = result_flx[i]
    my_list_nak = {'na': [fekm_na[-1], fikm_na[-1], fiks_na[-1], fike_na[-1], feks_na[-1]],
                   'k': [fekm_k[-1], fikm_k[-1], fiks_k[-1], fike_k[-1], feks_k[-1]],
                   'cl': [fekm_cl[-1], fikm_cl[-1], fiks_cl[-1], fike_cl[-1], feks_cl[-1]],
                   'gluc': [fekm_gluc[-1], fikm_gluc[-1], fiks_gluc[-1], fike_gluc[-1], feks_gluc[-1]],
                   'fekm': [fekm_na[-1], fekm_k[-1], fekm_cl[-1], fekm_gluc[-1]],
                   'fikm': [fikm_na[-1], fikm_k[-1], fikm_cl[-1], fikm_gluc[-1]],
                   'fiks': [fiks_na[-1], fiks_k[-1], fiks_cl[-1], fiks_gluc[-1]],
                   'fike': [fike_na[-1], fike_k[-1], fike_cl[-1], fike_gluc[-1]],
                   'feks': [feks_na[-1], feks_k[-1], feks_cl[-1], feks_gluc[-1]]}

    slt_con_nak = {'Na': ci_na[-1], 'K': ci_k[-1], 'Cl': ci_cl[-1], 'Gluc': ci_gluc[-1]}
    flx_na_mem_nak = [feks_na[-1], fike_na[-1], fiks_na[-1], fekm_na[-1], fikm_na[-1], flux_na[-1]]
    epithelial_flx_variation_nak = [fekm_na_csf[-1] + fikm_na_csf[-1], fekm_na_goldman[-1] + fikm_na_goldman[-1],
                                    jnhe3_na[-1] + na_mi_nagluc[-1] + na_mi_nah2po4[-1]]
    epithelial_flx_elct_chmcl_nak = [jnhe3_na[-1],  na_mi_nagluc[-1],  na_mi_nah2po4[-1]]

    t = 0
    guess = [var[t] for var in vars]
    while True:
        t += 1
        if t == T:
            break
        else:
            param_sodium_pumps = 1
            param_kcl = 0
            result = scipy.optimize.root(eQs, np.array(guess), 1)
            guess = [var[t] for var in vars]
            result_flx = eQs(guess, 0)
            for i in range(len(result_flx)):
                flx[i][t] = result_flx[i]
    my_list_kcl = {'na': [fekm_na[-1], fikm_na[-1], fiks_na[-1], fike_na[-1], feks_na[-1]],
                   'k': [fekm_k[-1], fikm_k[-1], fiks_k[-1], fike_k[-1], feks_k[-1]],
                   'cl': [fekm_cl[-1], fikm_cl[-1], fiks_cl[-1], fike_cl[-1], feks_cl[-1]],
                   'gluc': [fekm_gluc[-1], fikm_gluc[-1], fiks_gluc[-1], fike_gluc[-1], feks_gluc[-1]],
                   'fekm': [fekm_na[-1], fekm_k[-1], fekm_cl[-1], fekm_gluc[-1]],
                   'fikm': [fikm_na[-1], fikm_k[-1], fikm_cl[-1],  fikm_gluc[-1]],
                   'fiks': [fiks_na[-1], fiks_k[-1], fiks_cl[-1], fiks_gluc[-1]],
                   'fike': [fike_na[-1], fike_k[-1], fike_cl[-1],  fike_gluc[-1]],
                   'feks': [feks_na[-1], feks_k[-1], feks_cl[-1], feks_gluc[-1]]}

    slt_con_kcl = {'Na': ci_na[-1], 'K': ci_k[-1], 'Cl': ci_cl[-1], 'Gluc': ci_gluc[-1]}
    t = 0

    t = 0
    guess = [ var [ t ] for var in vars ]
    while True:
        t += 1
        if t == T:
            break
        elif t < T/4:
            param_sodium_pumps = 1
            param_kcl = 1
            param_nahco3 = 1
            result = scipy.optimize.root(eQs, np.array(guess), 1)
            guess = [var[t] for var in vars]
            result_flx = eQs(guess, 0)
            for i in range(len(result_flx)):
                flx[i][t] = result_flx[i]
        else:
            param_sodium_pumps = 1
            param_kcl = 1
            param_nahco3 = 0
            result = scipy.optimize.root(eQs, np.array(guess), 1)
            guess = [var[t] for var in vars]
            result_flx = eQs(guess, 0)
            for i in range(len(result_flx)):
                flx[i][t] = result_flx[i]
    my_list_nahco3 = {'na':[fekm_na[-1], fikm_na[-1], fiks_na[-1], fike_na[-1], feks_na[-1]],
                      'k':[fekm_k[-1], fikm_k[-1], fiks_k[-1], fike_k[-1], feks_k[-1]],
                      'cl':[fekm_cl[-1], fikm_cl[-1], fiks_cl[-1], fike_cl[-1], feks_cl[-1]],
                      'gluc':[fekm_gluc[-1], fikm_gluc[-1], fiks_gluc[-1], fike_gluc[-1], feks_gluc[-1]],
                      'fekm':[fekm_na[-1], fekm_k[-1], fekm_cl[-1], fekm_gluc[-1]],
                      'fikm':[fikm_na[-1], fikm_k[-1], fikm_cl[-1], fikm_gluc[-1]],
                      'fiks':[fiks_na[-1], fiks_k[-1], fiks_cl[-1], fiks_gluc[-1]],
                      'fike':[fike_na[-1], fike_k[-1], fike_cl[-1], fike_gluc[-1]],
                      'feks':[feks_na[-1], feks_k[-1], feks_cl[-1], feks_gluc[-1]]}

    slt_con_nahco3 = {'Na': ci_na[-1],  'K': ci_k[-1], 'Cl':ci_cl[-1], 'Gluc': ci_gluc[-1]}
    # print('ci_k[-1]',ci_k[-1])
    from collections import defaultdict
    my_list_Figure6A = defaultdict(list)
    for d in (my_list_default, my_list_nak, my_list_kcl, my_list_nahco3):
        for key, value in d.items():
            my_list_Figure6A[key].append(value)
    pickled_list = pickle.dumps(my_list_Figure6A)
    of = open('Data_Figure_6A_v3.py', 'wb')
    of.write(pickled_list)
    of.close()
    my_list_Figure6B = defaultdict(list)
    for d in (slt_con_default, slt_con_nak, slt_con_kcl, slt_con_nahco3):
        for key, value in d.items():
            my_list_Figure6B[key].append(value)
    pickled_list = pickle.dumps(my_list_Figure6B)
    of = open('Data_Figure_6B_v3.py', 'wb')
    of.write(pickled_list)
    of.close()
    t = 0
    guess = [var[t] for var in vars]
    while True:
        t += 1
        if t == T:
            break
        else:
            param_sodium_pumps = 1
            param_kcl = 1
            param_nahco3 = 1
            param_nhe3 = 0
            result = scipy.optimize.root(eQs, np.array(guess), 1)
            guess = [var[t] for var in vars]
            result_flx = eQs(guess, 0)
            for i in range(len(result_flx)):
                flx[i][t] = result_flx[i]
    my_list_nhe3 = {'na':[fekm_na[-1], fikm_na[-1], fiks_na[-1], fike_na[-1], feks_na[-1]],
                    'k':[fekm_k[-1], fikm_k[-1], fiks_k[-1], fike_k[-1], feks_k[-1]],
                    'cl':[fekm_cl[-1], fikm_cl[-1], fiks_cl[-1], fike_cl[-1], feks_cl[-1]],
                    'gluc':[fekm_gluc[-1], fikm_gluc[-1], fiks_gluc[-1], fike_gluc[-1], feks_gluc[-1]],
                    'fekm':[fekm_na[-1], fekm_k[-1], fekm_cl[-1], fekm_gluc[-1]],
                    'fikm':[fikm_na[-1], fikm_k[-1], fikm_cl[-1], fikm_gluc[-1]],
                    'fiks':[fiks_na[-1], fiks_k[-1], fiks_cl[-1], fiks_gluc[-1]],
                    'fike':[fike_na[-1], fike_k[-1], fike_cl[-1], fike_gluc[-1]],
                    'feks':[feks_na[-1], feks_k[-1], feks_cl[-1], feks_gluc[-1]]}

    slt_con_nhe3 = {'Na': ci_na[-1], 'K': ci_k[-1], 'Cl': ci_cl[-1], 'Gluc': ci_gluc[-1]}
    flx_na_mem_nhe3 = [feks_na[-1], fike_na[-1], fiks_na[-1], fekm_na[-1], fikm_na[-1], fekm_na[-1] + fikm_na[-1]]
    epithelial_flx_variation_nhe3 = [fekm_na_csf[-1] + fikm_na_csf[-1], fekm_na_goldman[-1] + fikm_na_goldman[-1],
                                     jnhe3_na[-1] + na_mi_nagluc[-1]+ na_mi_nah2po4[-1]]
    epithelial_flx_elct_chmcl_nhe3 = [jnhe3_na[-1],  na_mi_nagluc[-1],  na_mi_nah2po4[-1]]
    t = 0
    guess = [var[t] for var in vars]
    while True:
        t += 1
        if t == T:
            break
        else:
            param_sodium_pumps = 1
            param_kcl = 1
            param_nahco3 = 1
            param_nhe3 = 1
            param_sglt = 0
            result = scipy.optimize.root(eQs, np.array(guess), 1)
            guess = [var[t] for var in vars]
            result_flx = eQs(guess, 0)
            for i in range(len(result_flx)):
                flx[i][t] = result_flx[i]
    my_list_nagluc = {'na':[fekm_na[-1], fikm_na[-1], fiks_na[-1], fike_na[-1], feks_na[-1]],
                      'k':[fekm_k[-1], fikm_k[-1], fiks_k[-1], fike_k[-1], feks_k[-1]],
                      'cl':[fekm_cl[-1], fikm_cl[-1], fiks_cl[-1], fike_cl[-1], feks_cl[-1]],
                      'gluc':[fekm_gluc[-1], fikm_gluc[-1], fiks_gluc[-1],fike_gluc[-1], feks_gluc[-1]],
                      'fekm':[fekm_na[-1], fekm_k[-1], fekm_cl[-1], fekm_gluc[-1]],
                      'fikm':[fikm_na[-1], fikm_k[-1], fikm_cl[-1],  fikm_gluc[-1]],
                      'fiks':[fiks_na[-1], fiks_k[-1], fiks_cl[-1], fiks_gluc[-1]],
                      'fike':[fike_na[-1], fike_k[-1], fike_cl[-1],  fike_gluc[-1]],
                      'feks':[feks_na[-1], feks_k[-1], feks_cl[-1], feks_gluc[-1]]}

    slt_con_nagluc = {'Na': ci_na[-1], 'K': ci_k[-1], 'Cl': ci_cl[-1],  'Gluc': ci_gluc[-1]}
    t = 0
    guess = [var[t] for var in vars]
    while True:
        t += 1
        if t == T:
            break
        else:
            param_sodium_pumps = 1
            param_kcl = 1
            param_nahco3 = 1
            param_nhe3 = 1
            param_sglt = 1
            param_h2po4 = 0
            result = scipy.optimize.root(eQs, np.array(guess), 1)
            guess = [var[t] for var in vars]
            result_flx = eQs(guess, 0)
            for i in range(len(result_flx)):
                flx[i][t] = result_flx[i]

    my_list_nah2po4 = {'na':[fekm_na[-1], fikm_na[-1], fiks_na[-1], fike_na[-1], feks_na[-1]],
                       'k':[fekm_k[-1], fikm_k[-1], fiks_k[-1], fike_k[-1], feks_k[-1]],
                       'cl':[fekm_cl[-1], fikm_cl[-1], fiks_cl[-1], fike_cl[-1], feks_cl[-1]],
                       'gluc':[fekm_gluc[-1], fikm_gluc[-1], fiks_gluc[-1], fike_gluc[-1],  feks_gluc[-1]],
                       'fekm':[fekm_na[-1], fekm_k[-1], fekm_cl[-1], fekm_gluc[-1]],
                       'fikm':[fikm_na[-1], fikm_k[-1], fikm_cl[-1],  fikm_gluc[-1]],
                       'fiks':[fiks_na[-1], fiks_k[-1], fiks_cl[-1], fiks_gluc[-1]],
                       'fike':[fike_na[-1], fike_k[-1], fike_cl[-1],  fike_gluc[-1]],
                       'feks':[feks_na[-1], feks_k[-1], feks_cl[-1], feks_gluc[-1]]}
    slt_con_nah2po4 = {'Na': ci_na[-1],  'K': ci_k[-1], 'Cl': ci_cl[-1],'Gluc': ci_gluc[-1]}

    from collections import defaultdict
    my_list_Figure7 = {'flx_na_mem_default': flx_na_mem_default, 'epithelial_flx_variation_default':
                       epithelial_flx_variation_default, 'epithelial_flx_elct_chmcl_default': epithelial_flx_elct_chmcl_default,
                       'flx_na_mem_nak': flx_na_mem_nak, 'epithelial_flx_variation_nak':epithelial_flx_variation_nak, 'epithelial_flx_elct_chmcl_nak': epithelial_flx_elct_chmcl_nak,
                       'flx_na_mem_nhe3': flx_na_mem_nhe3, 'epithelial_flx_variation_nhe3':epithelial_flx_variation_nhe3, 'epithelial_flx_elct_chmcl_nhe3': epithelial_flx_elct_chmcl_nhe3}
    pickled_list = pickle.dumps(my_list_Figure7)
    of = open('Data_Figure_7_v3.py', 'wb')
    of.write(pickled_list)
    of.close()
    my_list_Figure8A = defaultdict(list)
    for d in (my_list_default, my_list_nhe3, my_list_nagluc, my_list_nah2po4):
        for key, value in d.items():
            my_list_Figure8A[key].append(value)
    pickled_list = pickle.dumps(my_list_Figure8A)
    of = open('Data_Figure_8A_v3.py', 'wb')
    of.write(pickled_list)
    of.close()
    my_list_Figure8B = defaultdict(list)
    for d in (slt_con_default, slt_con_nhe3, slt_con_nagluc, slt_con_nah2po4):
        for key, value in d.items():
            my_list_Figure8B[key].append(value)
    pickled_list = pickle.dumps(my_list_Figure8B)
    of = open('Data_Figure_8B_v3.py', 'wb')
    of.write(pickled_list)
    of.close()
    of = open('Data_Figure_6A_v3.py', 'rb')
    read_file = of.read()
    my_loaded_list = pickle.loads(read_file)
    of.close()

    flux_na = np.array(my_loaded_list['na'])
    flux_k = np.array(my_loaded_list['k'])
    flux_cl = np.array(my_loaded_list['cl'])
    flux_gluc = np.array(my_loaded_list['gluc'])
    df1 = pd.DataFrame(flux_na,
                       index= [" A Original Setup", "B NaK = 0", "C KCl = 0", "D NaHCO3 = 0"],
                       columns= ["EM", "IM", "IS", "IE", "ES"])
    df2 = pd.DataFrame(flux_k,
                       index= [" A Original Setup", "B NaK = 0", "C KCl = 0", "D NaHCO3 = 0"],
                       columns= ["EM", "IM", "IS", "IE", "ES"])
    df3 = pd.DataFrame(flux_cl,
                       index= [" A Original Setup", "B NaK = 0", "C KCl = 0", "D NaHCO3 = 0"],
                       columns= ["EM", "IM", "IS", "IE", "ES"])
    df4 = pd.DataFrame(flux_gluc,
                       index= [" A Original Setup", "B NaK = 0", "C KCl = 0", "D NaHCO3 = 0"],
                       columns= ["EM", "IM", "IS", "IE", "ES"])

    plt.figure(12, figsize=(12, 12))
    plot_clustered_stacked([df1, df2, df3, df4],["Na", "K", "Cl", "Gluc"], title= "Figure_6A, Membrane Fluxes[mmol/s.cm2]")
#    plt.savefig('Figure_6A' + '.png')
#     of = open('Data_Figure_8A_v3.py', 'rb')
#     read_file = of.read()
#     my_loaded_list = pickle.loads(read_file)
#     of.close()
#
#     flux_na = np.array( my_loaded_list['na'])
#     flux_k = np.array(my_loaded_list['k'])
#     flux_cl = np.array(my_loaded_list['cl'])
#     flux_gluc = np.array(my_loaded_list['gluc'])
#     df1 = pd.DataFrame(flux_na,
#                        index= ["A Original Setup", " B NHE3 = 0", "C SGLT = 0", "D NaH2PO4 = 0"],
#                        columns= ["EM", "IM", "IS", "IE", "ES"])
#     df2 = pd.DataFrame(flux_k,
#                        index= ["A Original Setup", " B NHE3 = 0", "C SGLT = 0", "D NaH2PO4 = 0"],
#                        columns= ["EM", "IM", "IS", "IE", "ES"])
#     df3 = pd.DataFrame(flux_cl,
#                        index= ["A Original Setup", " B NHE3 = 0", "C SGLT = 0", "D NaH2PO4 = 0"],
#                        columns= ["EM", "IM", "IS", "IE", "ES"])
#     df4 = pd.DataFrame(flux_gluc,
#                        index= ["A Original Setup", " B NHE3 = 0", "C SGLT = 0", "D NaH2PO4 = 0"],
#                        columns= ["EM", "IM", "IS", "IE", "ES"])
#
#     plt.figure(18, figsize=(12, 12))
#     plot_clustered_stacked([df1, df2, df3,df4],["Na", "K", "Cl", "Gluc"], title="Figure_8A, Membrane Fluxes[mmol/s.cm2]")
    # plt.savefig('Figure_8A' + '.png')
    of = open('Data_Figure_6B_v3.py', 'rb')
    read_file = of.read()
    my_loaded_list1 = pickle.loads(read_file)
    of.close()

    Na = np.array(my_loaded_list1['Na'])
    K = np.array(my_loaded_list1['K'])
    Cl = np.array(my_loaded_list1['Cl'])
    Gluc = np.array(my_loaded_list1['Gluc'])
    x = ["Original Setup", "NaK = 0", "KCl = 0", "NaHCO3 = 0"]

    bar1 = np.arange(len(x))
    w = 0.1
    bar2 = [i + 2.0 * w for i in bar1]
    bar3 = [i + 2.0 * w for i in bar2]
    bar4 = [i + 2.0 * w for i in bar3]
    bar5 = [i + w for i in bar2]
    plt.figure(19)
    plt.bar(bar1, Na, 0.2, color="#2B9F78", label="Na")
    plt.bar(bar2, K, 0.2, color="#007282", label="K")
    plt.bar(bar3, Cl, 0.2, color="black", label="Cl")
    plt.bar(bar4, Gluc, 0.2, color="#D55E00", label="Gluc")
    plt.title('Figure_6B', fontsize=12)
    plt.ylabel('Cellular Concentration[mol/l]', fontsize=12)
    plt.xlabel("Solutes", fontsize=12)
    plt.xticks(bar5, x, fontsize=11)
    plt.legend(frameon=False)
    # plt.savefig('Figure_8B' + '.png')
    of = open('Data_Figure_8B_v3.py', 'rb')
    read_file = of.read()
    my_loaded_list1 = pickle.loads(read_file)
    of.close()

    Na = np.array(my_loaded_list1['Na'])
    K = np.array(my_loaded_list1['K'])
    Cl = np.array(my_loaded_list1['Cl'])
    Gluc = np.array(my_loaded_list1['Gluc'])
    x = ["Original Setup", "NHE3 = 0", "SGLT = 0", "NaH2PO4 = 0"]

    bar1 = np.arange(len(x))
    bar2 = [i + 2.0 * w for i in bar1]
    bar3 = [i + 2.0 * w for i in bar2]
    bar4 = [i + 2.0 * w for i in bar3]
    bar5 = [i + w for i in bar2]
    plt.figure(20)
    plt.bar(bar1, Na, 0.2, color="#2B9F78", label="Na")
    plt.bar(bar2, K, 0.2, color="#007282", label="K")
    plt.bar(bar3, Cl, 0.2, color="black", label="Cl")
    plt.bar(bar4, Gluc, 0.2, color="#D55E00", label="Gluc")
    plt.title('Figure_8B', fontsize=12)
    plt.ylabel('Cellular Concentration[mol/l]', fontsize=12)
    plt.xlabel("Solutes", fontsize=12)
    plt.xticks(bar5, x, fontsize=11)
    plt.legend(frameon=False)
    # plt.savefig('Figure_6B' + '.png')
    of = open('Data_Figure_7_v3.py', 'rb')
    read_file = of.read()
    my_loaded_list = pickle.loads(read_file)
    of.close()

    fig = plt.figure(num=21, figsize=(12, 12))
    x1= ['ES', 'IE', 'IS', 'ME', 'MI','EPTL']
    ax1 = fig.add_subplot(3, 3, 1)
    ax1.bar(x1, np.array(my_loaded_list['flx_na_mem_default']), color= ['black', '#b66dff', '#24ff24', '#007282', '#2B9F78', '#920000'])
    # ax1.set_ylim(-0.2,7.5)
    plt.setp(plt.gca(), xticklabels=[])
    x2 = ['CSF', 'GLD', 'ECHMCL']
    ax2 = fig.add_subplot(3, 3, 2)
    Epithelial = my_loaded_list['epithelial_flx_variation_default']
    ax2.bar(x2, Epithelial, color= ['black', '#b66dff','#2B9F78'])
    # ax2.set_ylim(-0.2,8.4)
    plt.setp(plt.gca(), xticklabels=[])
    x3 = ['NHE3', 'SGLT','NaH2PO4']
    ax3 = fig.add_subplot(3, 3, 3)
    ax3.bar(x3, my_loaded_list['epithelial_flx_elct_chmcl_default'], color= ['#db6d00', '#b66dff', '#2B9F78'])
    # ax3.set_ylim(-0.2,3.8)
    plt.setp(plt.gca(), xticklabels=[])
    ax4 = fig.add_subplot(3, 3, 4)
    labels = ['Interspace Basement', 'Cell Lateral', 'Basal Cell', 'Tight Junction', 'Cell Apical', 'Epithelial']
    colors = ['black', '#b66dff', '#24ff24', '#007282', '#2B9F78', '#920000']
    for tmpX, tmpY, color, label in zip(x1, my_loaded_list['flx_na_mem_nak'], colors, labels):
        ax4.bar(tmpX, tmpY, color=color, label=label)
    # ax4.set_ylim(-0.2, 7.5)
    plt.setp(plt.gca(), xticklabels=[])
    ax4.legend(shadow=False, fancybox=False, frameon=False)
    ax5 = fig.add_subplot(3, 3, 5)
    labels = ['Covective ', 'Passive', 'Electochemical']
    colors = ['black', '#b66dff', '#2B9F78']
    for tmpX, tmpY, color, label in zip(x2, my_loaded_list['epithelial_flx_variation_nak'], colors, labels):
        ax5.bar(tmpX, tmpY, color=color, label=label)
    plt.setp(plt.gca(), xticklabels=[])
    # ax5.set_ylim(-0.2, 8.4)
    ax5.legend(shadow=False, fancybox=False, frameon=False)
    ax6 = fig.add_subplot(3, 3, 6)
    colors = ['#db6d00', '#b66dff','#2B9F78']
    labels = ['NHE3','SGLT','NaH2PO4']
    for tmpX, tmpY, color, label in zip(x3, my_loaded_list['epithelial_flx_elct_chmcl_nak'], colors, labels):
        ax6.bar(tmpX, tmpY, color=color, label=label)
    plt.setp(plt.gca(), xticklabels=[])
    # ax6.set_ylim(-0.1,3.8)
    ax6.legend(shadow=False, fancybox=False, frameon=False)
    ax7=fig.add_subplot(3, 3, 7)
    ax7.bar(x1, my_loaded_list['flx_na_mem_nhe3'], color= ['black', '#b66dff','#24ff24', '#007282','#2B9F78','#920000'])
    ax7.set_xlabel('Membrane Flux_Na[mmol/s.cm2]', color='black')
    # ax7.set_ylim(-0.2,7.5)
    plt.setp(plt.gca(), xticklabels=[])
    ax8=fig.add_subplot(3, 3, 8)
    ax8.bar(x2,my_loaded_list['epithelial_flx_variation_nhe3'], color= ['black', '#b66dff','#2B9F78'])
    ax8.set_xlabel('Epithelial Flux_Na[mmol/s.cm2]', color='black')
    plt.setp(plt.gca(), xticklabels=[])
    # ax8.set_ylim(-0.2,8.4)
    ax9 = fig.add_subplot(3, 3, 9)
    ax9.bar(x3, my_loaded_list['epithelial_flx_elct_chmcl_nhe3'], color= ['#db6d00', '#b66dff', '#2B9F78'])
    plt.setp(plt.gca(), xticklabels=[])
    ax9.set_xlabel('Electrochemical Flux_Na[mmol/s.cm2]', color='black')
    # ax9.set_ylim(-0.1,3.8)

    for ax, color in zip([ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9], ['black', '#920000', '#2B9F78', 'black', '#920000', '#2B9F78', 'black', '#920000','#2B9F78']):
        for ticks in ax.xaxis.get_ticklines() + ax.yaxis.get_ticklines():
            ticks.set_color(color)
        for pos in['top', 'bottom', 'right', 'left']:
            ax.spines[pos].set_edgecolor(color)
    plt.show()