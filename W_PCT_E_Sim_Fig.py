import numpy as np
import scipy
import scipy.optimize
from scipy import optimize
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
# from sympy import *
import math
from PCT_GLOB_new import *
import pandas as pd
import altair as alt
import matplotlib.cm as cm
import numpy as np
import pickle
# ELECTROLYTE  TRANSPORT   ACROSS  A SIMPLE EPITHELIUM STEADY - STATE AND  TRANSIENT ANALYSIS

def Extract(lst):
    return [ item [ -1 ] for item in lst ]

def reset(lst):
    return [reset(e) if isinstance(e, list) else 0 for e in lst]



def clear(seq):
    for i, item in enumerate(seq):
        try:
            clear(item)
        except TypeError:
            try:
                seq[i] = 0
            except IndexError and KeyError:
                pass

def csf(ca, cb, flux, s, param_csf):
    # Convective Solute Fluxes, see Eqs: (37, 38)
    if param_csf == 0:
        return 0

    def lmmsc(ca, cb):
        # Definition of the logarithmic mean membrane solute concentration
        if ca > 0 and cb > 0 and ca - cb != 0 and cb != 0 and np.log10(ca / cb) != 0:
            return (ca - cb) / (np.log10(ca / cb))
        else:
            return cb
    return flux * (1.00 - s) * lmmsc(ca, cb)


def AAPC(ca, cb, z, a, hab, va, vb, param_AAPC):
    # AREA    ADJUSTED    PARTIAL    CONDUCTANCE,  see Eqs: (39), (41)
    zab = (z * (va - vb)*f * 1.e-6) / rte
    if param_AAPC == 0:
        return 0
    elif zab == 0 or va == vb and ca > 0 and cb > 0:
        return hab * a * (ca - cb)
    elif zab > 0:
        return hab * a * zab * (ca - cb * math.exp(-zab)) / (1 - math.exp(-zab))
    else:
        return hab * a * zab * (ca * math.exp(zab) - cb) / (math.exp(zab) - 1)


def goldman(hab, a, z, va, vb, ca, cb, param_goldman):
    # see  Eqs: (39), (41)
    zab = (z * f * (va - vb) * 1.e-6) / rte
    if param_goldman == 0:
        return [ 0 ]
    elif zab == 0 or va == vb and ca > 0 and cb > 0:
        return hab * a * (ca - cb)
    elif zab > 0:
        return hab * a * zab * (ca - cb * np.exp(-zab)) / (1 - np.exp(-zab))
    else:
        return hab * a * zab * (ca * np.exp(zab) - cb) / (np.exp(zab) - 1)


def k_cl(ca_k, cb_k, ca_cl, cb_cl, z_k, z_cl, va, vb, a, l_kcl, param_kcl):
    # k-cl simple cotransporter with 1:1 stoichiometry, located on Peritubular Membrane which
    # includes both Cell-Lateral Membrane (ie) /Cell-Basal (is) Membrane.
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
    return [k_kcl, cl_kcl]


def sglt_mi(cm_na, ci_na, cm_gluc, ci_gluc, z_na, z_gluc, vm, vi, ami, lmi_nagluc, param_sglt_mi):
    # Na-glucose simple cotransporter with 1:1 stoichiometry, located on  Apical  Membrane
    # return is the transported flux for each solute due to the existance of Na-glucose simple cotransporter.
    # f_eps(): Modular function to calculate the electrochemical potential of species:(RT)*ln(c) + z*f*V
    # see Eq: (43)
    # checking the unit consistency for Electrodiffusive_flux Equation:
    # [mV] = 1.e-3[V]----> [mV] = 1.e-3 [Joule/Coul]
    # [xm]---->[Joule/mmol] = [Joule/mmol][1] + [1][Coul/mol][mV]
    # [xm]---->[Joule/mmol] = [Joule/mmol][1] + [1][Coul/mol]*1.e-3 [ Joule / Coul ]
    # [xm]---->[Joule/mmol] = [Joule/mmol][1] + [1][Coul/mmol]* 1.e-6 [ Joule / Coul ]
    # [Electrodiffusive_flux] ----> [mmol/s] = [mmol2/Joule.s]*[Joule/mmol]
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
    # Na-h2po4 simple cotransporter with 1:1 stoichiometry, located on  Apical  Membrane
    # see Eq: (44)
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
    # Weintein introduced two simple exchangers  na/h and na/nh4 in their work mentioned above and
    # it seems these two simple exchangers are equivalent to NHE transporter that they introduced in 1995.
    # They located on Apical membrane.
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
    # Weintein introduced two simple exchangers  na/h and na/nh4 in their work mentioned above and
    # it seems these two simple exchangers are equivalent to NHE transporter that they introduced in 1995.
    # They located on Apical membrane.
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
    # cl/hco2 simple exchanger with 1:-1 stoichiometry, located on Apical  Membrane.
    # Eq: (47)
    if param_clhco2_mi == 0:
        return [ 0, 0 ]
    else:
        xm_cl = f_eps(cm_cl, z_cl, vm)
        xm_hco2 = f_eps(cm_hco2, z_hco2, vm)
        xi_cl = f_eps(ci_cl, z_cl, vi)
        xi_hco2 = f_eps(ci_hco2, z_hco2, vi)
        cl_mi_clhco2 = lmi_clhco2 * ami * (xm_cl - xi_cl - xm_hco2 + xi_hco2)
        hco2_mi_clhco2 = - lmi_clhco2 * ami * (xm_cl - xi_cl - xm_hco2 + xi_hco2)
        return [ cl_mi_clhco2, hco2_mi_clhco2 ]


def clhco3_mi(cm_cl, ci_cl, cm_hco3, ci_hco3, z_cl, z_hco3, vm, vi, ami, lmi_clhco3, param_clhco3_mi):
    # cl/hco3 simple exchanger with 1:-1 stoichiometry, located on Apical  Membrane.
    # See Eq: (48)
    if param_clhco3_mi == 0:
        return [ 0, 0 ]
    else:
        xm_cl = f_eps(cm_cl, z_cl, vm)
        xm_hco3 = f_eps(cm_hco3, z_hco3, vm)
        xi_cl = f_eps(ci_cl, z_cl, vi)
        xi_hco3 = f_eps(ci_hco3, z_hco3, vi)
        cl_mi_clhco3 = lmi_clhco3 * ami * (xm_cl - xi_cl - xm_hco3 + xi_hco3)
        hco3_mi_clhco3 = - lmi_clhco3 * ami * (xm_cl - xi_cl - xm_hco3 + xi_hco3)
    return [ cl_mi_clhco3, hco3_mi_clhco3 ]


def na_hco3(ci_na, cs_na, ci_hco3, cs_hco3, z_na, z_hco3, vi, vs, ais, lis_nahco3, param_na_hco3):
    # See Eq: (49)
    if param_na_hco3 == 0:
        return[0, 0]
    else:
        xi_na = f_eps(ci_na, z_na, vi)
        xi_hco3 = f_eps(ci_hco3, z_hco3, vi)
        xs_na = f_eps(cs_na, z_na, vs)
        xs_hco3 = f_eps(cs_hco3, z_hco3, vs)
        na_is_nahco3 = lis_nahco3 * ais * (xi_na - xs_na + 3 * (xi_hco3 - xs_hco3))
        hco3_is_nahco3 = 3 * lis_nahco3 * ais * (xi_na - xs_na + 3 * (xi_hco3 - xs_hco3))
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
    # Proton pumps area special kind  of transporter  that push
    # hydrogen ions  from areas of low concentration
    # with high concentration. Ions moving down a gradient release energy,
    # but when they move up a gradient, it takes energy.Diffusion can then use this
    # gradient to capture energy again, as the ions move downhill.
    # See Eq: (53)
    xm_h = f_eps(cm_h, z_h, vm)
    xi_h = f_eps(ci_h, z_h, vi)
    gamma = -(xihp * (xm_h - xi_h - xhp))
    if parama_mi_h == 0:
        return 0
    elif gamma < 0:
        return -lhp * ami * (1 - 1 / (1 + np.exp(gamma)))
    return -lhp * ami * (1 / (1 + np.exp(-gamma)))


def nak_atp( ci_na, ce_na,ci_k,  ce_k, ce_nh4, param_nak_atp):
    # Powered by ATP, the pump moves
    # sodium and potassium  ions in opposite
    # directions, each against its concentration
    # gradient.In  a single cycle of the pump, three sodium ions are
    # extruded from and two  potassium ions are imported into  the cell.
    # See Eqs: (54, 55, 56, 57, 58)
    if param_nak_atp == 0:
        return[0, 0, 0]
    else:
        # sodium affinity
        knpn = 0.0002 * (1.0 + ci_k / .00833)
        # knh4 = 0.0001 * (1.0 + ce_na / .0185)
        knpk = 0.0001 * (1.0 + ce_na / .0185)
        # atpase transporter flux in is membrane
        a = (ci_na / (knpn + ci_na)) ** 3
        b = ((ce_k + ce_nh4) / (knpk + ce_k + ce_nh4)) ** 2
        atp_na = n_p * a * b
        # alloW for competition betWeen k+ and nh4+
        atp_k = -atp_na * (2/3) * ce_k / (ce_k + ce_nh4 / knh4)
        atp_nh4 = -atp_na * (2/3) * ce_nh4 / (knh4 * ce_k + ce_nh4)
        return[atp_na, atp_k, atp_nh4]


def lch(ca, cb):
    if ca > 0 and cb > 0 and (ca / cb) != 0 and cb != 0:
        return np.log10(ca / cb)
    else:
        return np.log10(abs(ca / cb))

def f_eps(c, z, v):
    # electrochemical potential of species
    # see Eq: (52)
    if c > 0 and c != 0:
        return rte * np.log(c) + z * f * v * 1.e-6
    else:
        print('uy')
        return rte * np.log(abs(c)) + z * f * v * 1.e-6


def ebuf(lch, pk, ca, cb, param_ebuf):
    # pH equilibria of four buffer pairs, see Eqs: (15), (16)
    if param_ebuf == 0:
        return 0
    if ca > 0 and cb > 0 and (ca / cb) != 0 and cb != 0:
        return lch - pk - np.log10(ca / cb)

    else:
        return lch - pk - np.log10(abs(ca / cb))


def buff_activation_co2_formate_phosphate_ammonia(q_h, q_nh4, q_hco3, q_h2co2, q_h2po4, q_co2, q_h2co3, q_cb, c_co2, c_h2co3, volume, scale,flow_dependent, Co2_Progressive_Activation_Param):
    if Co2_Progressive_Activation_Param == 1:
        # print('Buffer activation at time =' + str(int(4000 * Nfac)))
        # print('Buffer activation for co2, formate, phosphate, and ammonia', 't=', str(t))
        # co2, formate, phosphate, and ammonia content:
        # print('q_hco3,  q_h2co3, q_co2', str(q_hco3), str(q_h2co3), str(q_co2))

        if flow_dependent == 1:
            q_hco3 = -1 * (q_h + q_nh4 - q_hco3 + q_h2co2 + q_h2po4 - q_cb)
            q_h2co3 = -1 * (- q_h2co3 - q_h2co3 + scale * volume * (khy * c_co2 - kdhy * c_h2co3))
            q_co2  = q_hco3 + q_h2co3 + q_co2
        else:
            # print('No Buffer Effect for co2, formate, phosphate, and ammonia')
            q_hco3 = -1 * ( q_h + q_nh4 - q_hco3 + q_h2co2 + q_h2po4 )
            # LIS: hydration and dhydration of co2
            # see label {co2_hyd_dhyd}
            q_h2co3 = -1 * (- q_h2co3 - q_h2co3 + scale * volume * (khy * c_co2 - kdhy * c_h2co3))
            # LIS: see label {conser_charge_co2}
            q_co2 = q_hco3 + q_h2co3 + q_co2
            #print('q_hco3,  q_h2co3, q_co2', str(q_hco3),  str(q_h2co3), str(q_co2))
    return q_hco3,  q_h2co3, q_co2


def buff_activation(q_h_vary, q_h2_vary, lc_h, pk, c_h_vary, c_h2_vary, qi_nh3_cell, Buffer_Activation_Param):
    if Buffer_Activation_Param == 1:
        #  pH equilibria of four buffer pairs
        # LIS: see label {phosphate}
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

# common definition


def zab(z, va, vb):
    return z * f * (va - vb) * 1.e-6 / rte


def phi_scale(phi, scale):
    return phi * scale


def nhe3(ci_h, ci_na, ci_nh4, cm_h, cm_na, cm_nh4, param_nah):
    # Sodium–proton  exchanger(Na + / H + exchanger) is a membrane protein
    # that transports Na + into the cell, and H + out  of  the
    # cell, which is located in the  cell  basal and cell  basement.They
    # are  found in the epithelial  cells  of  the  proximal
    # convoluted tubule.Na + / H + exchangers  are  thought
    # to  be  implicated in other   disorders  such as hypertension.Defects in Na + / H + antiporters
    # may  result in kidney failure.
    if param_nah == 0:
        return[0, 0, 0]
    else:
        # the nhe3 model parameters
        cxt = 1.00

        # pnah_na = 0.8000e+04
        # pnah_h = 0.2400e+04
        # pnah_nh4 = 0.8000e+04
        # pnah_m = 0.2000e+01
        # pnah_mm = 0.0000

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

        # translate concentrations to the nhe3 model
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
    # update variables
    for i in range(len(guess)):
        vars [ i ] [ t ] = guess [ i ]

    # Weinstein (2007)
    # Flow-dependent transport in a mathematical model of rat proximal tubule
    ae [ t ] = ae0 * (1 + mua * (pe [ t ] - pm))
    ae [ t ] = ae0 if ae [ t ] < ae0 else ae [ t ]

    chvl [ t ] = chvl0 * (1.0 + muv * (pe [ t ] - pm))
    chvl [ t ] = chvl0 if chvl [ t ] < chvl0 else chvl [ t ]
    l [ t ] = chvl [ t ]

    # LOG Conc. of Hydrogen Interspace
    lche = pkc + lch(ce_hco3 [ t ], ce_h2co3 [ t ])
    ce_h [ t ] = 10 ** (-lche)

    # ELECTROLYTE TRANSPORT ACROSS A SIMPLE EPITHELIUM of the rat proximal tubule (1992)
    # The cell volume is presented as a function of the inverse cell impermeant species concentration
    clvl [ t ] = clvl0 * imp0 / imp [ t ]
    clvl [ t ] = clvl0 if imp [ t ] == 0 else clvl [ t ]
    clvl_imp = clvl0 * imp0 / imp [ t ]
    # The cell height (summation of Extracellular channel volume
    # and Intracellular compartment volume[cm3/cm2 epithelium])
    l [ t ] = l [ t ] + clvl [ t ]
    p_i = pm

    # LOG Conc. of Hydrogen Cell
    lchi = pkc + lch(ci_hco3 [ t ], ci_h2co3 [ t ])
    ci_h [ t ] = 10 ** (-lchi)

    # ELECTROLYTE TRANSPORT ACROSS A SIMPLE EPITHELIUM 1979
    # Applying the Kedem and Katchalsky relations:
    # The membrane transport for nonelectrolyte solutions generated by the hydrostatic
    # pressure difference and the osmotic pressure difference
    # The equations below represent the volume flow or the convective volume flow in between the different membrane
    fevm = ame * lpme * (pm - pe[t] - rt * impm) / rt
    fevs = ae[t] * lpes * (rt * imps + pe[t] - ps) / rt
    fivm = ami * lpmi * (rt * imp[t] - rt * impm + pm - p_i) / rt
    fivs = ais * lpis * (p_i - ps + rt * imps - rt * imp[t]) / rt
    jv =  aie * lpie *  (p_i - pe[t] - rt * imp[t]) / rt
    # Updating the volume flow
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
                    ce_h [t]- cm_h) + sme_hco2 * (
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

    jv = jv + lpie * aie * (
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

    Cellular_water_fluxes = fivm - fivs - jv


    # Convective Intraepithelial Solute Fluxes
    param_csf = 1
    fekm_na = csf(ce_na [ t ], cm_na, fevm, sme_na, param_csf)
    fekm_na_csf = csf(ce_na [ t ], cm_na, fevm, sme_na, param_csf)
    fekm_k = csf(ce_k [ t ], cm_k, fevm, sme_k, param_csf)
    fekm_cl = csf(ce_cl [ t ], cm_cl, fevm, sme_cl, param_csf)
    fekm_hco3 = csf(ce_hco3 [ t ], cm_hco3, fevm, sme_hco3, param_csf)
    fekm_h2co3 = csf(ce_h2co3 [ t ], cm_h2co3, fevm, sme_h2co3, param_csf)
    fekm_co2 = csf(ce_co2 [ t ], cm_co2, fevm, sme_co2, param_csf)
    fekm_hpo4 = csf(ce_hpo4 [ t ], cm_hpo4, fevm, sme_hpo4, param_csf)
    fekm_h2po4 = csf(ce_h2po4 [ t ], cm_h2po4, fevm, sme_h2po4, param_csf)
    fekm_urea = csf(ce_urea [ t ], cm_urea, fevm, sme_urea, param_csf)
    fekm_nh3 = csf(ce_nh3 [ t ], cm_nh3, fevm, sme_nh3, param_csf)
    fekm_nh4 = csf(ce_nh4 [ t ], cm_nh4, fevm, sme_nh4, param_csf)
    fekm_h = csf(ce_h [ t ], cm_h, fevm, sme_h, param_csf)
    fekm_hco2 = csf(ce_hco2 [ t ], cm_hco2, fevm, sme_hco2, param_csf)
    fekm_h2co2 = csf(ce_h2co2 [ t ], cm_h2co2, fevm, sme_h2co2, param_csf)
    fekm_gluc = csf(ce_gluc [ t ], cm_gluc, fevm, sme_gluc, param_csf)

    feks_na = csf(ce_na [ t ], cs_na, fevs, ses_na, param_csf)
    feks_na_csf = csf(ce_na [ t ], cs_na, fevs, ses_na, param_csf)
    feks_k = csf(ce_k [ t ], cs_k, fevs, ses_k, param_csf)
    feks_cl = csf(ce_cl [ t ], cs_cl, fevs, ses_cl, param_csf)
    feks_hco3 = csf(ce_hco3 [ t ], cs_hco3, fevs, ses_hco3, param_csf)
    feks_h2co3 = csf(ce_h2co3 [ t ], cs_h2co3, fevs, ses_h2co3, param_csf)
    feks_co2 = csf(ce_co2 [ t ], cs_co2, fevs, ses_co2, param_csf)
    feks_hpo4 = csf(ce_hpo4 [ t ], cs_hpo4, fevs, ses_hpo4, param_csf)
    feks_h2po4 = csf(ce_h2po4 [ t ], cs_h2po4, fevs, ses_h2po4, param_csf)
    feks_urea = csf(ce_urea [ t ], cs_urea, fevs, ses_urea, param_csf)
    feks_nh3 = csf(ce_nh3 [ t ], cs_nh3, fevs, ses_nh3, param_csf)
    feks_nh4 = csf(ce_nh4 [ t ], cs_nh4, fevs, ses_nh4, param_csf)
    feks_h = csf(ce_h [ t ], cs_h, fevs, ses_h, param_csf)
    feks_hco2 = csf(ce_hco2 [ t ], cs_hco2, fevs, ses_hco2, param_csf)
    feks_h2co2 = csf(ce_h2co2 [ t ], cs_h2co2, fevs, ses_h2co2, param_csf)
    feks_gluc = csf(ce_gluc [ t ], cs_gluc, fevs, ses_gluc, param_csf)

    fikm_na = csf(ci_na [ t ], cm_na, fivm, smi_na, param_csf)
    fikm_na_csf = csf(ci_na [ t ], cm_na, fivm, smi_na, param_csf)
    fikm_k = csf(ci_k [ t ], cm_k, fivm, smi_k, param_csf)
    fikm_cl = csf(ci_cl [ t ], cm_cl, fivm, smi_cl, param_csf)
    fikm_hco3 = csf(ci_hco3 [ t ], cm_hco3, fivm, smi_hco3, param_csf)
    fikm_h2co3 = csf(ci_h2co3 [ t ], cm_h2co3, fivm, smi_h2co3, param_csf)
    fikm_co2 = csf(ci_co2 [ t ], cm_co2, fivm, smi_co2, param_csf)
    fikm_hpo4 = csf(ci_hpo4 [ t ], cm_hpo4, fivm, smi_hpo4, param_csf)
    fikm_h2po4 = csf(ci_h2po4 [ t ], cm_h2po4, fivm, smi_h2po4, param_csf)
    fikm_urea = csf(ci_urea [ t ], cm_urea, fivm, smi_urea, param_csf)
    fikm_nh3 = csf(ci_nh3 [ t ], cm_nh3, fivm, smi_nh3, param_csf)
    fikm_nh4 = csf(ci_nh4 [ t ], cm_nh4, fivm, smi_nh4, param_csf)
    fikm_h = csf(ci_h [ t ], cm_h, fivm, smi_h, param_csf)
    fikm_hco2 = csf(ci_hco2 [ t ], cm_hco2, fivm, smi_hco2, param_csf)
    fikm_h2co2 = csf(ci_h2co2 [ t ], cm_h2co2, fivm, smi_h2co2, param_csf)
    fikm_gluc = csf(ci_gluc [ t ], cm_gluc, fivm, smi_gluc, param_csf)

    fiks_na = csf(ci_na [ t ], cs_na, fivs, sis_na, param_csf)
    fiks_na_csf = csf(ci_na [ t ], cs_na, fivs, sis_na, param_csf)
    fiks_k = csf(ci_k [ t ], cs_k, fivs, sis_k, param_csf)
    fiks_cl = csf(ci_cl [ t ], cs_cl, fivs, sis_cl, param_csf)
    fiks_hco3 = csf(ci_hco3 [ t ], cs_hco3, fivs, sis_hco3, param_csf)
    fiks_h2co3 = csf(ci_h2co3 [ t ], cs_h2co3, fivs, sis_h2co3, param_csf)
    fiks_co2 = csf(ci_co2 [ t ], cs_co2, fivs, sis_co2, param_csf)
    fiks_hpo4 = csf(ci_hpo4 [ t ], cs_hpo4, fivs, sis_hpo4, param_csf)
    fiks_h2po4 = csf(ci_h2po4 [ t ], cs_h2po4, fivs, sis_h2po4, param_csf)
    fiks_urea = csf(ci_urea [ t ], cs_urea, fivs, sis_urea, param_csf)
    fiks_nh3 = csf(ci_nh3 [ t ], cs_nh3, fivs, sis_nh3, param_csf)
    fiks_nh4 = csf(ci_nh4 [ t ], cs_nh4, fivs, sis_nh4, param_csf)
    fiks_h = csf(ci_h [ t ], cs_h, fivs, sis_h, param_csf)
    fiks_hco2 = csf(ci_hco2 [ t ], cs_hco2, fivs, sis_hco2, param_csf)
    fiks_h2co2 = csf(ci_h2co2 [ t ], cs_h2co2, fivs, sis_h2co2, param_csf)
    fiks_gluc = csf(ci_gluc [ t ], cs_gluc, fivs, sis_gluc, param_csf)

    jk_na = csf(ci_na [ t ], ce_na [ t ], jv, sis_na, param_csf)
    jk_na_csf = csf(ci_na [ t ], ce_na [ t ], jv, sis_na, param_csf)
    jk_k = csf(ci_k [ t ], ce_k [ t ], jv, sis_k, param_csf)
    jk_cl = csf(ci_cl [ t ], ce_cl [ t ], jv, sis_cl, param_csf)
    jk_hco3 = csf(ci_hco3 [ t ], ce_hco3 [ t ], jv, sis_hco3, param_csf)
    jk_h2co3 = csf(ci_h2co3 [ t ], ce_h2co3 [ t ], jv, sis_h2co3, param_csf)
    jk_co2 = csf(ci_co2 [ t ], ce_co2 [ t ], jv, sis_co2, param_csf)
    jk_hpo4 = csf(ci_hpo4 [ t ], ce_hpo4 [ t ], jv, sis_hpo4, param_csf)
    jk_h2po4 = csf(ci_h2po4 [ t ], ce_h2po4 [ t ], jv, sis_h2po4, param_csf)
    jk_urea = csf(ci_urea [ t ], ce_urea [ t ], jv, sis_urea, param_csf)
    jk_nh3 = csf(ci_nh3 [ t ], ce_nh3 [ t ], jv, sis_nh3, param_csf)
    jk_nh4 = csf(ci_nh4 [ t ], ce_nh4 [ t ], jv, sis_nh4, param_csf)
    jk_h = csf(ci_h [ t ], ce_h [ t ], jv, sis_h, param_csf)
    jk_hco2 = csf(ci_hco2 [ t ], ce_hco2 [ t ], jv, sis_hco2, param_csf)
    jk_h2co2 = csf(ci_h2co2 [ t ], ce_h2co2 [ t ], jv, sis_h2co2, param_csf)
    jk_gluc = csf(ci_gluc [ t ], ce_gluc [ t ], jv, sis_gluc, param_csf)

    # Goldman Fluxes: fluxes due to diffusion and electrical potential difference for
    # all of the ions that are permeant through the membrane
    param_goldman = 1
    fekm_na = fekm_na + goldman(hme_na, ame, z_na, vm [ t ], ve [ t ], cm_na,
                                ce_na [ t ], param_goldman)
    fekm_na_goldman = goldman(hme_na, ame, z_na, vm [ t ], ve [ t ], cm_na,
                      ce_na [ t ], param_goldman)
    fekm_k = fekm_k + goldman(hme_k, ame, z_k, vm [ t ], ve [ t ], cm_k, ce_k [ t ], param_goldman)
    fekm_cl = fekm_cl + goldman(hme_cl, ame, z_cl, vm [ t ], ve [ t ], cm_cl, ce_cl [ t ], param_goldman)
    fekm_hco3 = fekm_hco3 + goldman(hme_hco3, ame, z_hco3, vm [ t ], ve [ t ], cm_hco3, ce_hco3 [ t ], param_goldman)
    fekm_h2co3 = fekm_h2co3 + goldman(hme_h2co3, ame, z_h2co3, vm [ t ], ve [ t ], cm_h2co3,
                                      ce_h2co3 [ t ], param_goldman)
    fekm_co2 = fekm_co2 + goldman(hme_co2, ame, z_co2, vm [ t ], ve [ t ], cm_co2, ce_co2 [ t ], param_goldman)
    fekm_hpo4 = fekm_hpo4 + goldman(hme_hpo4, ame, z_hpo4, vm [ t ], ve [ t ], cm_hpo4,
                                    ce_hpo4 [ t ], param_goldman)
    fekm_h2po4 = fekm_h2po4 + goldman(hme_h2po4, ame, z_h2po4, vm [ t ], ve [ t ], cm_h2po4,
                                      ce_h2po4 [ t ], param_goldman)
    fekm_urea = fekm_urea + goldman(hme_urea, ame, z_urea, vm [ t ], ve [ t ], cm_urea, ce_urea [ t ], param_goldman)
    fekm_nh3 = fekm_nh3 + goldman(hme_nh3, ame, z_nh3, vm [ t ], ve [ t ], cm_nh3, ce_nh3 [ t ], param_goldman)
    fekm_nh4 = fekm_nh4 + goldman(hme_nh4, ame, z_nh4, vm [ t ], ve [ t ], cm_nh4, ce_nh4 [ t ], param_goldman)
    fekm_h = fekm_h + goldman(hme_h, ame, z_h, vm [ t ], ve [ t ], cm_h, ce_h [ t ], param_goldman)
    fekm_hco2 = fekm_hco2 + goldman(hme_hco2, ame, z_hco2, vm [ t ], ve [ t ], cm_hco2, ce_hco2 [ t ], param_goldman)
    fekm_h2co2 = fekm_h2co2 + goldman(hme_h2co2, ame, z_h2co2, vm [ t ], ve [ t ], cm_h2co2, ce_h2co2 [ t ],
                                      param_goldman)
    fekm_gluc = fekm_gluc + goldman(hme_gluc, ame, z_gluc, vm [ t ], ve [ t ], cm_gluc, ce_gluc [ t ], param_goldman)

    feks_na = feks_na + goldman(hes_na, ae [ t ], z_na, ve [ t ], vs, ce_na [ t ],
                                cs_na, param_goldman)
    feks_na_goldman = goldman(hes_na, ae [ t ], z_na, ve [ t ], vs, ce_na [ t ],
                                cs_na, param_goldman)
    feks_k = feks_k + goldman(hes_k, ae [ t ], z_k, ve [ t ], vs, ce_k [ t ], cs_k, param_goldman)
    feks_cl = feks_cl + goldman(hes_cl, ae [ t ], z_cl, ve [ t ], vs, ce_cl [ t ], cs_cl, param_goldman)
    feks_hco3 = feks_hco3 + goldman(hes_hco3, ae [ t ], z_hco3, ve [ t ], vs, ce_hco3 [ t ], cs_hco3, param_goldman)
    feks_h2co3 = feks_h2co3 + goldman(hes_h2co3, ae [ t ], z_h2co3, ve [ t ], vs, ce_h2co3 [ t ],
                                      cs_h2co3, param_goldman)
    feks_co2 = feks_co2 + goldman(hes_co2, ae [ t ], z_co2, ve [ t ], vs, ce_co2 [ t ], cs_co2, param_goldman)
    feks_hpo4 = feks_hpo4 + goldman(hes_hpo4, ae [ t ], z_hpo4, ve [ t ], vs, ce_hpo4 [ t ],
                                    cs_hpo4, param_goldman)
    feks_h2po4 = feks_h2po4 + goldman(hes_h2po4, ae [ t ], z_h2po4, ve [ t ], vs, ce_h2po4 [ t ],
                                      cs_h2po4, param_goldman)
    feks_urea = feks_urea + goldman(hes_urea, ae [ t ], z_urea, ve [ t ], vs, ce_urea [ t ], cs_urea, param_goldman)
    feks_nh3 = feks_nh3 + goldman(hes_nh3, ae [ t ], z_nh3, ve [ t ], vs, ce_nh3 [ t ], cs_nh3, param_goldman)
    feks_nh4 = feks_nh4 + goldman(hes_nh4, ae [ t ], z_nh4, ve [ t ], vs, ce_nh4 [ t ], cs_nh4, param_goldman)
    feks_h = feks_h + goldman(hes_h, ae [ t ], z_h, ve [ t ], vs, ce_h [ t ], cs_h, param_goldman)
    feks_hco2 = feks_hco2 + goldman(hes_hco2, ae [ t ], z_hco2, ve [ t ], vs, ce_hco2 [ t ], cs_hco2, param_goldman)

    feks_h2co2 = feks_h2co2 + goldman(hme_h2co2, ame, z_h2co2, ve [ t ], vs, ce_h2co2 [ t ], cs_h2co2, param_goldman)
    feks_gluc = feks_gluc + goldman(hes_gluc, ae [ t ], z_gluc, ve [ t ], vs, ce_gluc [ t ], cs_gluc, param_goldman)

    fikm_na = fikm_na + goldman(hmi_na, ami, z_na, vm [ t ], vi [ t ], cm_na,
                                ci_na [ t ], param_goldman)
    fikm_na_goldman = goldman(hmi_na, ami, z_na, vm [ t ], vi [ t ], cm_na,
                      ci_na [ t ], param_goldman)
    fikm_k = fikm_k + goldman(hmi_k, ami, z_k, vm [ t ], vi [ t ], cm_k, ci_k [ t ], param_goldman)
    fikm_cl = fikm_cl + goldman(hmi_cl, ami, z_cl, vm [ t ], vi [ t ], cm_cl, ci_cl [ t ], param_goldman)
    fikm_hco3 = fikm_hco3 + goldman(hmi_hco3, ami, z_hco3, vm [ t ], vi [ t ], cm_hco3, ci_hco3 [ t ], param_goldman)
    fikm_h2co3 = fikm_h2co3 + goldman(hmi_h2co3, ami, z_h2co3, vm [ t ], vi [ t ], cm_h2co3,
                                      ci_h2co3 [ t ], param_goldman)
    fikm_co2 = fikm_co2 + goldman(hmi_co2, ami, z_co2, vm [ t ], vi [ t ], cm_co2, ci_co2 [ t ], param_goldman)
    fikm_hpo4 = fikm_hpo4 + goldman(hmi_hpo4, ami, z_hpo4, vm [ t ], vi [ t ], cm_hpo4,
                                    ci_hpo4 [ t ], param_goldman)
    fikm_h2po4 = fikm_h2po4 + goldman(hmi_h2po4, ami, z_h2po4, vm [ t ], vi [ t ], cm_h2po4,
                                      ci_h2po4 [ t ], param_goldman)
    fikm_urea = fikm_urea + goldman(hmi_urea, ami, z_urea, vm [ t ], vi [ t ], cm_urea, ci_urea [ t ], param_goldman)
    fikm_nh3 = fikm_nh3 + goldman(hmi_nh3, ami, z_nh3, vm [ t ], vi [ t ], cm_nh3, ci_nh3 [ t ], param_goldman)
    fikm_nh4 = fikm_nh4 + goldman(hmi_nh4, ami, z_nh4, vm [ t ], vi [ t ], cm_nh4, ci_nh4 [ t ], param_goldman)
    fikm_h = fikm_h + goldman(hmi_h, ami, z_h, vm [ t ], vi [ t ], cm_h, ci_h [ t ], param_goldman)
    fikm_hco2 = fikm_hco2 + goldman(hmi_hco2, ami, z_hco2, vm [ t ], vi [ t ], cm_hco2, ci_hco2 [ t ], param_goldman)
    fikm_h2co2 = fikm_h2co2 + goldman(hmi_h2co2, ami, z_h2co2, vm [ t ], vi [ t ], cm_h2co2,
                                      ci_h2co2 [ t ], param_goldman)
    fikm_gluc = fikm_gluc + goldman(hmi_gluc, ami, z_gluc, vm [ t ], vi [ t ], cm_gluc, ci_gluc [ t ], param_goldman)

    jk_na = jk_na + goldman(his_na, aie, z_na, vi [ t ], ve [ t ], ci_na [ t ], ce_na [ t ], param_goldman)
    jk_na_goldman = goldman(his_na, aie, z_na, vi [ t ], ve [ t ], ci_na [ t ], ce_na [ t ], param_goldman)
    jk_k = jk_k + goldman(his_k, aie, z_k, vi [ t ], ve [ t ], ci_k [ t ], ce_k [ t ], param_goldman)
    jk_cl = jk_cl + goldman(his_cl, aie, z_cl, vi [ t ], ve [ t ], ci_cl [ t ], ce_cl [ t ], param_goldman)
    jk_hco3 = jk_hco3 + goldman(his_hco3, aie, z_hco3, vi [ t ], ve [ t ], ci_hco3 [ t ], ce_hco3 [ t ], param_goldman)
    jk_h2co3 = jk_h2co3 + goldman(his_h2co3, aie, z_h2co3, vi [ t ], ve [ t ], ci_h2co3 [ t ],
                                  ce_h2co3 [ t ], param_goldman)
    jk_co2 = jk_co2 + goldman(his_co2, aie, z_co2, vi [ t ], ve [ t ], ci_co2 [ t ], ce_co2 [ t ], param_goldman)
    jk_hpo4 = jk_hpo4 + goldman(his_hpo4, aie, z_hpo4, vi [ t ], ve [ t ], ci_hpo4 [ t ],
                                ce_hpo4 [ t ], param_goldman)
    jk_h2po4 = jk_h2po4 + goldman(his_h2po4, aie, z_h2po4, vi [ t ], ve [ t ], ci_h2po4 [ t ],
                                  ce_h2po4 [ t ], param_goldman)
    jk_urea = jk_urea + goldman(his_urea, aie, z_urea, vi [ t ], ve [ t ], ci_urea [ t ], ce_urea [ t ], param_goldman)
    jk_nh3 = jk_nh3 + goldman(his_nh3, aie, z_nh3, vi [ t ], ve [ t ], ci_nh3 [ t ], ce_nh3 [ t ], param_goldman)
    jk_nh4 = jk_nh4 + goldman(his_nh4, aie, z_nh4, vi [ t ], ve [ t ], ci_nh4 [ t ], ce_nh4 [ t ], param_goldman)
    jk_h = jk_h + goldman(his_h, aie, z_h, vi [ t ], ve [ t ], ci_h [ t ], ce_h [ t ], param_goldman)
    jk_hco2 = jk_hco2 + goldman(his_hco2, aie, z_hco2, vi [ t ], ve [ t ], ci_hco2 [ t ], ce_hco2 [ t ], param_goldman)
    jk_h2co2 = jk_h2co2 + goldman(his_h2co2, aie, z_h2co2, vi [ t ], ve [ t ], ci_h2co2 [ t ], ce_h2co2 [ t ],
                                  param_goldman)
    jk_gluc = jk_gluc + goldman(his_gluc, aie, z_gluc, vi [ t ], ve [ t ], ci_gluc [ t ], ce_gluc [ t ], param_goldman)

    fiks_na = fiks_na + goldman(his_na, ais, z_na, vi [ t ], vs, ci_na [ t ], cs_na, param_goldman)
    fiks_na_goldman = goldman(his_na, ais, z_na, vi [ t ], vs, ci_na [ t ], cs_na, param_goldman)
    fiks_k = fiks_k + goldman(his_k, ais, z_k, vi [ t ], vs, ci_k [ t ], cs_k, param_goldman)
    fiks_cl = fiks_cl + goldman(his_cl, ais, z_cl, vi [ t ], vs, ci_cl [ t ], cs_cl, param_goldman)
    fiks_hco3 = fiks_hco3 + goldman(his_hco3, ais, z_hco3, vi [ t ], vs, ci_hco3 [ t ], cs_hco3, param_goldman)
    fiks_h2co3 = fiks_h2co3 + goldman(his_h2co3, ais, z_h2co3, vi [ t ], vs, ci_h2co3 [ t ],
                                      cs_h2co3, param_goldman)
    fiks_co2 = fiks_co2 + goldman(his_co2, ais, z_co2, vi [ t ], vs, ci_co2 [ t ], cs_co2, param_goldman)
    fiks_hpo4 = fiks_hpo4 + goldman(his_hpo4, ais, z_hpo4, vi [ t ], vs, ci_hpo4 [ t ],
                                    cs_hpo4, param_goldman)
    fiks_h2po4 = fiks_h2po4 + goldman(his_h2po4, ais, z_h2po4, vi [ t ], vs, ci_h2po4 [ t ],
                                      cs_h2po4, param_goldman)
    fiks_urea = fiks_urea + goldman(his_urea, ais, z_urea, vi [ t ], vs, ci_urea [ t ], cs_urea, param_goldman)
    fiks_nh3 = fiks_nh3 + goldman(his_nh3, ais, z_nh3, vi [ t ], vs, ci_nh3 [ t ], cs_nh3, param_goldman)
    fiks_nh4 = fiks_nh4 + goldman(his_nh4, ais, z_nh4, vi [ t ], vs, ci_nh4 [ t ], cs_nh4, param_goldman)
    fiks_h = fiks_h + goldman(his_h, ais, z_h, vi [ t ], vs, ci_h [ t ], cs_h, param_goldman)
    fiks_hco2 = fiks_hco2 + goldman(his_hco2, ais, z_hco2, vi [ t ], vs, ci_hco2 [ t ], cs_hco2, param_goldman)
    fiks_h2co2 = fiks_h2co2 + goldman(his_h2co2, ais, z_h2co2, vi [ t ], vs, ci_h2co2 [ t ],
                                      cs_h2co2, param_goldman)
    fiks_gluc = fiks_gluc + goldman(his_gluc, ais, z_gluc, vi [ t ], vs, ci_gluc [ t ], cs_gluc, param_goldman)

    # net transporters on mi bourder
    
    sglt = sglt_mi(cm_na, ci_na [ t ], cm_gluc, ci_gluc [ t ], z_na, z_gluc, vm [ t ], vi [ t ], ami, lmi_nagluc,
                   param_sglt)
    na_mi_nagluc = sglt [ 0 ]
    gluc_mi_nagluc = sglt [ 1 ]

    nah2po4 = nah2po4_mi(cm_na, ci_na [ t ], cm_h2po4, ci_h2po4 [ t ], z_na, z_h2po4, vm [ t ], vi [ t ], ami,
                         lmi_nah2po4, param_h2po4)
    na_mi_nah2po4 = nah2po4 [ 0 ]
    h2po4_mi_nah2po4 = nah2po4 [ 1 ]
    
    clhco3 = clhco3_mi(cm_cl, ci_cl [ t ], cm_hco3, ci_hco3 [ t ], z_cl, z_hco3, vm [ t ], vi [ t ], ami,
                       lmi_clhco3, param_clhco3)
    cl_mi_clhco3 = clhco3 [ 0 ]
    hco3_mi_clhco3 = clhco3 [ 1 ]

    clhco2 = clhco2_mi(cm_cl, ci_cl [ t ], cm_hco2, ci_hco2 [ t ], z_cl, z_hco2, vm [ t ], vi [ t ], ami,
                       lmi_clhco2, param_clhco2)
    cl_mi_clhco2 = clhco2 [ 0 ]
    hco2_mi_clhco2 = clhco2 [ 1 ]
    # the nah exchanger translate concentrations to the nah model on  mi bourder

    nah_mi = nah(cm_na, ci_na[t], cm_h, ci_h[t], z_na, z_h, vm[t], vi[t], ami, lmi_nah, param_nah_mi)
    na_mi_nah = nah_mi[0]
    h_mi_nah = nah_mi[1]

    
    nanh4_mi = nanh4(cm_na, ci_na[t], cm_nh4, ci_nh4[t], z_na, z_nh4, vm[t], vi[t], ami, lmi_nanh4, param_nanh4_mi)
    na_mi_nanh4 = nanh4_mi[0]
    nh4_mi_nanh4 = nanh4_mi[1]

    # the nhe3 exchanger translate concentrations to the nhe3 model on  mi bourder
    
    mynah = nhe3(ci_h [ t ], ci_na [ t ], ci_nh4 [ t ], cm_h, cm_na, cm_nh4, param_nhe3)
    jnah_na = mynah [ 0 ]
    jnah_h = mynah [ 1 ]
    jnah_nh4 = mynah [ 2 ]
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
    # net cotransporters on is bourder
    
    nahco3 = na_hco3(ci_na [ t ], cs_na, ci_hco3 [ t ], cs_hco3, z_na, z_hco3, vi [ t ], vs, ais,
                     lis_nahco3, param_nahco3)
    na_is_nahco3 = nahco3 [ 0 ]
    hco3_is_nahco3 = nahco3 [ 1 ]
    kcl = k_cl(ci_k [ t ], cs_k, ci_cl [ t ], cs_cl, z_k, z_cl, vi [ t ], vs, ais, lis_kcl, param_kcl)
    k_is_kcl = kcl [ 0 ]
    cl_is_kcl = kcl [ 1 ]    
    na_clhco3 = na_cl_hco3(ci_na [ t ], cs_na, ci_cl [ t ], cs_cl, ci_hco3 [ t ], cs_hco3, z_na, z_cl, z_hco3,
                           vi [ t ], vs, ais, lis_na_clhco3, param_na_cl_hco3)
    na_is_na_clhco3 = na_clhco3 [ 0 ]
    cl_is_na_clhco3 = na_clhco3 [ 1 ]
    hco3_is_na_clhco3 = na_clhco3 [ 2 ]

    fiks_na = fiks_na + na_is_nahco3 + na_is_na_clhco3
    fiks_k = fiks_k + k_is_kcl
    fiks_cl = fiks_cl + cl_is_kcl + cl_is_na_clhco3
    fiks_hco3 = fiks_hco3 + hco3_is_na_clhco3 + hco3_is_nahco3

    nahco3_aie = na_hco3(ci_na [ t ], ce_na [ t ], ci_hco3 [ t ], ce_hco3 [ t ], z_na, z_hco3, vi [ t ], ve [ t ],
                           aie,
                           lis_nahco3, param_nahco3)
    na_aie_nahco3 = nahco3_aie [ 0 ]
    hco3_aie_nahco3 = nahco3_aie [ 1 ]

    kcl_aie = k_cl(ci_k [ t ], ce_k [ t ], ci_cl [ t ], ce_cl [ t ], z_k, z_cl, vi [ t ], ve [ t ], aie, lis_kcl, param_kcl)
    k_aie_kcl = kcl_aie [ 0 ]
    cl_aie_kcl = kcl_aie [ 1 ]

    na_clhco3_aie = na_cl_hco3(ci_na [ t ], ce_na [ t ], ci_cl [ t ], ce_cl [ t ], ci_hco3 [ t ], ce_hco3 [ t ], z_na,
                                 z_cl, z_hco3,
                                 vi [ t ], ve [ t ], aie, lis_na_clhco3, param_na_cl_hco3)
    na_aie_na_clhco3 = na_clhco3_aie [ 0 ]
    cl_aie_na_clhco3 = na_clhco3_aie [ 1 ]
    hco3_aie_na_clhco3 = na_clhco3_aie [ 2 ]

    jk_na = jk_na + na_aie_nahco3 + na_aie_na_clhco3
    jk_k = jk_k + k_aie_kcl
    jk_cl = jk_cl + cl_aie_kcl + cl_aie_na_clhco3
    jk_hco3 = jk_hco3 + hco3_aie_na_clhco3 + hco3_aie_nahco3

    nak = nak_atp(ci_na [ t ], cs_na, ci_k [ t ], cs_k,  cs_nh4, param_sodium_pumps)
    atis_na = nak [ 0 ]
    atis_k = nak [ 1 ]
    atis_nh4 = nak [ 2 ]

    fiks_na = fiks_na + ais * atis_na
    fiks_k = fiks_k + ais * atis_k
    fiks_nh4 = fiks_nh4 + ais * atis_nh4
    # sodium pumps on ie bourder
    nak = nak_atp(ci_na [ t ], ce_na[t], ci_k [ t ], ce_k[t],  ce_nh4[t], param_sodium_pumps)
    atie_na = nak [ 0 ]
    atie_k = nak [ 1 ]
    atie_nh4 = nak [ 2 ]

    jk_na = jk_na + aie * atie_na
    jk_k = jk_k + aie * atie_k
    jk_nh4 = jk_nh4 + aie * atie_nh4

    # proton pumps
    
    atmi_h = h_atp_mi(cm_h, ci_h [ t ], vm [ t ], vi [ t ], z_h, ami, param_H_pumps)

    # proton pumps flux update
    fikm_h = fikm_h + atmi_h
    # establish the error vectors, the "phi" array.
    # first for the interspace electroneutrality
    phie_en = 0
    phie_en = phie_en + z_na * ce_na [ t ] + z_k * ce_k [ t ] + z_cl * ce_cl [
        t ] + z_hco3 * ce_hco3 [ t ] + z_h2co3 * ce_h2co3 [ t ] + z_co2 * ce_co2 [ t ] + z_hpo4 * ce_hpo4 [
                  t ] + z_h2po4 * ce_h2po4 [ t ] + z_urea * ce_urea [ t ] + z_nh3 * ce_nh3 [ t ] + z_nh4 * ce_nh4 [
                  t ] + z_h * ce_h [ t ] + z_hco2 * ce_hco2 [ t ] + z_h2co2 * ce_h2co2 [ t ] + z_gluc * ce_gluc [ t ]

    phii_en = imp [ t ] * zimp
    phii_en = phii_en - cbuf [ t ] + z_na * ci_na [ t ] + z_k * ci_k [ t ] \
              + z_cl * ci_cl [ t ] + z_hco3 * ci_hco3 [ t ] + z_h2co3 * ci_h2co3 [ t ] + z_co2 * ci_co2 [ t ] + z_hpo4 * \
              ci_hpo4 [ t ] + z_h2po4 * ci_h2po4 [ t ] + z_urea * ci_urea [ t ] + z_nh3 * ci_nh3 [ t ] + z_nh4 * \
              ci_nh4 [ t ] + z_h * ci_h [ t ] + z_hco2 * ci_hco2 [ t ] + z_h2co2 * ci_h2co2 [ t ] + z_gluc * ci_gluc [
                  t ]

    # mass conservation in the Time-dependent case
    phie_vlm = fevs - fevm - jv + rtau * (chvl [ t ] - chvl [ t - 1 ])
    qe_na = feks_na - fekm_na - jk_na + rtau * (ce_na [ t ] * chvl [ t ] - ce_na [ t - 1 ] * chvl [ t - 1 ])
    qe_k = feks_k - fekm_k - jk_k + rtau * (ce_k [ t ] * chvl [ t ] - ce_k [ t - 1 ] * chvl [ t - 1 ])
    qe_cl = feks_cl - fekm_cl - jk_cl + rtau * (ce_cl [ t ] * chvl [ t ] - ce_cl [ t - 1 ] * chvl [ t - 1 ])
    qe_hco3 = feks_hco3 - fekm_hco3 - jk_hco3 + rtau * (ce_hco3 [ t ] * chvl [ t ] - ce_hco3 [ t - 1 ] * chvl [ t - 1 ])
    qe_h2co3 = feks_h2co3 - fekm_h2co3 - jk_h2co3 + rtau * (
            ce_h2co3 [ t ] * chvl [ t ] - ce_h2co3 [ t - 1 ] * chvl [ t - 1 ])
    qe_co2 = feks_co2 - fekm_co2 - jk_co2 + rtau * (ce_co2 [ t ] * chvl [ t ] - ce_co2 [ t - 1 ] * chvl [ t - 1 ])
    qe_hpo4 = feks_hpo4 - fekm_hpo4 - jk_hpo4 + rtau * (ce_hpo4 [ t ] * chvl [ t ] - ce_hpo4 [ t - 1 ] * chvl [ t - 1 ])
    qe_h2po4 = feks_h2po4 - fekm_h2po4 - jk_h2po4 + rtau * (
            ce_h2po4 [ t ] * chvl [ t ] - ce_h2po4 [ t - 1 ] * chvl [ t - 1 ])
    qe_urea = feks_urea - fekm_urea - jk_urea + rtau * (
            ce_urea [ t ] * chvl [ t ] - ce_urea [ t - 1 ] * chvl [ t - 1 ])
    qe_nh3 = feks_nh3 - fekm_nh3 - jk_nh3 + rtau * (ce_nh3 [ t ] * chvl [ t ] - ce_nh3 [ t - 1 ] * chvl [ t - 1 ])
    qe_nh4 = feks_nh4 - fekm_nh4 - jk_nh4 + rtau * (ce_nh4 [ t ] * chvl [ t ] - ce_nh4 [ t - 1 ] * chvl [ t - 1 ])
    qe_h = feks_h - fekm_h - jk_h + rtau * (ce_h [ t ] * chvl [ t ] - ce_h [ t - 1 ] * chvl [ t - 1 ])
    qe_hco2 = feks_hco2 - fekm_hco2 - jk_hco2 + rtau * (ce_hco2 [ t ] * chvl [ t ] - ce_hco2 [ t - 1 ] * chvl [ t - 1 ])
    qe_h2co2 = feks_h2co2 - fekm_h2co2 - jk_h2co2 + rtau * (
            ce_h2co2 [ t ] * chvl [ t ] - ce_h2co2 [ t - 1 ] * chvl [ t - 1 ])

    qe_gluc = feks_gluc - fekm_gluc - jk_gluc + rtau * (
            ce_gluc [ t ] * chvl [ t ] - ce_gluc [ t - 1 ] * chvl [ t - 1 ])

    # mass conservation in the Time - dependent case
    phii_vlm = fivs - fivm + jv + rtau * (clvl [ t ] - clvl [ t - 1 ])

    qi_na = fiks_na - fikm_na + jk_na + rtau * (
            ci_na [ t ] * clvl [ t ] - ci_na [ t - 1 ] * clvl [ t - 1 ])

    qi_k = fiks_k - fikm_k + jk_k + rtau * (
            ci_k [ t ] * clvl [ t ] - ci_k [ t - 1 ] * clvl [ t - 1 ])
    qi_cl = fiks_cl - fikm_cl + jk_cl + rtau * (
            ci_cl [ t ] * clvl [ t ] - ci_cl [ t - 1 ] * clvl [ t - 1 ])

    qi_hco3 = fiks_hco3 - fikm_hco3 + jk_hco3 + rtau * (
            ci_hco3 [ t ] * clvl [ t ] - ci_hco3 [ t - 1 ] * clvl [ t - 1 ])
    qi_h2co3 = fiks_h2co3 - fikm_h2co3 + jk_h2co3 + rtau * (
            ci_h2co3 [ t ] * clvl [ t ] - ci_h2co3 [ t - 1 ] * clvl [ t - 1 ])
    qi_co2 = fiks_co2 - fikm_co2 + jk_co2 + rtau * (
            ci_co2 [ t ] * clvl [ t ] - ci_co2 [ t - 1 ] * clvl [ t - 1 ])
    qi_hpo4 = fiks_hpo4 - fikm_hpo4 + jk_hpo4 + rtau * (
            ci_hpo4 [ t ] * clvl [ t ] - ci_hpo4 [ t - 1 ] * clvl [ t - 1 ])
    qi_h2po4 = fiks_h2po4 - fikm_h2po4 + jk_h2po4 + rtau * (
            ci_h2po4 [ t ] * clvl [ t ] - ci_h2po4 [ t - 1 ] * clvl [ t - 1 ])

    qi_urea = fiks_urea - fikm_urea + jk_urea + rtau * (
            ci_urea [ t ] * clvl [ t ] - ci_urea [ t - 1 ] * clvl [ t - 1 ])
    qi_nh3 = fiks_nh3 - fikm_nh3 + jk_nh3 + rtau * (
            ci_nh3 [ t ] * clvl [ t ] - ci_nh3 [ t - 1 ] * clvl [ t - 1 ])
    qi_nh4 = fiks_nh4 - fikm_nh4 + jk_nh4 + rtau * (
            ci_nh4 [ t ] * clvl [ t ] - ci_nh4 [ t - 1 ] * clvl [ t - 1 ])
    qi_h = fiks_h - fikm_h + jk_h + rtau * (
            ci_h [ t ] * clvl [ t ] - ci_h [ t - 1 ] * clvl [ t - 1 ])
    qi_hco2 = fiks_hco2 - fikm_hco2 + jk_hco2 + rtau * (
            ci_hco2 [ t ] * clvl [ t ] - ci_hco2 [ t - 1 ] * clvl [ t - 1 ])
    qi_h2co2 = fiks_h2co2 - fikm_h2co2 + jk_h2co2 + rtau * (
            ci_h2co2 [ t ] * clvl [ t ] - ci_h2co2 [ t - 1 ] * clvl [ t - 1 ])
    qi_gluc = fiks_gluc - fikm_gluc + jk_gluc + rtau * (
            ci_gluc [ t ] * clvl [ t ] - ci_gluc [ t - 1 ] * clvl [ t - 1 ])

    # the proton flux must include the cellular buffers
    qi_h = qi_h - rtau * (cbuf [ t ] * clvl [ t ] - cbuf [ t - 1 ] * clvl [ t - 1 ])
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

    # now set the phis in terms of the solute generation rates
    # first the non-reactive species
    # solver = 1
    if solver == 1:
        phi [ 0 ] = phie_en
        phi [ 1 ] = phie_vlm
        phie_na = qe_na
        phi [ 2 ] = phie_na
        phie_k = qe_k
        phi [ 3 ] = phie_k
        phie_cl = qe_cl
        phi [ 4 ] = phie_cl
        phi [ 10 ] = qe_urea
        phie_gluc = qe_gluc
        phi [ 15 ] = phie_gluc
        phi [ 16 ] = phii_en
        phi [ 17 ] = phii_vlm
        phii_na = qi_na
        phi [ 18 ] = phii_na
        phii_k = qi_k
        phi [ 19 ] = phii_k
        phii_cl = qi_cl
        phi [ 20 ] = phii_cl
        phii_urea = qi_urea
        phi [ 26 ] = phii_urea
        phii_gluc = qi_gluc
        phi [ 31 ] = phii_gluc

        Buffer_Co2_LIS_Param = 1
        Buffer_Co2_Cell_Param = 1
        Buffer_HPO4_LIS_Param = 1
        Buffer_HPO4_CELL_Param = 1
        Buffer_NH3_LIS_Param = 1
        Buffer_NH3_CELL_Param = 1
        Buffer_HCO2_LIS_Param = 1
        Buffer_HCO2_CELL_Param = 1
        # LIS & CELL: see label {phosphate}
        # LIS & CELL: see labela {pH_equilibria}, used as paired equations

        phi [ 8 ], phi [ 9 ] = buff_activation(qe_hpo4, qe_h2po4, lche, pkp, ce_hpo4 [ t ], ce_h2po4 [ t ],
                                               0, Buffer_HPO4_LIS_Param)
        # print('phi [ 8 ], phi [ 9 ]', phi [ 8 ], phi [ 9 ])

        phi [ 24 ], phi [ 25 ] = buff_activation(qi_hpo4, qi_h2po4, lchi, pkp, ci_hpo4 [ t ], ci_h2po4 [ t ],
                                                 0, Buffer_HPO4_CELL_Param)

        # LIS & CELL: see label {ammonia}
        # LIS & CELL: see labela {pH_equilibria}, used as paired equations
        phi [ 11 ], phi [ 12 ] = buff_activation(qe_nh3, qe_nh4, lche, pkn, ce_nh3 [ t ], ce_nh4 [ t ],
                                                 0, Buffer_NH3_LIS_Param)

        phi [ 27 ], phi [ 28 ] = buff_activation(qi_nh3, qi_nh4, lchi, pkn, ci_nh3 [ t ], ci_nh4 [ t ],
                                                 1, Buffer_NH3_CELL_Param)

        # LIS & CELL: see label {formate}
        # LIS & CELL: see labela {pH_equilibria}, used as paired equations
        phi [ 13 ], phi [ 14 ] = buff_activation(qe_hco2, qe_h2co2, lche, pkf, ce_hco2 [ t ], ce_h2co2 [ t ],
                                                 0, Buffer_HCO2_LIS_Param)

        phi [ 29 ], phi [ 30 ] = buff_activation(qi_hco2, qi_h2co2, lchi, pkf, ci_hco2 [ t ], ci_h2co2 [ t ],
                                                 0, Buffer_HCO2_CELL_Param)

        # cell buffer content and ph equilibrium
        qi_cb = cbuf [ t ] + hcbuf [ t ] - (tbuf * clvl0 / clvl [ t ])
        phi [ 32 ] = qi_cb
        c_eQ_h2po4 = np.log10(abs(
            (ci_hpo4 [ t ] * hcbuf [ t ]) / (ci_h2po4 [ t ] * cbuf [ t ]))) if ci_h2po4 [ t ] * cbuf [ t ] == 0 or (
                (ci_hpo4 [ t ] * hcbuf [ t ]) / (ci_h2po4 [ t ] * cbuf [ t ]) <= 0) else np.log10(
            (ci_hpo4 [ t ] * hcbuf [ t ]) / (ci_h2po4 [ t ] * cbuf [ t ]))
        phi_ph_eQ = pkb - pkp - c_eQ_h2po4

        phi [ 33 ] = phi_ph_eQ
        # co2, formate, phosphate, and ammonia content:
        phi [ 5 ], phi [ 6 ], phi [ 7 ] = buff_activation_co2_formate_phosphate_ammonia(qe_h, qe_nh4, qe_hco3, qe_h2co2,
                                                                                        qe_h2po4, qe_co2, qe_h2co3, 0,
                                                                                        ce_co2 [ t ],
                                                                                        ce_h2co3 [ t ], chvl [ t ],
                                                                                        scale, 0,
                                                                                        Buffer_Co2_LIS_Param)

        # print('phi [5 ], phi [ 6 ]', phi [ 5 ], phi [ 6 ])
        phi [ 21 ], phi [ 22 ], phi [ 23 ] = buff_activation_co2_formate_phosphate_ammonia(qi_h, qi_nh4, qi_hco3,
                                                                                           qi_h2co2,
                                                                                           qi_h2po4, qi_co2, qi_h2co3,
                                                                                           qi_cb,
                                                                                           ci_co2 [ t ],
                                                                                           ci_h2co3 [ t ], clvl [ t ],
                                                                                           scale,
                                                                                           1,
                                                                                           Buffer_Co2_Cell_Param)
        # cell buffer content and ph equilibrium
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
        phi [ 34 ] = phie_cur
        return phi
    else:
        # flux
        flux = 1
        if flux:
            flux_na = fikm_na + fekm_na
            na_nak =  aie * atie_na
            flux_na_goldman = fekm_na_goldman + fikm_na_goldman


    return [fekm_na, fekm_cl, fekm_k, fekm_gluc,
            fikm_na, fikm_cl, fikm_k, fikm_gluc,
            fiks_na, fiks_cl, fiks_k, fiks_gluc,
            jk_na, jk_cl, jk_k, jk_gluc,
            feks_na, feks_cl, feks_k, feks_gluc, fekm_na_csf, fekm_na_goldman, fikm_na, fikm_na_csf,
            fikm_na_goldman, jnhe3_na, na_mi_nagluc, na_mi_nah2po4, feks_na_csf, feks_na_goldman, flux_na,
            fekm_na_csf + fikm_na_csf, flux_na_goldman,
            jnhe3_na + na_mi_nagluc + na_mi_nah2po4, gluc_mi_nagluc, na_nak, hco3_clhco3, Cellular_water_fluxes,
            clvl_imp ]
t0 = 0
tf = 100
T = 1000
dt = (tf - t0) / (T - 1)
rtau = 1 / dt
ve = matrix(-0.89432938258185330771e-02, T)
pe = matrix(-0.23093764341580683919e+02, T)
ce_na = matrix(0.14040045563695688347e+00, T)
ce_k = matrix(0.46537228854932385230e-02, T)
ce_cl = matrix(0.11200865813019536543e+00, T)
ce_hco3 = matrix(0.25607874636347928432e-01, T)
ce_h2co3 = matrix(0.43754640651881161056e-05, T)
ce_co2 = matrix(0.14987447089135489450e-02, T)
ce_hpo4 = matrix(0.29852338056241254673e-02, T)
ce_h2po4 = matrix(0.86622197071977683099e-03, T)
ce_urea = matrix(0.49046959853563188228e-02, T)
ce_nh3 = matrix(0.26914919565666300082e-05, T)
ce_nh4 = matrix(0.17484103066624158852e-03, T)
ce_hco2 = matrix(0.77584319325420768570e-03, T)
ce_h2co2 = matrix(0.20531685246992843752e-06, T)
ce_gluc = matrix(0.77236905064622532815e-02, T)
vi = matrix(-0.54945785940474827669e+02, T)
p_i = matrix(0.0e-01, T)
ci_na = matrix(0.20443615264909884011e-01, T)
ci_k = matrix(0.13742335290954643678e+00, T)
ci_cl = matrix(0.16901004032978079322e-01, T)
ci_hco3 = matrix(0.25090576334173893269e-01, T)
ci_h2co3 = matrix(0.43750529566888027647e-05, T)
ci_co2 = matrix(0.14988418688838011684e-02, T)
ci_hpo4 = matrix(0.94254445085816557920e-02, T)
ci_h2po4 = matrix(0.27910960064784369992e-02, T)
ci_urea = matrix(0.49546890743651676378e-02, T)
ci_nh3 = matrix(0.35291888877349245763e-05, T)
ci_nh4 = matrix(0.23396304517990653771e-03, T)
ci_hco2 = matrix(0.53552249401043496273e-03, T)
ci_h2co2 = matrix(0.93379252815431763174e-07, T)
ci_gluc = matrix(0.14876174248492728819e-01, T)
hcbuf = matrix(0.26959905796615793450e-01, T)
cbuf = matrix(0.40012078673998335843e-01, T)
vm = matrix(-0.18573662594736270459e+00, T)
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
# flux_hco3 = matrix(0, T)
tot_flx_im_na = matrix(0, T)
tot_flx_im_k = matrix(0, T)
tot_flx_im_cl = matrix(0, T)
tot_flx_im_hco2 = matrix(0, T)
hco3_clhco3 = matrix(0, T)
cl_clhco3 = matrix(0, T)
cl_kcl = matrix(0, T)
cl_clhco2 = matrix(0, T)
cl_na_clhco3 = matrix(0, T)
na_nak = matrix(0, T)
k_nak = matrix(0, T)
tot_flx_im_h2po4 = matrix(0, T)
tot_flx_im_nh4 = matrix(0, T)
tot_flx_im_h = matrix(0, T)
na_nah = matrix(0, T)
tot_flx_im_gluc = matrix(0, T)
tot_flx_is_na = matrix(0, T)
tot_flx_is_k = matrix(0, T)
tot_flx_is_cl = matrix(0, T)
tot_flx_is_hco3 = matrix(0, T)
cjk_na = matrix(0, T)
cjk_k = matrix(0, T)
cjk_cl = matrix(0, T)
cjk_hco3 = matrix(0, T)
fikm_na_csf = matrix(0, T)
fikm_na_goldman =matrix(0, T)
fekm_na_goldman =matrix(0, T)
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

jk_na = matrix(0, T)
jk_k = matrix(0, T)
jk_cl = matrix(0, T)
jk_gluc = matrix(0, T)

feks_na = matrix(0, T)
feks_k = matrix(0, T)
feks_cl = matrix(0, T)
feks_gluc = matrix(0, T)


na_mi_nagluc = matrix(0, T)
gluc_mi_nagluc = matrix(0, T)
na_mi_nah2po4 = matrix(0, T)
jnhe3_na = matrix(0, T)
cl_mi_clhco2 = matrix(0, T)
cl_mi_clhco3 = matrix(0, T)
vars = [ ve, pe, ce_na, ce_k,
         ce_cl, ce_hco3, ce_h2co3, ce_co2, ce_hpo4, ce_h2po4, ce_urea, ce_nh3, ce_nh4, ce_hco2, ce_h2co2, ce_gluc,
         vi, imp, ci_na, ci_k, ci_cl, ci_hco3, ci_h2co3, ci_co2, ci_hpo4, ci_h2po4, ci_urea, ci_nh3, ci_nh4, ci_hco2,
         ci_h2co2, ci_gluc, cbuf, hcbuf, vm]

Figure_4a = 0
Figure_4b = 0
Figure_4c = 0
Figure_5 = 0
Figure_9_10 = 0
Figure_6_7_8 = 1
Figure_7 = 0
Figure_6A = 1
Figure_6B = 0
Figure_8A = 0
Figure_8B = 0
flx = [fekm_na, fekm_cl, fekm_k, fekm_gluc, fikm_na, fikm_cl, fikm_k, fikm_gluc,
      fiks_na, fiks_cl, fiks_k, fiks_gluc, jk_na, jk_cl, jk_k, jk_gluc,
      feks_na, feks_cl, feks_k, feks_gluc, fekm_na_csf, fekm_na_goldman, fikm_na, fikm_na_csf,
      fikm_na_goldman,jnhe3_na, na_mi_nagluc, na_mi_nah2po4, feks_na_csf, feks_na_goldman, flux_na,
      fekm_na_csf + fikm_na_csf, flux_na_goldman,
      jnhe3_na + na_mi_nagluc + na_mi_nah2po4, gluc_mi_nagluc, na_nak, hco3_clhco3, Cellular_water_fluxes,
      clvl_imp]
if Figure_9_10:
    t = 0
    guess = [ var [ t ] for var in vars ]
    from scipy import optimize
    mylist = [ 0.0 for i in range(T) ]
    cmna = [ 0.1400, 0.1410, 0.140, 0.14, 0.140 ]
    n_nhe3 = [0.27500e-08, 0.2000e-08, 0.1500e-08, 0.1000e-08,  0.01000e-08]
    while True:
        t += 1
        if t == T:
            break
        elif 0 < t <= T / len(cmna):
            nnhe3 = n_nhe3 [ 0 ]
            mylist [ t ] = n_nhe3 [ 0 ]
            result = scipy.optimize.root(eQs, np.array(guess), 1)
            guess = [ var [ t ] for var in vars ]
            result_flx = eQs(guess, 0)
            for i in range(len(result_flx)):
                flx [ i ] [ t ] = result_flx [ i ]
        elif 1 * T / len(cmna) < t <= 2 * T / len(cmna):
            nnhe3 = n_nhe3 [ 1 ]
            mylist [ t ] = n_nhe3 [ 1 ]
            result = scipy.optimize.root(eQs, np.array(guess), 1)
            guess = [ var [ t ] for var in vars ]
            result_flx = eQs(guess, 0)
            for i in range(len(result_flx)):
                flx [ i ] [ t ] = result_flx [ i ]
        elif 2 * T / len(cmna) < t <= 3 * T / len(cmna):
            nnhe3 = n_nhe3 [ 2 ]
            mylist [ t ] = n_nhe3 [ 2 ]
            result = scipy.optimize.root(eQs, np.array(guess), 1)
            guess = [ var [ t ] for var in vars ]
            result_flx = eQs(guess, 0)
            for i in range(len(result_flx)):
                flx [ i ] [ t ] = result_flx [ i ]
        elif 3 * T / len(cmna) < t <= 4 * T / len(cmna):
            nnhe3 = n_nhe3 [ 3 ]
            mylist [ t ] = n_nhe3 [ 3 ]
            result = scipy.optimize.root(eQs, np.array(guess), 1)
            guess = [ var [ t ] for var in vars ]
            result_flx = eQs(guess, 0)
            for i in range(len(result_flx)):
                flx [ i ] [ t ] = result_flx [ i ]
        elif 4 * T / len(cmna) < t <= 5 * T / len(cmna):
            nnhe3 = n_nhe3 [ 4 ]
            mylist [ t ] = n_nhe3 [ 4 ]
            result = scipy.optimize.root(eQs, np.array(guess), 1)
            guess = [ var [ t ] for var in vars ]
            result_flx = eQs(guess, 0)
            for i in range(len(result_flx)):
                flx [ i ] [ t ] = result_flx [ i ]

    t = np.linspace(t0, tf, T)
    t_t = np.transpose(t)
    mylist = np.asarray(mylist)
    my_list_Figure_9_10 = {'t': t,'L_NHE3': mylist, 'flux_na_nah': jnhe3_na, 'flux_hco3_clcho3': hco3_clhco3, 'flux_Na_nak': na_nak, 'ci_na': ci_na , 'ci_cl': ci_cl, 'ci_hco2': ci_hco2, 'ci_hco3': ci_hco3 ,'clvl_imp': clvl_imp}
    pickled_list = pickle.dumps(my_list_Figure_9_10)
    of = open('Data_Figure_9_10.py', 'wb')
    of.write(pickled_list)
    of.close()
    fig8 = plt.figure(constrained_layout=True)
    spec2 = gridspec.GridSpec(ncols=1, nrows=4, figure=fig8)
    ax1 = fig8.add_subplot(spec2 [ 0, 0 ])
    ax1.plot(t [ 0:-2 ], mylist [ 0:-2 ], 'r-')
    ax1.set_ylabel('L_NHE3', color='blue')
    ax4 = fig8.add_subplot(spec2 [ 1, 0 ])
    ax4.plot(t, jnhe3_na, 'b-')
    ax4.set_ylabel('flux_na_nah', color='blue')
    ax2 = fig8.add_subplot(spec2 [ 2, 0 ])
    ax2.plot(t, hco3_clhco3, 'b-')
    ax2.set_ylabel(' flux_hco3_clcho3 ', color='blue')
    ax3 = fig8.add_subplot(spec2 [ 3, 0 ])
    ax3.plot(t, na_nak, 'b-')
    ax3.set_ylabel('flux_Na_nak ', color='blue')
    fig3 = plt.figure(constrained_layout=True)
    spec2 = gridspec.GridSpec(ncols=1, nrows =  6, figure=fig3)
    ax1 = fig3.add_subplot(spec2 [ 0, 0 ])
    ax1.plot(t [ 0:-2 ], mylist [ 0:-2 ], 'r-')
    ax1.set_ylabel('L_NHE3', color='blue')
    ax2 = fig3.add_subplot(spec2 [ 1, 0 ])
    ax2.plot(t, ci_na, 'b-')
    ax2.set_ylabel(' ci_na ', color='blue')
    ax3 = fig3.add_subplot(spec2 [ 2, 0 ])
    ax3.plot(t, ci_cl, 'b-')
    ax3.set_ylabel('ci_cl ', color='blue')
    ax4 = fig3.add_subplot(spec2 [ 3, 0 ])
    ax4.plot(t, ci_hco3, 'b-')
    ax4.set_ylabel('ci_hco3 ', color='blue')
    ax5 = fig3.add_subplot(spec2 [ 4, 0 ])
    ax5.plot(t, ci_hco2, 'b-')
    ax5.set_ylabel('ci_hco2 ', color='blue')
    ax6 = fig3.add_subplot(spec2 [ 5, 0 ])
    ax6.plot(t, clvl_imp, 'b-')
    ax6.set_ylabel('clvl ', color='blue')
    ax6.set_xlabel('Time (s) ', color='blue')
    plt.show
if Figure_5:
    t = 0
    guess = [ var [ t ] for var in vars ]
    from scipy import optimize
    NaCl = 1
    scale_per = 0.12
    cmna = [ 0.14 + NaCl * j * 0.14 * scale_per for j in range(4)]
    cmna.append(0.14)
    cmcl = [ 0.1132 + NaCl * j * 0.1132 * scale_per for j in range(4)]
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
            cs_na = cmna [ 0]
            cs_cl = cmcl [ 0 ]
            mylist_na[t] = cmna[0]
            mylist_cl [ t ] = cmcl [ 0 ]
            result = scipy.optimize.root(eQs, np.array(guess), 1)
            guess =[var[t] for var in vars]
            result_flx=eQs(guess, 0)
            for i in range(len(result_flx)):
                flx[i][t]=result_flx[i]
        elif 1 * T / len(cmna) < t <= 2 * T / len(cmna):
            mylist_na[t]=cmna[1]
            mylist_cl [ t ]=cmcl [ 1]
            cm_na=cmna [1 ]
            cm_cl=cmcl [ 1 ]
            cs_na=cmna [ 1 ]
            cs_cl=cmcl [ 1]
            result=scipy.optimize.root(eQs, np.array(guess), 1)
            guess =[var[t] for var in vars]
            result_flx=eQs(guess, 0)
            for i in range(len(result_flx)):
                flx[i][t]=result_flx[i]
        elif 2 * T / len(cmna) < t <= 3 * T / len(cmna):
            cm_na=cmna [ 2 ]
            cm_cl=cmcl [ 2 ]
            cs_na=cmna [ 2 ]
            cs_cl=cmcl [ 2 ]
            mylist_na[t]=cmna[2]
            mylist_cl [ t ]=cmcl [2]
            result=scipy.optimize.root(eQs, np.array(guess), 1)
            guess =[var[t] for var in vars]
            result_flx=eQs(guess, 0)
            for i in range(len(result_flx)):
                flx[i][t]=result_flx[i]
        elif 3 * T / len(cmna) < t <= 4 * T / len(cmna):
            cm_na=cmna [3]
            cm_cl=cmcl [3]
            cs_na=cmna [3]
            cs_cl=cmcl [3]
            mylist_na[t]=cmna[3]
            mylist_cl[ t ]=cmcl [ 3 ]
            result=scipy.optimize.root(eQs, np.array(guess), 1)
            guess =[var[t] for var in vars]
            result_flx=eQs(guess, 0)
            for i in range(len(result_flx)):
                flx[i][t]=result_flx[i]
        elif 4 * T / len(cmna) < t <= 5 * T / len(cmna):
            cm_na=cmna [ 4 ]
            cm_cl=cmcl [4 ]
            cs_na=cmna [4 ]
            cs_cl=cmcl [ 4 ]
            mylist_na[t]=cmna[4]
            mylist_cl[ t ]=cmcl [4]
            result=scipy.optimize.root(eQs, np.array(guess), 1)
            guess =[var[t] for var in vars]
            result_flx=eQs(guess, 0)
            for i in range(len(result_flx)):
                flx[i][t]=result_flx[i]
        else:
            print('t is not in a defined range.')
    t=np.linspace(t0, tf, T)
    my_list_Figure_5 = {'t': t, 'mylist_na': mylist_na, 'mylist_cl': mylist_cl, 'vm': vm,
                            'ci_na': ci_na, 'ci_cl': ci_cl, 'ci_k': ci_k }
    pickled_list = pickle.dumps(my_list_Figure_5)
    of = open('Data_Figure_5.py', 'wb')
    of.write(pickled_list)
    of.close()
    fig3=plt.figure(constrained_layout=True)
    spec2=gridspec.GridSpec(ncols=1, nrows=5, figure=fig3)
    ax1=fig3.add_subplot(spec2 [ 0, 0 ])
    ax1.plot(t [ 0:-2 ], mylist_na [ 0:-2 ], color='#007282', label="Na ")
    ax1.plot(t [ 0:-2 ], mylist_cl [ 0:-2 ], color='#b66dff', label="Cl ")
    ax1.set_ylabel('Cm, Cs', color='black')
    ax1.set_xlim(-5, tf)
    ax1.legend(frameon=False)
    ax0 = fig3.add_subplot(spec2 [1, 0 ])
    ax0.set_ylabel('Vm', color='black')
    ax0.plot(t [ 0:-2 ], vm[ 0:-2 ], 'black')
    ax0.set_xlim(-5, tf)
    ax0.set_ylim(0.0, 0.15)
    ax2 = fig3.add_subplot(spec2 [2, 0 ])
    ax2.set_ylabel('ci_na', color='black')
    ax2.plot(t [ 0:-2 ], ci_na[ 0:-2 ], 'black')
    ax2.set_xlim(-5, tf)
    ax3 = fig3.add_subplot(spec2 [3, 0 ])
    ax3.plot(t [ 0:-2 ], ci_k[ 0:-2 ], 'black')
    ax3.set_ylabel('ci_k', color='black')
    ax3.set_xlim(-5, tf)
    ax4 = fig3.add_subplot(spec2 [4, 0 ])
    ax4.plot(t [ 0:-2 ], ci_cl[ 0:-2 ], 'black')
    ax4.set_ylabel('ci_cl', color='black')
    ax4.set_xlim(-5, tf)
    plt.show()
if Figure_4a:
    mylist_na = [ 0.1 for i in range(T) ]
    mylist_cl = [ 0.1 for i in range(T) ]
    mylist_hco3 = [ 0.1 for i in range(T) ]
    cm_hco3_0 = 0.004
    cm_hco3_f = 0.05
    slope = (cm_hco3_f - cm_hco3_0) / (tf - t0)
    mylist_hco3 = [ cm_hco3_0 + 0.1 * slope * i for i in range(T) ]
    t = 0
    guess = [ var [ t ] for var in vars ]

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
            cm_hco3 = mylist_hco3 [ t ]
            result = scipy.optimize.root(eQs, np.array(guess), 1)
            guess = [ var [ t ] for var in vars ]
            result_flx = eQs(guess, 0)
            for i in range(len(result_flx)):
                flx [ i ] [ t ] = result_flx [ i ]

    my_list = {'fekm_cl': fekm_cl, 'fikm_cl': fikm_cl, 'flux_cl': flux_cl, 'mylist_hco3': mylist_hco3}
    pickled_list = pickle.dumps(my_list)
    of = open('Data_Figure_4a.py', 'wb')
    of.write(pickled_list)
    of.close()
    fig4 = plt.figure(constrained_layout=True)
    spec2 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig4)
    ax3 = fig4.add_subplot(spec2 [ 0, 0 ])
    ax3.set_ylabel('flux_cl [mmol/sec]', color='black')
    x1 = ax3.plot(mylist_hco3 [ 30:-2 ], flux_cl [ 30:-2 ], 'black')
    x2 = ax3.plot(mylist_hco3 [ 30:-2 ], fikm_cl [ 30:-2 ], 'g')
    x3 = ax3.plot(mylist_hco3 [ 30:-2 ], fekm_cl [ 30:-2 ], 'blue')
    ax3.legend([ x1 [ 0 ], x2 [ 0 ], x3 [ 0 ] ], [ 'Total Epithelial', 'Cellular', 'Junctional' ])
if Figure_4b:
    mylist_hco2 = [0.1 for i in range(T)]
    cm_hco2_0 = 0.004
    cm_hco2_f = 0.05
    slope = (cm_hco2_f - cm_hco2_0) / (tf - t0)
    mylist_hco2 = [ cm_hco2_0 + 0.1 * slope * i for i in range(T)]
    t = 0
    guess = [ var [ t ] for var in vars ]

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
            cm_hco2 = mylist_hco2 [ t ]
            cs_hco2 = mylist_hco2 [ t ]
            cm_hco3 = 0.004
            result = scipy.optimize.root(eQs, np.array(guess), 1)
            guess = [ var [ t ] for var in vars ]
            result_flx = eQs(guess, 0)
            for i in range(len(result_flx)):
                flx [ i ] [ t ] = result_flx [ i ]
    fig4b = plt.figure(constrained_layout=True)
    spec2 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig4b)
    ax1 = fig4b.add_subplot(spec2 [ 0, 0 ])
    ax1.set_ylabel('flux_cl [mmol/sec]', color='black')
    x1 = ax1.plot(mylist_hco2 [ 30:-2 ], flux_cl [ 30:-2 ], 'black')
    x2 = ax1.plot(mylist_hco2 [ 30:-2 ], fikm_cl [ 30:-2 ], 'g')
    x3 = ax1.plot(mylist_hco2 [ 30:-2 ], fekm_cl [ 30:-2 ], 'blue')
    ax1.legend([ x1 [ 0 ], x2 [ 0 ], x3 [ 0 ] ], [ 'Total Epithelial', 'Cellular', 'Junctional' ])
    plt.show()
    my_list = {'fekm_cl': fekm_cl, 'fikm_cl': fikm_cl, 'flux_cl': flux_cl, 'mylist_hco2': mylist_hco2}
    pickled_list = pickle.dumps(my_list)
    of = open('Data_Figure_4b.py', 'wb')
    of.write(pickled_list)
    of.close()
if Figure_4c:
    t = 0
    guess = [ var [ t ] for var in vars ]
    from scipy import optimize
    cm_hco3 = 0.004
    cm_hco2_0 = 0.004
    cm_hco2_f = 0.05
    slope = (cm_hco2_f - cm_hco2_0) / (tf - t0)
    mylist_hco2 = [ cm_hco2_0 + 0.1 * slope * i for i in range(T) ]
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
            cm_hco2 = mylist_hco2 [ t ]
            cs_hco2 = mylist_hco2 [ t ]
            result = scipy.optimize.root(eQs, np.array(guess), 1)
            guess = [ var [ t ] for var in vars ]
            result_flx = eQs(guess, 0)
            for i in range(len(result_flx)):
                flx [ i ] [ t ] = result_flx [ i ]
    my_list = {'fekm_cl': fekm_cl, 'fikm_cl': fikm_cl, 'flux_cl': flux_cl, 'mylist_hco2': mylist_hco2}
    pickled_list = pickle.dumps(my_list)
    of = open('Data_Figure_4c.py', 'wb')
    of.write(pickled_list)
    of.close()
    fig4 = plt.figure(constrained_layout=True)
    spec2 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig4)
    ax1 = fig4.add_subplot(spec2 [0, 0])
    ax1.set_ylabel('flux_cl [mmol/sec]', color='black')
    x1 = ax1.plot(mylist_hco2 [ 30:-2 ], flux_cl [ 30:-2 ], 'black')
    x2 = ax1.plot(mylist_hco2 [ 30:-2 ], fikm_cl [ 30:-2 ], 'g')
    x3 = ax1.plot(mylist_hco2 [ 30:-2 ], fekm_cl [ 30:-2 ], 'blue')
    ax1.legend([ x1 [ 0 ], x2 [ 0 ], x3 [ 0 ] ], [ 'Total Epithelial', 'Cellular', 'Junctional' ])
if Figure_6_7_8:
    t = 0
    guess = [ var [ t ] for var in vars ]
    while True:
        t += 1
        if t == T:
            break
        else:
            result = scipy.optimize.root(eQs, np.array(guess), 1)
            guess = [ var [ t ] for var in vars ]
            result_flx = eQs(guess, 0)
            for i in range(len(result_flx)):
                flx [ i ] [ t ] = result_flx [ i ]
    my_list_default = {'fekm': [ fekm_na [ -1 ], fekm_k [ -1 ], fekm_cl [ -1 ],  fekm_gluc [ -1 ] ],
                       'fikm': [ fikm_na [ -1 ], fikm_k [ -1 ], fikm_cl [ -1 ],  fikm_gluc [ -1 ] ],
                       'fiks': [ fiks_na [ -1 ],  fiks_k [ -1 ], fiks_cl [ -1 ], fiks_gluc [ -1 ] ],
                       'jk': [ jk_na [ -1 ], jk_k [ -1 ], jk_cl [ -1 ], jk_gluc [ -1 ] ],
                       'feks': [ feks_na [ -1 ], feks_k [ -1 ], feks_cl [ -1 ], feks_gluc [ -1 ] ],
                       }
    slt_con_default = {'Na': ci_na [ -1 ], 'K': ci_k [ -1 ], 'Cl': ci_cl [ -1 ],  'Gluc': ci_gluc [ -1 ]}
    flx_na_mem_default = [ feks_na [ -1 ], jk_na[-1], fiks_na [ -1 ], fekm_na [ -1 ], fikm_na [ -1 ],
                       flux_na [ -1 ] ]

    epithelial_flx_variation_default = [ fekm_na_csf [ -1 ] + fikm_na_csf [ -1 ],
                                     fekm_na_goldman [ -1 ] + fikm_na_goldman [ -1 ],
                                     jnhe3_na [ -1 ] + na_mi_nagluc [ -1 ] + na_mi_nah2po4 [ -1 ] ]
    epithelial_flx_elct_chmcl_default = [ jnhe3_na [ -1 ], na_mi_nagluc [ -1 ], na_mi_nah2po4 [ -1 ] ]

    t = 0
    guess = [ var [ t ] for var in vars ]
    while True:
        t += 1
        if t == T:
            break
        else:
            param_sodium_pumps = 0
            result = scipy.optimize.root(eQs, np.array(guess), 1)
            guess = [ var [ t ] for var in vars ]
            result_flx = eQs(guess, 0)
            for i in range(len(result_flx)):
                flx [ i ] [ t ] = result_flx [ i ]
    my_list_nak = {'fekm':[fekm_na[-1], fekm_k[-1], fekm_cl[-1],  fekm_gluc[-1]],
            'fikm':[fikm_na[-1],  fikm_k[-1], fikm_cl[-1], fikm_gluc[-1]],
            'fiks':[fiks_na[-1], fiks_k[-1], fiks_cl[-1],  fiks_gluc[-1]],
            'jk': [jk_na[-1], jk_k[-1], jk_cl[-1], jk_gluc[-1]],
            'feks':[feks_na[-1], feks_k[-1], feks_cl[-1], feks_gluc[-1] ]}
    slt_con_nak = {'Na': ci_na [ -1 ], 'K': ci_k [ -1 ], 'Cl': ci_cl [ -1 ], 'Gluc': ci_gluc [ -1 ]}
    flx_na_mem_nak = [feks_na[-1], jk_na[-1], fiks_na[-1], fekm_na[-1], fikm_na[-1], flux_na[-1]]
    epithelial_flx_variation_nak = [ fekm_na_csf[-1] + fikm_na_csf[-1], fekm_na_goldman[-1] + fikm_na_goldman[-1], jnhe3_na[-1] + na_mi_nagluc [-1] + na_mi_nah2po4[-1]]
    epithelial_flx_elct_chmcl_nak = [ jnhe3_na[-1],  na_mi_nagluc [-1],  na_mi_nah2po4[-1]]

    t = 0
    guess = [ var [ t ] for var in vars ]
    while True:
        t += 1
        if t == T:
            break
        else:
            param_sodium_pumps = 1
            param_kcl = 0
            result = scipy.optimize.root(eQs, np.array(guess), 1)
            guess = [ var [ t ] for var in vars ]
            result_flx = eQs(guess, 0)
            for i in range(len(result_flx)):
                flx [ i ] [ t ] = result_flx [ i ]
    my_list_kcl = {'fekm':[fekm_na[-1], fekm_k[-1], fekm_cl[-1], fekm_gluc[-1]],
            'fikm':[fikm_na[-1], fikm_k[-1], fikm_cl[-1],  fikm_gluc[-1]],
            'fiks':[fiks_na[-1], fiks_k[-1], fiks_cl[-1], fiks_gluc[-1]],
            'jk':[jk_na[-1], jk_k[-1], jk_cl[-1],  jk_gluc[-1]],
            'feks':[feks_na[-1], feks_k[-1], feks_cl[-1], feks_gluc[-1] ]}
    slt_con_kcl = {'Na': ci_na [ -1 ], 'K': ci_k [ -1 ], 'Cl': ci_cl [ -1 ], 'Gluc': ci_gluc [ -1 ]}
    t = 0
    guess = [ var [ t ] for var in vars ]
    while True:
        t += 1
        if t == T:
            break
        else:
            param_sodium_pumps = 1
            param_kcl = 1
            param_nahco3 = 1.e-12
            result = scipy.optimize.root(eQs, np.array(guess), 1)
            guess = [ var [ t ] for var in vars ]
            result_flx = eQs(guess, 0)
            for i in range(len(result_flx)):
                flx [ i ] [ t ] = result_flx [ i ]
    my_list_nahco3 = {'fekm':[fekm_na[-1],  fekm_k[-1], fekm_cl[-1], fekm_gluc[-1]],
            'fikm':[fikm_na[-1],  fikm_k[-1], fikm_cl[-1], fikm_gluc[-1]],
            'fiks':[fiks_na[-1], fiks_k[-1], fiks_cl[-1], fiks_gluc[-1]],
            'jk':[jk_na[-1], jk_k[-1], jk_cl[-1], jk_gluc[-1]],
            'feks':[feks_na[-1], feks_k[-1], feks_cl[-1], feks_gluc[-1] ]}
    slt_con_nahco3 = {'Na': ci_na [ -1 ],  'K': ci_k [ -1 ], 'Cl':ci_cl [ -1 ], 'Gluc': ci_gluc [ -1 ]}
    from collections import defaultdict
    my_list_Figure6A = defaultdict(list)
    for d in (my_list_default, my_list_nak, my_list_kcl, my_list_nahco3):  # you can list as many input dicts as you want here
        for key, value in d.items():
            my_list_Figure6A [ key ].append(value)
    my_list_Figure6B = defaultdict(list)
    for d in (slt_con_default, slt_con_nak, slt_con_kcl, slt_con_nahco3):  # you can list as many input dicts as you want here
        for key, value in d.items():
            my_list_Figure6B [ key ].append(value)
    pickled_list = pickle.dumps(my_list_Figure6A)
    of = open('Data_Figure_6A.py', 'wb')
    of.write(pickled_list)
    of.close()
    pickled_list = pickle.dumps(my_list_Figure6B)
    of = open('Data_Figure_6B.py', 'wb')
    of.write(pickled_list)
    of.close()
    t = 0
    guess = [ var [ t ] for var in vars ]
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
            guess = [ var [ t ] for var in vars ]
            result_flx = eQs(guess, 0)
            for i in range(len(result_flx)):
                flx [ i ] [ t ] = result_flx [ i ]
    my_list_nhe3 = {'fekm':[fekm_na[-1],  fekm_k[-1], fekm_cl[-1], fekm_gluc[-1]],
            'fikm':[fikm_na[-1], fikm_k[-1], fikm_cl[-1],  fikm_gluc[-1]],
            'fiks':[fiks_na[-1], fiks_k[-1], fiks_cl[-1],  fiks_gluc[-1]],
            'jk':[jk_na[-1], jk_k[-1], jk_cl[-1], jk_gluc[-1]],
            'feks':[feks_na[-1], feks_k[-1], feks_cl[-1], feks_gluc[-1] ]}
    slt_con_nhe3 = {'Na': ci_na [ -1 ], 'K': ci_k [ -1 ], 'Cl': ci_cl [ -1 ], 'Gluc':ci_gluc [ -1 ]}
    flx_na_mem_nhe3 = [feks_na[-1], jk_na[-1], fiks_na[-1], fekm_na[-1], fikm_na[-1], fekm_na[-1] + fikm_na[-1]]
    epithelial_flx_variation_nhe3 = [ fekm_na_csf[-1] + fikm_na_csf[-1], fekm_na_goldman[-1] + fikm_na_goldman[-1], jnhe3_na[-1] + na_mi_nagluc [-1]+ na_mi_nah2po4[-1]]
    epithelial_flx_elct_chmcl_nhe3 = [ jnhe3_na[-1],  na_mi_nagluc [-1],  na_mi_nah2po4[-1]]
    t = 0
    guess = [var [ t ] for var in vars]
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
            guess = [ var [ t ] for var in vars ]
            result_flx = eQs(guess, 0)
            for i in range(len(result_flx)):
                flx [ i ] [ t ] = result_flx [ i ]
    my_list_nagluc = {'fekm':[fekm_na[-1], fekm_k[-1], fekm_cl[-1],  fekm_gluc[-1]],
            'fikm':[fikm_na[-1], fikm_k[-1], fikm_cl[-1],  fikm_gluc[-1]],
            'fiks':[fiks_na[-1], fiks_k[-1],  fiks_cl[-1], fiks_gluc[-1]],
            'jk':[jk_na[-1], jk_k[-1], jk_cl[-1],  jk_gluc[-1]],
            'feks':[feks_na[-1],  feks_k[-1], feks_cl[-1], feks_gluc[-1]]}
    slt_con_nagluc = {'Na': ci_na [ -1 ], 'K': ci_k [ -1 ], 'Cl': ci_cl [ -1 ],  'Gluc': ci_gluc [ -1 ]}
    t = 0
    guess = [ var [ t ] for var in vars ]
    while True:
        t += 1
        if t == T:
            break
        else:
            param_h2po4 = 0
            result = scipy.optimize.root(eQs, np.array(guess), 1)
            guess = [ var [ t ] for var in vars ]
            result_flx = eQs(guess, 0)
            for i in range(len(result_flx)):
                flx [ i ] [ t ] = result_flx [ i ]

    my_list_nah2po4 =  {'fekm':[fekm_na[-1],fekm_k[-1], fekm_cl[-1],  fekm_gluc[-1]],
            'fikm':[fikm_na[-1],  fikm_k[-1],fikm_cl[-1], fikm_gluc[-1]],
            'fiks':[fiks_na[-1], fiks_k[-1],fiks_cl[-1],  fiks_gluc[-1]],
            'jk':[jk_na[-1], jk_k[-1], jk_cl[-1],  jk_gluc[-1]],
            'feks':[feks_na[-1], feks_k[-1], feks_cl[-1], feks_gluc[-1] ]}
    slt_con_nah2po4 = {'Na': ci_na [ -1 ],  'K': ci_k [ -1 ], 'Cl': ci_cl [ -1 ],'Gluc': ci_gluc [ -1 ]}

    from collections import defaultdict
    my_list_Figure7 = {'flx_na_mem_default': flx_na_mem_default , 'epithelial_flx_variation_default':epithelial_flx_variation_default, 'epithelial_flx_elct_chmcl_default': epithelial_flx_elct_chmcl_default,
                       'flx_na_mem_nak': flx_na_mem_nak , 'epithelial_flx_variation_nak':epithelial_flx_variation_nak, 'epithelial_flx_elct_chmcl_nak': epithelial_flx_elct_chmcl_nak,
                       'flx_na_mem_nhe3': flx_na_mem_nhe3 , 'epithelial_flx_variation_nhe3':epithelial_flx_variation_nhe3, 'epithelial_flx_elct_chmcl_nhe3': epithelial_flx_elct_chmcl_nhe3}
    my_list_Figure8A = defaultdict(list)
    my_list_Figure8B = defaultdict(list)
    for d in (my_list_default, my_list_nhe3, my_list_nagluc, my_list_nah2po4):
        for key, value in d.items():
            my_list_Figure8A [ key ].append(value)
    for d in (slt_con_default, slt_con_nhe3, slt_con_nagluc, slt_con_nah2po4):
        for key, value in d.items():
            my_list_Figure8B [ key ].append(value)

    pickled_list = pickle.dumps(my_list_Figure7)
    of = open('Data_Figure_7.py', 'wb')
    of.write(pickled_list)
    of.close()
    pickled_list = pickle.dumps(my_list_Figure8A)
    of = open('Data_Figure_8A.py', 'wb')
    of.write(pickled_list)
    of.close()
    pickled_list = pickle.dumps( my_list_Figure8B)
    of = open('Data_Figure_8B.py', 'wb')
    of.write(pickled_list)
    of.close()
    if Figure_6A:
        of = open('Data_Figure_6A.py', 'rb')
        read_file = of.read()
        my_loaded_list = pickle.loads(read_file)
        of.close()
        scale_factor = 1
        flux_me = np.array( my_loaded_list [ 'fekm' ])
        flux_mi = np.array( my_loaded_list [ 'fikm' ])
        flux_is = np.array( my_loaded_list [ 'fiks' ])
        flux_ie = np.array(my_loaded_list [ 'jk' ])
        flux_es = np.array(my_loaded_list [ 'feks' ])
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
                         range=[ 'black', '#D55E00', '#F0E442', '#007282', '#2B9F78' ] ),
                      )) \
            .configure_view().configure_axis(grid=False )
        chart.show()
        chart.save('Figure_6A.svg')

        of = open('Data_Figure_6B.py', 'rb')
        read_file = of.read()
        my_loaded_list1 = pickle.loads(read_file)
        of.close()
        w =0.1
        Na = np.array(my_loaded_list1['Na'])
        K = np.array(my_loaded_list1 [ 'K'])
        Cl = np.array(my_loaded_list1 ['Cl'])
        Gluc = np.array(my_loaded_list1 ['Gluc'])
        x = ["Na  K  Cl Gluc", "Na  K  Cl Gluc", "Na  K  Cl Gluc","Na  K  Cl Gluc" ]

        bar1 = np.arange(len(x))
        bar2 = [ i + 2.0 * w for i in bar1 ]
        bar3 = [ i + 2.0 * w for i in bar2 ]
        bar4 = [ i + 2.0 * w for i in bar3 ]
        bar5 = [ i + w for i in bar2 ]
        plt.figure(1)
        plt.bar(bar1, Na, 0.2, color="#2B9F78", label="Original Setup")
        plt.bar(bar2, K, 0.2, color="#007282", label="NaK = 0")
        plt.bar(bar3, Cl, 0.2, color="black", label="KCl = 0")
        plt.bar(bar4, Gluc, 0.2, color="#D55E00", label="NaHCO3 = 0")
        plt.ylabel('Cellular Concentration', fontsize=12)
        plt.xlabel("Solutes", fontsize=12)
        plt.xticks(bar5, x, fontsize=11)
        plt.legend(frameon=False)
        plt.savefig(pltfolder + '/Figure_6B' + '.png')
    if Figure_8A:
        of = open('Data_Figure_8A.py', 'rb')
        read_file = f.read()
        my_loaded_list = pickle.loads(read_file)

        of.close()
        print('my_loaded_list', my_loaded_list)

        scale_factor = 1

        flux_me = np.array(
            my_loaded_list [ 'fekm' ])
        print(flux_me)
        flux_mi = np.array(
            my_loaded_list [ 'fikm' ])
        flux_is = np.array(
            my_loaded_list [ 'fiks' ])

        flux_ie = np.array(my_loaded_list [ 'jk' ])
        print(flux_ie)

        flux_es = np.array(my_loaded_list [ 'feks' ]
                       )
        df1 = pd.DataFrame( flux_me,
                           index=[ "Original Setup", "NHE3 = 0", "SGLT = 0", "NaH2PO4 = 0" ],
                           columns=[ "Na", "K", "Cl", "Gluc" ])
        df2 = pd.DataFrame( flux_is,
                           index=[ "Original Setup", "NHE3 = 0", "SGLT = 0", "NaH2PO4 = 0" ],
                           columns=[ "Na", "K", "Cl", "Gluc" ])
        df3 = pd.DataFrame( flux_ie,
                           index=[ "Original Setup", "NHE3 = 0", "SGLT = 0", "NaH2PO4 = 0" ],
                           columns=[ "Na", "K", "Cl", "Gluc" ])
        df4 = pd.DataFrame( flux_es,
                           index=[ "Original Setup", "NHE3 = 0", "SGLT = 0", "NaH2PO4 = 0" ],
                           columns=[ "Na", "K", "Cl", "Gluc" ])
        df5 = pd.DataFrame( flux_mi,
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

            # tell Altair which field to group columns on
            alt.X('c2', title=None, sort=[ "Na", "K", "Cl", "Gluc" ]),

            # tell Altair which field to use as Y values and how to calculate
            alt.Y('sum(values)',
                  axis=alt.Axis(
                      grid=False,
                      title=None)),

            # tell Altair which field to use to use as the set of columns to be  represented in each group
            alt.Column('c1', title=None, sort=[ "Original Setup", "NHE3 = 0", "SGLT = 0", "NaH2PO4 = 0" ]),

            # tell Altair which field to use for color segmentation
            alt.Color('Fluxes',
                      scale=alt.Scale(
                          # make it look pretty with an enjoyable color pallet
                          # COLOR VISION DEFICIENT :CVD 'black', orange:'#E69F00', blue:'#007282', vermilon:'#D55E00',
                          # BLUISHGREEN:'#2B9F78' yellow:'#F0E442'
                          range=[ 'black', '#D55E00', '#F0E442', '#007282', '#2B9F78' ]
                          #
                          # range=['bluishgreen', 'blue','yellow', 'black','vermillion']
                      ),
                      )) \
            .configure_view(
            # remove grid lines around column clusters
            # strokeOpacity=0
        ).configure_axis(
            grid=False
        )
        chart.show()
        # chart.save('chart.svg')
        chart.save('Figure_8A.png')

    of = open('Data_Figure_8B.py', 'rb')
    read_file = f.read()
    my_loaded_list1 = pickle.loads(read_file)
    of.close()
    print('my_loaded_list1', my_loaded_list1)
    Na = np.array(my_loaded_list1 [ 'Na' ])
    K = np.array(my_loaded_list1 [ 'K' ])
    Cl = np.array(my_loaded_list1 [ 'Cl' ])
    Gluc = np.array(my_loaded_list1 [ 'Gluc' ])
    w = 0.1
    Na = np.array(my_loaded_list1 [ 'Na' ])
    K = np.array(my_loaded_list1 [ 'K' ])
    Cl = np.array(my_loaded_list1 [ 'Cl' ])
    Gluc = np.array(my_loaded_list1 [ 'Gluc' ])
    x = [ "Na  K  Cl Gluc", "Na  K  Cl Gluc", "Na  K  Cl Gluc", "Na  K  Cl Gluc" ]

    bar1 = np.arange(len(x))
    bar2 = [ i + 2.0 * w for i in bar1 ]
    bar3 = [ i + 2.0 * w for i in bar2 ]
    bar4 = [ i + 2.0 * w for i in bar3 ]
    bar5 = [ i + w for i in bar2 ]
    plt.figure(1)
    plt.bar(bar1, Na, 0.2, color="#2B9F78", label="Original Setup")
    plt.bar(bar2, K, 0.2, color="#007282", label="NHE3 = 0")
    plt.bar(bar3, Cl, 0.2, color="black", label="SGLT = 0")
    plt.bar(bar4, Gluc, 0.2, color="#D55E00", label="NaH2PO4 = 0")
    plt.ylabel('Cellular Concentration', fontsize=12)
    plt.xlabel("Solutes", fontsize=12)
    plt.xticks(bar5, x, fontsize=11)
    plt.legend(frameon=False)
    plt.savefig(pltfolder + '/Figure_8B' + '.png')
    if Figure_7:
        of = open('Data_Figure_7.py', 'rb')
        read_file = f.read()
        my_loaded_list = pickle.loads(read_file)
        of.close()
        print('my_loaded_list', my_loaded_list['epithelial_flx_variation_default'])
        fig = plt.figure(num=None, figsize=(12, 12))
        x1=['ES', 'IE', 'IS', 'ME', 'MI','EPTL']
        ax1 = fig.add_subplot(3, 3, 1)
        ax1.bar(x1, np.array(my_loaded_list['flx_na_mem_default']), color=[ 'black', '#b66dff', '#24ff24', '#007282', '#2B9F78', '#920000'])
        # ax1.set_ylim(-0.2,7.5)
        plt.setp(plt.gca(), xticklabels=[])

        # Defines values for the second plot
        x2 = ['CSF', 'GLD', 'ECHMCL']
        print(type(my_loaded_list['epithelial_flx_variation_default']),my_loaded_list['epithelial_flx_variation_default'])
        ax2 = fig.add_subplot(3, 3, 2)

        Epithelial = my_loaded_list['epithelial_flx_variation_default']
        ax2.bar(x2, Epithelial, color=['black', '#b66dff','#2B9F78'])
        # ax2.set_ylim(-0.2,8.4)
        plt.setp(plt.gca(), xticklabels=[])
        x3 = ['NHE3', 'SGLT','NaH2PO4']
        ax3 = fig.add_subplot(3, 3, 3)
        ax3.bar(x3, my_loaded_list['epithelial_flx_elct_chmcl_default'], color=['#db6d00', '#b66dff', '#2B9F78'])
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
        ax6=fig.add_subplot(3, 3, 6)
        colors=['#db6d00', '#b66dff','#2B9F78']
        labels=['NHE3','SGLT','NaH2PO4']
        for tmpX, tmpY, color, label in zip(x3, my_loaded_list['epithelial_flx_elct_chmcl_nak'], colors, labels):
            ax6.bar(tmpX, tmpY, color=color, label=label)
        plt.setp(plt.gca(), xticklabels=[])
        # ax6.set_ylim(-0.1,3.8)
        ax6.legend(shadow=False, fancybox=False, frameon=False)
        ax7=fig.add_subplot(3, 3, 7)
        ax7.bar(x1, my_loaded_list['flx_na_mem_nhe3'], color=['black', '#b66dff','#24ff24', '#007282','#2B9F78','#920000'])
        ax7.set_xlabel('Membrane Flux_Na [mmol/s.cm2]', color='black')
        # ax7.set_ylim(-0.2,7.5)

        plt.setp(plt.gca(), xticklabels=[])
        # ax7.legend(shadow=False, fancybox=False, frameon=False)
        ax8=fig.add_subplot(3, 3, 8)
        ax8.bar(x2,my_loaded_list['epithelial_flx_variation_nhe3'], color=['black', '#b66dff','#2B9F78'])
        ax8.set_xlabel('Epithelial Flux_Na [mmol/s.cm2]', color='black')
        plt.setp(plt.gca(), xticklabels=[])
        # ax8.set_ylim(-0.2,8.4)
        # ax8.legend(shadow=False, fancybox=False, frameon=False)
        # label=['ES', 'IE', 'IS', 'ME', 'MI','EPTL']
        ax9 = fig.add_subplot(3, 3, 9)
        ax9.bar(x3, my_loaded_list['epithelial_flx_elct_chmcl_nhe3'], color=['#db6d00', '#b66dff', '#2B9F78'])
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