# A mathematical model of the proximal convoluted tubule base on Weinstein  model 2007.
# Leyla Noroozbabaee 11/1/2020
import os
import numpy as np
import scipy
import matplotlib.pyplot as plt
import math
from scipy import optimize
from PCT_GLOB_new import *


def sglt_mi(cm_na, ci_na, cm_gluc, ci_gluc, z_na, z_gluc, vm, vi, ami, lmi_nagluc, param_sglt_mi):
    # Na-glucose simple cotransporter with 1:1 stoichiometry, located on  Apical  Membrane
    # f_eps(): Modular function to calculate the electrochemical potential of species
    # see label: {eq:NAGluc}
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
    # see label {eq: NAH2PO4}
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


def nah_mi(cm_na, ci_na, cm_h, ci_h, z_na, z_h, vm, vi, ami, lmi_nah, param_nah_mi):
    # Weintein introduced two simple exchangers  na/h and na/nh4 in their work mentioned above and
    # it seems these two simple exchangers are equivalent to NHE transporter that they introduced in 1995.
    # They located on Apical membrane.
    # see label {eq:NAH}
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


def nanh4_mi(cm_na, ci_na, cm_nh4, ci_nh4 ,z_na, z_nh4, vm, vi, ami, lmi_nanh4, param_nanh4_mi):
    # Weintein introduced two simple exchangers  na/h and na/nh4 in their work mentioned above and
    # it seems these two simple exchangers are equivalent to NHE transporter that they introduced in 1995.
    # They located on Apical membrane.
    # see label {eq:NANH4}
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
    # see label {eq:CLHCO2}
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
    # cl/hco3 simple exchanger with 1:-1 stoichiometry, located on Apical  Membrane.
    # see label {eq: CLHCO3}
    if param_clhco3_mi == 0:
        return[0, 0]
    else:
        xm_cl = f_eps(cm_cl, z_cl, vm)
        xm_hco3 = f_eps(cm_hco3, z_hco3, vm)
        xi_cl = f_eps(ci_cl, z_cl, vi)
        xi_hco3 = f_eps(ci_hco3, z_hco3, vi)
        cl_mi_clhco3 = lmi_clhco3 * ami * (xm_cl - xi_cl - xm_hco3 + xi_hco3)
        hco3_mi_clhco3 = - lmi_clhco3 * ami * (xm_cl - xi_cl - xm_hco3 + xi_hco3)
    return [cl_mi_clhco3, hco3_mi_clhco3]


def kcl(ca_k, cb_k, ca_cl, cb_cl, z_k, z_cl, va, vb, a, l_kcl, param_kcl):
    # k-cl simple cotransporter with 1:1 stoichiometry, located on Peritubular Membrane which
    # includes both Cell-Lateral Membrane (ie) /Cell-Basal (is) Membrane.
    # see label {eq:KCL}
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


def nahco3(ca_na, cb_na, ca_hco3, cb_hco3, z_na, z_hco3, va, vb, a, l_nahco3, param_nahco3):
    # Na-hco3 simple cotransporter with 1:3 stoichiometry, located on  Peritubular Membrane which
    # includes both Cell-Lateral Membrane (ie) /Cell-Basal (is) Membrane.
    # see label {eq:Na-HCO3}
    if param_nahco3 == 0:
        return[0, 0]
    else:
        xa_na = f_eps(ca_na, z_na, va)
        xa_hco3 = f_eps(ca_hco3, z_hco3, va)
        xb_na = f_eps(cb_na, z_na, vb)
        xb_hco3 = f_eps(cb_hco3, z_hco3, vb)
        na_nahco3 = l_nahco3 * a * (xa_na - xb_na + 3 * (xa_hco3 - xb_hco3))
        hco3_nahco3 = 3 * l_nahco3 * a * (xa_na - xb_na + 3 * (xa_hco3 - xb_hco3))
        return [ na_nahco3, hco3_nahco3]


def na_clhco3(ca_na, cb_na, ca_cl, cb_cl, ca_hco3, cb_hco3, z_na, z_cl, z_hco3, va, vb, a, l_na_clhco3,
              param_na_clhco3):
    # na_clhco3 is a complex transporter 1:-1:2 stoichiometry, located on  Peritubular Membrane which
    # includes both Cell-Lateral Membrane (ie) /Cell-Basal (is) Membrane.
    # see label {eq:Na-HCO3-cl}
    if param_na_clhco3 == 0:
        return [0, 0, 0]
    else:
        xa_na = f_eps(ca_na, z_na, va)
        xa_cl = f_eps(ca_cl, z_cl, va)
        xa_hco3 = f_eps(ca_hco3, z_hco3, va)
        xb_na = f_eps(cb_na, z_na, vb)
        xb_cl = f_eps(cb_cl, z_cl, vb)
        xb_hco3 = f_eps(cb_hco3, z_hco3, vb)
        na_na_clhco3 = + a * l_na_clhco3 * (xa_na - xb_na - xa_cl + xb_cl + 2 * (xa_hco3 - xb_hco3))
        cl_na_clhco3 = - a * l_na_clhco3 * (xa_na - xb_na - xa_cl + xb_cl + 2 * (xa_hco3 - xb_hco3))
        hco3_na_clhco3 = 2 * a * l_na_clhco3 * (xa_na - xb_na - xa_cl + xb_cl + 2 * (xa_hco3 - xb_hco3))
        return [na_na_clhco3, cl_na_clhco3, hco3_na_clhco3]


def h_atp_mi(cm_h, ci_h, vm, vi, z_h, ami, parama_mi_h):
    # proton pumps located on Apical membrane.
    # f_eps(): Modular function to calculate the electrochemical potential of species
    # see label {h_ATP}
    if parama_mi_h == 0:
        return 0
    else:
        xm_h = f_eps(cm_h, z_h, vm)
        xi_h = f_eps(ci_h, z_h, vi)
        gamma = xihp * (xm_h - xi_h - xhp)
        return lhp * ami * (1 + 1 / (1 + math.exp(gamma)))


def nak_atp(ca_k, cb_k, ca_na, cb_na, cb_nh4, a, nak_atp_param):
    # Na-K pumps located on Peritubular Membrane which
    # includes both Cell-Lateral Membrane (ie) /Cell-Basal (is) Membrane.
    # see label {nak_ATP}
    if nak_atp_param == 0:
        return[0, 0, 0]
    else:
        # sodium affinity
        knpn = 0.0002 * (1.0 + ca_k / .00833)
        # potesium affinity
        knpk = 0.0001 * (1.0 + cb_na / .0185)
        # atpase transporter flux in is membrane
        na_atp = n_p * a * (ca_na / (ca_na + knpn )) ** 3 * (cb_k / (cb_k + knpk)) ** 2
        # allow for competition betWeen k+ and nh4+
        k_atp = -na_atp * a * 0.667 * cb_k / (cb_k + cb_nh4 / knh4)
        nh4_atp = -na_atp * a * 0.667 * cb_nh4 / (knh4 * cb_k + cb_nh4)
        return[na_atp, k_atp, nh4_atp]


def nah_nh4(ci_h, ci_na, ci_nh4, cm_h, cm_na, cm_nh4, param_nah):
    # the Na/H exchanger (NHE3) kinetic model by Weinstein in 1995.
    # However, Na/H exchanger introduced by Weinstein et al. 1992 using
    # two modular functions: mnah_mi() and nanh4_mi()

    if param_nah == 0:
        return[0, 0, 0]
    else:
        # the nah model parameters
        cxt = 1.00
        knah_na = 0.3000e-01
        knah_h = 0.7200e-07
        knah_nh4 = 0.2700e-01
        knah_i = 0.1000e-05

        pnah_na = 0.8000e+04
        pnah_h = 0.2400e+04
        pnah_nh4 = 0.8000e+04
        pnah_m = 0.2000e+01
        pnah_mm = 0.0000
        # translate concentrations to the nah model
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
        return[jnah_na, jnah_h, jnah_nh4]


def goldman(hab, a, z, va, vb, ca, cb, param_goldman):
    # Goldman Fluxes (passive fluxes) describes the ionic flux across a cell membrane
    # as a function of the transmembrane potential and the concentrations of the ion inside and outside of the cell.
    # see label {normalized_PD}
    zab = (z * f * (va - vb) * 1.e-6) / rte
    if param_goldman == 0:
        return[0]
    elif zab == 0 or va == vb and ca > 0 and cb > 0:
        # see label {Fick_passive}
        return hab * a * (ca - cb)
    else:
        # see label {Goldman_passive}
        return hab * a * zab * (ca - cb * math.exp(-zab)) / (1 - math.exp(-zab))


def csf(ca, cb, flux, s, param_csf):
    # convective solute fluxes
    if param_csf == 0:
        return 0

    def lmmsc(ca, cb):
        # logarithmic mean membrane solute concentration
        # see label {mean_concen}
        import math
        if ca > 0 and cb > 0 and ca - cb != 0:
            return (ca - cb) / (math.log10(ca / cb))
        else:
            print('ups')
            return cb
    # see label {convective_flux}
    # the input flux is the water fluxes, see label {water_flux}
    return flux * (1.00 - s) * lmmsc(ca, cb)


def matrix(init, h, tag):
    if tag == 'clever':
        a = 0.0e-4
        b = 1.0e0
    elif tag == 'arbitrary':
        a = 1.0e-4
        b = 0.0e0
    return[a+b*init for x in range(h)]


def f_eps(c, z, v):
    # electrochemical potential of species
    # see label {Electrochemical_P}
    if c > 0 and c != 0:
        return rte * math.log(c) + z * f * v * 1.e-6
    else:
        print('uy')

        return rte * math.log(abs(c)) + z * f * v * 1.e-6


def lch(ca, cb):
    import math
    if ca > 0 and cb > 0 and (ca / cb) != 0 and cb != 0:
        return math.log10(ca / cb)
    else:
        return math.log10(abs(ca / cb))


def ebuf(lch, pk, ca, cb, param_ebuf):
    # pH equilibria of four buffer pairs
    if param_ebuf == 0:
        return 0
    if ca > 0 and cb > 0 and (ca / cb) != 0 and cb != 0:
        return lch - pk - math.log10(ca / cb)

    else:
        return lch - pk - math.log10(abs(ca / cb))


# common definition
def zab(z, va, vb):
    return z * f * (va - vb) * 1.e-6 / rte


def phi_scale(phi, scale):
    return phi * scale


solver = 1

def eQs(guess, solver):
    # update variables
    for i in range(len(guess)):
        vars[i][t] = guess[i]
    # Weinstein (2007)
    # see label {ES_area}
    aes[t] = ae0 * (1.0 + mua * (pe[t] - pm))
    aes[t] = ae0 if aes[t] < ae0 else aes[t]
    # see label {LIS_volume}
    chvl[t] = chvl0 * (1.0 + muv * (pe[t] - pm))
    chvl[t] = chvl0 if chvl[t] < chvl0 else chvl[t]
    l[t] = chvl[t]

    # LOG Conc. of Hydrogen Interspace
    lche = pkc + lch(ce_hco3[t], ce_h2co3[t])
    ce_h[t] = 10 ** (-lche)

    # ELECTROLYTE TRANSPORT ACROSS A SIMPLE EPITHELIUM of the rat proximal tubule (1992)
    # The cell volume is presented as a function of the inverse cell impermeant species concentration
    # see label {cell_volume}
    clvl[t] = clvl0 * imp0 / imp[t]
    clvl[t] = clvl0 if imp[t] == 0 else clvl[t]

    # The cell height (summation of Extracellular channel volume
    # and Intracellular compartment volume[cm3/cm2 epithelium])
    # see label {epithelial_thickness}
    l[t] = l[t] + clvl[t]
    p_i = pm

    # LOG Conc. of Hydrogen Cell
    lchi = pkc + lch(ci_hco3[t], ci_h2co3[t])
    ci_h[t] = 10 ** (-lchi)
    # ELECTROLYTE TRANSPORT ACROSS A SIMPLE EPITHELIUM 1979
    # Applying the Kedem and Katchalsky relations:
    # The membrane transport for nonelectrolyte solutions generated by the hydrostatic
    # pressure difference and the osmotic pressure difference
    # The equations below represent the volume flow or the convective
    # volume flow in between the different membrane
    # The order of the two characters indicates the positive direction for mass flow.
    # see label {water_flux} for the first two terms, where π = rt*C_imp
    fevm = ame * lpme * (pm - pe[t] - rt * impm) / rt
    fevs = aes[t] * lpes * (pe[t] - ps + rt * imps) / rt
    fivm = ami * lpmi * (pm - p_i + rt * imp[t] - rt * impm) / rt
    fivs = ais * lpis * (p_i - ps + rt * imps - rt * imp[t]) / rt
    five = aie * lpie * (p_i - pe[t] - rt * imp[t]) / rt  #jv
    # Updating the volume flow with the last term
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
                    ce_nh4[t] - cm_nh4) + sme_hco2 * (
                    ce_hco2[t] - cm_hco2) + sme_h2co2 * (
                    ce_h2co2[t] - cm_h2co2) + sme_gluc * (
                    ce_gluc[t] - cm_gluc))

    fevs = fevs + aes[t] * lpes * (
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

    five = five + lpie * aie * (sie_na * (ce_na[t] - ci_na[t]) + sie_k * (
                    ce_k[t] - ci_k[t]) + sie_cl * (
                    ce_cl[t] - ci_cl[t]) + sie_hco3 * (
                    ce_hco3[t] - ci_hco3[t]) + sie_h2co3 * (
                    ce_h2co3[t] - ci_h2co3[t]) + sie_co2 * (
                    ce_co2[t] - ci_co2[t]) + sie_hpo4 * (
                    ce_hpo4[t] - ci_hpo4[t]) + sie_h2po4 * (
                    ce_h2po4[t] - ci_h2po4[t]) + sie_urea * (
                    ce_urea[t] - ci_urea[t]) + sie_nh3 * (
                    ce_nh3[t] - ci_nh3[t]) + sie_nh4 * (
                    ce_nh4[t] - ci_nh4[t]) + sie_h * (
                    ce_h[t] - ci_h[t]) + sie_hco2 * (
                    ce_hco2[t] - ci_hco2[t]) + sie_h2co2 * (
                    ce_h2co2[t] - ci_h2co2[t]) + sie_gluc * (
                    ce_gluc[t] - ci_gluc[t]))

    # Convective Intraepithelial Solute Fluxes
    # see label {convective_flux}
    param_csf = 1
    fekm_na = csf(ce_na[t], cm_na, fevm, sme_na, param_csf)
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

    fike_na = csf(ci_na[t], ce_na[t], five, sie_na, param_csf)
    fike_k = csf(ci_k[t], ce_k[t], five, sie_k, param_csf)
    fike_cl = csf(ci_cl[t], ce_cl[t], five, sie_cl, param_csf)
    fike_hco3 = csf(ci_hco3[t], ce_hco3[t], five, sie_hco3, param_csf)
    fike_h2co3 = csf(ci_h2co3[t], ce_h2co3[t], five, sie_h2co3, param_csf)
    fike_co2 = csf(ci_co2[t], ce_co2[t], five, sie_co2, param_csf)
    fike_hpo4 = csf(ci_hpo4[t], ce_hpo4[t], five, sie_hpo4, param_csf)
    fike_h2po4 = csf(ci_h2po4[t], ce_h2po4[t], five, sie_h2po4, param_csf)
    fike_urea = csf(ci_urea[t], ce_urea[t], five, sie_urea, param_csf)
    fike_nh3 = csf(ci_nh3[t], ce_nh3[t], five, sie_nh3, param_csf)
    fike_nh4 = csf(ci_nh4[t], ce_nh4[t], five, sie_nh4, param_csf)
    fike_h = csf(ci_h[t], ce_h[t], five, sie_h, param_csf)
    fike_hco2 = csf(ci_hco2[t], ce_hco2[t], five, sie_hco2, param_csf)
    fike_h2co2 = csf(ci_h2co2[t], ce_h2co2[t], five, sie_h2co2, param_csf)
    fike_gluc = csf(ci_gluc[t], ce_gluc[t], five, sie_gluc, param_csf)

    # Goldman Fluxes: fluxes due to diffusion and electrical potential difference for
    # all of the ions that are permeable through the membrane
    # see label {Goldman_passive}
    goldman_param = 1
    fekm_na = fekm_na + goldman(hme_na, ame, z_na, vm[t], ve[t], cm_na,
                                ce_na[t], goldman_param)
    fekm_k = fekm_k + goldman(hme_k, ame, z_k, vm[t], ve[t], cm_k, ce_k[t], goldman_param)
    fekm_cl = fekm_cl + goldman(hme_cl, ame, z_cl, vm[t], ve[t], cm_cl, ce_cl[t], goldman_param)
    fekm_hco3 = fekm_hco3 + goldman(hme_hco3, ame, z_hco3, vm[t], ve[t], cm_hco3, ce_hco3[t], goldman_param)
    fekm_h2co3 = fekm_h2co3 + goldman(hme_h2co3, ame, z_h2co3, vm[t], ve[t], cm_h2co3,
                                      ce_h2co3[t], goldman_param)
    fekm_co2 = fekm_co2 + goldman(hme_co2, ame, z_co2, vm[t], ve[t], cm_co2, ce_co2[t], goldman_param)
    fekm_hpo4 = fekm_hpo4 + goldman(hme_hpo4, ame, z_hpo4, vm[t], ve[t], cm_hpo4,
                                    ce_hpo4[t], goldman_param)
    fekm_h2po4 = fekm_h2po4 + goldman(hme_h2po4, ame, z_h2po4, vm[t], ve[t], cm_h2po4,
                                      ce_h2po4[t], goldman_param)
    fekm_urea = fekm_urea + goldman(hme_urea, ame, z_urea, vm[t], ve[t], cm_urea, ce_urea[t], goldman_param)
    fekm_nh3 = fekm_nh3 + goldman(hme_nh3, ame, z_nh3, vm[t], ve[t], cm_nh3, ce_nh3[t], goldman_param)
    fekm_nh4 = fekm_nh4 + goldman(hme_nh4, ame, z_nh4, vm[t], ve[t], cm_nh4, ce_nh4[t], goldman_param)
    fekm_h = fekm_h + goldman(hme_h, ame, z_h, vm[t], ve[t], cm_h, ce_h[t], goldman_param)
    fekm_hco2 = fekm_hco2 + goldman(hme_hco2, ame, z_hco2, vm[t], ve[t], cm_hco2, ce_hco2[t], goldman_param)
    fekm_h2co2 = fekm_h2co2 + goldman(hme_h2co2, ame, z_h2co2, vm[t], ve[t], cm_h2co2, ce_h2co2[t],
                                      goldman_param)
    fekm_gluc = fekm_gluc + goldman(hme_gluc, ame, z_gluc, vm[t], ve[t], cm_gluc, ce_gluc[t], goldman_param)

    feks_na = feks_na + goldman(hes_na, aes[t], z_na, ve[t], vs, ce_na[t],
                                cs_na, goldman_param)
    feks_k = feks_k + goldman(hes_k, aes[t], z_k, ve[t], vs, ce_k[t], cs_k, goldman_param)
    feks_cl = feks_cl + goldman(hes_cl, aes[t], z_cl, ve[t], vs, ce_cl[t], cs_cl, goldman_param)
    feks_hco3 = feks_hco3 + goldman(hes_hco3, aes[t], z_hco3, ve[t], vs, ce_hco3[t], cs_hco3, goldman_param)
    feks_h2co3 = feks_h2co3 + goldman(hes_h2co3, aes[t], z_h2co3, ve[t], vs, ce_h2co3[t],
                                      cs_h2co3, goldman_param)
    feks_co2 = feks_co2 + goldman(hes_co2, aes[t], z_co2, ve[t], vs, ce_co2[t], cs_co2, goldman_param)
    feks_hpo4 = feks_hpo4 + goldman(hes_hpo4, aes[t], z_hpo4, ve[t], vs, ce_hpo4[t],
                                    cs_hpo4, goldman_param)
    feks_h2po4 = feks_h2po4 + goldman(hes_h2po4, aes[t], z_h2po4, ve[t], vs, ce_h2po4[t],
                                      cs_h2po4, goldman_param)
    feks_urea = feks_urea + goldman(hes_urea, aes[t], z_urea, ve[t], vs, ce_urea[t], cs_urea, goldman_param)
    feks_nh3 = feks_nh3 + goldman(hes_nh3, aes[t], z_nh3, ve[t], vs, ce_nh3[t], cs_nh3, goldman_param)
    feks_nh4 = feks_nh4 + goldman(hes_nh4, aes[t], z_nh4, ve[t], vs, ce_nh4[t], cs_nh4, goldman_param)
    feks_h = feks_h + goldman(hes_h, aes[t], z_h, ve[t], vs, ce_h[t], cs_h, goldman_param)
    feks_hco2 = feks_hco2 + goldman(hes_hco2, aes[t], z_hco2, ve[t], vs, ce_hco2[t], cs_hco2, goldman_param)
    feks_h2co2 = feks_h2co2 + goldman(hes_h2co2, aes[t], z_h2co2, ve[t], vs, ce_h2co2[t], cs_h2co2, goldman_param)
    feks_gluc = feks_gluc + goldman(hes_gluc, aes[t], z_gluc, ve[t], vs, ce_gluc[t], cs_gluc, goldman_param)

    fikm_na = fikm_na + goldman(hmi_na, ami, z_na, vm[t], vi[t], cm_na,
                                ci_na[t], goldman_param)
    fikm_k = fikm_k + goldman(hmi_k, ami, z_k, vm[t], vi[t], cm_k, ci_k[t], goldman_param)
    fikm_cl = fikm_cl + goldman(hmi_cl, ami, z_cl, vm[t], vi[t], cm_cl, ci_cl[t], goldman_param)
    fikm_hco3 = fikm_hco3 + goldman(hmi_hco3, ami, z_hco3, vm[t], vi[t], cm_hco3, ci_hco3[t], goldman_param)
    fikm_h2co3 = fikm_h2co3 + goldman(hmi_h2co3, ami, z_h2co3, vm[t], vi[t], cm_h2co3,
                                      ci_h2co3[t], goldman_param)
    fikm_co2 = fikm_co2 + goldman(hmi_co2, ami, z_co2, vm[t], vi[t], cm_co2, ci_co2[t], goldman_param)
    fikm_hpo4 = fikm_hpo4 + goldman(hmi_hpo4, ami, z_hpo4, vm[t], vi[t], cm_hpo4,
                                    ci_hpo4[t], goldman_param)
    fikm_h2po4 = fikm_h2po4 + goldman(hmi_h2po4, ami, z_h2po4, vm[t], vi[t], cm_h2po4,
                                      ci_h2po4[t], goldman_param)
    fikm_urea = fikm_urea + goldman(hmi_urea, ami, z_urea, vm[t], vi[t], cm_urea, ci_urea[t], goldman_param)
    fikm_nh3 = fikm_nh3 + goldman(hmi_nh3, ami, z_nh3, vm[t], vi[t], cm_nh3, ci_nh3[t], goldman_param)
    fikm_nh4 = fikm_nh4 + goldman(hmi_nh4, ami, z_nh4, vm[t], vi[t], cm_nh4, ci_nh4[t], goldman_param)
    fikm_h = fikm_h + goldman(hmi_h, ami, z_h, vm[t], vi[t], cm_h, ci_h[t], goldman_param)
    fikm_hco2 = fikm_hco2 + goldman(hmi_hco2, ami, z_hco2, vm[t], vi[t], cm_hco2, ci_hco2[t], goldman_param)
    fikm_h2co2 = fikm_h2co2 + goldman(hmi_h2co2, ami, z_h2co2, vm[t], vi[t], cm_h2co2,
                                      ci_h2co2[t], goldman_param)
    fikm_gluc = fikm_gluc + goldman(hmi_gluc, ami, z_gluc, vm[t], vi[t], cm_gluc, ci_gluc[t], goldman_param)

    fike_na = fike_na + goldman(hie_na, aie, z_na, vi[t], ve[t], ci_na[t], ce_na[t], goldman_param)
    fike_k = fike_k + goldman(hie_k, aie, z_k, vi[t], ve[t], ci_k[t], ce_k[t], goldman_param)
    fike_cl = fike_cl + goldman(hie_cl, aie, z_cl, vi[t], ve[t], ci_cl[t], ce_cl[t], goldman_param)
    fike_hco3 = fike_hco3 + goldman(hie_hco3, aie, z_hco3, vi[t], ve[t], ci_hco3[t], ce_hco3[t], goldman_param)
    fike_h2co3 = fike_h2co3 + goldman(hie_h2co3, aie, z_h2co3, vi[t], ve[t], ci_h2co3[t], ce_h2co3[t], goldman_param)
    fike_co2 = fike_co2 + goldman(hie_co2, aie, z_co2, vi[t], ve[t], ci_co2[t], ce_co2[t], goldman_param)
    fike_hpo4 = fike_hpo4 + goldman(hie_hpo4, aie, z_hpo4, vi[t], ve[t], ci_hpo4[t], ce_hpo4[t], goldman_param)
    fike_h2po4 = fike_h2po4 + goldman(hie_h2po4, aie, z_h2po4, vi[t], ve[t], ci_h2po4[t], ce_h2po4[t], goldman_param)
    fike_urea = fike_urea + goldman(hie_urea, aie, z_urea, vi[t], ve[t], ci_urea[t], ce_urea[t], goldman_param)
    fike_nh3 = fike_nh3 + goldman(hie_nh3, aie, z_nh3, vi[t], ve[t], ci_nh3[t], ce_nh3[t], goldman_param)
    fike_nh4 = fike_nh4 + goldman(hie_nh4, aie, z_nh4, vi[t], ve[t], ci_nh4[t], ce_nh4[t], goldman_param)
    fike_h = fike_h + goldman(hie_h, aie, z_h, vi[t], ve[t], ci_h[t], ce_h[t], goldman_param)
    fike_hco2 = fike_hco2 + goldman(hie_hco2, aie, z_hco2, vi[t], ve[t], ci_hco2[t], ce_hco2[t], goldman_param)
    fike_h2co2 = fike_h2co2 + goldman(hie_h2co2, aie, z_h2co2, vi[t], ve[t], ci_h2co2[t], ce_h2co2[t], goldman_param)
    fike_gluc = fike_gluc + goldman(hie_gluc, aie, z_gluc, vi[t], ve[t], ci_gluc[t], ce_gluc[t], goldman_param)

    fiks_na = fiks_na + goldman(his_na, ais, z_na, vi[t], vs, ci_na[t], cs_na, goldman_param)
    fiks_k = fiks_k + goldman(his_k, ais, z_k, vi[t], vs, ci_k[t], cs_k, goldman_param)
    fiks_cl = fiks_cl + goldman(his_cl, ais, z_cl, vi[t], vs, ci_cl[t], cs_cl, goldman_param)
    fiks_hco3 = fiks_hco3 + goldman(his_hco3, ais, z_hco3, vi[t], vs, ci_hco3[t], cs_hco3, goldman_param)
    fiks_h2co3 = fiks_h2co3 + goldman(his_h2co3, ais, z_h2co3, vi[t], vs, ci_h2co3[t],
                                      cs_h2co3, goldman_param)
    fiks_co2 = fiks_co2 + goldman(his_co2, ais, z_co2, vi[t], vs, ci_co2[t], cs_co2, goldman_param)
    fiks_hpo4 = fiks_hpo4 + goldman(his_hpo4, ais, z_hpo4, vi[t], vs, ci_hpo4[t],
                                    cs_hpo4, goldman_param)
    fiks_h2po4 = fiks_h2po4 + goldman(his_h2po4, ais, z_h2po4, vi[t], vs, ci_h2po4[t],
                                      cs_h2po4, goldman_param)
    fiks_urea = fiks_urea + goldman(his_urea, ais, z_urea, vi[t], vs, ci_urea[t], cs_urea, goldman_param)
    fiks_nh3 = fiks_nh3 + goldman(his_nh3, ais, z_nh3, vi[t], vs, ci_nh3[t], cs_nh3, goldman_param)
    fiks_nh4 = fiks_nh4 + goldman(his_nh4, ais, z_nh4, vi[t], vs, ci_nh4[t], cs_nh4, goldman_param)
    fiks_h = fiks_h + goldman(his_h, ais, z_h, vi[t], vs, ci_h[t], cs_h, goldman_param)
    fiks_hco2 = fiks_hco2 + goldman(his_hco2, ais, z_hco2, vi[t], vs, ci_hco2[t], cs_hco2, goldman_param)
    fiks_h2co2 = fiks_h2co2 + goldman(his_h2co2, ais, z_h2co2, vi[t], vs, ci_h2co2[t],
                                      cs_h2co2, goldman_param)
    fiks_gluc = fiks_gluc + goldman(his_gluc, ais, z_gluc, vi[t], vs, ci_gluc[t], cs_gluc, goldman_param)

    # Updating flux equations with electodiffusive fluxes (transporters).
    # Electodiffusive fluxes  are proportional to the differences in electrochemical driving forces
    # progressive activation of Transporters along time.
    if 0 < t <= int(4000*Nfac):
        sglt_mi_param = 0
        nah2po4_mi_param = 0
        clhco3_mi_param = 0
        clhco2_mi_param = 0
        nahco3_is_param = 0
        nahco3_ie_param = 0
        kcl_is_param = 0
        kcl_ie_param = 0
        na_clhco3_is_param = 0
        na_clhco3_ie_param = 0
        nah_param = 0
        nak_atp_param = 0
        h_mi_atp_param = 0

    elif int(4000*Nfac) < t <= int(8000*Nfac):
        sglt_mi_param = 1
        nah2po4_mi_param = 0
        clhco3_mi_param = 0
        clhco2_mi_param = 0
        nahco3_is_param = 0
        nahco3_ie_param = 0
        kcl_is_param = 0
        kcl_ie_param = 0
        na_clhco3_is_param = 0
        na_clhco3_ie_param = 0
        nah_param = 0
        nak_atp_param = 0
        h_mi_atp_param = 0

    elif int(8000*Nfac) < t <= int(12000*Nfac):
        sglt_mi_param = 1
        nah2po4_mi_param = 1
        clhco3_mi_param = 0
        clhco2_mi_param = 0
        nahco3_is_param = 0
        nahco3_ie_param = 0
        kcl_is_param = 0
        kcl_ie_param = 0
        na_clhco3_is_param = 0
        na_clhco3_ie_param = 0
        nah_param = 0
        nak_atp_param = 0
        h_mi_atp_param = 0

    elif int(12000*Nfac) < t <= int(16000*Nfac):
        sglt_mi_param = 1
        nah2po4_mi_param = 1
        clhco3_mi_param = 1
        clhco2_mi_param = 0
        nahco3_is_param = 0
        nahco3_ie_param = 0
        kcl_is_param = 0
        kcl_ie_param = 0
        na_clhco3_is_param = 0
        na_clhco3_ie_param = 0
        nah_param = 0
        nak_atp_param = 0
        h_mi_atp_param = 0

    elif int(16000*Nfac) < t <= int(20000*Nfac):
        sglt_mi_param = 1
        nah2po4_mi_param = 1
        clhco3_mi_param = 1
        clhco2_mi_param = 1
        nahco3_is_param = 0
        nahco3_ie_param = 0
        kcl_is_param = 0
        kcl_ie_param = 0
        na_clhco3_is_param = 0
        na_clhco3_ie_param = 0
        nah_param = 0
        nak_atp_param = 0
        h_mi_atp_param = 0

    elif int(20000*Nfac) < t <= int(24000*Nfac):
        sglt_mi_param = 1
        nah2po4_mi_param = 1
        clhco3_mi_param = 1
        clhco2_mi_param = 1
        nahco3_is_param = 1
        nahco3_ie_param = 1
        kcl_is_param = 0
        kcl_ie_param = 0
        na_clhco3_is_param = 0
        na_clhco3_ie_param = 0
        nah_param = 0
        nak_atp_param = 0
        h_mi_atp_param = 0

    elif int(24000*Nfac) < t <= int(28000*Nfac):
        sglt_mi_param = 1
        nah2po4_mi_param = 1
        clhco3_mi_param = 1
        clhco2_mi_param = 1
        nahco3_is_param = 1
        nahco3_ie_param = 1
        kcl_is_param = 1
        kcl_ie_param = 1
        na_clhco3_is_param = 0
        na_clhco3_ie_param = 0
        nah_param = 0
        nak_atp_param = 0
        h_mi_atp_param = 0

    elif int(28000*Nfac) < t <= int(32000*Nfac):
        sglt_mi_param = 1
        nah2po4_mi_param = 1
        clhco3_mi_param = 1
        clhco2_mi_param = 1
        nahco3_is_param = 1
        nahco3_ie_param = 1
        kcl_is_param = 1
        kcl_ie_param = 1
        na_clhco3_is_param = 1
        na_clhco3_ie_param = 1
        nah_param = 0
        nak_atp_param = 0
        h_mi_atp_param = 0

    elif int(32000*Nfac) < t <= int(36000*Nfac):
        sglt_mi_param = 1
        nah2po4_mi_param = 1
        clhco3_mi_param = 1
        clhco2_mi_param = 1
        nahco3_is_param = 1
        nahco3_ie_param = 1
        kcl_is_param = 1
        kcl_ie_param = 1
        na_clhco3_is_param = 1
        na_clhco3_ie_param = 1
        nah_param = 1
        nak_atp_param = 0
        h_mi_atp_param = 0

    elif int(36000*Nfac) < t <= int(40000*Nfac):
        sglt_mi_param = 1
        nah2po4_mi_param = 1
        clhco3_mi_param = 1
        clhco2_mi_param = 1
        nahco3_is_param = 1
        nahco3_ie_param = 1
        kcl_is_param = 1
        kcl_ie_param = 1
        na_clhco3_is_param = 1
        na_clhco3_ie_param = 1
        nah_param = 1
        nak_atp_param = 1
        h_mi_atp_param = 0

    else:
        sglt_mi_param = 1
        nah2po4_mi_param = 1
        clhco3_mi_param = 1
        clhco2_mi_param = 1
        nahco3_is_param = 1
        nahco3_ie_param = 1
        kcl_is_param = 1
        kcl_ie_param = 1
        na_clhco3_is_param = 1
        na_clhco3_ie_param = 1
        nah_param = 1
        nak_atp_param = 1
        h_mi_atp_param = 1

    # Net electodiffusive fluxes on the mi border
    # see label {eq: NAGLUC}
    sglt = sglt_mi(cm_na, ci_na[t], cm_gluc, ci_gluc[t], z_na, z_gluc, vm[t], vi[t], ami, lmi_nagluc,
                   sglt_mi_param)
    na_mi_nagluc, gluc_mi_nagluc = sglt
    # see label {eq: NAH2PO4}
    nah2po4 = nah2po4_mi(cm_na, ci_na[t], cm_h2po4, ci_h2po4[t], z_na, z_h2po4, vm[t], vi[t], ami,
                         lmi_nah2po4, nah2po4_mi_param)
    na_mi_nah2po4, h2po4_mi_nah2po4 = nah2po4
    # see label {eq: CLHCO3}
    clhco3 = clhco3_mi(cm_cl, ci_cl[t], cm_hco3, ci_hco3[t], z_cl, z_hco3, vm[t], vi[t], ami,
                       lmi_clhco3, clhco3_mi_param)
    cl_mi_clhco3, hco3_mi_clhco3 = clhco3
    # see label {eq: CLHCO2}
    clhco2 = clhco2_mi(cm_cl, ci_cl[t], cm_hco2, ci_hco2[t], z_cl, z_hco2, vm[t], vi[t], ami,
                       lmi_clhco2, clhco2_mi_param)
    cl_mi_clhco2, hco2_mi_clhco2 = clhco2
    # the nah exchanger translate concentrations to the nah kinetic model on the mi border
    mynah = nah_nh4(ci_h[t], ci_na[t], ci_nh4[t], cm_h, cm_na, cm_nh4, nah_param)
    jnah_na = mynah[0]
    jnah_h = mynah[1]
    jnah_nh4 = mynah[2]
    na_mi_nhe = nnhe3 * ami * jnah_na
    h_mi_nhe = nnhe3 * ami * jnah_h
    nh4_mi_nhe = nnhe3 * ami * jnah_nh4
    # Updating apical fluxes by including the electodiffusive fluxes
    fikm_na = fikm_na + na_mi_nagluc + na_mi_nah2po4 + na_mi_nhe
    fikm_cl = fikm_cl + cl_mi_clhco2 + cl_mi_clhco3
    fikm_hco3 = fikm_hco3 + hco3_mi_clhco3
    fikm_h2po4 = fikm_h2po4 + h2po4_mi_nah2po4
    fikm_hco2 = fikm_hco2 + hco2_mi_clhco2
    fikm_gluc = fikm_gluc + gluc_mi_nagluc
    fikm_h = fikm_h + h_mi_nhe
    fikm_nh4 = fikm_nh4 + nh4_mi_nhe

    # Net electodiffusive fluxes on the is border
    # see label {eq: Na_HCO3}
    nahco3_is = nahco3(ci_na[t], cs_na, ci_hco3[t], cs_hco3, z_na, z_hco3, vi[t], vs, ais,
                       lis_nahco3, nahco3_is_param)
    na_is_nahco3, hco3_is_nahco3 = nahco3_is
    # see label {eq: KCL}
    kcl_is = kcl(ci_k[t], cs_k, ci_cl[t], cs_cl, z_k, z_cl, vi[t], vs, ais, lis_kcl, kcl_is_param)
    k_is_kcl, cl_is_kcl = kcl_is
    # see label {eq: Na_HCO3_CL}
    na_clhco3_is = na_clhco3(ci_na[t], cs_na, ci_cl[t], cs_cl, ci_hco3[t], cs_hco3, z_na, z_cl, z_hco3,
                             vi[t], vs, ais, lis_na_clhco3, na_clhco3_is_param)
    na_is_na_clhco3, cl_is_na_clhco3, hco3_is_na_clhco3 = na_clhco3_is
    # Updating basal fluxes by including the electodiffusive fluxes
    fiks_na = fiks_na + na_is_nahco3 + na_is_na_clhco3
    fiks_k = fiks_k + k_is_kcl
    fiks_cl = fiks_cl + cl_is_kcl + cl_is_na_clhco3
    fiks_hco3 = fiks_hco3 + hco3_is_na_clhco3 + hco3_is_nahco3

    # Net electodiffusive fluxes on the ie border
    # see label {eq: Na_HCO3}
    nahco3_ie = nahco3(ci_na[t], ce_na[t], ci_hco3[t], ce_hco3[t], z_na, z_hco3, vi[t], ve[t],
                           aie,
                           lie_nahco3, nahco3_ie_param)
    na_ie_nahco3, hco3_ie_nahco3 = nahco3_ie
    # see label {eq: KCL}
    kcl_ie = kcl(ci_k[t], ce_k[t], ci_cl[t], ce_cl[t], z_k, z_cl, vi[t], ve[t], aie, lie_kcl, kcl_ie_param)
    k_ie_kcl , cl_ie_kcl  = kcl_ie
    # see label {eq: Na_HCO3_CL}
    na_clhco3_ie = na_clhco3(ci_na[t], ce_na[t], ci_cl[t], ce_cl[t], ci_hco3[t], ce_hco3[t], z_na,
                                 z_cl, z_hco3,
                                 vi[t], ve[t], aie, lie_na_clhco3, na_clhco3_ie_param)
    # see label {eq: Na_HCO3_CL}
    na_ie_na_clhco3, cl_ie_na_clhco3, hco3_ie_na_clhco3 = na_clhco3_ie
    # Updating lateral fluxes by including the electodiffusive fluxes
    fike_na = fike_na + na_ie_nahco3 + na_ie_na_clhco3
    fike_k = fike_k + k_ie_kcl
    fike_cl = fike_cl + cl_ie_kcl + cl_ie_na_clhco3
    fike_hco3 = fike_hco3 + hco3_ie_na_clhco3 + hco3_ie_nahco3
    # NaK-ATPase (active transporter or sodium pumps on the is border
    # see label {eq: NaK_ATP}
    nak = nak_atp(ci_k[t], cs_k, ci_na[t], cs_na, cs_nh4, ais, nak_atp_param)
    na_is_atp, k_is_atp, nh4_is_atp = nak
    fiks_na = fiks_na + na_is_atp
    fiks_k = fiks_k + k_is_atp
    fiks_nh4 = fiks_nh4 + nh4_is_atp
    # NaK-ATPase (active transporter or sodium pumps on the ie border
    nak = nak_atp(ci_k[t], ce_k[t], ci_na[t], ce_na[t], ce_nh4[t], aie, nak_atp_param)
    na_ie_atp, k_ie_atp, nh4_ie_atp = nak
    fike_na = fike_na + na_ie_atp
    fike_k = fike_k + k_ie_atp
    fike_nh4 = fike_nh4 + nh4_ie_atp

    # H-ATPase (active transporter or proton pumps on the mi border# proton pumps
    # see label {eq: H_ATP}
    atmi_h = h_atp_mi(cm_h, ci_h[t], vm[t], vi[t], z_h, ami, h_mi_atp_param)
    # updating the related flux on the apical membrane
    fikm_h = fikm_h + atmi_h
    # establish the error vectors, the "phi" array.
    # first for the interspace electroneutrality
    # see label {conser_charge_LIS}
    phie_en = 0
    phie_en = phie_en + z_na * ce_na[t] + z_k * ce_k[t] + z_cl * ce_cl[
        t] + z_hco3 * ce_hco3[t] + z_h2co3 * ce_h2co3[t] + z_co2 * ce_co2[t] + z_hpo4 * ce_hpo4[
                  t] + z_h2po4 * ce_h2po4[t] + z_urea * ce_urea[t] + z_nh3 * ce_nh3[t] + z_nh4 * ce_nh4[
                  t] + z_h * ce_h[t] + z_hco2 * ce_hco2[t] + z_h2co2 * ce_h2co2[t] + z_gluc * ce_gluc[t]
    # see label {conser_charge_cell}
    phii_en = imp[t] * zimp
    phii_en = phii_en - cbuf[t] + z_na * ci_na[t] + z_k * ci_k[t] + z_cl * ci_cl[t] + z_hco3 *\
        ci_hco3[t] + z_h2co3 * ci_h2co3[t] + z_co2 * ci_co2[t] + z_hpo4 * \
        ci_hpo4[t] + z_h2po4 * ci_h2po4[t] + z_urea * ci_urea[t] + z_nh3 * ci_nh3[t] + z_nh4 * \
        ci_nh4[t] + z_h * ci_h[t] + z_hco2 * ci_hco2[t] + z_h2co2 * ci_h2co2[t] + z_gluc * ci_gluc[t]

    # mass conservation in the time-dependent case
    # see label {volume_gen_LIS}
    phie_vlm = fevs - fevm - five + rtau * (chvl[t] - chvl[t - 1])
    # see label {solute_gen_LIS}
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

    # mass conservation in the time - dependent case
    # see label {volume_gen_cell}
    phii_vlm = fivs - fivm + five + rtau * (clvl[t] - clvl[t - 1])
    # see label {solute_gen_cell}
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

    # the proton flux must include the cellular buffers
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

    # now set the phis in terms of the solute generation rates
    # first the non-reactive species
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

        # co2, formate, phosphate, and ammonia content:
        if int(40000*Nfac) < t:
            # LIS: the error term for bicarbonate generation is replaced by conservation
            # of charge in the buffer reactions.
            # see label {conser_charge_buf}
            qe_hco3 = + qe_h + qe_nh4 - qe_hco3 + qe_h2co2 + qe_h2po4
            phi[5] = qe_hco3
            # LIS: hydration and dhydration of co2
            # see label {co2_hyd_dhyd}
            qe_co2 = qe_co2 + 1.e6 * chvl[t] * (khy_4 * ce_co2[t] - kdhy_4 * ce_h2co3[t])
            phi[7] = qe_co2
            # LIS: see label {conser_charge_co2}
            qe_h2co3 = qe_hco3 + qe_h2co3 + qe_co2
            phi[6] = qe_h2co3
            # cell: see label {conser_charge_buf}
            qi_hco3 = - qi_hco3 + qi_h2po4 + qi_nh4 + qi_h + qi_h2co2
            phi[21] = qi_hco3
            # cell: hydration and dhydration of co2
            # see label {co2_hyd_dhyd}
            qi_co2 = qi_co2 + 1.e6 * clvl[t] * (khy * ci_co2[t] - kdhy * ci_h2co3[t])
            phi[23] = qi_co2
            # cell: see label {conser_charge_co2}
            qi_h2co3 = qi_hco3 + qi_h2co3 + qi_co2
            phi[22] = qi_h2co3
        else:

            phi[6] = qe_h2co3
            phi[7] = qe_co2
            phi[5] = qe_hco3
            phi[21] = qi_hco3
            phi[23] = qi_co2
            phi[22] = qi_h2co3
        if int(44000*Nfac) < t:
            # LIS: see label {phosphate}
            qe_hpo4 = qe_hpo4 + qe_h2po4
            phi[8] = qe_hpo4
            # LIS: see labela {pH_equilibria}, used as paired equations
            dihydrogen_phosphate_param_e = 1
            qe_h2po4 = ebuf(lche, pkp, ce_hpo4[t], ce_h2po4[t], dihydrogen_phosphate_param_e)
            phi[9] = qe_h2po4
            # cell: see label {phosphate}
            qi_hpo4 = qi_hpo4 + qi_h2po4
            phi[24] = qi_hpo4
            # cell: see labela {pH_equilibria}, used as paired equations
            dihydrogen_phosphate_param_i = 1
            qi_h2po4 = ebuf(lchi, pkp, ci_hpo4[t], ci_h2po4[t], dihydrogen_phosphate_param_i)
            phi[25] = qi_h2po4

        else:

            phi[8] = qe_hpo4
            phi[9] = qe_h2po4
            phi[24] = qi_hpo4
            phi[25] = qi_h2po4
    if int(48000*Nfac) < t:
        # LIS: see label {ammonia}
        qe_nh3 = qe_nh3 + qe_nh4
        phi[11] = qe_nh3
        # LIS: see labela {pH_equilibria}, used as paired equations
        ammonium_param_e = 1
        qe_nh4 = ebuf(lche, pkn, ce_nh3[t], ce_nh4[t], ammonium_param_e)
        phi[12] = qe_nh4
        # cell: see label {ammonia_Q}
        qi_nh3 = qi_nh3 + qi_nh4 - 1e6 * qiamm
        phi[27] = qi_nh3
        # cell: see labela {pH_equilibria}, used as paired equations for {ammonia_Q}
        ammonium_param_i = 1
        qi_nh4 = ebuf(lchi, pkn, ci_nh3[t], ci_nh4[t], ammonium_param_i)
        phi[28] = qi_nh4
    else:

        phi[11] = qe_nh3
        phi[12] = qe_nh4
        phi[27] = qi_nh3
        phi[28] = qi_nh4

    if int(52000*Nfac) < t:
        # LIS: see label {formate}
        qe_hco2 = qe_hco2 + qe_h2co2
        phi[13] = qe_hco2
        # LIS: see labela {pH_equilibria}, used as paired equations for {formate}
        dihydroxymethylidene_param_e = 1
        qe_h2co2 = ebuf(lche, pkf, ce_hco2[t], ce_h2co2[t], dihydroxymethylidene_param_e)
        phi[14] = qe_h2co2
        # cell: see label {formate}
        qi_hco2 = qi_hco2 + qi_h2co2
        phi[29] = qi_hco2
        # LIS: see labela {pH_equilibria}, used as paired equations for {formate}
        dihydroxymethylidene_param_i = 1
        qi_h2co2 = ebuf(lchi, pkf, ci_hco2[t], ci_h2co2[t], dihydroxymethylidene_param_i)
        phi[30] = qi_h2co2

    else:
        phi[13] = qe_hco2
        phi[14] = qe_h2co2
        phi[29] = qi_hco2
        phi[30] = qi_h2co2

    # cell buffer content and ph equilibrium
    # cell: see label {tbuf_volume}
    qi_cb = cbuf[t] + hcbuf[t] - (tbuf * clvl0 / clvl[t])
    phi[32] = qi_cb
    # cell: see label {ph_H2PO4_Hbuf}
    c_eQ_h2po4 = math.log10(abs(
        (ci_hpo4[t] * hcbuf[t]) / (ci_h2po4[t] * cbuf[t]))) if ci_h2po4[t] * cbuf[t] == 0 or (
            (ci_hpo4[t] * hcbuf[t]) / (ci_h2po4[t] * cbuf[t]) <= 0) else math.log10(
        (ci_hpo4[t] * hcbuf[t]) / (ci_h2po4[t] * cbuf[t]))
    phi_ph_eQ = pkb - pkp - c_eQ_h2po4
    phi[33] = phi_ph_eQ

    # cell buffer content and ph equilibrium
    # see label {electrical_Iin}
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


# # sol(1)='  na '
# # sol(2)='  k  '
# # sol(3)='  cl '
# # sol(4)=' hco3'
# # sol(5)='h2co3'
# # sol(6)=' co2 '
# # sol(7)=' hpo4'
# # sol(8)='h2po4'
# # sol(9)=' urea'
# # sol(10)=' nh3 '
# # sol(11)=' nh4 '
# # sol(12)='  h  '
# # sol(13)=' hco2'
# # sol(14)='h2co2'
# # sol(15)=' gluc'
Nfac = 1
t0 = 0
tf = 5600
T = int(56000*Nfac)
dt = float(tf - t0) / float(T)
print(dt)

rtau = 1./dt
tag = 'clever' # clever/arbitrary
ve = matrix(-0.894e-02, T, tag)
pe = matrix(-0.23e+02, T, tag)
ce_na = matrix(0.14e+00, T, tag)
ce_k = matrix(0.465e-02, T, tag)
ce_cl = matrix(0.112e+00, T, tag)
ce_hco3 = matrix(0.256e-01, T, tag)
ce_h2co3 = matrix(0.437e-05, T, tag)
ce_co2 = matrix(0.149e-02, T, tag)
ce_hpo4 = matrix(0.298e-02, T, tag)
ce_h2po4 = matrix(0.866e-03, T, tag)
ce_urea = matrix(0.490e-02, T, tag)
ce_nh3 = matrix(0.269e-05, T, tag)
ce_nh4 = matrix(0.174e-03, T, tag)
ce_hco2 = matrix(0.775e-03, T, tag)
ce_h2co2 = matrix(0.205e-06, T, tag)
ce_gluc = matrix(0.772e-02, T, tag)
vi = matrix(-0.549e+02, T, tag)
p_i = matrix(0.669e-01, T, tag)
ci_na = matrix(0.204e-01, T, tag)
ci_k = matrix(0.137e+00, T, tag)
ci_cl = matrix(0.169e-01, T, tag)
ci_hco3 = matrix(0.250e-01, T, tag)
ci_h2co3 = matrix(0.437e-05, T, tag)
ci_co2 = matrix(0.149e-02, T, tag)
ci_hpo4 = matrix(0.942e-02, T, tag)
ci_h2po4 = matrix(0.279e-02, T, tag)
ci_urea = matrix(0.495e-02, T, tag)
ci_nh3 = matrix(0.352e-05, T, tag)
ci_nh4 = matrix(0.233e-03, T, tag)
ci_hco2 = matrix(0.535e-03, T, tag)
ci_h2co2 = matrix(0.933e-07, T, tag)
ci_gluc = matrix(0.148e-01, T, tag)
hcbuf = matrix(0.269e-01, T, tag)
cbuf = matrix(0.40e-01, T, tag)
vm = matrix(-0.185e+00, T, tag)
ce_h = matrix(4.59e-11, T, tag)
ci_h = matrix(4.69e-11, T, tag)
aes = matrix(0.02000, T, tag)
chvl = matrix(0.7000e-04, T, tag)
clvl = matrix(0.1000e-02, T, tag)
l = matrix(0.1000e-02, T, tag)
imp = matrix(0.6000e-01, T, tag)
rm = matrix(0.1250e-02, T, tag)
am = matrix(0, T, tag)
phi = matrix(0, 35, tag)
lchm = pkc + math.log10(cm_hco3 / cm_h2co3)
lchs = pkc + math.log10(cs_hco3 / cs_h2co3)
cm_h = 10. ** (-lchm)
cs_h = 10. ** (-lchs)
ps = 9.00
vs = 0.00


vars =[ve, pe, ce_na, ce_k,ce_cl, ce_hco3, ce_h2co3, ce_co2, ce_hpo4, ce_h2po4, ce_urea, ce_nh3, ce_nh4, ce_hco2,
       ce_h2co2, ce_gluc,
       vi, imp, ci_na, ci_k, ci_cl, ci_hco3, ci_h2co3, ci_co2, ci_hpo4, ci_h2po4, ci_urea, ci_nh3, ci_nh4, ci_hco2,
       ci_h2co2, ci_gluc, cbuf, hcbuf, vm]

t = 0
guess =[var[t] for var in vars]
# The finite difference method was applied, stephenson (1978)
# Analysis of the transient behavior of the kidney model
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

print(np.array(guess))

while True:
    t += 1
    print('t = % 7.1f' % (t))
    if t == T:
        break
    else:
        cm_na = 0.11
        result = scipy.optimize.root(eQs, np.array(guess), 1)
        guess =[var[t] for var in vars]
#        print(t)
#        if t == 100:
#            print(result.x)
#            break
#            quit()
sizeX = 30
sizeY = 15
plot_concen=1
t = np.linspace(t0, tf, T)
pltfolder = 'New_plots_dt' + str(Nfac)
if not os.path.exists(pltfolder):
    os.makedirs(pltfolder)
if plot_concen:
    figA, ax = plt.subplots(7, 1, figsize=(sizeX, sizeY), constrained_layout=True)
    ax [ 0 ].plot(t, ci_na, 'b-')
    ax [ 0 ].set_ylabel(' ci_na ', color='blue')
    ax [ 1 ].plot(t, ci_k, 'b-')
    ax [ 1 ].set_ylabel('ci_k ', color='blue')
    ax [ 2 ].plot(t, ci_cl, 'b-')
    ax [ 2 ].set_ylabel('ci_cl ', color='blue')
    ax [ 3 ].plot(t, ci_urea, 'b-')
    ax [ 3 ].set_ylabel('ci_urea ', color='blue')
    ax [ 4 ].plot(t, ci_nh3, 'b-')
    ax [ 4 ].set_ylabel('ci_nh3', color='blue')
    ax [ 5 ].plot(t, ci_nh4, 'b-')
    ax [ 5 ].set_ylabel('ci_nh4', color='blue')
    ax [ 6 ].plot(t, ci_h, 'b-')
    ax [ 6 ].set_ylabel('ci_h ', color='blue')
    ax [ 6 ].set_xlabel('time (s) ', color='blue')
    figA.savefig(pltfolder + '/figA' + tag + '.png')
    plt.close(figA)

    figB, ax = plt.subplots(7, 1, figsize=(sizeX, sizeY), constrained_layout=True)
    ax [ 0 ].plot(t, ci_gluc, 'b-')
    ax [ 0 ].set_ylabel('ci_gluc ', color='blue')
    ax [ 0 ].set_title('Solute Concentration OF THE CELL[mmol] ')
    ax [ 1 ].plot(t, ci_hco3, 'b-')
    ax [ 1 ].set_ylabel('ci_hco3 ', color='blue')
    ax [ 2 ].plot(t, ci_h2co3, 'b-')
    ax [ 2 ].set_ylabel('ci_h2co3 ', color='blue')
    ax [ 3 ].plot(t, ci_hpo4, 'b-')
    ax [ 3 ].set_ylabel('ci_hpo4 ', color='blue')
    ax [ 4 ].plot(t, ci_h2po4, 'b-')
    ax [ 4 ].set_ylabel('ci_h2po4', color='blue')
    ax [ 5 ].plot(t, ci_hco2, 'b-')
    ax [ 5 ].set_ylabel('ci_hco2', color='blue')
    ax [ 6 ].plot(t, ci_h2co2, 'b-')
    ax [ 6 ].set_ylabel('ci_h2co2 ', color='blue')
    ax [ 6 ].set_xlabel('time (s) ', color='blue')
    figB.savefig(pltfolder + '/figB' + tag + '.png')
    plt.close(figB)

    figC, ax = plt.subplots(7, 1, figsize=(sizeX, sizeY), constrained_layout=True)
    ax [ 0 ].plot(t, ce_na, 'b-')
    ax [ 0 ].set_title('Solute Concentration OF THE CHANNEL[mmol]')
    ax [ 0 ].set_ylabel(' ce_na ', color='blue')
    ax [ 1 ].plot(t, ce_k, 'b-')
    ax [ 1 ].set_ylabel('ce_k ', color='blue')
    ax [ 2 ].plot(t, ce_cl, 'b-')
    ax [ 2 ].set_ylabel('ce_cl ', color='blue')
    ax [ 3 ].plot(t, ce_urea, 'b-')
    ax [ 3 ].set_ylabel('ce_urea ', color='blue')
    ax [ 4 ].plot(t, ce_nh3, 'b-')
    ax [ 4 ].set_ylabel('ce_nh3', color='blue')
    ax [ 5 ].plot(t, ce_nh4, 'b-')
    ax [ 5 ].set_ylabel('ce_nh4', color='blue')
    ax [ 6 ].plot(t, ce_h, 'b-')
    ax [ 6 ].set_ylabel('ce_h ', color='blue')
    ax [ 6 ].set_xlabel('time (s) ', color='blue')
    figC.savefig(pltfolder + '/figC' + tag + '.png')
    # plt.close(figC)
    plt.show()
    figD, ax = plt.subplots(7, 1, figsize=(sizeX, sizeY), constrained_layout=True)
    ax [ 0 ].plot(t, ce_gluc, 'b-')
    ax [ 0 ].set_ylabel('ce_gluc ', color='blue')
    ax [ 0 ].set_title('Solute Concentration OF THE CELL[mmol] ')
    ax [ 1 ].plot(t, ce_hco3, 'b-')
    ax [ 1 ].set_ylabel('ce_hco3 ', color='blue')
    ax [ 2 ].plot(t, ce_h2co3, 'b-')
    ax [ 2 ].set_ylabel('ce_h2co3 ', color='blue')
    ax [ 3 ].plot(t, ce_hpo4, 'b-')
    ax [ 3 ].set_ylabel('ce_hpo4 ', color='blue')
    ax [ 4 ].plot(t, ce_h2po4, 'b-')
    ax [ 4 ].set_ylabel('ce_h2po4', color='blue')
    ax [ 5 ].plot(t, ce_hco2, 'b-')
    ax [ 5 ].set_ylabel('ce_hco2', color='blue')
    ax [ 6 ].plot(t, ce_h2co2, 'b-')
    ax [ 6 ].set_ylabel('ce_h2co2 ', color='blue')
    ax [ 6 ].set_xlabel('time (s) ', color='blue')
    figD.savefig(pltfolder + '/figD' + tag + '.png')
    # plt.close(figD)
    plt.show()
    figE, ax = plt.subplots(6, 1, figsize=(sizeX, sizeY), constrained_layout=True)
    ax [ 0 ].plot(t, imp, 'r-', linewidth=2.0)
    ax [ 0 ].set_xlim(-5, tf)
    ax [ 0 ].set_ylabel('imp(mm) ', color='blue')
    ax [ 2 ].plot(t, pe, 'b-', linewidth=2.0)
    ax [ 2 ].set_ylabel('pressure(mm.Hg) ', color='blue')
    ax [ 2 ].set_xlabel('time (s) ', color='blue')
    ax [ 2 ].set_xlim(-5, tf)
    ax [ 1 ].set_title('MEMBRANE POTENTIAL[mVolt]')
    vem = np.asarray(ve) - np.asarray(vm)
    v_em = np.asarray(vem)
    ax [ 1 ].plot(t, v_em, 'b-', linewidth=2.0)
    ax [ 1 ].set_ylabel('voltage_em (mv) ', color='blue')
    ax [ 1 ].set_xlabel('time (s) ', color='blue')
    ax [ 1 ].set_xlim(-5, tf)
    vie = np.asarray(vi) - np.asarray(ve)
    v_ie = np.asarray(vie)
    ax [ 3 ].plot(t, v_ie, 'b-', linewidth=2.0)
    ax [ 3 ].set_ylabel('voltage_ie (mv) ', color='blue')
    ax [ 3 ].set_xlabel('time (s) ', color='blue')
    ax [ 3 ].set_xlim(-5, tf)
    vim = np.asarray(vi) - np.asarray(vm)
    v_im = np.asarray(vim)
    ax [ 4 ].plot(t, v_im, 'b-', linewidth=2.0)
    ax [ 4 ].set_ylabel('voltage_im (mv) ', color='blue')
    ax [ 4 ].set_xlabel('time (s) ', color='blue')
    ax [ 4 ].set_xlim(-5, tf)
    vis = np.asarray(vi) - np.asarray(vs)
    v_is = np.asarray(vis)
    ax [ 5 ].plot(t, v_is, 'b-', linewidth=2.0)
    ax [ 5 ].set_ylabel('voltage_is (mv) ', color='blue')
    ax [ 5 ].set_xlabel('time (s) ', color='blue')
    ax [ 5 ].set_xlim(-5, tf)
    figE.savefig(pltfolder + '/figE' + tag + '.png')
    # plt.close(figE)
    plt.show()