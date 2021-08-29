import numpy as np
import scipy
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
# from sympy import *
import math
from PCT_GLOB_new import *

# ELECTROLYTE  TRANSPORT   ACROSS  A SIMPLE EPITHELIUM STEADY - STATE AND  TRANSIENT ANALYSIS EQ.8
# AREA    ADJUSTED    PARTIAL    CONDUCTANCE(MHO / CM2)
# must be chacked, this equation is not consistant with other articles
def AAPC(ca, cb, z, a, hab, va, vb, param_AAPC):
    zab = (z * (va - vb)*f * 1.e-6) / rte
    if param_AAPC == 0:
        return 0
    # def lmmsc(ca, cb):
    #     import math
    #     print('ca',ca)
    #     print('cb', cb)
    #     if ca > 0 and cb > 0 and ca - cb != 0 and cb != 0 and math.log10(ca / cb) != 0:
    #         return (ca - cb) / (math.log10(ca / cb))
    #     else:
    #         return cb
    #def goldman_conduct(hab, a, z, va, vb, ca, cb, param_goldman):
    elif zab == 0 or va == vb and ca > 0 and cb > 0:
        return hab * a * (ca - cb)
    elif zab > 0:
        return hab * a * zab * (ca - cb * math.exp(-zab)) / (1 - math.exp(-zab))
    else:
        return hab * a * zab * (ca * math.exp(zab) - cb) / (math.exp(zab) - 1)

def sglt_mi(cm_na, ci_na, cm_gluc, ci_gluc, z_na, z_gluc, vm, vi, ami, lmi_nagluc, param_sglt_mi):
    # Na-glucose simple cotransporter with 1:1 stoichiometry, located on  Apical  Membrane
    # return is the transported flux for each solute due to the existance of Na-glucose simple cotransporter.
    # f_eps(): Modular function to calculate the electrochemical potential of species:(RT)*ln(c) + z*f*V
    # see label: {eq:NAGluc}
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


def clhco3_mi(cm_cl, ci_cl, cm_hco3, ci_hco3, z_cl, z_hco3, vm, vi, ami, lmi_clhco3, param_clhco3_mi):
    # cl/hco3 simple exchanger with 1:-1 stoichiometry, located on Apical  Membrane.
    # see label {eq: CLHCO3}
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


def clhco2_mi(cm_cl, ci_cl, cm_hco2, ci_hco2, z_cl, z_hco2, vm, vi, ami, lmi_clhco2, param_clhco2_mi):
    # cl/hco2 simple exchanger with 1:-1 stoichiometry, located on Apical  Membrane.
    # see label {eq:CLHCO2}
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

# proton pumps #checked
def h_atp_mi(cm_h, ci_h, vm, vi, z_h, ami, parama_mi_h):
    # Proton pumps area special kind  of transporter  that push
    # hydrogen ions  from areas of low concentration
    # with high concentration. Ions moving down a gradient release energy,
    # but when they move up a gradient, it takes energy.Diffusion can then use this
    # gradient to capture energy again, as the ions move downhill.
    xm_h = f_eps(cm_h, z_h, vm)
    xi_h = f_eps(ci_h, z_h, vi)
    gamma = -(xihp * (xm_h - xi_h - xhp))
    if parama_mi_h == 0:
        return 0
    elif gamma < 0:
        return -lhp * ami * (1 - 1 / (1 + np.exp(gamma)))
    return -lhp * ami * (1 / (1 + np.exp(-gamma)))


def na_hco3(ci_na, cs_na, ci_hco3, cs_hco3, z_na, z_hco3, vi, vs, ais, lis_nahco3, param_na_hco3):
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


def k_cl(ca_k, cb_k, ca_cl, cb_cl, z_k, z_cl, va, vb, a, l_kcl, param_kcl):
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

def nah(ci_h, ci_na, ci_nh4, cm_h, cm_na, cm_nh4, param_nah):
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
        # the nah model parameters
        cxt = 1.00
        # knah_na = 0.3000e-01
        # knah_h = 0.7200e-07
        # knah_nh4 = 0.2700e-01
        # knah_i = 0.1000e-05
        #
        # pnah_na = 0.8000e+04
        # pnah_h = 0.2400e+04
        # pnah_nh4 = 0.8000e+04
        # pnah_m = 0.2000e+01
        # pnah_mm = 0.0000

        knah_na = 0.3000e-01
        knah_h = 0.7200e-07
        knah_nh4 = 0.2700e-01
        knah_i = 0.1000e-05
        #
        pnah_na = 0.792000e+04
        pnah_h = 0.23800e+04
        pnah_nh4 = 0.792000e+04
        pnah_mm = 0.2000e+01
        pnah_m = 0.0000

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
        # jnah_na_max = cxt * psnah_na * psnah_h / (psnah_na + psnah_h)
        return[jnah_na, jnah_h, jnah_nh4]



def nak_atp( ci_na, ce_na,ci_k,  ce_k, ce_nh4, param_nak_atp):
    #powered by ATP, the pump moves
    # sodium and potassium  ions in opposite
    # directions, each against its concentration
    # gradient.In  a single cycle of the pump, three sodium ions are
    # extruded from and two  potassium ions are imported into  the cell.

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

# goldman fluxes (passive fluxes) describes the ionic flux across a cell membrane
# as a function of the transmembrane potential and the concentrations of the ion inside and outside of the cell.
# checked
# goldman fluxes (passive fluxes) describes the ionic flux across a cell membrane
# as a function of the transmembrane potential and the concentrations of the ion inside and outside of the cell.
# checked


def goldman(hab, a, z, va, vb, ca, cb, param_goldman):
    zab = (z * f * (va - vb) * 1.e-6) / rte
    if param_goldman == 0:
        return [ 0 ]
    elif zab == 0 or va == vb and ca > 0 and cb > 0:
        return hab * a * (ca - cb)
    elif zab > 0:
        return hab * a * zab * (ca - cb * np.exp(-zab)) / (1 - np.exp(-zab))
    else:
        return hab * a * zab * (ca * np.exp(zab) - cb) / (np.exp(zab) - 1)


# Convective Solute Fluxes
def csf(ca, cb, flux, s, param_csf):
    if param_csf == 0:
        return 0

    # Definition of the mean membrane solute concentration #checked
    def lmmsc(ca, cb):

        import math
        if ca > 0 and cb > 0 and ca - cb != 0 and cb != 0 and np.log10(ca / cb) != 0:
            return (ca - cb) / (np.log10(ca / cb))
        else:
            return cb

    return flux * (1.00 - s) * lmmsc(ca, cb)

def lch(ca, cb):
    import math
    if ca > 0 and cb > 0 and (ca / cb) != 0 and cb != 0:
        return np.log10(ca / cb)
    else:
        return np.log10(abs(ca / cb))


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
        return rte * np.log(c) + z * f * v * 1.e-6
    else:
        print('uy')

        return rte * np.log(abs(c)) + z * f * v * 1.e-6





# pH equilibria of four buffer pairs
def ebuf(lch, pk, ca, cb, param_ebuf):
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

           #q_h2_vary = lc_h - pk - np.log10(c_h_vary / c_h2_vary)
        else:
           q_h_vary = q_h_vary + q_h2_vary
           # print('c_h_vary', c_h_vary)
           # print('c_h2_vary', c_h2_vary)
           q_h2_vary = lc_h - pk - lch(c_h_vary, c_h2_vary)
    else:
        #print('No Buffer Effect for phosphate')
        pass
    return q_h_vary, q_h2_vary

# common definition
def zab(z, va, vb):
    return z * f * (va - vb) * 1.e-6 / rte


def phi_scale(phi, scale):
    return phi * scale


solver = 1


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

    # The cell height (summation of Extracellular channel volume
    # and Intracellular compartment volume[cm3/cm2 epithelium])
    l [ t ] = l [ t ] + clvl [ t ]
    p_i = pm

    # LOG Conc. of Hydrogen Cell
    lchi = pkc + lch(ci_hco3 [ t ], ci_h2co3 [ t ])
    ci_h [ t ] = 10 ** (-lchi)
    # LOG Conc. of Hydrogen Cell
    # lchi = pkp + lch(ci_hpo4 [ t ], ci_h2po4 [ t ])
    # ci_h [ t ] = 10 ** (-lchi)

    # lchi = pkp + lch(ci_hpo4 [ t ], ci_h2po4 [ t ])
    # ci_h [ t ] = 10 ** (-lchi)

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

    # Convective Intraepithelial Solute Fluxes
    param_csf = 1
    fekm_na = csf(ce_na [ t ], cm_na, fevm, sme_na, param_csf)
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

    # goldman  fluxes
    param_goldman = 1
    fekm_na = fekm_na + goldman(hme_na, ame, z_na, vm [ t ], ve [ t ], cm_na,
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
    param_sglt = 1
    sglt = sglt_mi(cm_na, ci_na [ t ], cm_gluc, ci_gluc [ t ], z_na, z_gluc, vm [ t ], vi [ t ], ami, lmi_nagluc,
                   param_sglt)
    na_mi_nagluc = sglt [ 0 ]
    gluc_mi_nagluc = sglt [ 1 ]
    param_h2po4 = 1
    nah2po4 = nah2po4_mi(cm_na, ci_na [ t ], cm_h2po4, ci_h2po4 [ t ], z_na, z_h2po4, vm [ t ], vi [ t ], ami,
                         lmi_nah2po4, param_h2po4)
    na_mi_nah2po4 = nah2po4 [ 0 ]
    h2po4_mi_nah2po4 = nah2po4 [ 1 ]
    param_clhco3 = 1
    clhco3 = clhco3_mi(cm_cl, ci_cl [ t ], cm_hco3, ci_hco3 [ t ], z_cl, z_hco3, vm [ t ], vi [ t ], ami,
                       lmi_clhco3, param_clhco3)
    cl_mi_clhco3 = clhco3 [ 0 ]
    hco3_mi_clhco3 = clhco3 [ 1 ]

    clhco2 = clhco2_mi(cm_cl, ci_cl [ t ], cm_hco2, ci_hco2 [ t ], z_cl, z_hco2, vm [ t ], vi [ t ], ami,
                       lmi_clhco2, 1)
    cl_mi_clhco2 = clhco2 [ 0 ]
    hco2_mi_clhco2 = clhco2 [ 1 ]

    # net cotransporters on is bourder
    param_nahco3 = 1
    nahco3 = na_hco3(ci_na [ t ], cs_na, ci_hco3 [ t ], cs_hco3, z_na, z_hco3, vi [ t ], vs, ais,
                       lis_nahco3, param_nahco3)
    na_is_nahco3 = nahco3 [ 0 ]
    hco3_is_nahco3 = nahco3 [ 1 ]
    param_kcl = 1
    kcl = k_cl(ci_k [ t ], cs_k, ci_cl [ t ], cs_cl, z_k, z_cl, vi [ t ], vs, ais, lis_kcl,  param_kcl)
    k_is_kcl = kcl [ 0 ]
    cl_is_kcl = kcl [ 1 ]

    na_clhco3 = na_cl_hco3(ci_na [ t ], cs_na, ci_cl [ t ], cs_cl, ci_hco3 [ t ], cs_hco3, z_na, z_cl, z_hco3,
                             vi [ t ], vs, ais, lis_na_clhco3, 1)
    na_is_na_clhco3 = na_clhco3 [ 0 ]
    cl_is_na_clhco3 = na_clhco3 [ 1 ]
    hco3_is_na_clhco3 = na_clhco3 [ 2 ]
    # the nah exchanger translate concentrations to the nah model on  mi bourder
    param_nhe3 = 1

    mynah = nah(ci_h [ t ], ci_na [ t ], ci_nh4 [ t ], cm_h, cm_na, cm_nh4, param_nhe3)
    jnah_na = mynah [ 0 ]
    jnah_h = mynah [ 1 ]
    jnah_nh4 = mynah [ 2 ]
    jnhe3_na = nnhe3 * ami * jnah_na
    jnhe3_h = nnhe3 * ami * jnah_h
    jnhe3_nh4 = nnhe3 * ami * jnah_nh4

    fikm_na = fikm_na + na_mi_nagluc + na_mi_nah2po4 + jnhe3_na
    fikm_cl = fikm_cl + cl_mi_clhco2 + cl_mi_clhco3
    fikm_hco3 = fikm_hco3 + hco3_mi_clhco3
    fikm_h2po4 = fikm_h2po4 + h2po4_mi_nah2po4
    fikm_hco2 = fikm_hco2 + hco2_mi_clhco2
    fikm_gluc = fikm_gluc + gluc_mi_nagluc
    fikm_h = fikm_h + jnhe3_h
    fikm_nh4 = fikm_nh4 + jnhe3_nh4

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
                                 vi [ t ], ve [ t ], aie, lis_na_clhco3, 1)
    na_aie_na_clhco3 = na_clhco3_aie [ 0 ]
    cl_aie_na_clhco3 = na_clhco3_aie [ 1 ]
    hco3_aie_na_clhco3 = na_clhco3_aie [ 2 ]

    jk_na = jk_na + na_aie_nahco3 + na_aie_na_clhco3
    jk_k = jk_k + k_aie_kcl
    jk_cl = jk_cl + cl_aie_kcl + cl_aie_na_clhco3
    jk_hco3 = jk_hco3 + hco3_aie_na_clhco3 + hco3_aie_nahco3

    # sodium pumps on is bourder
    param_sodium_pumps = 1
    nak = nak_atp(ci_na [ t ], cs_na, ci_k [ t ], cs_k,  cs_nh4, param_sodium_pumps)
    atis_na = nak [ 0 ]
    atis_k = nak [ 1 ]
    atis_nh4 = nak [ 2 ]
    # sodium pumps on ie bourder

    nak = nak_atp(ci_na [ t ], ce_na[t], ci_k [ t ], ce_k[t],  ce_nh4[t], param_sodium_pumps)
    atie_na = nak [ 0 ]
    atie_k = nak [ 1 ]
    atie_nh4 = nak [ 2 ]
    # proton pumps
    # proton pumps
    param_H_pumps = 1
    atmi_h = h_atp_mi(cm_h, ci_h [ t ], vm [ t ], vi [ t ], z_h, ami, param_H_pumps)

    jk_na = jk_na + aie * atie_na
    jk_k = jk_k + aie * atie_k
    jk_nh4 = jk_nh4 + aie * atie_nh4

    fiks_na = fiks_na + ais * atis_na
    fiks_k = fiks_k + ais * atis_k
    fiks_nh4 = fiks_nh4 + ais * atis_nh4
    # proton pumps
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

        # co2, formate, phosphate, and ammonia content:
        BuffPairs_progressive_activation_alongtime = 1
        # dihydrogen_phosphate_param_e = 1
        # dihydrogen_phosphate_param_i = 1
        # dihydrogen_ammonium_param_e = 1
        # dihydrogen_ammonium_param_i = 1
        # dihydroxymethylidene_param_e = 1
        # dihydroxymethylidene_param_i = 1
        n = 3
        # print('buffer', BuffPairs_Progres_Activ(BuffPairs_progressive_activation_alongtime, t, T, T_window, n))
        # Buffer_Co2_LIS_Param, Buffer_Co2_Cell_Param, Buffer_HPO4_LIS_Param, \
        # Buffer_HPO4_CELL_Param, Buffer_NH3_LIS_Param, \
        # Buffer_NH3_CELL_Param, Buffer_HCO2_LIS_Param, \
        # Buffer_HCO2_CELL_Param = BuffPairs_Progres_Activ(BuffPairs_progressive_activation_alongtime, t, T, T_window, n)
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
        # c_eQ_h2co3 = np.log10(abs(
        #     (ci_hco3 [ t ] * hcbuf [ t ]) / (ci_h2co3 [ t ] * cbuf [ t ]))) if ci_h2co3 [ t ] * cbuf [ t ] == 0 or (
        #         (ci_hco3 [ t ] * hcbuf [ t ]) / (ci_h2co3 [ t ] * cbuf [ t ]) <= 0) else np.log10(
        #     (ci_hco3 [ t ] * hcbuf [ t ]) / (ci_h2co3 [ t ] * cbuf [ t ]))
        # phi_ph_eQ = pkb - pkp - c_eQ_h2co3

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

        cures = f * (z_na * feks_na + z_k * feks_k + z_cl * feks_cl + z_hco3 * feks_hco3 + z_h2co3 * feks_h2co3 + z_co2
                     * feks_co2 + z_hpo4 * feks_hpo4 + z_h2po4 * feks_h2po4 + z_urea * feks_urea + z_nh3 * feks_nh3 + z_nh4
                     * feks_nh4 + z_h * feks_h + z_hco2 * feks_hco2 + z_h2co2 * feks_h2co2 + z_gluc * feks_gluc)

        curie = f * (z_na * jk_na + z_k * jk_k + z_cl * jk_cl + z_hco3 * jk_hco3 + z_h2co3 * jk_h2co3 + z_co2 * jk_co2 +
                     z_hpo4 * jk_hpo4 + z_h2po4 * jk_h2po4 + z_urea * jk_urea + z_nh3 * jk_nh3 + z_nh4 * jk_nh4 +
                     z_h * jk_h + z_hco2 * jk_hco2 + z_h2co2 * jk_h2co2 + z_gluc * jk_gluc)
        curis = f * (
                z_na * fiks_na + z_k * fiks_k + z_cl * fiks_cl + z_hco3 * fiks_hco3 + z_h2co3 * fiks_h2co3 + z_co2 * fiks_co2 + \
                z_hpo4 * fiks_hpo4 + z_h2po4 * fiks_h2po4 + z_urea * fiks_urea + z_nh3 * fiks_nh3 + z_nh4 * fiks_nh4 + \
                z_h * fiks_h + z_hco2 * fiks_hco2 + z_h2co2 * fiks_h2co2 + z_gluc * fiks_gluc)

        # flux - TOTAL EPITHELIAL SOLUTE FLUX
        epithelial_flux = 1
        if epithelial_flux:
            flux_na = fikm_na + fekm_na
            flux_k = fikm_k + fekm_k
            flux_cl = fikm_cl + fekm_cl
            flux_hco3 = fikm_hco3 + fekm_hco3
            flux_h2co3 = fikm_h2co3 + fekm_h2co3
            flux_co2 = fikm_co2 + fekm_co2
            flux_hpo4 = fikm_hpo4 + fekm_hpo4
            flux_h2po4 = fikm_h2po4 + fekm_h2po4
            flux_urea = fikm_urea + fekm_urea
            flux_nh3 = fikm_nh3 + fekm_nh3
            flux_nh4 = fikm_nh4 + fekm_nh4
            flux_h = fikm_h + fekm_h
            flux_hco2 = fikm_hco2 + fekm_hco2
            flux_h2co2 = fikm_h2co2 + fekm_h2co2
            flux_gluc = fikm_gluc + fekm_gluc
            flux_v = fivm + fevm
            # flux_isj = fivs + jv

            osmilarity = 1
            if osmilarity:
                osme = ce_na [ t ] + ce_k [ t ] + ce_cl [ t ] + ce_hco3 [ t ] + ce_h2co3 [ t ] + ce_co2 [ t ] + \
                       ce_hpo4 [ t ] + \
                       ce_h2po4 [ t ] + ce_urea [ t ] + ce_nh3 [ t ] + ce_nh4 [ t ] + ce_h [ t ] + ce_hco2 [ t ] + \
                       ce_h2co2 [ t ] + ce_gluc [ t ]
                osmi = imp [ t ] + ci_na [ t ] + ci_k [ t ] + ci_cl [ t ] + ci_hco3 [ t ] + ci_h2co3 [ t ] + ci_co2 [
                    t ] + ci_hpo4 [ t ] + \
                       ci_h2po4 [ t ] + ci_urea [ t ] + ci_nh3 [ t ] + ci_nh4 [ t ] + ci_h [ t ] + ci_hco2 [ t ] + \
                       ci_h2co2 [ t ] + ci_gluc [ t ]
                osmm = impm + cm_na + cm_k + cm_cl + cm_hco3 + cm_h2co3 + cm_co2 + cm_hpo4 + cm_h2po4 + cm_urea + cm_nh3 + cm_nh4 \
                       + cm_h + cm_hco2 + cm_h2co2 + cm_gluc
                osms = imps + cs_na + cs_k + cs_cl + cs_hco3 + cs_h2co3 + cs_co2 + cs_hpo4 + cs_h2po4 + cs_urea + cs_nh3 + cs_nh4 \
                       + cs_h + cs_hco2 + cs_h2co2 + cs_gluc

            transporterson = 1
            if transporterson:
                #  net transporters on mi bourder:[sglt_mi, nah2po4_mi,clhco3_mi]
                #  net cotransporters on is bourder:[na_hco3,k_cl,na_cl_hco3]
                # the nah exchanger translate concentrations to the nah model on  mi bourder

                # cikm-COTRANSPORT OF SOLUTE ACROSS APICAL MEMBRANE
                cikm_na = na_mi_nagluc + na_mi_nah2po4 + jnhe3_na
                cikm_cl = cl_mi_clhco2 + cl_mi_clhco3
                cikm_hco3 = hco3_mi_clhco3
                cikm_h2po4 = h2po4_mi_nah2po4
                cikm_hco2 = hco2_mi_clhco2
                cikm_gluc = gluc_mi_nagluc
                cikm_h = jnhe3_h
                cikm_nh4 = jnhe3_nh4
                #  ciks-COTRANSPORT OF SOLUTE ACROSS BASAL MEMBRANE
                ciks_na = na_is_nahco3 + na_is_na_clhco3
                ciks_k = k_is_kcl
                ciks_cl = cl_is_kcl + cl_is_na_clhco3
                ciks_hco3 = hco3_is_na_clhco3 + hco3_is_nahco3

                # jk-COTRANSPORT OF SOLUTE ACROSS LATERAL MEMBRANE=ike
                cjk_na = na_aie_nahco3 + na_aie_na_clhco3
                cjk_k = k_aie_kcl
                cjk_cl = cl_aie_kcl + cl_aie_na_clhco3
                cjk_hco3 = hco3_aie_na_clhco3 + hco3_aie_nahco3

                # pjk-ACTIVE TRANSPORT OF SOLUTE ACROSS LATERAL MEMBRANE
                # sodium pumps on mi bourder:nak_atp
                # nak = (ci_k [ t ], cs_k, ci_na [ t ], cs_na, ce_k [ t ], ce_nh4 [ t ], 1)
                # proton pumps:h_atp_mi
                atmi_h = h_atp_mi(cm_h, ci_h [ t ], vm [ t ], vi [ t ], z_h, ami, 1)
                # pikm-ACTIVE TRANSPORT OF SOLUTE ACROSS APICAL MEMBRANE
                pikm_h = + atmi_h
                #	piks- ACTIVE TRANSPORT OF SOLUTE ACROSS BASAL MEMBRANE
                # sodium pumps on is bourder:nak_atp
                piks_na = + ais * atis_na
                piks_k = + ais * atis_k
                piks_nh4 = ais * atis_nh4
                #	pjk- ACTIVE TRANSPORT OF SOLUTE ACROSS LATERAL MEMBRANE
                pjk_na = + aie * atie_na
                pjk_k = + aie * atie_k
                pjk_nh4 = + aie * atie_nh4

            conductance = 1
            if conductance:
                param_AAPC = 1
                gcem_na = AAPC(cm_na, ce_na [ t ], z_na, ame, hme_na, vm [ t ], ve [ t ], param_AAPC)

                gcem_k = AAPC(cm_k, ce_k [ t ], z_k, ame, hme_k, vm [ t ], ve [ t ], param_AAPC)
                gcem_cl = AAPC(cm_cl, ce_cl [ t ], z_cl, ame, hme_cl, vm [ t ], ve [ t ], param_AAPC)
                gcem_hco3 = AAPC(cm_hco3, ce_hco3 [ t ], z_hco3, ame, hme_hco3, vm [ t ], ve [ t ], param_AAPC)
                gcem_h2co3 = AAPC(cm_h2co3, ce_h2co3 [ t ], z_h2co3, ame, hme_h2co3, vm [ t ], ve [ t ],
                                  param_AAPC)
                gcem_co2 = AAPC(cm_co2, ce_co2 [ t ], z_co2, ame, hme_co2, vm [ t ], ve [ t ], param_AAPC)
                gcem_hpo4 = AAPC(cm_hpo4, ce_hpo4 [ t ], z_hpo4, ame, hme_hpo4, vm [ t ], ve [ t ], param_AAPC)
                gcem_h2po4 = AAPC(cm_h2po4, ce_h2po4 [ t ], z_h2po4, ame, hme_h2po4, vm [ t ], ve [ t ], param_AAPC)
                gcem_urea = AAPC(cm_urea, ce_urea [ t ], z_urea, ame, hme_urea, vm [ t ], ve [ t ], param_AAPC)
                gcem_nh3 = AAPC(cm_nh3, ce_nh3 [ t ], z_nh3, ame, hme_nh3, vm [ t ], ve [ t ], param_AAPC)
                gcem_nh4 = AAPC(cm_nh4, ce_nh4 [ t ], z_nh4, ame, hme_nh4, vm [ t ], ve [ t ], param_AAPC)
                gcem_h = AAPC(cm_h, ce_h [ t ], z_h, ame, hme_h, vm [ t ], ve [ t ], param_AAPC)
                gcem_hco2 = AAPC(cm_hco2, ce_hco2 [ t ], z_hco2, ame, hme_hco2, vm [ t ], ve [ t ], param_AAPC)
                gcem_h2co2 = AAPC(cm_h2co2, ce_h2co2 [ t ], z_h2co2, ame, hme_h2co2, vm [ t ], ve [ t ],
                                  param_AAPC)
                gcem_gluc = AAPC(cm_gluc, ce_gluc [ t ], z_gluc, ame, hme_gluc, vm [ t ], ve [ t ], param_AAPC)

                tgcem = gcem_na + gcem_k + gcem_cl + gcem_hco3 + gcem_h2co3 + gcem_co2 + gcem_hpo4 + gcem_h2po4 \
                        + gcem_gluc + gcem_urea + gcem_nh3 + gcem_nh4 + gcem_h + gcem_hco2 + gcem_h2co2
                tgcem = abs(tgcem / (ve [ t ] - vm [ t ]))
                gces_na = AAPC(ce_na [ t ], cs_na, z_na, ae [ t ], hes_na, ve [ t ], vs, param_AAPC)
                gces_k = AAPC(ce_k [ t ], cs_k, z_k, ae [ t ], hes_k, ve [ t ], vs, param_AAPC)
                gces_cl = AAPC(ce_cl [ t ], cs_cl, z_cl, ae [ t ], hes_cl, ve [ t ], vs, param_AAPC)
                gces_hco3 = AAPC(ce_hco3 [ t ], cs_hco3, z_hco3, ae [ t ], hes_hco3, ve [ t ], vs, param_AAPC)
                gces_h2co3 = AAPC(ce_h2co3 [ t ], cs_h2co3, z_h2co3, ae [ t ], hes_h2co3, ve [ t ], vs, param_AAPC)
                gces_co2 = AAPC(ce_co2 [ t ], cs_co2, z_co2, ae [ t ], hes_co2, ve [ t ], vs, param_AAPC)
                gces_hpo4 = AAPC(ce_hpo4 [ t ], cs_hpo4, z_hpo4, ae [ t ], hes_hpo4, ve [ t ], vs, param_AAPC)
                gces_h2po4 = AAPC(ce_h2po4 [ t ], cs_h2po4, z_h2po4, ae [ t ], hes_h2po4, ve [ t ], vs, param_AAPC)
                gces_urea = AAPC(ce_urea [ t ], cs_urea, z_urea, ae [ t ], hes_urea, ve [ t ], vs, param_AAPC)
                gces_nh3 = AAPC(ce_nh3 [ t ], cs_nh3, z_nh3, ae [ t ], hes_nh3, ve [ t ], vs, param_AAPC)
                gces_nh4 = AAPC(ce_nh4 [ t ], cs_nh4, z_nh4, ae [ t ], hes_nh4, ve [ t ], vs, param_AAPC)
                gces_h = AAPC(ce_h [ t ], cs_h, z_h, ae [ t ], hes_h, ve [ t ], vs, param_AAPC)
                gces_hco2 = AAPC(ce_hco2 [ t ], cs_hco2, z_hco2, ae [ t ], hes_hco2, ve [ t ], vs, param_AAPC)
                gces_h2co2 = AAPC(ce_h2co2 [ t ], cs_h2co2, z_h2co2, ae [ t ], hes_h2co2, ve [ t ], vs, param_AAPC)
                gces_gluc = AAPC(ce_gluc [ t ], cs_gluc, z_gluc, ae [ t ], hes_gluc, ve [ t ], vs, param_AAPC)

                tgces = gces_na + gces_k + gces_cl + gces_hco3 + gces_h2co3 + gces_co2 + gces_hpo4 + gces_h2po4 + \
                        gces_urea + gces_nh3 + gces_nh4 + gces_h + gces_hco2 + gces_h2co2 + gces_gluc
                tgces = tgces / (ve [ t ])

                gcmi_na = AAPC(cm_na, ci_na [ t ], z_na, ae [ t ], hmi_na, vm [ t ], vi [ t ], param_AAPC)
                gcmi_k = AAPC(cm_k, ci_k [ t ], z_k, ae [ t ], hmi_k, vm [ t ], vi [ t ], param_AAPC)
                gcmi_cl = AAPC(cm_cl, ci_cl [ t ], z_cl, ae [ t ], hmi_cl, vm [ t ], vi [ t ], param_AAPC)
                gcmi_hco3 = AAPC(cm_hco3, ci_hco3 [ t ], z_hco3, ae [ t ], hmi_hco3, vm [ t ], vi [ t ], param_AAPC)
                gcmi_h2co3 = AAPC(cm_h2co3, ci_h2co3 [ t ], z_h2co3, ae [ t ], hmi_h2co3, vm [ t ], vi [ t ],
                                  param_AAPC)
                gcmi_co2 = AAPC(cm_co2, ci_co2 [ t ], z_co2, ae [ t ], hmi_co2, vm [ t ], vi [ t ], param_AAPC)
                gcmi_hpo4 = AAPC(cm_hpo4, ci_hpo4 [ t ], z_hpo4, ae [ t ], hmi_hpo4, vm [ t ], vi [ t ], param_AAPC)
                gcmi_h2po4 = AAPC(cm_h2po4, ci_h2po4 [ t ], z_h2po4, ae [ t ], hmi_h2po4, vm [ t ], vi [ t ],
                                  param_AAPC)
                gcmi_urea = AAPC(cm_urea, ci_urea [ t ], z_urea, ae [ t ], hmi_urea, vm [ t ], vi [ t ], param_AAPC)
                gcmi_nh3 = AAPC(cm_nh3, ci_nh3 [ t ], z_nh3, ae [ t ], hmi_nh3, vm [ t ], vi [ t ], param_AAPC)
                gcmi_nh4 = AAPC(cm_nh4, ci_nh4 [ t ], z_nh4, ae [ t ], hmi_nh4, vm [ t ], vi [ t ], param_AAPC)
                gcmi_h = AAPC(cm_h, ci_h [ t ], z_h, ae [ t ], hmi_h, vm [ t ], vi [ t ], param_AAPC)
                gcmi_hco2 = AAPC(cm_hco2, ci_hco2 [ t ], z_hco2, ae [ t ], hmi_hco2, vm [ t ], vi [ t ], param_AAPC)
                gcmi_h2co2 = AAPC(cm_h2co2, ci_h2co2 [ t ], z_h2co2, ae [ t ], hmi_h2co2, vm [ t ], vi [ t ],
                                  param_AAPC)
                gcmi_gluc = AAPC(cm_gluc, ci_gluc [ t ], z_gluc, ae [ t ], hmi_gluc, vm [ t ], vi [ t ], param_AAPC)

                tgcmi = gcmi_na + gcmi_k + gcmi_cl + gcmi_hco3 + gcmi_h2co3 + gcmi_co2 + gcmi_hpo4 + gcmi_h2po4 + \
                        gcmi_urea + gcmi_nh3 + gcmi_nh4 + gcmi_h + gcmi_hco2 + gcmi_h2co2 + gcmi_gluc

                tgcmi = tgcmi / (vm [ t ] - vi [ t ])
                gcie_na = AAPC(ci_na [ t ], ce_na [ t ], z_na, aie, his_na, vi [ t ], ve [ t ], param_AAPC)
                gcie_k = AAPC(ci_k [ t ], ce_k [ t ], z_k, aie, his_k, vi [ t ], ve [ t ], param_AAPC)
                gcie_cl = AAPC(ci_cl [ t ], ce_cl [ t ], z_cl, aie, his_cl, vi [ t ], ve [ t ], param_AAPC)
                gcie_hco3 = AAPC(ci_hco3 [ t ], ce_hco3 [ t ], z_hco3, aie, his_hco3, vi [ t ], ve [ t ], param_AAPC)
                gcie_h2co3 = AAPC(ci_h2co3 [ t ], ce_h2co3 [ t ], z_h2co3, aie, his_h2co3, vi [ t ], ve [ t ],
                                  param_AAPC)
                gcie_co2 = AAPC(ci_co2 [ t ], ce_co2 [ t ], z_co2, aie, his_co2, vi [ t ], ve [ t ], param_AAPC)
                gcie_hpo4 = AAPC(ci_hpo4 [ t ], ce_hpo4 [ t ], z_hpo4, aie, his_hpo4, vi [ t ], ve [ t ], param_AAPC)
                gcie_h2po4 = AAPC(ci_h2po4 [ t ], ce_h2po4 [ t ], z_h2po4, aie, his_h2po4, vi [ t ], ve [ t ],
                                  param_AAPC)
                gcie_urea = AAPC(ci_urea [ t ], ce_urea [ t ], z_urea, aie, his_urea, vi [ t ], ve [ t ], param_AAPC)
                gcie_nh3 = AAPC(ci_nh3 [ t ], ce_nh3 [ t ], z_nh3, aie, his_nh3, vi [ t ], ve [ t ], param_AAPC)
                gcie_nh4 = AAPC(ci_nh4 [ t ], ce_nh4 [ t ], z_nh4, aie, his_nh4, vi [ t ], ve [ t ], param_AAPC)
                gcie_h = AAPC(ci_h [ t ], ce_h [ t ], z_h, aie, his_h, vi [ t ], ve [ t ], param_AAPC)
                gcie_hco2 = AAPC(ci_hco2 [ t ], ce_hco2 [ t ], z_hco2, aie, his_hco2, vi [ t ], ve [ t ], param_AAPC)
                gcie_h2co2 = AAPC(ci_h2co2 [ t ], ce_h2co2 [ t ], z_h2co2, aie, his_h2co2, vi [ t ], ve [ t ],
                                  param_AAPC)
                gcie_gluc = AAPC(ci_gluc [ t ], ce_gluc [ t ], z_gluc, aie, his_gluc, vi [ t ], ve [ t ], param_AAPC)
                # print('gcie_gluc ', gcie_gluc)
                tgcie = gcie_na + gcie_k + gcie_cl + gcie_hco3 + gcie_h2co3 + gcie_co2 + gcie_hpo4 + gcie_h2po4 + gcie_urea \
                        + gcie_nh3 + gcie_nh4 + gcie_h + gcie_hco2 + gcie_h2co2 + gcie_gluc
                tgcie = abs(tgcie / (vi [ t ] - ve [ t ]))

                gcis_na = AAPC(ci_na [ t ], cs_na, z_na, ais, his_na, vi [ t ], vs, param_AAPC)
                gcis_k = AAPC(ci_k [ t ], cs_k, z_k, ais, his_k, vi [ t ], vs, param_AAPC)
                gcis_cl = AAPC(ci_cl [ t ], cs_cl, z_cl, ais, his_cl, vi [ t ], vs, param_AAPC)
                gcis_hco3 = AAPC(ci_hco3 [ t ], cs_hco3, z_hco3, ais, his_hco3, vi [ t ], vs, param_AAPC)
                gcis_h2co3 = AAPC(ci_h2co3 [ t ], cs_h2co3, z_h2co3, ais, his_h2co3, vi [ t ], vs, param_AAPC)
                gcis_co2 = AAPC(ci_co2 [ t ], cs_co2, z_co2, ais, his_co2, vi [ t ], vs, param_AAPC)
                gcis_hpo4 = AAPC(ci_hpo4 [ t ], cs_hpo4, z_hpo4, ais, his_hpo4, vi [ t ], vs, param_AAPC)
                gcis_h2po4 = AAPC(ci_h2po4 [ t ], cs_h2po4, z_h2po4, ais, his_h2po4, vi [ t ], vs, param_AAPC)
                gcis_urea = AAPC(ci_urea [ t ], cs_urea, z_urea, ais, his_urea, vi [ t ], vs, param_AAPC)
                gcis_nh3 = AAPC(ci_nh3 [ t ], cs_nh3, z_nh3, ais, his_nh3, vi [ t ], vs, param_AAPC)
                gcis_nh4 = AAPC(ci_nh4 [ t ], cs_nh4, z_nh4, ais, his_nh4, vi [ t ], vs, param_AAPC)
                gcis_h = AAPC(ci_h [ t ], cs_h, z_h, ais, his_h, vi [ t ], vs, param_AAPC)
                gcis_hco2 = AAPC(ci_hco2 [ t ], cs_hco2, z_hco2, ais, his_hco2, vi [ t ], vs, param_AAPC)
                gcis_h2co2 = AAPC(ci_h2co2 [ t ], cs_h2co2, z_h2co2, ais, his_h2co2, vi [ t ], vs, param_AAPC)
                gcis_gluc = AAPC(ci_gluc [ t ], cs_gluc, z_gluc, ais, his_gluc, vi [ t ], vs, param_AAPC)

                tgcis = gcis_na + gcis_k + gcis_cl + gcis_hco3 + gcis_h2co3 + gcis_co2 + gcis_hpo4 + gcis_h2po4 + \
                        gcis_urea + gcis_nh3 + gcis_nh4 + gcis_h + gcis_hco2 + gcis_h2co2 + gcis_gluc

                tgcis = abs(tgcis / (vi [ t ]))

    return [ fekm_na, fekm_cl, fekm_k, fekm_gluc,
              fikm_na, fikm_cl,fikm_k, fikm_gluc,
              fiks_na, fiks_cl, fiks_k, fiks_gluc,
              jk_na, jk_cl, jk_k, jk_gluc,
              feks_na, feks_cl, feks_k, feks_gluc]




        # [ p_i, gluc_mi_nagluc, na_mi_nagluc, na_mi_nah2po4, jnhe3_na, cl_mi_clhco2, cl_mi_clhco3, fekm_na, fekm_cl, fekm_k,
        #      fekm_gluc, fikm_na, fikm_cl,fikm_k,
        #      fikm_gluc, flux_na, flux_k, flux_cl, flux_hco3, flux_h2co3, flux_co2, flux_hpo4,
        #      flux_h2po4, flux_urea, flux_nh3,
        #      flux_nh4, flux_h, flux_hco2, flux_h2co2, flux_gluc, flux_v, fivm, fevm, cure, curi, cures, curis, curie,
        #      osme, osmi, osmm, osms,
        #      cikm_na, cikm_cl, cikm_hco3, cikm_h2po4, cikm_hco2, cikm_gluc, cikm_h, cikm_nh4, ciks_na, ciks_k, ciks_cl,
        #      ciks_hco3, cjk_na, cjk_k, cjk_cl, cjk_hco3, pikm_h, piks_na, piks_k, piks_nh4, pjk_na, pjk_k, pjk_nh4,
        #      tgcem,
        #      tgcis, tgcmi, tgcie, tgces ]


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
tag = 'clever'
t0 = 0
tf = 500
T = 1500
dt = (tf - t0) / (T - 1)
rtau = 1 / dt

ve = matrix(-0.89432938258185330771e-02, T, tag)
pe = matrix(-0.23093764341580683919e+02, T, tag)
ce_na = matrix(0.14040045563695688347e+00, T, tag)
ce_k = matrix(0.46537228854932385230e-02, T, tag)
ce_cl = matrix(0.11200865813019536543e+00, T, tag)
ce_hco3 = matrix(0.25607874636347928432e-01, T, tag)
ce_h2co3 = matrix(0.43754640651881161056e-05, T, tag)
ce_co2 = matrix(0.14987447089135489450e-02, T, tag)
ce_hpo4 = matrix(0.29852338056241254673e-02, T, tag)
ce_h2po4 = matrix(0.86622197071977683099e-03, T, tag)
ce_urea = matrix(0.49046959853563188228e-02, T, tag)
ce_nh3 = matrix(0.26914919565666300082e-05, T, tag)
ce_nh4 = matrix(0.17484103066624158852e-03, T, tag)
ce_hco2 = matrix(0.77584319325420768570e-03, T, tag)
ce_h2co2 = matrix(0.20531685246992843752e-06, T, tag)
ce_gluc = matrix(0.77236905064622532815e-02, T, tag)
vi = matrix(-0.54945785940474827669e+02, T, tag)
p_i = matrix(0.0e-01, T, tag)
ci_na = matrix(0.20443615264909884011e-01, T, tag)
ci_k = matrix(0.13742335290954643678e+00, T, tag)
ci_cl = matrix(0.16901004032978079322e-01, T, tag)
ci_hco3 = matrix(0.25090576334173893269e-01, T, tag)
ci_h2co3 = matrix(0.43750529566888027647e-05, T, tag)
ci_co2 = matrix(0.14988418688838011684e-02, T, tag)
ci_hpo4 = matrix(0.94254445085816557920e-02, T, tag)
ci_h2po4 = matrix(0.27910960064784369992e-02, T, tag)
ci_urea = matrix(0.49546890743651676378e-02, T, tag)
ci_nh3 = matrix(0.35291888877349245763e-05, T, tag)
ci_nh4 = matrix(0.23396304517990653771e-03, T, tag)
ci_hco2 = matrix(0.53552249401043496273e-03, T, tag)
ci_h2co2 = matrix(0.93379252815431763174e-07, T, tag)
ci_gluc = matrix(0.14876174248492728819e-01, T, tag)
hcbuf = matrix(0.26959905796615793450e-01, T, tag)
cbuf = matrix(0.40012078673998335843e-01, T, tag)
vm = matrix(-0.18573662594736270459e+00, T, tag)
ce_h = matrix(4.59e-8, T, tag)
ci_h = matrix(4.69e-8, T, tag)
ae = matrix(0.02000, T, tag)
chvl = matrix(0.7000e-04, T, tag)
clvl = matrix(0.1000e-02, T, tag)
l = matrix(0.1000e-02, T, tag)
imp = matrix(0.6000e-01, T, tag)
rm = matrix(0.1250e-02, T, tag)
am = matrix(0, T, tag)
phi = matrix(0, 35, tag)
# lchm = pkc + math.log10(cm_hco3 / cm_h2co3)
# lchs = pkc + math.log10(cs_hco3 / cs_h2co3)
# cm_h = 10. ** (-lchm)
# cs_h = 10. ** (-lchs)
# cm_h = 4.95e-8
# cs_h = 4.95e-8

lchm = pkp + np.log10(cm_hpo4 / cm_h2po4)
lchs = pkp + np.log10(cs_hpo4 / cs_h2po4)
cm_h = 10. ** (-lchm)
cs_h = 10. ** (-lchs)
ps = 9.00
vs = 0.00
flux_na = matrix(0, T, tag)
flux_k = matrix(0, T, tag)
flux_cl = matrix(0, T, tag)
flux_hco3 = matrix(0, T, tag)
flux_h2co3 = matrix(0, T, tag)
flux_co2 = matrix(0, T, tag)
flux_hpo4 = matrix(0, T, tag)
flux_h2po4 = matrix(0, T, tag)
flux_urea = matrix(0, T, tag)
flux_nh3 = matrix(0, T, tag)
flux_nh4 = matrix(0, T, tag)
flux_h = matrix(0, T, tag)
flux_hco2 = matrix(0, T, tag)
flux_h2co2 = matrix(0, T, tag)
flux_gluc = matrix(0, T, tag)
flux_v = matrix(0, T, tag)
fivm = matrix(0, T, tag)
fevm = matrix(0, T, tag)
fluxvisj = matrix(0, T, tag)
cure = matrix(0, T, tag)
curi = matrix(0, T, tag)
cures = matrix(0, T, tag)
curie = matrix(0, T, tag)
curis = matrix(0, T, tag)
osme = matrix(0, T, tag)
osmi = matrix(0, T, tag)
osmm = matrix(0, T, tag)
osms = matrix(0, T, tag)
cikm_na = matrix(0, T, tag)
cikm_k = matrix(0, T, tag)
cikm_cl = matrix(0, T, tag)
cikm_hco3 = matrix(0, T, tag)
cikm_h2po4 = matrix(0, T, tag)
cikm_nh4 = matrix(0, T, tag)
cikm_h = matrix(0, T, tag)
cikm_hco2 = matrix(0, T, tag)
cikm_gluc = matrix(0, T, tag)
ciks_na = matrix(0, T, tag)
ciks_k = matrix(0, T, tag)
ciks_cl = matrix(0, T, tag)
ciks_hco3 = matrix(0, T, tag)
cjk_na = matrix(0, T, tag)
cjk_k = matrix(0, T, tag)
cjk_cl = matrix(0, T, tag)
cjk_hco3 = matrix(0, T, tag)
pikm_h = matrix(0, T, tag)
piks_na = matrix(0, T, tag)
piks_k = matrix(0, T, tag)
piks_nh4 = matrix(0, T, tag)
pjk_na = matrix(0, T, tag)
pjk_k = matrix(0, T, tag)
pjk_nh4 = matrix(0, T, tag)

fikm_na = matrix(0, T, tag)
fikm_k = matrix(0, T, tag)
fikm_cl = matrix(0, T, tag)
fikm_gluc = matrix(0, T, tag)

fiks_na = matrix(0, T, tag)
fiks_k = matrix(0, T, tag)
fiks_cl = matrix(0, T, tag)
fiks_gluc = matrix(0, T, tag)

fekm_na = matrix(0, T, tag)
fekm_k = matrix(0, T, tag)
fekm_cl = matrix(0, T, tag)
fekm_gluc = matrix(0, T, tag)

jk_na = matrix(0, T, tag)
jk_k = matrix(0, T, tag)
jk_cl = matrix(0, T, tag)
jk_gluc = matrix(0, T, tag)

feks_na = matrix(0, T, tag)
feks_k = matrix(0, T, tag)
feks_cl = matrix(0, T, tag)
feks_gluc = matrix(0, T, tag)


na_mi_nagluc = matrix(0, T, tag)
gluc_mi_nagluc = matrix(0, T, tag)
na_mi_nah2po4 = matrix(0, T, tag)
jnhe3_na = matrix(0, T, tag)
cl_mi_clhco2 = matrix(0, T, tag)
cl_mi_clhco3 = matrix(0, T, tag)
tgcem = matrix(0, T, tag)
tgcis = matrix(0, T, tag)
tgcmi = matrix(0, T, tag)
tgcie = matrix(0, T, tag)
tgces = matrix(0, T, tag)

# sol(1)='  na '
# sol(2)='  k  '
# sol(3)='  cl '
# sol(4)=' hco3'
# sol(5)='h2co3'
# sol(6)=' co2 '
# sol(7)=' hpo4'
# sol(8)='h2po4'
# sol(9)=' urea'
# sol(10)=' nh3 '
# sol(11)=' nh4 '
# sol(12)='  h  '
# sol(13)=' hco2'
# sol(14)='h2co2'
# sol(15)=' gluc'
vars = [ ve, pe, ce_na, ce_k,
         ce_cl, ce_hco3, ce_h2co3, ce_co2, ce_hpo4, ce_h2po4, ce_urea, ce_nh3, ce_nh4, ce_hco2, ce_h2co2, ce_gluc,
         vi, imp, ci_na, ci_k, ci_cl, ci_hco3, ci_h2co3, ci_co2, ci_hpo4, ci_h2po4, ci_urea, ci_nh3, ci_nh4, ci_hco2,
         ci_h2co2, ci_gluc, cbuf, hcbuf, vm ]
# flx = [ p_i, gluc_mi_nagluc, na_mi_nagluc, na_mi_nah2po4, jnhe3_na, cl_mi_clhco2, cl_mi_clhco3, fekm_na, fekm_cl,
#         fekm_gluc, fikm_na, fikm_cl,
#         fikm_gluc, flux_na, flux_k, flux_cl, flux_hco3, flux_h2co3, flux_co2, flux_hpo4,
#         flux_h2po4, flux_urea, flux_nh3,
#         flux_nh4, flux_h, flux_hco2, flux_h2co2, flux_gluc, flux_v, fivm, fevm, cure, curi, cures, curis, curie, osme,
#         osmi, osmm, osms,
#         cikm_na, cikm_cl, cikm_hco3, cikm_h2po4, cikm_hco2, cikm_gluc, cikm_h, cikm_nh4, ciks_na, ciks_k, ciks_cl,
#         ciks_hco3, cjk_na, cjk_k, cjk_cl, cjk_hco3, pikm_h, piks_na, piks_k, piks_nh4, pjk_na, pjk_k, pjk_nh4, tgcem,
#         tgcis, tgcmi, tgcie, tgces ]
flx = [ fekm_na, fekm_cl, fekm_k, fekm_gluc,
              fikm_na, fikm_cl,fikm_k, fikm_gluc,
              fiks_na, fiks_cl, fiks_k, fiks_gluc,
              jk_na, jk_cl, jk_k, jk_gluc,
              feks_na, feks_cl, feks_k, feks_gluc]
Pressureimpact = 0
sodiumimpact = 0
glucoseimpact = 0
chlorideimpact = 0
plot_concen = 1
if plot_concen:

    t = 0
    guess = [ var [ t ] for var in vars ]

    from scipy import optimize

    # changing the cm_na concentration in step_wise
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mylist = [ 0.1 for i in range(T) ]
    pm = 15.000
    cm_cl = 0.11322499
    cm_gluc = 0.005
    stepwise = 0
    if stepwise:
        cmna = [ 0.100, 0.110, 0.120, 0.13, 0.140 ]
        while True:
            t += 1
            if t == T:
                break
            elif 0 < t <= T / len(cmna):
                cm_na = cmna [ 0 ]
                mylist [ t ] = cmna [ 0 ]
                result = scipy.optimize.root(eQs, np.array(guess), 1)
                guess = [ var [ t ] for var in vars ]
                # print('result=',result.x)
                result_flx = eQs(guess, 0)
                for i in range(len(result_flx)):
                    flx [ i ] [ t ] = result_flx [ i ]
            elif 1 * T / len(cmna) < t <= 2 * T / len(cmna):
                mylist [ t ] = cmna [ 1 ]
                cm_na = cmna [ 1 ]
                result = scipy.optimize.root(eQs, np.array(guess), 1)
                guess = [ var [ t ] for var in vars ]
                result_flx = eQs(guess, 0)
                for i in range(len(result_flx)):
                    flx [ i ] [ t ] = result_flx [ i ]
            elif 2 * T / len(cmna) < t <= 3 * T / len(cmna):
                cm_na = cmna [ 2 ]
                mylist [ t ] = cmna [ 2 ]
                result = scipy.optimize.root(eQs, np.array(guess), 1)
                guess = [ var [ t ] for var in vars ]
                result_flx = eQs(guess, 0)
                for i in range(len(result_flx)):
                    flx [ i ] [ t ] = result_flx [ i ]
            elif 3 * T / len(cmna) < t <= 4 * T / len(cmna):
                cm_na = cmna [ 3 ]
                mylist [ t ] = cmna [ 3 ]
                result = scipy.optimize.root(eQs, np.array(guess), 1)
                guess = [ var [ t ] for var in vars ]
                result_flx = eQs(guess, 0)
                for i in range(len(result_flx)):
                    flx [ i ] [ t ] = result_flx [ i ]
            elif 4 * T / len(cmna) < t <= 5 * T / len(cmna):
                cm_na = cmna [ 4 ]
                mylist [ t ] = cmna [ 4 ]
                result = scipy.optimize.root(eQs, np.array(guess), 1)
                guess = [ var [ t ] for var in vars ]
                result_flx = eQs(guess, 0)
                for i in range(len(result_flx)):
                    flx [ i ] [ t ] = result_flx [ i ]
    else:
        cm_na = 0.14000
        while True:
            t += 1
            if t == T:
                break
            else:
                cm_na = 0.14
                result = scipy.optimize.root(eQs, np.array(guess), 1)
                guess = [ var [ t ] for var in vars ]
                result_flx = eQs(guess, 0)
                for i in range(len(result_flx)):
                    flx [ i ] [ t ] = result_flx [ i ]

#    print('param-nhe3', param_nhe3)
    print('guess', guess)
    print('result_flx', result_flx)
    t = np.linspace(t0, tf, T)
    t_t = np.transpose(t)
    mylist = np.asarray(mylist)
    fig2 = plt.figure(constrained_layout=True)
    spec2 = gridspec.GridSpec(ncols=1, nrows=8, figure=fig2)
    ax1 = fig2.add_subplot(spec2 [ 0, 0 ])
    ax1.plot(t [ 0:-2 ], mylist [ 0:-2 ], 'r-')
    ax1.set_title('Solute Concentration OF THE CELL [mmol] ')
    ax1.set_ylabel('cm_na (mm)', color='blue')
    ax2 = fig2.add_subplot(spec2 [ 1, 0 ])
    ax2.plot(t, ci_na, 'b-')
    ax2.set_ylabel(' ci_na ', color='blue')
    ax3 = fig2.add_subplot(spec2 [ 2, 0 ])
    ax3.plot(t, ci_k, 'b-')
    ax3.set_ylabel('ci_k ', color='blue')
    ax4 = fig2.add_subplot(spec2 [ 3, 0 ])
    ax4.plot(t, ci_cl, 'b-')
    ax4.set_ylabel('ci_cl ', color='blue')
    ax5 = fig2.add_subplot(spec2 [ 4, 0 ])
    ax5.plot(t, ci_urea, 'b-')
    ax5.set_ylabel('ci_urea ', color='blue')
    ax6 = fig2.add_subplot(spec2 [ 5, 0 ])
    ax6.plot(t, ci_nh3, 'b-')
    ax6.set_ylabel('ci_nh3', color='blue')
    ax7 = fig2.add_subplot(spec2 [ 6, 0 ])
    ax7.plot(t, ci_nh4, 'b-')
    ax7.set_ylabel('ci_nh4', color='blue')
    ax8 = fig2.add_subplot(spec2 [ 7, 0 ])
    ax8.plot(t, ci_h, 'b-')
    ax8.set_ylabel('ci_h ', color='blue')
    ax8.set_xlabel('Time (s) ', color='blue')

    fig3 = plt.figure(constrained_layout=True)
    spec2 = gridspec.GridSpec(ncols=1, nrows=8, figure=fig3)
    ax1 = fig3.add_subplot(spec2 [ 0, 0 ])
    ax1.plot(t [ 0:-2 ], mylist [ 0:-2 ], 'r-')
    ax1.set_title('Solute Concentration OF THE CELL [mmol] ')
    ax1.set_ylabel('cm_na (mm)', color='blue')
    ax2 = fig3.add_subplot(spec2 [ 1, 0 ])
    ax2.plot(t, ci_gluc, 'b-')
    ax2.set_ylabel('ci_gluc ', color='blue')
    ax3 = fig3.add_subplot(spec2 [ 2, 0 ])
    ax3.plot(t, ci_hco3, 'b-')
    ax3.set_ylabel('ci_hco3 ', color='blue')
    ax4 = fig3.add_subplot(spec2 [ 3, 0 ])
    ax4.plot(t, ci_h2co3, 'b-')
    ax4.set_ylabel('ci_h2co3 ', color='blue')
    ax5 = fig3.add_subplot(spec2 [ 4, 0 ])
    ax5.plot(t, ci_hpo4, 'b-')
    ax5.set_ylabel('ci_hpo4 ', color='blue')
    ax6 = fig3.add_subplot(spec2 [ 5, 0 ])
    ax6.plot(t, ci_h2po4, 'b-')
    ax6.set_ylabel('ci_h2po4', color='blue')
    ax7 = fig3.add_subplot(spec2 [ 6, 0 ])
    ax7.plot(t, ci_hco2, 'b-')
    ax7.set_ylabel('ci_hco2', color='blue')
    ax8 = fig3.add_subplot(spec2 [ 7, 0 ])
    ax8.plot(t, ci_h2co2, 'b-')
    ax8.set_ylabel('ci_h2co2 ', color='blue')
    ax8.set_xlabel('Time (s) ', color='blue')

    fig5 = plt.figure(constrained_layout=True)
    spec2 = gridspec.GridSpec(ncols=1, nrows=8, figure=fig5)
    ax1 = fig5.add_subplot(spec2 [ 0, 0 ])
    ax1.plot(t [ 0:-2 ], mylist [ 0:-2 ], 'r-')
    ax1.set_title('Solute Concentration OF THE CHANNEL [mmol]')
    ax1.set_ylabel('cm_na (mm)', color='blue')
    ax2 = fig5.add_subplot(spec2 [ 1, 0 ])
    ax2.plot(t, ce_na, 'b-')
    ax2.set_ylabel(' ce_na ', color='blue')
    ax3 = fig5.add_subplot(spec2 [ 2, 0 ])
    ax3.plot(t, ce_k, 'b-')
    ax3.set_ylabel('ce_k ', color='blue')
    ax4 = fig5.add_subplot(spec2 [ 3, 0 ])
    ax4.plot(t, ce_cl, 'b-')
    ax4.set_ylabel('ce_cl ', color='blue')
    ax5 = fig5.add_subplot(spec2 [ 4, 0 ])
    ax5.plot(t, ce_urea, 'b-')
    ax5.set_ylabel('ce_urea ', color='blue')
    ax6 = fig5.add_subplot(spec2 [ 5, 0 ])
    ax6.plot(t, ce_nh3, 'b-')
    ax6.set_ylabel('ce_nh3', color='blue')
    ax7 = fig5.add_subplot(spec2 [ 6, 0 ])
    ax7.plot(t, ce_nh4, 'b-')
    ax7.set_ylabel('ce_nh4', color='blue')
    ax8 = fig5.add_subplot(spec2 [ 7, 0 ])
    ax8.plot(t, ce_h, 'b-')
    ax8.set_ylabel('ce_h ', color='blue')
    ax8.set_xlabel('Time (s) ', color='blue')
    fig3 = plt.figure(constrained_layout=True)
    spec2 = gridspec.GridSpec(ncols=1, nrows=8, figure=fig3)
    ax1 = fig3.add_subplot(spec2 [ 0, 0 ])
    ax1.plot(t [ 0:-2 ], mylist [ 0:-2 ], 'r-')
    ax1.set_title('Solute Concentration OF THE CELL [mmol] ')
    ax1.set_ylabel('cm_na (mm)', color='blue')
    ax2 = fig3.add_subplot(spec2 [ 1, 0 ])
    ax2.plot(t, ce_gluc, 'b-')
    ax2.set_ylabel('ce_gluc ', color='blue')
    ax3 = fig3.add_subplot(spec2 [ 2, 0 ])
    ax3.plot(t, ce_hco3, 'b-')
    ax3.set_ylabel('ce_hco3 ', color='blue')
    ax4 = fig3.add_subplot(spec2 [ 3, 0 ])
    ax4.plot(t, ce_h2co3, 'b-')
    ax4.set_ylabel('ce_h2co3 ', color='blue')
    ax5 = fig3.add_subplot(spec2 [ 4, 0 ])
    ax5.plot(t, ce_hpo4, 'b-')
    ax5.set_ylabel('ce_hpo4 ', color='blue')
    ax6 = fig3.add_subplot(spec2 [ 5, 0 ])
    ax6.plot(t, ce_h2po4, 'b-')
    ax6.set_ylabel('ce_h2po4', color='blue')
    ax7 = fig3.add_subplot(spec2 [ 6, 0 ])
    ax7.plot(t, ce_hco2, 'b-')
    ax7.set_ylabel('ce_hco2', color='blue')
    ax8 = fig3.add_subplot(spec2 [ 7, 0 ])
    ax8.plot(t, ce_h2co2, 'b-')
    ax8.set_ylabel('ce_h2co2 ', color='blue')
    ax8.set_xlabel('Time (s) ', color='blue')

    fig6 = plt.figure(constrained_layout=True)
    spec2 = gridspec.GridSpec(ncols=1, nrows=4, figure=fig6)
    ax5 = fig6.add_subplot(spec2 [ 0, 0 ])
    ax5.plot(t [ 0:-2 ], mylist [ 0:-2 ], 'r-', linewidth=2.0)
    ax5.set_xlim(-5, tf)
    ax5.set_title('MEMBRANE POTENTIAL [mVolt]')
    ax5.set_ylabel('cm_na (mm) ', color='blue')
    ax6 = fig6.add_subplot(spec2 [ 1, 0 ])
    vem = np.asarray(ve) - np.asarray(vm)
    v_em = np.asarray(vem)
    ax6.plot(t, v_em, 'b-', linewidth=2.0)
    ax6.set_ylabel('voltage_em (mv) ', color='blue')
    ax6.set_xlabel('Time (s) ', color='blue')
    ax6.set_xlim(-5, tf)
    ax7 = fig6.add_subplot(spec2 [ 2, 0 ])
    vie = np.asarray(vi) - np.asarray(ve)
    v_ie = np.asarray(vie)
    ax7.plot(t, v_ie, 'b-', linewidth=2.0)
    ax7.set_ylabel('voltage_ie (mv) ', color='blue')
    ax7.set_xlabel('Time (s) ', color='blue')
    ax7.set_xlim(-5, tf)
    ax7 = fig6.add_subplot(spec2 [ 3, 0 ])
    vmi = np.asarray(vm) - np.asarray(vi)
    v_mi = np.asarray(vie)
    ax7.plot(t, v_mi, 'b-', linewidth=2.0)
    ax7.set_ylabel('voltage_mi (mv) ', color='blue')
    ax7.set_xlabel('Time (s) ', color='blue')
    ax7.set_xlim(-5, tf)
    # plt.show()
    # epithelial_flux=1,
    # if epithelial_flux:
    fig5 = plt.figure(constrained_layout=True)
    spec2 = gridspec.GridSpec(ncols=1, nrows=8, figure=fig5)
    ax1 = fig5.add_subplot(spec2 [ 0, 0 ])
    ax1.plot(t [ 0:-2 ], mylist [ 0:-2 ], 'r-')
    ax1.set_title('TOTAL EPITHELIAL SOLUTE FLUX [mmol/sec]')
    ax1.set_ylabel('cm_na (mm)', color='blue')
    ax2 = fig5.add_subplot(spec2 [ 1, 0 ])
    ax2.plot(t, flux_na, 'b-')
    ax2.set_ylabel(' flux_na ', color='blue')
    ax3 = fig5.add_subplot(spec2 [ 2, 0 ])
    ax3.plot(t, flux_k, 'b-')
    ax3.set_ylabel('flux_k ', color='blue')
    ax4 = fig5.add_subplot(spec2 [ 3, 0 ])
    ax4.plot(t, flux_cl, 'b-')
    ax4.set_ylabel('flux_cl ', color='blue')
    ax5 = fig5.add_subplot(spec2 [ 4, 0 ])
    ax5.plot(t, flux_urea, 'b-')
    ax5.set_ylabel('flux_urea ', color='blue')
    ax6 = fig5.add_subplot(spec2 [ 5, 0 ])
    ax6.plot(t, flux_nh3, 'b-')
    ax6.set_ylabel('flux_nh3', color='blue')
    ax7 = fig5.add_subplot(spec2 [ 6, 0 ])
    ax7.plot(t, flux_nh4, 'b-')
    ax7.set_ylabel('flux_nh4', color='blue')
    ax8 = fig5.add_subplot(spec2 [ 7, 0 ])
    ax8.plot(t, flux_h, 'b-')
    ax8.set_ylabel('flux_h ', color='blue')
    fig8 = plt.figure(constrained_layout=True)
    spec2 = gridspec.GridSpec(ncols=1, nrows=8, figure=fig8)
    ax1 = fig8.add_subplot(spec2 [ 0, 0 ])
    ax1.plot(t [ 0:-2 ], mylist [ 0:-2 ], 'r-')
    ax1.set_title('TOTAL EPITHELIAL SOLUTE FLUX [mmol/sec]')
    ax1.set_ylabel('cm_na (mm)', color='blue')
    ax2 = fig8.add_subplot(spec2 [ 1, 0 ])
    ax2.plot(t, flux_gluc, 'b-')
    ax2.set_ylabel('flux_gluc ', color='blue')
    ax3 = fig8.add_subplot(spec2 [ 2, 0 ])
    ax3.plot(t, flux_hco3, 'b-')
    ax3.set_ylabel('flux_hco3 ', color='blue')
    ax4 = fig8.add_subplot(spec2 [ 3, 0 ])
    ax4.plot(t, flux_h2co3, 'b-')
    ax4.set_ylabel('flux_h2co3 ', color='blue')
    ax5 = fig8.add_subplot(spec2 [ 4, 0 ])
    ax5.plot(t, flux_hpo4, 'b-')
    ax5.set_ylabel('flux_hpo4 ', color='blue')
    ax6 = fig8.add_subplot(spec2 [ 5, 0 ])
    ax6.plot(t, flux_h2po4, 'b-')
    ax6.set_ylabel('flux_h2po4', color='blue')
    ax7 = fig8.add_subplot(spec2 [ 6, 0 ])
    ax7.plot(t, flux_hco2, 'b-')
    ax7.set_ylabel('flux_hco2', color='blue')
    ax8 = fig8.add_subplot(spec2 [ 7, 0 ])
    ax8.plot(t, flux_h2co2, 'b-')
    ax8.set_ylabel('flux_h2co2 ', color='blue')
    ax8.set_xlabel('Time (s) ', color='blue')
    # plt.show()
    # osmilarity = 1
    # if osmilarity:
    fig8 = plt.figure(constrained_layout=True)
    spec2 = gridspec.GridSpec(ncols=1, nrows=4, figure=fig8)
    ax1 = fig8.add_subplot(spec2 [ 0, 0 ])
    ax1.plot(t, osmi, 'r-')
    ax1.set_title('OSMOLALITY OF INDICATED COMPARTMENT. = E, I, M, S [mmol]')
    ax1.set_ylabel('osm_i', color='blue')
    ax2 = fig8.add_subplot(spec2 [ 1, 0 ])
    ax2.plot(t, osme, 'b-')
    ax2.set_ylabel(' osm_e', color='blue')
    ax3 = fig8.add_subplot(spec2 [ 2, 0 ])
    ax3.plot(t, osmm, 'b-')
    ax3.set_ylabel('osm_m ', color='blue')
    ax4 = fig8.add_subplot(spec2 [ 3, 0 ])
    ax4.plot(t, osms, 'b-')
    ax4.set_ylabel('osm_s ', color='blue')
    # plt.show()
    # flux - TOTAL EPITHELIAL SOLUTE FLUX
    # epithelial_flux = 1
    # if  epithelial_flux:
    fig2 = plt.figure(constrained_layout=True)
    spec2 = gridspec.GridSpec(ncols=1, nrows=5, figure=fig2)
    ax1 = fig2.add_subplot(spec2 [ 0, 0 ])
    ax1.plot(t [ 0:-2 ], mylist [ 0:-2 ], 'r-')
    ax1.set_title('COTRANSPORT  ACROSS APICAL M [mmol/sec]')
    ax1.set_ylabel('cm_na (mm)', color='blue')
    ax2 = fig2.add_subplot(spec2 [ 1, 0 ])
    ax2.plot(t, cikm_na, 'b-')
    ax2.set_ylabel(' flux_na ', color='blue')
    ax3 = fig2.add_subplot(spec2 [ 2, 0 ])
    ax3.plot(t, cikm_cl, 'b-')
    ax3.set_ylabel('flux_cl ', color='blue')
    ax4 = fig2.add_subplot(spec2 [ 3, 0 ])
    ax4.plot(t, cikm_hco3, 'b-')
    ax4.set_ylabel('flux_hco3 ', color='blue')
    ax5 = fig2.add_subplot(spec2 [ 4, 0 ])
    ax5.plot(t, cikm_h2po4, 'b-')
    ax5.set_ylabel('flux_h2po4 ', color='blue')

    fig3 = plt.figure(constrained_layout=True)
    spec2 = gridspec.GridSpec(ncols=1, nrows=5, figure=fig3)
    ax1 = fig3.add_subplot(spec2 [ 0, 0 ])
    ax1.plot(t [ 0:-2 ], mylist [ 0:-2 ], 'r-')
    ax1.set_title('COTRANSPORT ACROSS APICAL M [mmol/sec]')
    ax1.set_ylabel('cm_na (mm)', color='blue')
    ax2 = fig3.add_subplot(spec2 [ 1, 0 ])
    ax2.plot(t, cikm_hco2, 'b-')
    ax2.set_ylabel(' flux_hco2 ', color='blue')
    ax3 = fig3.add_subplot(spec2 [ 2, 0 ])
    ax3.plot(t, cikm_gluc, 'b-')
    ax3.set_ylabel('flux_gluc ', color='blue')
    ax4 = fig3.add_subplot(spec2 [ 3, 0 ])
    ax4.plot(t, cikm_h, 'b-')
    ax4.set_ylabel('flux_h ', color='blue')
    ax5 = fig3.add_subplot(spec2 [ 4, 0 ])
    ax5.plot(t, cikm_nh4, 'b-')
    ax5.set_ylabel('flux_nh4 ', color='blue')

    fig4 = plt.figure(constrained_layout=True)
    spec2 = gridspec.GridSpec(ncols=1, nrows=5, figure=fig4)
    ax1 = fig4.add_subplot(spec2 [ 0, 0 ])
    ax1.plot(t [ 0:-2 ], mylist [ 0:-2 ], 'r-')
    ax1.set_title('COTRANSPORT ACROSS BASAL M [mmol/sec]')
    ax1.set_ylabel('cm_na (mm)', color='blue')
    ax2 = fig4.add_subplot(spec2 [ 1, 0 ])
    ax2.plot(t, ciks_na, 'b-')
    ax2.set_ylabel('flux_na', color='blue')
    ax3 = fig4.add_subplot(spec2 [ 2, 0 ])
    ax3.plot(t, ciks_k, 'b-')
    ax3.set_ylabel('flux_k', color='blue')
    ax4 = fig4.add_subplot(spec2 [ 3, 0 ])
    ax4.plot(t, ciks_cl, 'b-')
    ax4.set_ylabel('flux_cl ', color='blue')
    ax5 = fig4.add_subplot(spec2 [ 4, 0 ])
    ax5.plot(t, ciks_hco3, 'b-')
    ax5.set_ylabel('flux_hco3 ', color='blue')

    fig5 = plt.figure(constrained_layout=True)
    spec2 = gridspec.GridSpec(ncols=1, nrows=5, figure=fig5)
    ax1 = fig5.add_subplot(spec2 [ 0, 0 ])
    ax1.plot(t [ 0:-2 ], mylist [ 0:-2 ], 'r-')
    ax1.set_title('COTRANSPORT ACROSS LATERAL M [mmol/sec]')
    ax1.set_ylabel('cm_na (mm)', color='blue')
    ax2 = fig5.add_subplot(spec2 [ 1, 0 ])
    ax2.plot(t, cjk_na, 'b-')
    ax2.set_ylabel('flux_na ', color='blue')
    ax3 = fig5.add_subplot(spec2 [ 2, 0 ])
    ax3.plot(t, cjk_k, 'b-')
    ax3.set_ylabel('flux_k ', color='blue')
    ax4 = fig5.add_subplot(spec2 [ 3, 0 ])
    ax4.plot(t, cjk_cl, 'b-')
    ax4.set_ylabel('flux_cl ', color='blue')
    ax5 = fig5.add_subplot(spec2 [ 4, 0 ])
    ax5.plot(t, cjk_hco3, 'b-')
    ax5.set_ylabel('flux_hco3 ', color='blue')

    fig6 = plt.figure(constrained_layout=True)
    spec2 = gridspec.GridSpec(ncols=1, nrows=4, figure=fig6)
    ax1 = fig6.add_subplot(spec2 [ 0, 0 ])
    ax1.plot(t [ 0:-2 ], mylist [ 0:-2 ], 'r-')
    ax1.set_title('ACTIVE TRANSPORT ACROSS LATERAL M, Na PUMP [mmol/sec]')
    ax1.set_ylabel('cm_na (mm)', color='blue')
    ax2 = fig6.add_subplot(spec2 [ 1, 0 ])
    ax2.plot(t, pjk_na, 'b-')
    ax2.set_ylabel(' flux_na ', color='blue')
    ax3 = fig6.add_subplot(spec2 [ 2, 0 ])
    ax3.plot(t, pjk_k, 'b-')
    ax3.set_ylabel('flux_k ', color='blue')
    ax4 = fig6.add_subplot(spec2 [ 3, 0 ])
    ax4.plot(t, pjk_nh4, 'b-')
    ax4.set_ylabel('flux_nh4 ', color='blue')

    fig7 = plt.figure(constrained_layout=True)
    spec2 = gridspec.GridSpec(ncols=1, nrows=4, figure=fig7)
    ax1 = fig7.add_subplot(spec2 [ 0, 0 ])
    ax1.plot(t [ 0:-2 ], mylist [ 0:-2 ], 'r-')
    ax1.set_title('ACTIVE TRANSPORT ACROSS BASAL M,Na PUMP [mmol/sec]')
    ax1.set_ylabel('cm_na (mm)', color='blue')
    ax2 = fig7.add_subplot(spec2 [ 1, 0 ])
    ax2.plot(t, piks_na, 'b-')
    ax2.set_ylabel(' flux_na ', color='blue')
    ax3 = fig7.add_subplot(spec2 [ 2, 0 ])
    ax3.plot(t, piks_k, 'b-')
    ax3.set_ylabel('flux_k ', color='blue')
    ax4 = fig7.add_subplot(spec2 [ 3, 0 ])
    ax4.plot(t, piks_nh4, 'b-')
    ax4.set_ylabel('flux_nh4 ', color='blue')

    fig8 = plt.figure(constrained_layout=True)
    spec2 = gridspec.GridSpec(ncols=1, nrows=2, figure=fig8)
    ax1 = fig8.add_subplot(spec2 [ 0, 0 ])
    ax1.plot(t [ 0:-2 ], mylist [ 0:-2 ], 'r-')
    ax1.set_title('ACTIVE TRANSPORT ACROSS APICAL M,H PUMP [mmol/sec]')
    ax1.set_ylabel('cm_na (mm)', color='blue')
    ax2 = fig8.add_subplot(spec2 [ 1, 0 ])
    ax2.plot(t, pikm_h, 'b-')
    ax2.set_ylabel(' flux_h ', color='blue')

    fig7 = plt.figure(constrained_layout=True)
    spec2 = gridspec.GridSpec(ncols=1, nrows=5, figure=fig7)
    ax1 = fig7.add_subplot(spec2 [ 0, 0 ])
    ax1.plot(t, tgcis, 'b-')
    ax1.set_title('AREA ADJUSTED TOTAL CONDUCTANCE [MHO/CM2]')
    ax1.set_ylabel('conduct_is', color='blue')
    ax2 = fig7.add_subplot(spec2 [ 1, 0 ])
    ax2.plot(t, tgcie, 'b-')
    ax2.set_ylabel('conduct_ie', color='blue')
    ax3 = fig7.add_subplot(spec2 [ 2, 0 ])
    ax3.plot(t, tgces, 'b-')
    ax3.set_ylabel('conduct_es', color='blue')
    ax4 = fig7.add_subplot(spec2 [ 3, 0 ])
    ax4.plot(t, tgcem, 'b-')
    ax4.set_ylabel('conduct_em', color='blue')
    ax5 = fig7.add_subplot(spec2 [ 4, 0 ])
    ax5.plot(t, tgcmi, 'b-')
    ax5.set_ylabel('conduct_mi', color='blue')

    fig7 = plt.figure(constrained_layout=True)
    spec2 = gridspec.GridSpec(ncols=1, nrows=5, figure=fig7)
    ax1 = fig7.add_subplot(spec2 [ 0, 0 ])
    ax1.plot(t, cure, 'blue')
    ax1.set_title('CURRENT [mAmp]')
    ax1.set_ylabel('I_em', color='blue')
    ax2 = fig7.add_subplot(spec2 [ 1, 0 ])
    ax2.plot(t, curi, 'b-')
    ax2.set_ylabel('I_im', color='blue')
    ax3 = fig7.add_subplot(spec2 [ 2, 0 ])
    ax3.plot(t, cures, 'b-')
    ax3.set_ylabel('I_es', color='blue')
    ax4 = fig7.add_subplot(spec2 [ 3, 0 ])
    ax4.plot(t, curis, 'b-')
    ax4.set_ylabel('I_is', color='blue')
    ax5 = fig7.add_subplot(spec2 [ 4, 0 ])
    ax5.plot(t, curie, 'b-')
    ax5.set_ylabel('I_ie', color='blue')
    plt.show()
elif sodiumimpact:
    t = 0
    guess = [ var [ t ] for var in vars ]
    from scipy import optimize

    # changing the cm_na concentration in step_wise
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mylist = [ 0.1 for i in range(T) ]
    pm = 15.000
    cm_cl = 0.11322499
    cm_gluc = 0.005
    stepwise = 1
    if stepwise:
        cmna = [ 0.100, 0.105, 0.110, 0.115, 0.120 ]
        while True:
            t += 1
            if t == T:
                break
            elif 0 < t <= T / len(cmna):
                cm_na = cmna [ 0 ]
                mylist [ t ] = cmna [ 0 ]
                result = scipy.optimize.root(eQs, np.array(guess), 1)
                guess = [ var [ t ] for var in vars ]
                # print('result=',result.x)
                result_flx = eQs(guess, 0)
                for i in range(len(result_flx)):
                    flx [ i ] [ t ] = result_flx [ i ]
            elif 1 * T / len(cmna) < t <= 2 * T / len(cmna):
                mylist [ t ] = cmna [ 1 ]
                cm_na = cmna [ 1 ]
                result = scipy.optimize.root(eQs, np.array(guess), 1)
                guess = [ var [ t ] for var in vars ]
                result_flx = eQs(guess, 0)
                for i in range(len(result_flx)):
                    flx [ i ] [ t ] = result_flx [ i ]
            elif 2 * T / len(cmna) < t <= 3 * T / len(cmna):
                cm_na = cmna [ 2 ]
                mylist [ t ] = cmna [ 2 ]
                result = scipy.optimize.root(eQs, np.array(guess), 1)
                guess = [ var [ t ] for var in vars ]
                result_flx = eQs(guess, 0)
                for i in range(len(result_flx)):
                    flx [ i ] [ t ] = result_flx [ i ]
            elif 3 * T / len(cmna) < t <= 4 * T / len(cmna):
                cm_na = cmna [ 3 ]
                mylist [ t ] = cmna [ 3 ]
                result = scipy.optimize.root(eQs, np.array(guess), 1)
                guess = [ var [ t ] for var in vars ]
                result_flx = eQs(guess, 0)
                for i in range(len(result_flx)):
                    flx [ i ] [ t ] = result_flx [ i ]
            elif 4 * T / len(cmna) < t <= 5 * T / len(cmna):
                cm_na = cmna [ 4 ]
                mylist [ t ] = cmna [ 4 ]
                result = scipy.optimize.root(eQs, np.array(guess), 1)
                guess = [ var [ t ] for var in vars ]
                result_flx = eQs(guess, 0)
                for i in range(len(result_flx)):
                    flx [ i ] [ t ] = result_flx [ i ]
    else:
        cm_na = 0.11000
        while True:
            t += 1
            if t == T:
                break
            else:
                cm_na = 0.11
                result = scipy.optimize.root(eQs, np.array(guess), 1)
                guess = [ var [ t ] for var in vars ]
                result_flx = eQs(guess, 0)
                for i in range(len(result_flx)):
                    flx [ i ] [ t ] = result_flx [ i ]

    t = np.linspace(t0, tf, T, tag)
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec


    def format_axes(fig18):
        for i, ax in enumerate(fig18.axes):
            ax.text(0.5, 0.5, "ax%d" % (i + 1), va="center", ha="center")
            ax.tick_params(labelbottom=False, labelleft=True, labelsize=10)


    fig18 = plt.figure(constrained_layout=False, figsize=(6, 9))
    spec2 = gridspec.GridSpec(ncols=1, nrows=7, hspace=0, figure=fig18)
    anno_opts = dict(xy=(0.8, 0.1), xycoords='axes fraction',
                     va='center', ha='center')

    ax1 = fig18.add_subplot(spec2 [ 0, 0 ])
    x1 = ax1.plot(t [ 0:-2 ], mylist [ 0:-2 ], 'r-')
    ax1.annotate('Luminal Concentration', **anno_opts)
    ax1.tick_params(labelbottom=False, labelleft=True, labelsize=10)
    ax2 = fig18.add_subplot(spec2 [ 1, 0 ])
    x2 = ax2.plot(t, ci_na, 'blue')
    ax2.tick_params(labelbottom=False, labelleft=True, labelsize=10)
    ax2.set_ylabel('c_na [mm]', color='blue', fontsize=12)
    ax2.annotate('Cellular Concentration', **anno_opts)
    ax3 = fig18.add_subplot(spec2 [ 2, 0 ])
    x2 = ax3.plot(t, ce_na, 'blue')
    # ax3.tick_params(labelbottom=False, labelleft=True, labelsize=14)
    ax3.annotate('Junctional Concentration', **anno_opts)
    ax3.set_ylim([ 0.12, 0.138 ])
    ax3.tick_params(labelbottom=False, labelleft=True, labelsize=10)
    ax4 = fig18.add_subplot(spec2 [ 3, 0 ])
    # ax4.set_ylabel('f_na [mmol/sec]', color='blue')
    fikm_nas = [ phi_scale(x, 1e5) for x in fikm_na ]
    fekm_nas = [ phi_scale(x, 1e5) for x in fekm_na ]
    x2 = ax4.plot(t, fikm_nas, 'r')
    x3 = ax4.plot(t, fekm_nas, 'blue')
    ax4.set_ylim([ 0.0, 1.25 ])
    ax4.legend([ x2 [ 0 ], x3 [ 0 ] ], [ 'Cellular Na Fluxes', 'Junctional Na Fluxes' ])
    ax4.tick_params(labelbottom=False, labelleft=True, labelsize=10)
    ax2 = fig18.add_subplot(spec2 [ 4, 0 ])
    ax2.set_ylabel('f$1e5$ [mmol/sec]', color='blue', fontsize=12)
    ax2.tick_params(labelbottom=False, labelleft=True, labelsize=10)
    fivms = [ phi_scale(x, 1e5) for x in fivm ]

    fevms = [ phi_scale(x, 1e5) for x in fevm ]
    x2 = ax2.plot(t, fivms, 'g')
    x3 = ax2.plot(t, fevms, 'blue')
    ax2.legend([ x2 [ 0 ], x3 [ 0 ] ], [ 'Cellular-Volume', 'Junctional-Volume' ])

    # ax4.annotate('Sodium Fluxes', **anno_opts)
    # ax2.set_title('Sodium Fluxes')
    ax5 = fig18.add_subplot(spec2 [ -2:, 0 ])
    # ax5.set_ylabel('            f_na$1e5$ [mmol/sec]', color='blue', fontsize=12)

    na_mi_naglucs = [ phi_scale(x, 1e5) for x in na_mi_nagluc ]
    na_mi_nah2po4s = [ phi_scale(x, 1e5) for x in na_mi_nah2po4 ]
    jnhe3_nas = [ phi_scale(x, 1e5) for x in jnhe3_na ]
    x1 = ax5.plot(t, na_mi_naglucs, 'r')
    x2 = ax5.plot(t, na_mi_nah2po4s, 'b')
    x3 = ax5.plot(t, jnhe3_nas, 'g')

    # ax5.annotate('Cotransporters', **anno_opts)
    ax5.legend([ x3 [ 0 ], x2 [ 0 ], x1 [ 0 ] ], [ 'nhe3', 'na/h2po4', 'na-gluc' ])
    ax5.set_xlabel('Time [sec]', color='blue', fontsize=12)
    # ax2.set_title('Sodium Cotransporters on Apical Membrane')
    ax5.tick_params(labelbottom=True, labelleft=True, labelsize=10)

    # format_axes(fig18)

    fig19 = plt.figure(constrained_layout=True)
    spec2 = gridspec.GridSpec(ncols=1, nrows=2, figure=fig19)
    ax1 = fig19.add_subplot(spec2 [ 0, 0 ])
    x1 = ax1.plot(t, fikm_na, 'g')
    x2 = ax1.plot(t, jnhe3_na, 'b')
    ax1.legend([ x1 [ 0 ], x2 [ 0 ] ], [ 'Cellular', 'nhe3' ])
    ax1.set_title('Sodium Fluxes')
    ax2 = fig19.add_subplot(spec2 [ 1, 0 ])
    ax2.set_ylabel('flux_na [mmol/sec]', color='blue')
    x1 = ax2.plot(t, na_mi_nagluc, 'r')
    x2 = ax2.plot(t, na_mi_nah2po4, 'b')
    x3 = ax2.plot(t, jnhe3_na, 'g')
    ax2.legend([ x3 [ 0 ], x2 [ 0 ], x1 [ 0 ] ], [ 'nhe3', 'na/h2po4', 'na-gluc' ])
    ax2.set_title('Sodium Cotransporters on Apical Membrane')

    fig20 = plt.figure(constrained_layout=True)
    spec2 = gridspec.GridSpec(ncols=1, nrows=2, figure=fig20)
    ax1 = fig20.add_subplot(spec2 [ 0, 0 ])
    x1 = ax1.plot(t [ 0:-2 ], mylist [ 0:-2 ], 'r')
    x2 = ax1.plot(t, ci_cl, 'blue')
    ax1.legend([ x1 [ 0 ], x2 [ 0 ] ], [ 'Luminal_na', 'Cellular_cl' ])
    ax1.set_ylabel('c (mm)', color='blue')
    ax1.set_title('Solute Concentration')
    ax2 = fig20.add_subplot(spec2 [ 1, 0 ])
    ax2.set_ylabel('flux_cl [mmol/sec]', color='blue')
    x1 = ax2.plot(t, flux_cl, 'r')
    x2 = ax2.plot(t, fikm_cl, 'g')
    x3 = ax2.plot(t, fekm_cl, 'blue')
    ax2.legend([ x1 [ 0 ], x2 [ 0 ], x3 [ 0 ] ], [ 'Total Epithelial', 'Cellular', 'Junctional' ])
    ax2.set_title('chloride fluxes')

    fig21 = plt.figure(constrained_layout=True)
    spec2 = gridspec.GridSpec(ncols=1, nrows=2, figure=fig21)
    ax1 = fig21.add_subplot(spec2 [ 0, 0 ])
    x1 = ax1.plot(t, cl_mi_clhco2, 'r')
    x2 = ax1.plot(t, fikm_cl, 'g')
    ax1.set_title('chloride fluxes')
    ax1.legend([ x1 [ 0 ], x2 [ 0 ] ], [ 'cl/hco2', 'Cellular' ])
    ax2 = fig21.add_subplot(spec2 [ 1, 0 ])
    ax2.set_ylabel('flux_cl [mmol/sec]', color='blue')
    x1 = ax2.plot(t, cl_mi_clhco2, 'r')
    x2 = ax2.plot(t, cl_mi_clhco3, 'b')
    ax2.legend([ x1 [ 0 ], x2 [ 0 ] ], [ 'cl/hco2', 'cl/hco3' ])
    ax2.set_title('chloride cotransporters on apical membrane')

    fig22 = plt.figure(constrained_layout=True)
    spec2 = gridspec.GridSpec(ncols=1, nrows=3, figure=fig22)
    ax1 = fig22.add_subplot(spec2 [ 0, 0 ])
    x1 = ax1.plot(t [ 0:-2 ], mylist [ 0:-2 ], 'r')
    x2 = ax1.plot(t, ci_gluc, 'blue')
    ax1.legend([ x1 [ 0 ], x2 [ 0 ] ], [ 'Luminal_na', 'Cellular_gluc' ])
    ax1.set_ylabel('c (mm)', color='blue')
    ax1.set_title('Solute Concentration')
    ax2 = fig22.add_subplot(spec2 [ 1, 0 ])
    ax2.set_ylabel('flux_gluc [mmol/sec]', color='blue')
    x1 = ax2.plot(t, flux_gluc, 'r')
    x2 = ax2.plot(t, fikm_gluc, 'g')
    x3 = ax2.plot(t, fekm_gluc, 'blue')
    ax2.legend([ x1 [ 0 ], x2 [ 0 ], x3 [ 0 ] ], [ 'Total Epithelial', 'Cellular', 'Junctional' ])

    ax3 = fig22.add_subplot(spec2 [ 2, 0 ])
    ax3.set_ylabel('flux_gluc [mmol/sec]', color='blue')
    x2 = ax3.plot(t, fikm_gluc, 'g')
    x1 = ax3.plot(t, gluc_mi_nagluc, 'r')
    ax3.legend([ x1 [ 0 ], x2 [ 0 ] ], [ 'na-gluc', 'Cellular' ])
    ax3.set_title('glucose fluxes')

    fig23 = plt.figure(constrained_layout=True)
    spec2 = gridspec.GridSpec(ncols=1, nrows=3, figure=fig23)
    ax1 = fig23.add_subplot(spec2 [ 0, 0 ])
    x1 = ax1.plot(t [ 0:-2 ], mylist [ 0:-2 ], 'r')
    ax1.set_ylabel('cm_na (mm)', color='blue')
    ax1.set_title('Luminal Sodium Concentration')
    ax2 = fig23.add_subplot(spec2 [ 1, 0 ])
    ax2.set_ylabel('flux_volume [mmol/sec]', color='blue')
    x1 = ax2.plot(t, flux_v, 'r')
    x2 = ax2.plot(t, fivm, 'g')
    x3 = ax2.plot(t, fevm, 'blue')
    ax2.legend([ x1 [ 0 ], x2 [ 0 ], x3 [ 0 ] ], [ 'Total Epithelial', 'Cellular', 'Junctional' ])

    ax3 = fig23.add_subplot(spec2 [ 2, 0 ])
    ax3.set_ylabel('p [mm.hg]', color='blue')
    p_m = [ 15 for i in range(T) ]
    x1 = ax3.plot(t, p_m, 'r')
    x2 = ax3.plot(t, p_i, 'b')
    x3 = ax3.plot(t, pe, 'g')
    ax3.legend([ x1 [ 0 ], x2 [ 0 ], x3 [ 0 ] ], [ 'Luminal ', 'Cellular', 'Junctional' ])
    ax3.set_title('Pressure')
    for ax in fig23.get_axes():
        ax.label_outer()
    plt.show()

elif glucoseimpact:
    t = 0
    guess = [ var [ t ] for var in vars ]

    from scipy import optimize

    # changing the cm_na concentration in step_wise
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mylist = [ 0.005 for i in range(T) ]
    pm = 15.000
    cm_cl = 0.11322499
    cm_na = 0.11
    stepwise = 1
    if stepwise:
        cmgluc = [ 0.0050, 0.0055, 0.006, 0.00650, 0.0070 ]
        while True:
            t += 1
            if t == T:
                break
            elif 0 < t <= T / len(cmgluc):
                cm_gluc = cmgluc [ 0 ]
                mylist [ t ] = cmgluc [ 0 ]
                result = scipy.optimize.root(eQs, np.array(guess), 1)
                guess = [ var [ t ] for var in vars ]
                # print('result=',result.x)
                result_flx = eQs(guess, 0)
                for i in range(len(result_flx)):
                    flx [ i ] [ t ] = result_flx [ i ]
            elif 1 * T / len(cmgluc) < t <= 2 * T / len(cmgluc):
                cm_gluc = cmgluc [ 1 ]
                mylist [ t ] = cmgluc [ 1 ]
                result = scipy.optimize.root(eQs, np.array(guess), 1)
                guess = [ var [ t ] for var in vars ]
                result_flx = eQs(guess, 0)
                for i in range(len(result_flx)):
                    flx [ i ] [ t ] = result_flx [ i ]
            elif 2 * T / len(cmgluc) < t <= 3 * T / len(cmgluc):
                cm_gluc = cmgluc [ 2 ]
                mylist [ t ] = cmgluc [ 2 ]
                result = scipy.optimize.root(eQs, np.array(guess), 1)
                guess = [ var [ t ] for var in vars ]
                result_flx = eQs(guess, 0)
                for i in range(len(result_flx)):
                    flx [ i ] [ t ] = result_flx [ i ]
            elif 3 * T / len(cmgluc) < t <= 4 * T / len(cmgluc):
                cm_gluc = cmgluc [ 3 ]
                mylist [ t ] = cmgluc [ 3 ]
                result = scipy.optimize.root(eQs, np.array(guess), 1)
                guess = [ var [ t ] for var in vars ]
                result_flx = eQs(guess, 0)
                for i in range(len(result_flx)):
                    flx [ i ] [ t ] = result_flx [ i ]
            elif 4 * T / len(cmgluc) < t <= 5 * T / len(cmgluc):
                cm_gluc = cmgluc [ 4 ]
                mylist [ t ] = cmgluc [ 4 ]
                result = scipy.optimize.root(eQs, np.array(guess), 1)
                guess = [ var [ t ] for var in vars ]
                result_flx = eQs(guess, 0)
                for i in range(len(result_flx)):
                    flx [ i ] [ t ] = result_flx [ i ]
        else:
            cm_na = 0.11000
            while True:
                t += 1
                if t == T:
                    break
                else:
                    cm_na = 0.11
                    result = scipy.optimize.root(eQs, np.array(guess), 1)
                    guess = [ var [ t ] for var in vars ]
                    result_flx = eQs(guess, 0)
                    for i in range(len(result_flx)):
                        flx [ i ] [ t ] = result_flx [ i ]

        t = np.linspace(t0, tf, T, tag)
        mylist = np.asarray(mylist)
        fig18 = plt.figure(constrained_layout=True)
        spec2 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig18)
        ax1 = fig18.add_subplot(spec2 [ 0, 0 ])
        x1 = ax1.plot(t [ 0:-2 ], mylist [ 0:-2 ], 'r-')
        x2 = ax1.plot(t, ci_na, 'blue')
        x3 = ax1.plot(t, ci_cl, 'g')
        x4 = ax1.plot(t, ci_gluc, 'k')
        ax1.legend([ x1 [ 0 ], x2 [ 0 ], x3 [ 0 ], x4 [ 0 ] ],
                   [ 'Luminal-gluc', 'Cellular-na', 'Cellular-cl', 'Cellular-gluc' ])
        ax1.set_ylabel('c (mm)', color='blue')
        ax1.set_title('Solute Concentration')
        ax1.set_xlabel('Time [sec]', color='blue')

        fig19 = plt.figure(constrained_layout=True)
        spec2 = gridspec.GridSpec(ncols=1, nrows=3, figure=fig19)
        ax1 = fig19.add_subplot(spec2 [ 2, 0 ])
        x1 = ax1.plot(t, fikm_na, 'g')
        x2 = ax1.plot(t, jnhe3_na, 'b')
        ax1.legend([ x1 [ 0 ], x2 [ 0 ] ], [ 'Cellular', 'nhe3' ])
        ax1.set_xlabel('Time [sec]', color='blue')
        ax1.set_title('Sodium Fluxes')
        ax2 = fig19.add_subplot(spec2 [ 1, 0 ])
        ax2.set_ylabel('flux_na [mmol/sec]', color='blue')
        x1 = ax2.plot(t, na_mi_nagluc, 'r')
        x2 = ax2.plot(t, na_mi_nah2po4, 'g')
        x3 = ax2.plot(t, jnhe3_na, 'b')
        ax2.legend([ x3 [ 0 ], x2 [ 0 ], x1 [ 0 ] ], [ 'nhe3', 'na/h2po4', 'na-gluc' ])
        ax2.set_title('sodium cotransporters on apical membrane')
        ax2 = fig19.add_subplot(spec2 [ 0, 0 ])
        # ax2.set_ylabel('flux_na [mmol/sec]', color='blue')
        x1 = ax2.plot(t, flux_na, 'r')
        x2 = ax2.plot(t, fikm_na, 'g')
        x3 = ax2.plot(t, fekm_na, 'blue')
        ax2.set_title('Sodium Fluxes')
        ax2.legend([ x1 [ 0 ], x2 [ 0 ], x3 [ 0 ] ], [ 'Total Epithelial', 'Cellular', 'Junctional' ])

        fig20 = plt.figure(constrained_layout=True)
        spec2 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig20)
        ax1 = fig20.add_subplot(spec2 [ 0, 0 ])
        x1 = ax1.plot(t [ 0:-2 ], mylist [ 0:-2 ], 'r')
        x2 = ax1.plot(t, ci_gluc, 'blue')
        ax1.legend([ x1 [ 0 ], x2 [ 0 ] ], [ 'Luminal_gluc', 'Cellular_gluc' ])
        ax1.set_ylabel('c (mm)', color='blue')
        ax1.set_title('Solute Concentration')

        fig21 = plt.figure(constrained_layout=True)
        spec2 = gridspec.GridSpec(ncols=1, nrows=3, figure=fig21)
        ax2 = fig21.add_subplot(spec2 [ 0, 0 ])
        # ax2.set_ylabel('flux_cl [mmol/sec]', color='blue')
        x1 = ax2.plot(t, flux_gluc, 'r')
        x2 = ax2.plot(t, fikm_gluc, 'g')
        x3 = ax2.plot(t, fekm_gluc, 'blue')
        ax2.legend([ x1 [ 0 ], x2 [ 0 ], x3 [ 0 ] ], [ 'Total Epithelial', 'Cellular', 'Junctional' ])
        ax2.set_title('glucose fluxes')
        ax1 = fig21.add_subplot(spec2 [ 2, 0 ])
        x1 = ax1.plot(t, gluc_mi_nagluc, 'r')
        x2 = ax1.plot(t, fikm_gluc, 'g')
        x3 = ax1.plot(t, na_mi_nagluc, 'b')
        ax1.set_xlabel('Time [sec]', color='blue')
        ax1.set_title('glucose fluxes')
        ax1.legend([ x1 [ 0 ], x2 [ 0 ], x3 [ 0 ] ], [ 'gluc-na-gluc', 'Cellular', 'na-gluc' ])
        ax2 = fig21.add_subplot(spec2 [ 1, 0 ])
        ax2.set_ylabel('flux_cl [mmol/sec]', color='blue')
        # ax2.set_xlabel('Time [sec]', color='blue')
        x1 = ax2.plot(t, cl_mi_clhco2, 'r')
        x2 = ax2.plot(t, cl_mi_clhco3, 'b')
        ax2.legend([ x1 [ 0 ], x2 [ 0 ] ], [ 'cl/hco2', 'cl/hco3' ])
        ax2.set_title('Chloride Cotransporters on Apical Membrane')

        fig22 = plt.figure(constrained_layout=True)
        spec2 = gridspec.GridSpec(ncols=1, nrows=3, figure=fig22)
        ax1 = fig22.add_subplot(spec2 [ 0, 0 ])
        x1 = ax1.plot(t [ 0:-2 ], mylist [ 0:-2 ], 'r')
        x2 = ax1.plot(t, ci_gluc, 'blue')
        ax1.legend([ x1 [ 0 ], x2 [ 0 ] ], [ 'Luminal_cl', 'Cellular_gluc' ])
        ax1.set_ylabel('c (mm)', color='blue')
        ax1.set_title('Solute Concentration')
        ax2 = fig22.add_subplot(spec2 [ 1, 0 ])
        ax2.set_ylabel('flux_gluc [mmol/sec]', color='blue')
        x1 = ax2.plot(t, flux_gluc, 'r')
        x2 = ax2.plot(t, fikm_gluc, 'g')
        x3 = ax2.plot(t, fekm_gluc, 'blue')
        ax2.legend([ x1 [ 0 ], x2 [ 0 ], x3 [ 0 ] ], [ 'Total Epithelial', 'Cellular', 'Junctional' ])

        ax3 = fig22.add_subplot(spec2 [ 2, 0 ])
        ax3.set_ylabel('flux_gluc [mmol/sec]', color='blue')
        x2 = ax3.plot(t, fikm_gluc, 'g')
        x1 = ax3.plot(t, gluc_mi_nagluc, 'r')
        ax3.legend([ x1 [ 0 ], x2 [ 0 ] ], [ 'na-gluc', 'Cellular' ])
        ax3.set_title('glucose fluxes')
        ax3.set_xlabel('Time [sec]', color='blue')
        fig23 = plt.figure(constrained_layout=True)
        spec2 = gridspec.GridSpec(ncols=1, nrows=3, figure=fig23)
        ax1 = fig23.add_subplot(spec2 [ 0, 0 ])
        x1 = ax1.plot(t [ 0:-2 ], mylist [ 0:-2 ], 'r')
        ax1.set_ylabel('cm_gluc (mm)', color='blue')
        ax1.set_title('Luminal glucose concentration')
        ax2 = fig23.add_subplot(spec2 [ 1, 0 ])
        ax2.set_ylabel('flux_volume [mmol/sec]', color='blue')
        x1 = ax2.plot(t, flux_v, 'r')
        x2 = ax2.plot(t, fivm, 'g')
        x3 = ax2.plot(t, fevm, 'blue')
        ax2.legend([ x1 [ 0 ], x2 [ 0 ], x3 [ 0 ] ], [ 'Total Epithelial', 'Cellular', 'Junctional' ])

        ax3 = fig23.add_subplot(spec2 [ 2, 0 ])
        ax3.set_ylabel('p [mm.hg]', color='blue')
        p_m = [ 15 for i in range(T) ]
        x1 = ax3.plot(t, p_m, 'r')
        x2 = ax3.plot(t, p_i, 'b')
        x3 = ax3.plot(t, pe, 'g')
        ax3.legend([ x1 [ 0 ], x2 [ 0 ], x3 [ 0 ] ], [ 'Luminal ', 'Cellular', 'Junctional' ])
        ax3.set_title('Pressure')
        ax3.set_xlabel('Time [sec]', color='blue')
        plt.show()

elif Pressureimpact:
    t = 0
    guess = [ var [ t ] for var in vars ]
    from scipy import optimize

    # changing the luminal Pressure in step_wise
    mylist = [ 15 for i in range(T) ]
    # pm = 15.000
    cm_cl = 0.11322499
    cm_gluc = 0.005
    cm_na = 0.1

    stepwise = 1
    if stepwise:
        cmgluc = [ 15, 20, 25, 30, 35 ]
        while True:
            t += 1
            if t == T:
                break
            elif 0 < t <= T / len(cmgluc):
                pm = cmgluc [ 0 ]
                mylist [ t ] = cmgluc [ 0 ]
                result = scipy.optimize.root(eQs, np.array(guess), 1)
                guess = [ var [ t ] for var in vars ]
                # print('result=',result.x)
                result_flx = eQs(guess, 0)
                for i in range(len(result_flx)):
                    flx [ i ] [ t ] = result_flx [ i ]
            elif 1 * T / len(cmgluc) < t <= 2 * T / len(cmgluc):
                pm = cmgluc [ 1 ]
                mylist [ t ] = cmgluc [ 1 ]
                result = scipy.optimize.root(eQs, np.array(guess), 1)
                guess = [ var [ t ] for var in vars ]
                result_flx = eQs(guess, 0)
                for i in range(len(result_flx)):
                    flx [ i ] [ t ] = result_flx [ i ]
            elif 2 * T / len(cmgluc) < t <= 3 * T / len(cmgluc):
                pm = cmgluc [ 2 ]
                mylist [ t ] = cmgluc [ 2 ]
                result = scipy.optimize.root(eQs, np.array(guess), 1)
                guess = [ var [ t ] for var in vars ]
                result_flx = eQs(guess, 0)
                for i in range(len(result_flx)):
                    flx [ i ] [ t ] = result_flx [ i ]
            elif 3 * T / len(cmgluc) < t <= 4 * T / len(cmgluc):
                pm = cmgluc [ 3 ]
                mylist [ t ] = cmgluc [ 3 ]
                result = scipy.optimize.root(eQs, np.array(guess), 1)
                guess = [ var [ t ] for var in vars ]
                result_flx = eQs(guess, 0)
                for i in range(len(result_flx)):
                    flx [ i ] [ t ] = result_flx [ i ]
            elif 4 * T / len(cmgluc) < t <= 5 * T / len(cmgluc):
                pm = cmgluc [ 4 ]
                mylist [ t ] = cmgluc [ 4 ]
                result = scipy.optimize.root(eQs, np.array(guess), 1)
                guess = [ var [ t ] for var in vars ]
                result_flx = eQs(guess, 0)
                for i in range(len(result_flx)):
                    flx [ i ] [ t ] = result_flx [ i ]
            # elif 0 < t <= T / len(cmgluc):
            #     pm = cmgluc [ 0 ]
            #     mylist [ t ] = cmgluc [ 0 ]
            #     result = scipy.optimize.root(eQs, np.array(guess), 1)
            #     guess = [ var [ t ] for var in vars ]
            #     # print('result=',result.x)
            #     result_flx = eQs(guess, 0)
            #     for i in range(len(result_flx)):
            #         flx [ i ] [ t ] = result_flx [ i ]
    else:
        cm_na = 0.14000
        while True:
            t += 1
            if t == T:
                break
            else:
                cm_na = 0.14
                result = scipy.optimize.root(eQs, np.array(guess), 1)
                guess = [ var [ t ] for var in vars ]
                result_flx = eQs(guess, 0)
                for i in range(len(result_flx)):
                    flx [ i ] [ t ] = result_flx [ i ]

    t = np.linspace(t0, tf, T, tag)
    mylist = np.asarray(mylist)
    fig18 = plt.figure(constrained_layout=True)
    spec2 = gridspec.GridSpec(ncols=1, nrows=4, figure=fig18)
    ax1 = fig18.add_subplot(spec2 [ 0, 0 ])
    x1 = ax1.plot(t [ 0:-2 ], mylist [ 0:-2 ], 'r-')
    ax1.set_title('Luminal Pressure')
    ax1.set_ylabel('pm (mm.Hg)', color='blue')

    ax2 = fig18.add_subplot(spec2 [ 1, 0 ])
    ax2.set_title('Junctional-Solute Concentration')
    x1 = ax2.plot(t, ce_na, 'b')
    ax2.legend([ x1 [ 0 ] ], [ 'ce_na' ])
    # ax2.set_ylim([ 0.1275, 0.129 ])

    ax3 = fig18.add_subplot(spec2 [ 2, 0 ])
    x2 = ax3.plot(t, ce_cl, 'b')
    ax3.set_ylabel('c (mm)', color='blue')
    ax3.legend([ x2 [ 0 ] ], [ 'ce_cl' ])
    # ax3.set_ylim([ 0.1025, 0.1035 ])
    ax4 = fig18.add_subplot(spec2 [ 3, 0 ])
    x3 = ax4.plot(t, ce_gluc, 'b')
    ax4.set_xlabel('Time [sec]', color='blue')
    ax4.legend([ x3 [ 0 ] ], [ 'ce_gluc' ])
    # ax4.set_ylim([ 0.0057, 0.006 ])

    fig18 = plt.figure(constrained_layout=True)
    spec2 = gridspec.GridSpec(ncols=1, nrows=4, figure=fig18)
    ax1 = fig18.add_subplot(spec2 [ 0, 0 ])
    x1 = ax1.plot(t [ 0:-2 ], mylist [ 0:-2 ], 'r-')
    ax1.set_title('Luminal Pressure')
    ax1.set_ylabel('pm (mm.Hg)', color='blue')
    ax2 = fig18.add_subplot(spec2 [ 1, 0 ])
    ax2.set_title('Cellular-Solute Concentration')
    x2 = ax2.plot(t, ci_na, 'blue')
    ax2.legend([ x2 [ 0 ] ], [ 'ci-na' ])
    # ax2.set_ylim([ 0.024, 0.025])
    ax3 = fig18.add_subplot(spec2 [ 2, 0 ])
    x3 = ax3.plot(t, ci_cl, 'b')
    ax3.set_ylabel('c (mm)', color='blue')
    ax3.legend([ x3 [ 0 ] ], [ 'ci_cl' ])
    # ax3.set_ylim([ 0.003, 0.004 ])
    ax4 = fig18.add_subplot(spec2 [ 3, 0 ])
    x4 = ax4.plot(t, ci_gluc, 'b')
    # ax4.set_ylim([ 0.0111, 0.01125])
    ax4.set_xlabel('Time [sec]', color='blue')
    ax4.legend([ x4 [ 0 ] ], [ 'ci_gluc' ])
    fig19 = plt.figure(constrained_layout=True)
    spec2 = gridspec.GridSpec(ncols=1, nrows=3, figure=fig19)
    ax1 = fig19.add_subplot(spec2 [ 2, 0 ])
    x1 = ax1.plot(t, fikm_na, 'g')
    x2 = ax1.plot(t, jnhe3_na, 'b')
    # ax1.set_ylim([ 0.0000045, 0.000007])
    ax1.legend([ x1 [ 0 ], x2 [ 0 ] ], [ 'Cellular', 'NHE3' ])
    ax1.set_xlabel('Time [sec]', color='blue')
    ax1.set_title('Sodium Fluxes')
    ax2 = fig19.add_subplot(spec2 [ 1, 0 ])
    ax2.set_ylabel('flux_na [mmol/sec]', color='blue')
    x1 = ax2.plot(t, na_mi_nagluc, 'r')
    x2 = ax2.plot(t, na_mi_nah2po4, 'g')
    x3 = ax2.plot(t, jnhe3_na, 'b')
    ax2.legend([ x3 [ 0 ], x2 [ 0 ], x1 [ 0 ] ], [ 'NHE3', 'na/h2po4', 'na-gluc' ])
    ax2.set_title('Sodium Cotransporters on Apical Membrane')
    ax2 = fig19.add_subplot(spec2 [ 0, 0 ])
    # ax2.set_ylabel('flux_na [mmol/sec]', color='blue')
    x1 = ax2.plot(t, flux_na, 'r')
    x2 = ax2.plot(t, fikm_na, 'g')
    x3 = ax2.plot(t, fekm_na, 'blue')
    ax2.set_title('Sodium Fluxes')
    ax2.legend([ x1 [ 0 ], x2 [ 0 ], x3 [ 0 ] ], [ 'Total Epithelial', 'Cellular', 'Junctional' ])

    # fig20 = plt.figure(constrained_layout=True)
    # spec2 = gridspec.GridSpec(ncols=1, nrows=1, figure=fig20)
    # ax1 = fig20.add_subplot(spec2 [ 0, 0 ])
    # ax1.plot(t [ 0:-2 ], mylist [ 0:-2 ], 'r')
    # #ax1.legend([ x1 [ 0 ], x2 [ 0 ] ], [ 'Luminal_gluc', 'Cellular_gluc' ])
    # ax1.set_ylabel('pm [mm.Hg]', color='blue')
    # ax1.set_title('Luminal Pressure')

    fig21 = plt.figure(constrained_layout=True)
    spec2 = gridspec.GridSpec(ncols=1, nrows=3, figure=fig21)
    ax1 = fig21.add_subplot(spec2 [ 0, 0 ])
    ax1.plot(t [ 0:-2 ], mylist [ 0:-2 ], 'r')
    # x2 = ax1.plot(t, ci_gluc, 'blue')
    # ax1.legend([ x1 [ 0 ], x2 [ 0 ] ], [ 'Luminal_gluc', 'Cellular_gluc' ])
    ax1.set_ylabel('pm [mm.Hg]', color='blue')
    ax1.set_title('Luminal Pressure')

    ax2 = fig21.add_subplot(spec2 [ 2, 0 ])
    x1 = ax2.plot(t, flux_cl, 'r')
    x2 = ax2.plot(t, fikm_cl, 'g')
    x3 = ax2.plot(t, fekm_cl, 'blue')
    ax2.set_title('Chloride Fluxes (due to changes of luminal Pressure)')
    ax2.legend([ x1 [ 0 ], x2 [ 0 ], x3 [ 0 ] ], [ 'Total Epithelial', 'Cellular', 'Junctional' ])
    ax2.set_xlabel('Time [sec]', color='blue')

    # ax1.set_title('chloride fluxes')
    # ax1.legend([ x1 [ 0 ], x2 [ 0 ], x3 [ 0 ] ], [ 'gluc-na-gluc', 'Cellular', 'na-gluc' ])
    ax3 = fig21.add_subplot(spec2 [ 1, 0 ])
    ax3.set_ylabel('flux_cl [mmol/sec]', color='blue')
    # ax2.set_xlabel('Time [sec]', color='blue')
    x1 = ax3.plot(t, cl_mi_clhco2, 'r')
    x2 = ax3.plot(t, cl_mi_clhco3, 'b')
    ax3.legend([ x1 [ 0 ], x2 [ 0 ] ], [ 'cl/hco2', 'cl/hco3' ])
    ax3.set_title('Chloride Cotransporters on Apical Membrane')

    fig22 = plt.figure(constrained_layout=True)
    spec2 = gridspec.GridSpec(ncols=1, nrows=3, figure=fig22)
    ax1 = fig22.add_subplot(spec2 [ 0, 0 ])
    ax1.plot(t [ 0:-2 ], mylist [ 0:-2 ], 'r')
    # x2 = ax1.plot(t, ci_gluc, 'blue')
    # ax1.legend([ x1 [ 0 ], x2 [ 0 ] ], [ 'Luminal_gluc', 'Cellular_gluc' ])
    ax1.set_ylabel('pm [mm.Hg]', color='blue')
    ax1.set_title('Luminal Pressure')
    ax2 = fig22.add_subplot(spec2 [ 1, 0 ])
    ax2.set_ylabel('flux_gluc [mmol/sec]', color='blue')
    x1 = ax2.plot(t, flux_gluc, 'r')
    x2 = ax2.plot(t, fikm_gluc, 'g')
    x3 = ax2.plot(t, fekm_gluc, 'blue')
    ax2.legend([ x1 [ 0 ], x2 [ 0 ], x3 [ 0 ] ], [ 'Total Epithelial', 'Cellular', 'Junctional' ])
    ax2.set_title('Glucose Fluxes')
    ax3 = fig22.add_subplot(spec2 [ 2, 0 ])
    ax3.set_ylabel('flux_gluc [mmol/sec]', color='blue')
    x2 = ax3.plot(t, fikm_gluc, 'g')
    x1 = ax3.plot(t, gluc_mi_nagluc, 'r')
    ax3.legend([ x1 [ 0 ], x2 [ 0 ] ], [ 'na-gluc', 'Cellular' ])
    ax3.set_title('glucose fluxes')
    ax3.set_xlabel('Time [sec]', color='blue')
    # ax3.set_ylim([0.00000144, 0.00000155])

    fig23 = plt.figure(constrained_layout=True)
    spec2 = gridspec.GridSpec(ncols=1, nrows=3, figure=fig23)
    ax1 = fig23.add_subplot(spec2 [ 0, 0 ])
    x1 = ax1.plot(t [ 0:-2 ], mylist [ 0:-2 ], 'r')
    ax1.set_ylabel('p [mm.Hg]', color='blue')
    ax1.set_title('Luminal and Cellular Pressure ')

    ax3 = fig23.add_subplot(spec2 [ 1, 0 ])
    ax3.set_ylabel('pe [mm.Hg]', color='blue')
    # ax3.set_ylim([ -18,-17])
    ax3.plot(t, pe, 'g')
    ax3.set_title('Junctional Pressure')

    ax2 = fig23.add_subplot(spec2 [ 2, 0 ])
    ax2.set_ylabel('Flux_v [mmol/sec]', color='blue')
    x1 = ax2.plot(t, flux_v, 'r')
    x2 = ax2.plot(t, fivm, 'g')
    x3 = ax2.plot(t, fevm, 'blue')
    ax2.set_title('Flux_Volume ')
    ax2.legend([ x1 [ 0 ], x2 [ 0 ], x3 [ 0 ] ], [ 'Total Epithelial', 'Cellular', 'Junctional' ])
    ax2.set_xlabel('Time [sec]', color='blue')


    # plt.show()
    def format_axes(fig18):
        for i, ax in enumerate(fig18.axes):
            ax.text(0.5, 0.5, "ax%d" % (i + 1), va="center", ha="center")
            ax.tick_params(labelbottom=False, labelleft=True, labelsize=10)


    fig18 = plt.figure(constrained_layout=False, figsize=(6, 6))
    spec2 = gridspec.GridSpec(ncols=1, nrows=5, hspace=0, figure=fig18)
    anno_opts = dict(xy=(0.78, 0.1), xycoords='axes fraction',
                     va='center', ha='center')

    ax1 = fig18.add_subplot(spec2 [ 0, 0 ])
    x1 = ax1.plot(t [ 0:-2 ], mylist [ 0:-2 ], 'r-')
    ax1.annotate('Luminal/ Cellular Pressure', **anno_opts)
    # ax1.set_ylabel('p [mm.Hg]', color='blue', fontsize=12)
    ax1.tick_params(labelbottom=True, labelleft=True, labelsize=10)

    ax2 = fig18.add_subplot(spec2 [ 1, 0 ])
    x2 = ax2.plot(t, ci_na, 'blue')
    ax2.tick_params(labelbottom=False, labelleft=True, labelsize=10)
    # ax2.set_ylabel('c_na [mm]', color='blue', fontsize=12)
    ax2.annotate('Cellular ', **anno_opts)
    ax3 = fig18.add_subplot(spec2 [ 2, 0 ])
    x2 = ax3.plot(t, ce_na, 'blue')
    # ax3.set_ylabel('c_na [mm]', color='blue', fontsize=12)
    ax3.tick_params(labelbottom=False, labelleft=True, labelsize=14)
    ax3.annotate('Junctional', **anno_opts)
    ax3.set_ylim([ 0.12, 0.138 ])
    ax3.tick_params(labelbottom=False, labelleft=True, labelsize=10)
    ax4 = fig18.add_subplot(spec2 [ 3, 0 ])
    # ax4.set_ylabel('f_na [mmol/sec]', color='blue')
    fikm_nas = [ phi_scale(x, 1e5) for x in fikm_na ]
    fekm_nas = [ phi_scale(x, 1e5) for x in fekm_na ]
    x2 = ax4.plot(t, fikm_nas, 'r')
    x3 = ax4.plot(t, fekm_nas, 'blue')
    ax4.set_ylim([ 0.0, 1.25 ])
    ax4.legend([ x2 [ 0 ], x3 [ 0 ] ], [ 'Cellular-Na Fluxes', 'Junctional-Na Fluxes' ])
    ax4.tick_params(labelbottom=True, labelleft=True, labelsize=10)
    ax2 = fig18.add_subplot(spec2 [ 4, 0 ])
    # ax2.set_ylabel('      flux$*1e5$ [mmol/sec]', color='blue', fontsize=12)
    ax2.tick_params(labelbottom=False, labelleft=True, labelsize=10)
    fivms = [ phi_scale(x, 1e5) for x in fivm ]
    fevms = [ phi_scale(x, 1e5) for x in fevm ]
    x2 = ax2.plot(t, fivms, 'g')
    x3 = ax2.plot(t, fevms, 'blue')
    ax2.legend([ x2 [ 0 ], x3 [ 0 ] ], [ 'Cellular-Volume', 'Junctional-Volume' ])

    # ax4.annotate('Sodium Fluxes', **anno_opts)
    # ax2.set_title('Sodium Fluxes')
    # ax5 = fig18.add_subplot(spec2 [ 5, 0 ])
    # ax5.set_ylabel('            f_na$1e5$ [mmol/sec]', color='blue', fontsize=12)

    # na_mi_naglucs = [ phi_scale(x, 1e5) for x in na_mi_nagluc ]
    # na_mi_nah2po4s = [ phi_scale(x, 1e5) for x in na_mi_nah2po4 ]
    # jnhe3_nas = [ phi_scale(x, 1e5) for x in jnhe3_na ]
    # x1 = ax5.plot(t, na_mi_naglucs, 'r')
    # x2 = ax5.plot(t, na_mi_nah2po4s, 'b')
    # x3 = ax5.plot(t, jnhe3_nas, 'g')
    # x1 = ax5.plot(t, vi, 'b')
    # x2 = ax5.plot(t, ve, 'b')
    # ax5.legend([ x1 [ 0 ], x2[ 0 ] ], [ 'Cellular-Voltage', 'Junctional-Voltage' ])
    # ax5.set_ylabel('V [mV]', color='blue', fontsize=12)

    # ax5.annotate('Cotransporters', **anno_opts)
    # ax5.legend([ x3 [ 0 ], x2 [ 0 ], x1 [ 0 ] ], [ 'nhe3', 'na/h2po4', 'na-gluc' ])
    # ax5.set_xlabel('Time [sec]', color='blue', fontsize=12)
    # ax5.annotate('Cellular-Voltage', **anno_opts)
    # ax2.set_title('Sodium Cotransporters on Apical Membrane')
    # ax5.tick_params(labelbottom=True, labelleft=True, labelsize=10)
    plt.show()



elif chlorideimpact:
    t = 0
    guess = [ var [ t ] for var in vars ]
    from scipy import optimize

    # changing the cm_cl concentration in step_wise
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mylist = [ 0.11322499 for i in range(T) ]
    pm = 15.000
    # cm_cl = 0.11322499
    cm_gluc = 0.005
    cm_na = 0.14
    t = 0
    stepwise = 1
    if stepwise:
        cmgluc = [ 0.11322499, 0.11822499, 0.12322499, 0.12822499, 0.13322499 ]
        while True:
            t += 1
            if t == T:
                break
            elif 0 < t <= T / len(cmgluc):
                cm_cl = cmgluc [ 0 ]
                mylist [ t ] = cmgluc [ 0 ]
                result = scipy.optimize.root(eQs, np.array(guess), 1)
                guess = [ var [ t ] for var in vars ]
                # print('result=',result.x)
                result_flx = eQs(guess, 0)
                for i in range(len(result_flx)):
                    flx [ i ] [ t ] = result_flx [ i ]
            elif 1 * T / len(cmgluc) < t <= 2 * T / len(cmgluc):
                cm_cl = cmgluc [ 1 ]
                mylist [ t ] = cmgluc [ 1 ]
                result = scipy.optimize.root(eQs, np.array(guess), 1)
                guess = [ var [ t ] for var in vars ]
                result_flx = eQs(guess, 0)
                for i in range(len(result_flx)):
                    flx [ i ] [ t ] = result_flx [ i ]
            elif 2 * T / len(cmgluc) < t <= 3 * T / len(cmgluc):
                cm_cl = cmgluc [ 2 ]
                mylist [ t ] = cmgluc [ 2 ]
                result = scipy.optimize.root(eQs, np.array(guess), 1)
                guess = [ var [ t ] for var in vars ]
                result_flx = eQs(guess, 0)
                for i in range(len(result_flx)):
                    flx [ i ] [ t ] = result_flx [ i ]
            elif 3 * T / len(cmgluc) < t <= 4 * T / len(cmgluc):
                cm_cl = cmgluc [ 3 ]
                mylist [ t ] = cmgluc [ 3 ]
                result = scipy.optimize.root(eQs, np.array(guess), 1)
                guess = [ var [ t ] for var in vars ]
                result_flx = eQs(guess, 0)
                for i in range(len(result_flx)):
                    flx [ i ] [ t ] = result_flx [ i ]
            elif 4 * T / len(cmgluc) < t <= 5 * T / len(cmgluc):
                cm_cl = cmgluc [ 4 ]
                mylist [ t ] = cmgluc [ 4 ]
                result = scipy.optimize.root(eQs, np.array(guess), 1)
                guess = [ var [ t ] for var in vars ]
                result_flx = eQs(guess, 0)
                for i in range(len(result_flx)):
                    flx [ i ] [ t ] = result_flx [ i ]
            # elif 5 * T / len(cmgluc) < t <= T / len(cmgluc):
            #     cm_cl = cmgluc [ 0 ]
            #     mylist [ t ] = cmgluc [ 0 ]
            #     result = scipy.optimize.root(eQs, np.array(guess), 1)
            #     guess = [ var [ t ] for var in vars ]
            #     # print('result=',result.x)
            #     result_flx = eQs(guess, 0)
            #     for i in range(len(result_flx)):
            #         flx [ i ] [ t ] = result_flx [ i ]
    else:
        cm_na = 0.14000
        while True:
            t += 1
            if t == T:
                break
            else:
                cm_cl = 0.11322499
                result = scipy.optimize.root(eQs, np.array(guess), 1)
                guess = [ var [ t ] for var in vars ]
                result_flx = eQs(guess, 0)
                for i in range(len(result_flx)):
                    flx [ i ] [ t ] = result_flx [ i ]
    #
    t = np.linspace(t0, tf, T, tag)
    mylist = np.asarray(mylist)
    fig18 = plt.figure(constrained_layout=True)
    spec2 = gridspec.GridSpec(ncols=1, nrows=4, figure=fig18)
    ax1 = fig18.add_subplot(spec2 [ 0, 0 ])
    x1 = ax1.plot(t [ 0:-2 ], mylist [ 0:-2 ], 'r-')
    ax1.set_title('Luminal-chloride concentration')
    ax1.set_ylabel('cm_cl [mm]', color='blue')
    ax3 = fig18.add_subplot(spec2 [ 1, 0 ])
    x2 = ax3.plot(t, ce_cl, 'b')
    ax3.set_title('Junctional-Solute Concentration')
    ax3.legend([ x2 [ 0 ] ], [ 'ce_cl' ])
    # ax3.set_ylim([ 0.0983, 0.1 ])
    ax2 = fig18.add_subplot(spec2 [ 2, 0 ])
    x1 = ax2.plot(t, ce_na, 'b')
    ax2.legend([ x1 [ 0 ] ], [ 'ce_na' ])
    ax2.set_ylabel('c (mm)', color='blue')
    # ax2.set_ylim([ 0.1245, 0.126 ])
    ax4 = fig18.add_subplot(spec2 [ 3, 0 ])
    x3 = ax4.plot(t, ce_gluc, 'b')
    ax4.set_xlabel('Time [sec]', color='blue')
    ax4.legend([ x3 [ 0 ] ], [ 'ce_gluc' ])
    # ax4.set_ylim([ 0.0048, 0.006 ])

    fig18 = plt.figure(constrained_layout=True)
    spec2 = gridspec.GridSpec(ncols=1, nrows=4, figure=fig18)
    ax1 = fig18.add_subplot(spec2 [ 0, 0 ])
    x1 = ax1.plot(t [ 0:-2 ], mylist [ 0:-2 ], 'r-')
    ax1.set_title('Luminal-Chloride Concentration')
    ax1.set_ylabel('cm_cl [mm]', color='blue')
    # ax2.set_ylim([ 0.027, 0.035])
    ax3 = fig18.add_subplot(spec2 [ 1, 0 ])
    ax3.set_title('Cellular-Solute Concentration')
    x3 = ax3.plot(t, ci_cl, 'b')
    ax3.legend([ x3 [ 0 ] ], [ 'ci_cl' ])
    # ax3.set_ylim([ 0.012, 0.02 ])
    ax2 = fig18.add_subplot(spec2 [ 2, 0 ])
    x2 = ax2.plot(t, ci_na, 'blue')
    ax2.legend([ x2 [ 0 ] ], [ 'ci-na' ])
    ax2.set_ylabel('c (mm)', color='blue')
    ax4 = fig18.add_subplot(spec2 [ 3, 0 ])
    x4 = ax4.plot(t, ci_gluc, 'b')
    # ax4.set_ylim([ 0.009, 0.012])
    ax4.set_xlabel('Time [sec]', color='blue')
    ax4.legend([ x4 [ 0 ] ], [ 'ci_gluc' ])

    fig19 = plt.figure(constrained_layout=True)
    spec2 = gridspec.GridSpec(ncols=1, nrows=3, figure=fig19)
    ax2 = fig19.add_subplot(spec2 [ 0, 0 ])
    # ax2.set_ylabel('flux_na [mmol/sec]', color='blue')
    x1 = ax2.plot(t, flux_na, 'r')
    x2 = ax2.plot(t, fikm_na, 'g')
    x3 = ax2.plot(t, fekm_na, 'blue')
    ax2.set_title('Sodium Fluxes')
    ax2.legend([ x1 [ 0 ], x2 [ 0 ], x3 [ 0 ] ], [ 'Total Epithelial', 'Cellular', 'Junctional' ])
    ax2 = fig19.add_subplot(spec2 [ 1, 0 ])
    ax2.set_ylabel('flux_na [mmol/sec]', color='blue')
    x1 = ax2.plot(t, na_mi_nagluc, 'r')
    x2 = ax2.plot(t, na_mi_nah2po4, 'g')
    x3 = ax2.plot(t, jnhe3_na, 'b')
    ax2.legend([ x3 [ 0 ], x2 [ 0 ], x1 [ 0 ] ], [ 'NHE3', 'na/h2po4', 'na-gluc' ])
    ax2.set_title('Sodium Cotransporters on Apical Membrane')
    ax1 = fig19.add_subplot(spec2 [ 2, 0 ])
    x1 = ax1.plot(t, fikm_na, 'g')
    x2 = ax1.plot(t, jnhe3_na, 'b')
    ax1.legend([ x1 [ 0 ], x2 [ 0 ] ], [ 'Cellular', 'NHE3' ])
    ax1.set_xlabel('Time [sec]', color='blue')
    ax1.set_title('Sodium Fluxes')

    fig21 = plt.figure(constrained_layout=True)
    spec2 = gridspec.GridSpec(ncols=1, nrows=4, figure=fig21)
    ax1 = fig21.add_subplot(spec2 [ 0, 0 ])
    ax1.plot(t [ 0:-2 ], mylist [ 0:-2 ], 'r')
    ax1.set_ylabel('p (mm.Hg)', color='blue')
    ax1.set_title('Luminal Pressure [due to changes of cm_cl]')
    ax2 = fig21.add_subplot(spec2 [ 1, 0 ])
    x1 = ax2.plot(t, flux_cl, 'r')
    x2 = ax2.plot(t, fikm_cl, 'g')
    x3 = ax2.plot(t, fekm_cl, 'blue')
    ax2.set_title('Chloride Fluxes')
    ax2.legend([ x1 [ 0 ], x2 [ 0 ], x3 [ 0 ] ], [ 'Total Epithelial', 'Cellular', 'Junctional' ])
    ax3 = fig21.add_subplot(spec2 [ 2, 0 ])
    ax3.set_ylabel('flux_cl [mmol/sec]', color='blue')
    # ax2.set_xlabel('Time [sec]', color='blue')
    x1 = ax3.plot(t, cl_mi_clhco2, 'r')
    x2 = ax3.plot(t, cl_mi_clhco3, 'b')
    ax3.legend([ x1 [ 0 ], x2 [ 0 ] ], [ 'cl/hco2', 'cl/hco3' ])
    ax3.set_title('Chloride Cotransporters on Apical Membrane')
    ax2 = fig21.add_subplot(spec2 [ 3, 0 ])
    x1 = ax2.plot(t, cl_mi_clhco2, 'r')
    x3 = ax2.plot(t, fekm_cl, 'blue')
    ax2.set_title('Chloride Fluxes')
    ax2.legend([ x1 [ 0 ], x2 [ 0 ] ], [ 'cl/hco2', 'Junctional' ])
    ax2.set_xlabel('Time [sec]', color='blue')

    fig22 = plt.figure(constrained_layout=True)
    spec2 = gridspec.GridSpec(ncols=1, nrows=3, figure=fig22)
    ax1 = fig22.add_subplot(spec2 [ 0, 0 ])
    ax1.plot(t [ 0:-2 ], mylist [ 0:-2 ], 'r')
    ax1.set_ylabel('p (mm.hg)', color='blue')
    ax1.set_title('Luminal Pressure [due to changes of cm_cl]')
    ax2 = fig22.add_subplot(spec2 [ 1, 0 ])
    ax2.set_ylabel('flux_gluc [mmol/sec]', color='blue')
    x1 = ax2.plot(t, flux_gluc, 'r')
    x2 = ax2.plot(t, fikm_gluc, 'g')
    x3 = ax2.plot(t, fekm_gluc, 'blue')
    ax2.legend([ x1 [ 0 ], x2 [ 0 ], x3 [ 0 ] ], [ 'Total Epithelial', 'Cellular', 'Junctional' ])
    ax2.set_title('glucose fluxes')
    ax3 = fig22.add_subplot(spec2 [ 2, 0 ])
    ax3.set_ylabel('flux_gluc [mmol/sec]', color='blue')
    x2 = ax3.plot(t, fikm_gluc, 'g')
    x1 = ax3.plot(t, gluc_mi_nagluc, 'r')
    ax3.legend([ x1 [ 0 ], x2 [ 0 ] ], [ 'na-gluc', 'Cellular' ])
    # ax3.set_title('glucose fluxes')
    ax3.set_xlabel('Time [sec]', color='blue')

    fig23 = plt.figure(constrained_layout=True)
    spec2 = gridspec.GridSpec(ncols=1, nrows=3, figure=fig23)
    ax1 = fig23.add_subplot(spec2 [ 0, 0 ])
    x1 = ax1.plot(t [ 0:-2 ], mylist [ 0:-2 ], color='r', label='Luminal Pressure')
    ax1.set_ylabel('p [mm.Hg]', color='blue')
    ax1.set_title('Luminal Pressure')
    ax1.set_Text(0.5, 1.0, 'A single plot')
    ax3 = fig23.add_subplot(spec2 [ 1, 0 ])
    ax3.set_ylabel('p [mm.Hg]', color='blue')
    ax3.plot(t, pe, 'g')
    ax3.set_title('Junctional Pressure [due to changes of cm_cl]')
    ax2 = fig23.add_subplot(spec2 [ 2, 0 ])
    ax2.set_ylabel('flux_v [mmol/sec]', color='blue')
    x1 = ax2.plot(t, flux_v, 'r')
    x2 = ax2.plot(t, fivm, 'g')
    x3 = ax2.plot(t, fevm, 'blue')
    ax2.set_title('Flux_Volume ')
    ax2.legend([ x1 [ 0 ], x2 [ 0 ], x3 [ 0 ] ], [ 'Total Epithelial', 'Cellular', 'Junctional' ])
    ax2.set_xlabel('Time [sec]', color='blue')
    plt.show()

    # def format_axes(fig18):
    #     for i, ax in enumerate(fig18.axes):
    #         ax.text(0.5, 0.5, "ax%d" % (i + 1), va="center", ha="center")
    #         ax.tick_params(labelbottom=False, labelleft=True, labelsize=10)

    # fig18 = plt.figure(constrained_layout=False, figsize=(6, 9))
    # spec2 = gridspec.GridSpec(ncols=1, nrows=7, hspace=0, figure=fig18)
    # anno_opts = dict(xy=(0.8, 0.1), xycoords='axes fraction',
    #                  va='center', ha='center')
    #
    # ax1 = fig18.add_subplot(spec2 [ 0, 0 ])
    # x1 = ax1.plot(t [ 0:-2 ], mylist [ 0:-2 ], 'r-')
    # ax1.annotate('Luminal Concentration', **anno_opts)
    # ax1.tick_params(labelbottom=False, labelleft=True, labelsize=10)
    # ax2 = fig18.add_subplot(spec2 [ 1, 0 ])
    # x2 = ax2.plot(t, ci_na, 'blue')
    # ax2.tick_params(labelbottom=False, labelleft=True, labelsize=10)
    # ax2.set_ylabel('c_na [mm]', color='blue', fontsize=12)
    # ax2.annotate('Cellular Concentration', **anno_opts)
    # ax3 = fig18.add_subplot(spec2 [ 2, 0 ])
    # x2 = ax3.plot(t, ce_na, 'blue')
    # # ax3.tick_params(labelbottom=False, labelleft=True, labelsize=14)
    # ax3.annotate('Junctional Concentration', **anno_opts)
    # ax3.set_ylim([ 0.12, 0.138 ])
    # ax3.tick_params(labelbottom=False, labelleft=True, labelsize=10)
    # ax4 = fig18.add_subplot(spec2 [ 3, 0 ])
    # # ax4.set_ylabel('f_na [mmol/sec]', color='blue')
    # fikm_nas = [ phi_scale(x, 1e5) for x in fikm_na ]
    # fekm_nas = [ phi_scale(x, 1e5) for x in fekm_na ]
    # x2 = ax4.plot(t, fikm_nas, 'r')
    # x3 = ax4.plot(t, fekm_nas, 'blue')
    # ax4.set_ylim([ 0.0, 1.25 ])
    # ax4.legend([ x2 [ 0 ], x3 [ 0 ] ], [ 'Cellular Na Fluxes', 'Junctional Na Fluxes' ])
    # ax4.tick_params(labelbottom=False, labelleft=True, labelsize=10)
    # ax2 = fig18.add_subplot(spec2 [ 4, 0 ])
    # ax2.set_ylabel('f$1e5$ [mmol/sec]', color='blue', fontsize=12)
    # ax2.tick_params(labelbottom=False, labelleft=True, labelsize=10)
    # fivms = [ phi_scale(x, 1e5) for x in fivm ]
    #
    # fevms = [ phi_scale(x, 1e5) for x in fevm ]
    # x2 = ax2.plot(t, fivms, 'g')
    # x3 = ax2.plot(t, fevms, 'blue')
    # ax2.legend([ x2 [ 0 ], x3 [ 0 ] ], [ 'Cellular-Volume', 'Junctional-Volume' ])
    #
    # # ax4.annotate('Sodium Fluxes', **anno_opts)
    # # ax2.set_title('Sodium Fluxes')
    # ax5 = fig18.add_subplot(spec2 [ -2:, 0 ])
    # # ax5.set_ylabel('            f_na$1e5$ [mmol/sec]', color='blue', fontsize=12)
    #
    # na_mi_naglucs = [ phi_scale(x, 1e5) for x in na_mi_nagluc ]
    # na_mi_nah2po4s = [ phi_scale(x, 1e5) for x in na_mi_nah2po4 ]
    # jnhe3_nas = [ phi_scale(x, 1e5) for x in jnhe3_na ]
    # x1 = ax5.plot(t, na_mi_naglucs, 'r')
    # x2 = ax5.plot(t, na_mi_nah2po4s, 'b')
    # x3 = ax5.plot(t, jnhe3_nas, 'g')
    #
    # # ax5.annotate('Cotransporters', **anno_opts)
    # ax5.legend([ x3 [ 0 ], x2 [ 0 ], x1 [ 0 ] ], [ 'nhe3', 'na/h2po4', 'na-gluc' ])
    # ax5.set_xlabel('Time [sec]', color='blue', fontsize=12)
    # # ax2.set_title('Sodium Cotransporters on Apical Membrane')
    # ax5.tick_params(labelbottom=True, labelleft=True, labelsize=10)


