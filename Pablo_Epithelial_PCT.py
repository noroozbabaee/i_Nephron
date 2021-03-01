import numpy as np
import scipy
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
#from sympy import *
import math
from PCT_GLOB import *
from scipy import optimize


def sglt_mi(cm_na, ci_na, cm_gluc, ci_gluc, z_gluc, z_na, vm, vi, ami, lmi_nagluc, param_sglt_mi):
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


# proton pumps #checked
def at_mi_h(cm_h, ci_h, vm, vi, z_h, ami, parama_mi_h):
    xm_h = f_eps(cm_h, z_h, vm)
    xi_h = f_eps(ci_h, z_h, vi)
    gamma = -(xihp * (xm_h - xi_h - xhp))
    if parama_mi_h == 0:
        return 0
    elif gamma < 0:
        return - lhp * ami * (1 - 1 / (1 + math.exp(gamma)))
    return - lhp * ami * (1 / (1 + math.exp(-gamma)))


# goldman fluxes (passive fluxes) describes the ionic flux across a cell membrane
# as a function of the transmembrane potential and the concentrations of the ion inside and outside of the cell.
# checked
def goldman(hab, a, z, va, vb, ca, cb, param_goldman):
    zab = (z * f * (va - vb) * 1.e-6) / rte
    if param_goldman == 0:
        return[0]
    elif zab == 0 or va == vb and ca > 0 and cb > 0:
        return hab * a * (ca - cb)
    elif zab > 0:
        return hab * a * zab * (ca - cb * math.exp(-zab)) / (1 - math.exp(-zab))
    else:
        return hab * a * zab * (ca * math.exp(zab) - cb) / (math.exp(zab) - 1)


# Convective Solute Fluxes
def CSF(ca, cb, flux, s, param_CSF):
    if param_CSF == 0:
        return 0

    # Definition of the mean membrane solute concentration #checked
    def lmmsc(ca, cb):
        
        import math
        if ca > 0 and cb > 0 and ca - cb != 0 and cb != 0 and math.log10(ca / cb) != 0:
            return (ca - cb) / (math.log10(ca / cb))
        else:
            return cb

    return flux * (1.00 - s) * lmmsc(ca, cb)


def matrix(init, h):
    return[init for x in range(h)]


def f_eps(c, z, v):
    if c > 0 and c != 0:
        return rte * math.log10(c) + z * f * v * 1.e-6
    else:
        return rte * math.log10(abs(c)) + z * f * v * 1.e-6


def lch(ca, cb):
    import math
    if ca > 0 and cb > 0 and (ca / cb) != 0 and cb != 0:
        return math.log10(ca / cb)
    else:
        return math.log10(abs(ca / cb))


def clhco3_mi(cm_cl, ci_cl, cm_hco3, ci_hco3, z_cl, z_hco3, vm, vi, ami, lmi_clhco3, param_clhco3_mi):
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


def nahco3_is(ci_na, cs_na, ci_hco3, cs_hco3, z_na, z_hco3, vi, vs, ais, lis_nahco3, param_nahco3_is):
    if param_nahco3_is == 0:
        return[0, 0]
    else:
        xi_na = f_eps(ci_na, z_na, vi)
        xi_hco3 = f_eps(ci_hco3, z_hco3, vi)
        xs_na = f_eps(cs_na, z_na, vs)
        xs_hco3 = f_eps(cs_hco3, z_hco3, vs)
        na_is_nahco3 = lis_nahco3 * ais * (xi_na - xs_na + 3 * (xi_hco3 - xs_hco3))
        hco3_is_nahco3 = 3 * lis_nahco3 * ais * (xi_na - xs_na + 3 * (xi_hco3 - xs_hco3))
        return[na_is_nahco3, hco3_is_nahco3]


def na_clhco3_is(ci_na, cs_na, ci_cl, cs_cl, ci_hco3, cs_hco3, z_na, z_cl, z_hco3, vi, vs, ais, lis_na_clhco3,
                 param_na_clhco3_is):
    if param_na_clhco3_is == 0:
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


def clhco2_mi(cm_cl, ci_cl, cm_hco2, ci_hco2, z_cl, z_hco2, vm, vi, ami, lmi_clhco2, param_clhco2_mi):
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


def nah2po4_mi(cm_na, ci_na, cm_h2po4, ci_h2po4, z_na, z_h2po4, vm, vi, ami, lmi_nah2po4, param_nah2po4_mi):
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


def nah(ci_h, ci_na, ci_nh4, cm_h, cm_na, cm_nh4, param_nah):
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
        # jnah_na_max = cxt * psnah_na * psnah_h / (psnah_na + psnah_h)
        return[jnah_na, jnah_h, jnah_nh4]


def kcl_is(ci_k, cs_k, ci_cl, cs_cl, z_k, z_cl, vi, vs, ais, lis_kcl, param_kcl_is):
    if param_kcl_is == 0:
        return[0, 0]
    else:
        xi_k = f_eps(ci_k, z_k, vi)
        xi_cl = f_eps(ci_cl, z_cl, vi)
        xs_k = f_eps(cs_k, z_k, vs)
        xs_cl = f_eps(cs_cl, z_cl, vs)
        k_is_kcl = lis_kcl * ais * (xi_k - xs_k + xi_cl - xs_cl)
        cl_is_kcl = lis_kcl * ais * (xi_k - xs_k + xi_cl - xs_cl)
    return[k_is_kcl, cl_is_kcl]





def nak_atp(ci_k, cs_k, ci_na, cs_na, ce_k, ce_nh4, param_nak_atp):
    if param_nak_atp == 0:
        return[0, 0, 0]
    else:
        knpn = 0.0002 * (1.0 + ci_k / .00833)
        # sodium affinity
        knpk = 0.0001 * (1.0 + cs_na / .0185)
        # atpase transporter flux in is membrane
        atis_na = n_p * (ci_na / (knpn + ci_na)) ** 3 * (cs_k / (knpk + cs_k)) ** 2
        # alloW for competition betWeen k+ and nh4+
        atis_k = -atis_na * 0.667 * ce_k / (ce_k + ce_nh4 / knh4)
        atis_nh4 = -atis_na * 0.667 * ce_nh4 / (knh4 * ce_k + ce_nh4)
        return[atis_na, atis_k, atis_nh4]


# pH equilibria of four buffer pairs
def ebuf(lch, pk, ca, cb, param_ebuf):
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
def Transp_Progres_Activ(Transporters_progressive_activation_alongtime, t, T, T_window):
    if Transporters_progressive_activation_alongtime == 1:
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
        Time_Durtn = int(T / T_window)
        n = 0 #int(2 * Time_Durtn)
        if n < t < n + Time_Durtn:
            print('No Transporter Yet at Time_Durtn' if t == Time_Durtn else '')
        elif n + Time_Durtn <= t < n + 2 * Time_Durtn :
            sglt_mi_param = 1
            print('sglt_mi_param_activation at time =' + str(Time_Durtn) if t == n + 2 * Time_Durtn else '')
        elif n + 2 * Time_Durtn <= t< n + 3* Time_Durtn :
            sglt_mi_param = 1
            nah2po4_mi_param = 1
        elif n + 3 * Time_Durtn <= t < n + 4 * Time_Durtn:
            sglt_mi_param = 1
            nah2po4_mi_param = 1
            clhco3_mi_param = 1
        elif n +4 * Time_Durtn <= t < n + 5 * Time_Durtn :
            sglt_mi_param = 1
            nah2po4_mi_param = 1
            clhco3_mi_param = 1
            clhco2_mi_param = 1
        elif n + 5 * Time_Durtn <= t< n + 6 * Time_Durtn :
            sglt_mi_param = 1
            nah2po4_mi_param = 1
            clhco3_mi_param = 1
            clhco2_mi_param = 1
            nahco3_is_param = 1
            nahco3_ie_param = 1
        elif n + 6 * Time_Durtn <= t< n + 7 * Time_Durtn:
            sglt_mi_param = 1
            nah2po4_mi_param = 1
            clhco3_mi_param = 1
            clhco2_mi_param = 1
            nahco3_is_param = 1
            nahco3_ie_param = 1
            kcl_is_param = 1
            kcl_ie_param = 1
        elif n + 7 * Time_Durtn <= t < n + 8 * Time_Durtn:
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
        elif n + 8 * Time_Durtn <= t < n + 9 * Time_Durtn:
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
        elif n + 9 * Time_Durtn <= t < n + 10 * Time_Durtn :
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
        # return sglt_mi_param, nah2po4_mi_param, clhco3_mi_param, clhco2_mi_param, nahco3_is_param, nahco3_ie_param, \
        #        kcl_is_param, kcl_ie_param, na_clhco3_is_param, na_clhco3_ie_param, nah_param, nak_atp_param, h_mi_atp_param
    else:
        if t >= int(T/4):
            print(
                'Transporters added to the system at time=' + str(t) if t == int(T / 4) else '')
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
        else:
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
    return sglt_mi_param, nah2po4_mi_param, clhco3_mi_param, clhco2_mi_param, nahco3_is_param, nahco3_ie_param, \
    kcl_is_param, kcl_ie_param, na_clhco3_is_param, na_clhco3_ie_param, nah_param, nak_atp_param, h_mi_atp_param
def Buff_Activ_co2_formate_phosphate_ammonia(q_h, q_nh4, q_hco3, q_h2co2, q_h2po4, q_co2, q_h2co3, c_co2, c_h2co3, volume, scale,flow_dependent, Co2_Progressive_Activation_Param):
    if Co2_Progressive_Activation_Param == 1:
        # print('Buffer activation at time =' + str(int(4000 * Nfac)))
        # print('Buffer activation for co2, formate, phosphate, and ammonia', 't=', str(t))
        # co2, formate, phosphate, and ammonia content:
        # print('q_hco3,  q_h2co3, q_co2', str(q_hco3), str(q_h2co3), str(q_co2))
         if flow_dependent == 1:
            q_hco3 = + q_h + q_nh4 - q_hco3 + q_h2co2 + q_h2po4
            q_h2co3 = q_co2 + scale * volume * (khy * c_co2 - kdhy * c_h2co3)
            q_co2  = q_hco3 + q_h2co3 + q_co2
         else:
        # print('No Buffer Effect for co2, formate, phosphate, and ammonia')
            q_hco3 = + q_h + q_nh4 - q_hco3 + q_h2co2 + q_h2po4

            # LIS: hydration and dhydration of co2
            # see label {co2_hyd_dhyd}
            q_h2co3 = q_co2 + scale * volume * (khy * c_co2 - kdhy * c_h2co3)
            # LIS: see label {conser_charge_co2}
            q_co2 = q_hco3 + q_h2co3 + q_co2
    #print('q_hco3,  q_h2co3, q_co2', str(q_hco3),  str(q_h2co3), str(q_co2))
    return q_hco3,  q_h2co3, q_co2


def Buff_Activ(q_h_vary, q_h2_vary, lc_h, pk, c_h_vary, c_h2_vary, qi_nh3_cell, ebuf_Param,  Buffer_Activation_Param):
    if Buffer_Activation_Param == 1:
        # print('buffer_effect_activate_over_time_for_phosphate', 't=', str(t))
        # def lch(ca, cb):
        #
        #     if ca > 0 and cb > 0 and (ca / cb) != 0 and cb != 0:
        #         return math.log10(ca / cb)
        #     else:
        #         return math.log10(abs(ca / cb))
        #
        # def ebuf(lc_h, pk, ca, cb, param_ebuf):
        #     # pH equilibria of four buffer pairs
        #     if param_ebuf == 0:
        #         return 0
        #     if ca > 0 and cb > 0 and (ca / cb) != 0 and cb != 0:
        #         return lc_h - pk - math.log10(ca / cb)
        #
        #     else:
        #         return lc_h - pk - math.log10(abs(ca / cb))

        # LIS: see label {phosphate}
        if qi_nh3_cell==1:
           q_h_vary = q_h_vary + q_h2_vary - 1.e6 * qiamm
           q_h2_vary = ebuf(lc_h, pk, c_h_vary, c_h2_vary, ebuf_Param)
        else:
           q_h_vary = q_h_vary + q_h2_vary
           print('Buffer Effect for phosphate')
           # LIS: see labela {pH_equilibria}, used as paired equations
           q_h2_vary = ebuf(lc_h, pk, c_h_vary, c_h2_vary, ebuf_Param)
    else:
        print('No Buffer Effect for phosphate')
        pass
    return q_h_vary, q_h2_vary

def eQs(guess, solver):
    # update variables
    for i in range(len(guess)):
        vars[i][t] = guess[i]
    # Weinstein (2007)
    # Flow-dependent transport in a mathematical model of rat proximal tubule
    ae[t] = ae0 * (1 + mua * (pe[t] - pm))
    ae[t] = ae0 if ae[t] < ae0 else ae[t]

    chvl[t] = chvl0 * (1.0 + muv * (pe[t] - pm))
    chvl[t] = chvl0 if chvl[t] < chvl0 else chvl[t]
    l[t] = chvl[t]

    # LOG Conc. of Hydrogen Interspace
    lche = pkc + lch(ce_hco3[t], ce_h2co3[t])
    ce_h[t] = 10 ** (-lche)

    # ELECTROLYTE TRANSPORT ACROSS A SIMPLE EPITHELIUM of the rat proximal tubule (1992)
    # The cell volume is presented as a function of the inverse cell impermeant species concentration
    clvl[t] = clvl0 * imp0 / imp[t]
    clvl[t] = clvl0 if imp[t] == 0 else clvl[t]

    # The cell height (summation of Extracellular channel volume
    # and Intracellular compartment volume[cm3/cm2 epithelium])
    l[t] = l[t] + clvl[t]
    p_i = pm

    # LOG Conc. of Hydrogen Cell
    lchi = pkc + lch(ci_hco3[t], ci_h2co3[t])
    ci_h[t] = 10 ** (-lchi)
    # ELECTROLYTE TRANSPORT ACROSS A SIMPLE EPITHELIUM 1979
    # Applying the Kedem and Katchalsky relations:
    # The membrane transport for nonelectrolyte solutions generated by the hydrostatic
    # pressure difference and the osmotic pressure difference
    # The equations below represent the volume flow or the convective volume flow in between the different membrane
    fevm = ame * lpme * (pm - pe[t] - rt * impm) / rt
    fevs = ae[t] * lpes * (rt * imps + pe[t] - ps) / rt
    fivm = ami * lpmi * (rt * imp[t] - rt * impm + pm - p_i) / rt
    fivs = ais * lpis * (p_i - ps + rt * imps - rt * imp[t]) / rt
    jv = lpis * aie * (p_i - pe[t] - rt * imp[t]) / rt
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

    jv = jv + lpis * aie * (
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
    param_CSF = 1
    fekm_na = CSF(ce_na[t], cm_na, fevm, sme_na, param_CSF)
    fekm_k = CSF(ce_k[t], cm_k, fevm, sme_k, param_CSF)
    fekm_cl = CSF(ce_cl[t], cm_cl, fevm, sme_cl, param_CSF)
    fekm_hco3 = CSF(ce_hco3[t], cm_hco3, fevm, sme_hco3, param_CSF)
    fekm_h2co3 = CSF(ce_h2co3[t], cm_h2co3, fevm, sme_h2co3, param_CSF)
    fekm_co2 = CSF(ce_co2[t], cm_co2, fevm, sme_co2, param_CSF)
    fekm_hpo4 = CSF(ce_hpo4[t], cm_hpo4, fevm, sme_hpo4, param_CSF)
    fekm_h2po4 = CSF(ce_h2po4[t], cm_h2po4, fevm, sme_h2po4, param_CSF)
    fekm_urea = CSF(ce_urea[t], cm_urea, fevm, sme_urea, param_CSF)
    fekm_nh3 = CSF(ce_nh3[t], cm_nh3, fevm, sme_nh3, param_CSF)
    fekm_nh4 = CSF(ce_nh4[t], cm_nh4, fevm, sme_nh4, param_CSF)
    fekm_h = CSF(ce_h[t], cm_h, fevm, sme_h, param_CSF)
    fekm_hco2 = CSF(ce_hco2[t], cm_hco2, fevm, sme_hco2, param_CSF)
    fekm_h2co2 = CSF(ce_h2co2[t], cm_h2co2, fevm, sme_h2co2, param_CSF)
    fekm_gluc = CSF(ce_gluc[t], cm_gluc, fevm, sme_gluc, param_CSF)

    feks_na = CSF(ce_na[t], cs_na, fevs, ses_na, param_CSF)
    feks_k = CSF(ce_k[t], cs_k, fevs, ses_k, param_CSF)
    feks_cl = CSF(ce_cl[t], cs_cl, fevs, ses_cl, param_CSF)
    feks_hco3 = CSF(ce_hco3[t], cs_hco3, fevs, ses_hco3, param_CSF)
    feks_h2co3 = CSF(ce_h2co3[t], cs_h2co3, fevs, ses_h2co3, param_CSF)
    feks_co2 = CSF(ce_co2[t], cs_co2, fevs, ses_co2, param_CSF)
    feks_hpo4 = CSF(ce_hpo4[t], cs_hpo4, fevs, ses_hpo4, param_CSF)
    feks_h2po4 = CSF(ce_h2po4[t], cs_h2po4, fevs, ses_h2po4, param_CSF)
    feks_urea = CSF(ce_urea[t], cs_urea, fevs, ses_urea, param_CSF)
    feks_nh3 = CSF(ce_nh3[t], cs_nh3, fevs, ses_nh3, param_CSF)
    feks_nh4 = CSF(ce_nh4[t], cs_nh4, fevs, ses_nh4, param_CSF)
    feks_h = CSF(ce_h[t], cs_h, fevs, ses_h, param_CSF)
    feks_hco2 = CSF(ce_hco2[t], cs_hco2, fevs, ses_hco2, param_CSF)
    feks_h2co2 = CSF(ce_h2co2[t], cs_h2co2, fevs, ses_h2co2, param_CSF)
    feks_gluc = CSF(ce_gluc[t], cs_gluc, fevs, ses_gluc, param_CSF)

    fikm_na = CSF(ci_na[t], cm_na, fivm, smi_na, param_CSF)
    fikm_k = CSF(ci_k[t], cm_k, fivm, smi_k, param_CSF)
    fikm_cl = CSF(ci_cl[t], cm_cl, fivm, smi_cl, param_CSF)
    fikm_hco3 = CSF(ci_hco3[t], cm_hco3, fivm, smi_hco3, param_CSF)
    fikm_h2co3 = CSF(ci_h2co3[t], cm_h2co3, fivm, smi_h2co3, param_CSF)
    fikm_co2 = CSF(ci_co2[t], cm_co2, fivm, smi_co2, param_CSF)
    fikm_hpo4 = CSF(ci_hpo4[t], cm_hpo4, fivm, smi_hpo4, param_CSF)
    fikm_h2po4 = CSF(ci_h2po4[t], cm_h2po4, fivm, smi_h2po4, param_CSF)
    fikm_urea = CSF(ci_urea[t], cm_urea, fivm, smi_urea, param_CSF)
    fikm_nh3 = CSF(ci_nh3[t], cm_nh3, fivm, smi_nh3, param_CSF)
    fikm_nh4 = CSF(ci_nh4[t], cm_nh4, fivm, smi_nh4, param_CSF)
    fikm_h = CSF(ci_h[t], cm_h, fivm, smi_h, param_CSF)
    fikm_hco2 = CSF(ci_hco2[t], cm_hco2, fivm, smi_hco2, param_CSF)
    fikm_h2co2 = CSF(ci_h2co2[t], cm_h2co2, fivm, smi_h2co2, param_CSF)
    fikm_gluc = CSF(ci_gluc[t], cm_gluc, fivm, smi_gluc, param_CSF)

    fiks_na = CSF(ci_na[t], cs_na, fivs, sis_na, param_CSF)
    fiks_k = CSF(ci_k[t], cs_k, fivs, sis_k, param_CSF)
    fiks_cl = CSF(ci_cl[t], cs_cl, fivs, sis_cl, param_CSF)
    fiks_hco3 = CSF(ci_hco3[t], cs_hco3, fivs, sis_hco3, param_CSF)
    fiks_h2co3 = CSF(ci_h2co3[t], cs_h2co3, fivs, sis_h2co3, param_CSF)
    fiks_co2 = CSF(ci_co2[t], cs_co2, fivs, sis_co2, param_CSF)
    fiks_hpo4 = CSF(ci_hpo4[t], cs_hpo4, fivs, sis_hpo4, param_CSF)
    fiks_h2po4 = CSF(ci_h2po4[t], cs_h2po4, fivs, sis_h2po4, param_CSF)
    fiks_urea = CSF(ci_urea[t], cs_urea, fivs, sis_urea, param_CSF)
    fiks_nh3 = CSF(ci_nh3[t], cs_nh3, fivs, sis_nh3, param_CSF)
    fiks_nh4 = CSF(ci_nh4[t], cs_nh4, fivs, sis_nh4, param_CSF)
    fiks_h = CSF(ci_h[t], cs_h, fivs, sis_h, param_CSF)
    fiks_hco2 = CSF(ci_hco2[t], cs_hco2, fivs, sis_hco2, param_CSF)
    fiks_h2co2 = CSF(ci_h2co2[t], cs_h2co2, fivs, sis_h2co2, param_CSF)
    fiks_gluc = CSF(ci_gluc[t], cs_gluc, fivs, sis_gluc, param_CSF)

    jk_na = CSF(ci_na[t], ce_na[t], jv, sis_na, param_CSF)
    jk_k = CSF(ci_k[t], ce_k[t], jv, sis_k, param_CSF)
    jk_cl = CSF(ci_cl[t], ce_cl[t], jv, sis_cl, param_CSF)
    jk_hco3 = CSF(ci_hco3[t], ce_hco3[t], jv, sis_hco3, param_CSF)
    jk_h2co3 = CSF(ci_h2co3[t], ce_h2co3[t], jv, sis_h2co3, param_CSF)
    jk_co2 = CSF(ci_co2[t], ce_co2[t], jv, sis_co2, param_CSF)
    jk_hpo4 = CSF(ci_hpo4[t], ce_hpo4[t], jv, sis_hpo4, param_CSF)
    jk_h2po4 = CSF(ci_h2po4[t], ce_h2po4[t], jv, sis_h2po4, param_CSF)
    jk_urea = CSF(ci_urea[t], ce_urea[t], jv, sis_urea, param_CSF)
    jk_nh3 = CSF(ci_nh3[t], ce_nh3[t], jv, sis_nh3, param_CSF)
    jk_nh4 = CSF(ci_nh4[t], ce_nh4[t], jv, sis_nh4, param_CSF)
    jk_h = CSF(ci_h[t], ce_h[t], jv, sis_h, param_CSF)
    jk_hco2 = CSF(ci_hco2[t], ce_hco2[t], jv, sis_hco2, param_CSF)
    jk_h2co2 = CSF(ci_h2co2[t], ce_h2co2[t], jv, sis_h2co2, param_CSF)
    jk_gluc = CSF(ci_gluc[t], ce_gluc[t], jv, sis_gluc, param_CSF)

    # Goldman Fluxes: fluxes due to diffusion and electrical potential difference for
    # all of the ions that are permeant through the membrane
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

    feks_na = feks_na + goldman(hes_na, ae[t], z_na, ve[t], vs, ce_na[t],
                                cs_na, goldman_param)
    feks_k = feks_k + goldman(hes_k, ae[t], z_k, ve[t], vs, ce_k[t], cs_k, goldman_param)
    feks_cl = feks_cl + goldman(hes_cl, ae[t], z_cl, ve[t], vs, ce_cl[t], cs_cl, goldman_param)
    feks_hco3 = feks_hco3 + goldman(hes_hco3, ae[t], z_hco3, ve[t], vs, ce_hco3[t], cs_hco3, goldman_param)
    feks_h2co3 = feks_h2co3 + goldman(hes_h2co3, ae[t], z_h2co3, ve[t], vs, ce_h2co3[t],
                                      cs_h2co3, goldman_param)
    feks_co2 = feks_co2 + goldman(hes_co2, ae[t], z_co2, ve[t], vs, ce_co2[t], cs_co2, goldman_param)
    feks_hpo4 = feks_hpo4 + goldman(hes_hpo4, ae[t], z_hpo4, ve[t], vs, ce_hpo4[t],
                                    cs_hpo4, goldman_param)
    feks_h2po4 = feks_h2po4 + goldman(hes_h2po4, ae[t], z_h2po4, ve[t], vs, ce_h2po4[t],
                                      cs_h2po4, goldman_param)
    feks_urea = feks_urea + goldman(hes_urea, ae[t], z_urea, ve[t], vs, ce_urea[t], cs_urea, goldman_param)
    feks_nh3 = feks_nh3 + goldman(hes_nh3, ae[t], z_nh3, ve[t], vs, ce_nh3[t], cs_nh3, goldman_param)
    feks_nh4 = feks_nh4 + goldman(hes_nh4, ae[t], z_nh4, ve[t], vs, ce_nh4[t], cs_nh4, goldman_param)
    feks_h = feks_h + goldman(hes_h, ae[t], z_h, ve[t], vs, ce_h[t], cs_h, goldman_param)
    feks_hco2 = feks_hco2 + goldman(hes_hco2, ae[t], z_hco2, ve[t], vs, ce_hco2[t], cs_hco2, goldman_param)

    feks_h2co2 = feks_h2co2 + goldman(hme_h2co2, ame, z_h2co2, ve[t], vs, ce_h2co2[t], cs_h2co2, goldman_param)
    feks_gluc = feks_gluc + goldman(hes_gluc, ae[t], z_gluc, ve[t], vs, ce_gluc[t], cs_gluc, goldman_param)

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

    jk_na = jk_na + goldman(his_na, aie, z_na, vi[t], ve[t], ci_na[t], ce_na[t], goldman_param)
    jk_k = jk_k + goldman(his_k, aie, z_k, vi[t], ve[t], ci_k[t], ce_k[t], goldman_param)
    jk_cl = jk_cl + goldman(his_cl, aie, z_cl, vi[t], ve[t], ci_cl[t], ce_cl[t], goldman_param)
    jk_hco3 = jk_hco3 + goldman(his_hco3, aie, z_hco3, vi[t], ve[t], ci_hco3[t], ce_hco3[t], goldman_param)
    jk_h2co3 = jk_h2co3 + goldman(his_h2co3, aie, z_h2co3, vi[t], ve[t], ci_h2co3[t],
                                  ce_h2co3[t], goldman_param)
    jk_co2 = jk_co2 + goldman(his_co2, aie, z_co2, vi[t], ve[t], ci_co2[t], ce_co2[t], goldman_param)
    jk_hpo4 = jk_hpo4 + goldman(his_hpo4, aie, z_hpo4, vi[t], ve[t], ci_hpo4[t],
                                ce_hpo4[t], goldman_param)
    jk_h2po4 = jk_h2po4 + goldman(his_h2po4, aie, z_h2po4, vi[t], ve[t], ci_h2po4[t],
                                  ce_h2po4[t], goldman_param)
    jk_urea = jk_urea + goldman(his_urea, aie, z_urea, vi[t], ve[t], ci_urea[t], ce_urea[t], goldman_param)
    jk_nh3 = jk_nh3 + goldman(his_nh3, aie, z_nh3, vi[t], ve[t], ci_nh3[t], ce_nh3[t], goldman_param)
    jk_nh4 = jk_nh4 + goldman(his_nh4, aie, z_nh4, vi[t], ve[t], ci_nh4[t], ce_nh4[t], goldman_param)
    jk_h = jk_h + goldman(his_h, aie, z_h, vi[t], ve[t], ci_h[t], ce_h[t], goldman_param)
    jk_hco2 = jk_hco2 + goldman(his_hco2, aie, z_hco2, vi[t], ve[t], ci_hco2[t], ce_hco2[t], goldman_param)
    jk_h2co2 = jk_h2co2 + goldman(his_h2co2, aie, z_h2co2, vi[t], ve[t], ci_h2co2[t], ce_h2co2[t],
                                  goldman_param)
    jk_gluc = jk_gluc + goldman(his_gluc, aie, z_gluc, vi[t], ve[t], ci_gluc[t], ce_gluc[t], goldman_param)

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
    # electodiffusive fluxes  are proportional to the differences
    # in electrochemical driving forces
    # Updating flux equations with electodiffusive fluxes (transporters).
    # Electodiffusive fluxes  are proportional to the differences in electrochemical driving forces
    T_window = 15
    Transporters_progressive_activation_alongtime = 1
    # sglt_mi_param, nah2po4_mi_param, clhco3_mi_param, clhco2_mi_param, nahco3_is_param, nahco3_ie_param, \
    # kcl_is_param, kcl_ie_param, na_clhco3_is_param, na_clhco3_ie_param, nah_param, nak_atp_param, h_mi_atp_param = \
    #     Transp_Progres_Activ(Transporters_progressive_activation_alongtime, t, T, T_window)
    # print('Transp_Progres_Activ', Transp_Progres_Activ(Transporters_progressive_activation_alongtime, t, T, T_window))
    if Transporters_progressive_activation_alongtime == 1:
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
        Time_Durtn = int(T / T_window)
        n = 0 #int(2 * Time_Durtn)
        if n < t < n + Time_Durtn:
            print('No Transporter Yet at Time_Durtn' if t == Time_Durtn else '')
        elif n + Time_Durtn <= t < n + 2 * Time_Durtn :
            sglt_mi_param = 1
            print('sglt_mi_param_activation at time =' + str(Time_Durtn) if t == n + 2 * Time_Durtn else '')
        elif n + 2 * Time_Durtn <= t< n + 3* Time_Durtn :
            sglt_mi_param = 1
            nah2po4_mi_param = 1
        elif n + 3 * Time_Durtn <= t < n + 4 * Time_Durtn:
            sglt_mi_param = 1
            nah2po4_mi_param = 1
            clhco3_mi_param = 1
        elif n +4 * Time_Durtn <= t < n + 5 * Time_Durtn :
            sglt_mi_param = 1
            nah2po4_mi_param = 1
            clhco3_mi_param = 1
            clhco2_mi_param = 1
        elif n + 5 * Time_Durtn <= t< n + 6 * Time_Durtn :
            sglt_mi_param = 1
            nah2po4_mi_param = 1
            clhco3_mi_param = 1
            clhco2_mi_param = 1
            nahco3_is_param = 1
            nahco3_ie_param = 1
        elif n + 6 * Time_Durtn <= t< n + 7 * Time_Durtn:
            sglt_mi_param = 1
            nah2po4_mi_param = 1
            clhco3_mi_param = 1
            clhco2_mi_param = 1
            nahco3_is_param = 1
            nahco3_ie_param = 1
            kcl_is_param = 1
            kcl_ie_param = 1
        elif n + 7 * Time_Durtn <= t < n + 8 * Time_Durtn:
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
        elif n + 8 * Time_Durtn <= t < n + 9 * Time_Durtn:
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
        elif n + 9 * Time_Durtn <= t < n + 10 * Time_Durtn :
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

    # net transporters on mi bourder
    sglt = sglt_mi(cm_na, ci_na[t], cm_gluc, ci_gluc[t], z_na, z_gluc, vm[t], vi[t], ami, lmi_nagluc,
                   sglt_mi_param)
    na_mi_nagluc = sglt[0]
    gluc_mi_nagluc = sglt[1]

    nah2po4 = nah2po4_mi(cm_na, ci_na[t], cm_h2po4, ci_h2po4[t], z_na, z_h2po4, vm[t], vi[t], ami,
                         lmi_nah2po4, nah2po4_mi_param)
    na_mi_nah2po4 = nah2po4[0]
    h2po4_mi_nah2po4 = nah2po4[1]

    clhco3 = clhco3_mi(cm_cl, ci_cl[t], cm_hco3, ci_hco3[t], z_cl, z_hco3, vm[t], vi[t], ami,
                       lmi_clhco3, clhco3_mi_param)
    cl_mi_clhco3 = clhco3[0]
    hco3_mi_clhco3 = clhco3[1]

    clhco2 = clhco2_mi(cm_cl, ci_cl[t], cm_hco2, ci_hco2[t], z_cl, z_hco2, vm[t], vi[t], ami,
                       lmi_clhco2, clhco2_mi_param)
    cl_mi_clhco2 = clhco2[0]
    hco2_mi_clhco2 = clhco2[1]

    # net cotransporters on is bourder
    nahco3 = nahco3_is(ci_na[t], cs_na, ci_hco3[t], cs_hco3, z_na, z_hco3, vi[t], vs, ais,
                       lis_nahco3, nahco3_is_param)
    na_is_nahco3 = nahco3[0]
    hco3_is_nahco3 = nahco3[1]

    kcl = kcl_is(ci_k[t], cs_k, ci_cl[t], cs_cl, z_k, z_cl, vi[t], vs, ais, lis_kcl, kcl_is_param)
    k_is_kcl = kcl[0]
    cl_is_kcl = kcl[1]

    na_clhco3 = na_clhco3_is(ci_na[t], cs_na, ci_cl[t], cs_cl, ci_hco3[t], cs_hco3, z_na, z_cl, z_hco3,
                             vi[t], vs, ais, lis_na_clhco3, na_clhco3_is_param)
    na_is_na_clhco3 = na_clhco3[0]
    cl_is_na_clhco3 = na_clhco3[ 1]
    hco3_is_na_clhco3 = na_clhco3[2]

    # the nah exchanger translate concentrations to the nah model on  mi bourder
    mynah = nah(ci_h[t], ci_na[t], ci_nh4[t], cm_h, cm_na, cm_nh4, nah_param)
    jnah_na = mynah[0]
    jnah_h = mynah[1]
    jnah_nh4 = mynah[2]
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

    nahco3_aie = nahco3_is(ci_na[t], ce_na[t], ci_hco3[t], ce_hco3[t], z_na, z_hco3, vi[t], ve[t],
                           aie,
                           lis_nahco3, nahco3_is_param)
    na_aie_nahco3 = nahco3_aie[0]
    hco3_aie_nahco3 = nahco3_aie[1]

    kcl_aie = kcl_is(ci_k[t], ce_k[t], ci_cl[t], ce_cl[t], z_k, z_cl, vi[t], ve[t], aie, lis_kcl, kcl_is_param)
    k_aie_kcl = kcl_aie[0]
    cl_aie_kcl = kcl_aie[1]

    na_clhco3_aie = na_clhco3_is(ci_na[t], ce_na[t], ci_cl[t], ce_cl[t], ci_hco3[t], ce_hco3[t], z_na,
                                 z_cl, z_hco3,
                                 vi[t], ve[t], aie, lis_na_clhco3, na_clhco3_is_param)
    na_aie_na_clhco3 = na_clhco3_aie[0]
    cl_aie_na_clhco3 = na_clhco3_aie[1]
    hco3_aie_na_clhco3 = na_clhco3_aie[2]

    jk_na = jk_na + na_aie_nahco3 + na_aie_na_clhco3
    jk_k = jk_k + k_aie_kcl
    jk_cl = jk_cl + cl_aie_kcl + cl_aie_na_clhco3
    jk_hco3 = jk_hco3 + hco3_aie_na_clhco3 + hco3_aie_nahco3

    # sodium pumps on mi bourder
    nak = nak_atp(ci_k[t], cs_k, ci_na[t], cs_na, ce_k[t], ce_nh4[t], nak_atp_param)
    atis_na = nak[0]
    atis_k = nak[1]
    atis_nh4 = nak[2]
    # proton pumps
    atmi_h = at_mi_h(cm_h, ci_h[t], vm[t], vi[t], z_h, ami,  h_mi_atp_param)

    jk_na = jk_na + aie * atis_na
    jk_k = jk_k + aie * atis_k
    jk_nh4 = jk_nh4 + aie * atis_nh4

    fiks_na = fiks_na + ais * atis_na
    fiks_k = fiks_k + ais * atis_k
    fiks_nh4 = fiks_nh4 + ais * atis_nh4
    # proton pumps
    fikm_h = fikm_h + atmi_h
    # establish the error vectors, the "phi" array.
    # first for the interspace electroneutrality
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

    # mass conservation in the time-dependent case
    phie_vlm = fevs - fevm - jv + rtau * (chvl[t] - chvl[t - 1])
    qe_na = feks_na - fekm_na - jk_na + rtau * (ce_na[t] * chvl[t] - ce_na[t - 1] * chvl[t - 1])
    qe_k = feks_k - fekm_k - jk_k + rtau * (ce_k[t] * chvl[t] - ce_k[t - 1] * chvl[t - 1])
    qe_cl = feks_cl - fekm_cl - jk_cl + rtau * (ce_cl[t] * chvl[t] - ce_cl[t - 1] * chvl[t - 1])
    qe_hco3 = feks_hco3 - fekm_hco3 - jk_hco3 + rtau * (ce_hco3[t] * chvl[t] - ce_hco3[t - 1] * chvl[t - 1])
    qe_h2co3 = feks_h2co3 - fekm_h2co3 - jk_h2co3 + rtau * (
            ce_h2co3[t] * chvl[t] - ce_h2co3[t - 1] * chvl[t - 1])
    qe_co2 = feks_co2 - fekm_co2 - jk_co2 + rtau * (ce_co2[t] * chvl[t] - ce_co2[t - 1] * chvl[t - 1])
    qe_hpo4 = feks_hpo4 - fekm_hpo4 - jk_hpo4 + rtau * (ce_hpo4[t] * chvl[t] - ce_hpo4[t - 1] * chvl[t - 1])
    qe_h2po4 = feks_h2po4 - fekm_h2po4 - jk_h2po4 + rtau * (
            ce_h2po4[t] * chvl[t] - ce_h2po4[t - 1] * chvl[t - 1])
    qe_urea = feks_urea - fekm_urea - jk_urea + rtau * (
            ce_urea[t] * chvl[t] - ce_urea[t - 1] * chvl[t - 1])
    qe_nh3 = feks_nh3 - fekm_nh3 - jk_nh3 + rtau * (ce_nh3[t] * chvl[t] - ce_nh3[t - 1] * chvl[t - 1])
    qe_nh4 = feks_nh4 - fekm_nh4 - jk_nh4 + rtau * (ce_nh4[t] * chvl[t] - ce_nh4[t - 1] * chvl[t - 1])
    qe_h = feks_h - fekm_h - jk_h + rtau * (ce_h[t] * chvl[t] - ce_h[t - 1] * chvl[t - 1])
    qe_hco2 = feks_hco2 - fekm_hco2 - jk_hco2 + rtau * (ce_hco2[t] * chvl[t] - ce_hco2[t - 1] * chvl[t - 1])
    qe_h2co2 = feks_h2co2 - fekm_h2co2 - jk_h2co2 + rtau * (
            ce_h2co2[t] * chvl[t] - ce_h2co2[t - 1] * chvl[t - 1])

    qe_gluc = feks_gluc - fekm_gluc - jk_gluc + rtau * (
            ce_gluc[t] * chvl[t] - ce_gluc[t - 1] * chvl[t - 1])

    # mass conservation in the time - dependent case
    phii_vlm = fivs - fivm + jv + rtau * (clvl[t] - clvl[t - 1])

    qi_na = fiks_na - fikm_na + jk_na + rtau * (
            ci_na[t] * clvl[t] - ci_na[t - 1] * clvl[t - 1])

    qi_k = fiks_k - fikm_k + jk_k + rtau * (
            ci_k[t] * clvl[t] - ci_k[t - 1] * clvl[t - 1])
    qi_cl = fiks_cl - fikm_cl + jk_cl + rtau * (
            ci_cl[t] * clvl[t] - ci_cl[t - 1] * clvl[t - 1])

    qi_hco3 = fiks_hco3 - fikm_hco3 + jk_hco3 + rtau * (
            ci_hco3[t] * clvl[t] - ci_hco3[t - 1] * clvl[t - 1])
    qi_h2co3 = fiks_h2co3 - fikm_h2co3 + jk_h2co3 + rtau * (
            ci_h2co3[t] * clvl[t] - ci_h2co3[t - 1] * clvl[t - 1])
    qi_co2 = fiks_co2 - fikm_co2 + jk_co2 + rtau * (
            ci_co2[t] * clvl[t] - ci_co2[t - 1] * clvl[t - 1])
    qi_hpo4 = fiks_hpo4 - fikm_hpo4 + jk_hpo4 + rtau * (
            ci_hpo4[t] * clvl[t] - ci_hpo4[t - 1] * clvl[t - 1])
    qi_h2po4 = fiks_h2po4 - fikm_h2po4 + jk_h2po4 + rtau * (
            ci_h2po4[t] * clvl[t] - ci_h2po4[t - 1] * clvl[t - 1])

    qi_urea = fiks_urea - fikm_urea + jk_urea + rtau * (
            ci_urea[t] * clvl[t] - ci_urea[t - 1] * clvl[t - 1])
    qi_nh3 = fiks_nh3 - fikm_nh3 + jk_nh3 + rtau * (
            ci_nh3[t] * clvl[t] - ci_nh3[t - 1] * clvl[t - 1])
    qi_nh4 = fiks_nh4 - fikm_nh4 + jk_nh4 + rtau * (
            ci_nh4[t] * clvl[t] - ci_nh4[t - 1] * clvl[t - 1])
    qi_h = fiks_h - fikm_h + jk_h + rtau * (
            ci_h[t] * clvl[t] - ci_h[t - 1] * clvl[t - 1])
    qi_hco2 = fiks_hco2 - fikm_hco2 + jk_hco2 + rtau * (
            ci_hco2[t] * clvl[t] - ci_hco2[t - 1] * clvl[t - 1])
    qi_h2co2 = fiks_h2co2 - fikm_h2co2 + jk_h2co2 + rtau * (
            ci_h2co2[t] * clvl[t] - ci_h2co2[t - 1] * clvl[t - 1])
    qi_gluc = fiks_gluc - fikm_gluc + jk_gluc + rtau * (
            ci_gluc[t] * clvl[t] - ci_gluc[t - 1] * clvl[t - 1])

    # the proton flux must include the cellular buffers
    qi_h = qi_h - rtau * (cbuf[t] * clvl[t] - cbuf[t - 1] * clvl[t - 1])
    scale = 1e6
    phie_vlm =phi_scale(phie_vlm, scale)
    qe_na =phi_scale(qe_na, scale)
    qe_k =phi_scale(qe_k, scale)
    qe_cl =phi_scale(qe_cl, scale)
    qe_hco3 =phi_scale(qe_hco3, scale)
    qe_h2co3 =phi_scale(qe_h2co3, scale)
    qe_co2 =phi_scale(qe_co2, scale)
    qe_hpo4 =phi_scale(qe_hpo4, scale)
    qe_h2po4 =phi_scale(qe_h2po4, scale)
    qe_urea =phi_scale(qe_urea, scale)
    qe_nh3 =phi_scale(qe_nh3, scale)
    qe_nh4 =phi_scale(qe_nh4, scale)
    qe_h =phi_scale(qe_h, scale)
    qe_hco2 =phi_scale(qe_hco2, scale)
    qe_h2co2 =phi_scale(qe_h2co2, scale)
    qe_gluc =phi_scale(qe_gluc, scale)

    phii_vlm =phi_scale(phii_vlm, scale)
    qi_na =phi_scale(qi_na, scale)
    qi_k =phi_scale(qi_k, scale)
    qi_cl =phi_scale(qi_cl, scale)
    qi_hco3 =phi_scale(qi_hco3, scale)
    qi_h2co3 =phi_scale(qi_h2co3, scale)
    qi_co2 =phi_scale(qi_co2, scale)
    qi_hpo4 =phi_scale(qi_hpo4, scale)
    qi_h2po4 =phi_scale(qi_h2po4, scale)
    qi_urea =phi_scale(qi_urea, scale)
    qi_nh3 =phi_scale(qi_nh3, scale)
    qi_nh4 =phi_scale(qi_nh4, scale)
    qi_h =phi_scale(qi_h, scale)
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
        Time_Durtn = int(T / T_window)
        n = int(2 * Time_Durtn)
        if n + 10 * Time_Durtn < t :
            # the error term for bicarbonate generation is replaced by conservation
            # of charge in the buffer reactions.
            qe_hco3 = - qe_hco3 + qe_h2po4 + qe_nh4 + qe_h + qe_h2co2
            phi[5] = qe_hco3
            # hydration and dhydration of co2
            qe_h2co3 = qe_co2 + 1.e6 * chvl[t] * (khy_4 * ce_co2[t] - kdhy_4 * ce_h2co3[t])
            phi[6] = qe_h2co3
            qe_co2 = qe_hco3 + qe_h2co3 + qe_co2
            phi[7] = qe_co2

            qi_hco3 = - qi_hco3 + qi_h2po4 + qi_nh4 + qi_h + qi_h2co2
            phi[21] = qi_hco3
            qi_h2co3 = qi_co2 + 1.e6 * clvl[t] * (khy * ci_co2[t] - kdhy * ci_h2co3[t])
            phi[22] = qi_h2co3
            qi_co2 = qi_hco3 + qi_h2co3 + qi_co2
            phi[23] = qi_co2
        else:

            phi[5] = qe_hco3
            phi[6] = qe_h2co3
            phi[7] = qe_co2
            phi[21] = qi_hco3
            phi[22] = qi_h2co3
            phi[23] = qi_co2
        if n + 11 * Time_Durtn < t :
            qe_hpo4 = qe_hpo4 + qe_h2po4
            phi[8] = qe_hpo4

            dihydrogen_phosphate_param_e = 1
            qe_h2po4 = ebuf(lche, pkp, ce_hpo4[t], ce_h2po4[t], dihydrogen_phosphate_param_e)
            phi[9] = qe_h2po4

            qi_hpo4 = qi_hpo4 + qi_h2po4
            phi[24] = qi_hpo4
            dihydrogen_phosphate_param_i = 1
            qi_h2po4 = ebuf(lchi, pkp, ci_hpo4[t], ci_h2po4[t], dihydrogen_phosphate_param_i)
            phi[25] = qi_h2po4

        else:

            phi[8] = qe_hpo4
            phi[9] = qe_h2po4
            phi[24] = qi_hpo4
            phi[25] = qi_h2po4
        if n + 12 * Time_Durtn < t :
            qe_nh3 = qe_nh3 + qe_nh4
            phi[11] = qe_nh3
            ammonium_param_e = 1
            qe_nh4 = ebuf(lche, pkn, ce_nh3[t], ce_nh4[t], ammonium_param_e)
            phi[12] = qe_nh4
            qi_nh3 = qi_nh3 + qi_nh4 - 1e6 * qiamm
            phi[27] = qi_nh3
            ammonium_param_i = 1
            qi_nh4 = ebuf(lchi, pkn, ci_nh3[t], ci_nh4[t], ammonium_param_i)
            phi[28] = qi_nh4
        else:

            phi[11] = qe_nh3
            phi[12] = qe_nh4
            phi[27] = qi_nh3
            phi[28] = qi_nh4

        if n + 13 * Time_Durtn < t :
            qe_hco2 = qe_hco2 + qe_h2co2
            phi[13] = qe_hco2
            dihydroxymethylidene_param_e = 1
            qe_h2co2 = ebuf(lche, pkf, ce_hco2[t], ce_h2co2[t], dihydroxymethylidene_param_e)
            phi[14] = qe_h2co2
            qi_hco2 = qi_hco2 + qi_h2co2
            phi[29] = qi_hco2
            dihydroxymethylidene_param_i = 1
            qi_h2co2 = ebuf(lchi, pkf, ci_hco2[t], ci_h2co2[t], dihydroxymethylidene_param_i)
            phi[30] = qi_h2co2

        else:
            phi[13] = qe_hco2
            phi[14] = qe_h2co2
            phi[29] = qi_hco2
            phi[30] = qi_h2co2

        # cell buffer content and ph equilibrium
    qi_cb = cbuf[t] + hcbuf[t] - (tbuf * clvl0 / clvl[t])
    phi[32] = qi_cb
    c_eQ_h2po4 = math.log10(abs(
        (ci_hpo4[t] * hcbuf[t]) / (ci_h2po4[t] * cbuf[t]))) if ci_h2po4[t] * cbuf[t] == 0 or (
            (ci_hpo4[t] * hcbuf[t]) / (ci_h2po4[t] * cbuf[t]) <= 0) else math.log10(
        (ci_hpo4[t] * hcbuf[t]) / (ci_h2po4[t] * cbuf[t]))
    phi_ph_eQ = pkb - pkp - c_eQ_h2po4
    phi[33] = phi_ph_eQ

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
Nfac = 10
t0 = 0
tf = 3000
T = int(tf * Nfac)
dt = float(tf - t0) / float(T)
print('dt', dt)
rtau = 1/dt
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
p_i = matrix(0.66971984470614129292e-01, T)
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
ce_h = matrix(4.59e-11, T)
ci_h = matrix(4.69e-11, T)
ae = matrix(0.02000, T)
chvl = matrix(0.7000e-04, T)
clvl = matrix(0.1000e-02, T)
l = matrix(0.1000e-02, T)
imp = matrix(0.6000e-01, T)
rm = matrix(0.1250e-02, T)
am = matrix(0, T)
phi = matrix(0, 35)
lchm = pkc + math.log10(cm_hco3 / cm_h2co3)
lchs = pkc + math.log10(cs_hco3 / cs_h2co3)
cm_h = 10. ** (-lchm)
cs_h = 10. ** (-lchs)
ps = 9.00
vs = 0.00


vars = [ve, pe, ce_na, ce_k,ce_cl, ce_hco3, ce_h2co3, ce_co2, ce_hpo4, ce_h2po4, ce_urea, ce_nh3, ce_nh4, ce_hco2,
       ce_h2co2, ce_gluc,
       vi, imp, ci_na, ci_k, ci_cl, ci_hco3, ci_h2co3, ci_co2, ci_hpo4, ci_h2po4, ci_urea, ci_nh3, ci_nh4, ci_hco2,
       ci_h2co2, ci_gluc, cbuf, hcbuf, vm]

t = 0
guess =[var[t] for var in vars]
# The finite difference method was applied, stephenson (1978)
# Analysis of the transient behavior of the kidney model
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while True:
    t += 1
    if t == T:
        break
    else:
        cm_na = 0.11
        result = scipy.optimize.root(eQs, np.array(guess), 1)
        guess =[var[t] for var in vars]

t = np.linspace(t0, tf, T)
plot_concen = 1
if plot_concen:
    fig2 = plt.figure(constrained_layout=True)
    spec2 = gridspec.GridSpec(ncols=1, nrows=8, figure=fig2)
    ax2 = fig2.add_subplot(spec2[0, 0])
    ax2.plot(t, ci_na, 'b-')
    ax2.set_title('Solute Concentration OF THE CELL[mmol] ')
    ax2.set_ylabel(' ci_na ', color='blue')
    ax3 = fig2.add_subplot(spec2[1, 0])
    ax3.plot(t, ci_k, 'b-')
    ax3.set_ylabel('ci_k ', color='blue')
    ax4 = fig2.add_subplot(spec2[2, 0])
    ax4.plot(t, ci_cl, 'b-')
    ax4.set_ylabel('ci_cl ', color='blue')
    ax5 = fig2.add_subplot(spec2[3, 0])
    ax5.plot(t, ci_urea, 'b-')
    ax5.set_ylabel('ci_urea ', color='blue')
    ax6 = fig2.add_subplot(spec2[4, 0])
    ax6.plot(t, ci_nh3, 'b-')
    ax6.set_ylabel('ci_nh3', color='blue')
    ax7 = fig2.add_subplot(spec2[5, 0])
    ax7.plot(t, ci_nh4, 'b-')
    ax7.set_ylabel('ci_nh4', color='blue')
    ax8 = fig2.add_subplot(spec2[6, 0])
    ax8.plot(t, ci_h, 'b-')
    ax8.set_ylabel('ci_h ', color='blue')
    ax8.set_xlabel('time (s) ', color='blue')

    fig3 = plt.figure(constrained_layout=True)
    spec2 = gridspec.GridSpec(ncols=1, nrows=7, figure=fig3)
    ax2 = fig3.add_subplot(spec2[0, 0])
    ax2.plot(t, ci_gluc, 'b-')
    ax2.set_ylabel('ci_gluc ', color='blue')
    ax2.set_title('Solute Concentration OF THE CELL[mmol] ')
    ax3 = fig3.add_subplot(spec2[1, 0])
    ax3.plot(t, ci_hco3, 'b-')
    ax3.set_ylabel('ci_hco3 ', color='blue')
    ax4 = fig3.add_subplot(spec2[2, 0])
    ax4.plot(t, ci_h2co3, 'b-')
    ax4.set_ylabel('ci_h2co3 ', color='blue')
    ax5 = fig3.add_subplot(spec2[3, 0])
    ax5.plot(t, ci_hpo4, 'b-')
    ax5.set_ylabel('ci_hpo4 ', color='blue')
    ax6 = fig3.add_subplot(spec2[4, 0])
    ax6.plot(t, ci_h2po4, 'b-')
    ax6.set_ylabel('ci_h2po4', color='blue')
    ax7 = fig3.add_subplot(spec2[5, 0])
    ax7.plot(t, ci_hco2, 'b-')
    ax7.set_ylabel('ci_hco2', color='blue')
    ax8 = fig3.add_subplot(spec2[6, 0])
    ax8.plot(t, ci_h2co2, 'b-')
    ax8.set_ylabel('ci_h2co2 ', color='blue')
    ax8.set_xlabel('time (s) ', color='blue')

    fig5 = plt.figure(constrained_layout=True)
    spec2 = gridspec.GridSpec(ncols=1, nrows=7, figure=fig5)
    ax2 = fig5.add_subplot(spec2[0, 0])
    ax2.plot(t, ce_na, 'b-')
    ax2.set_title('Solute Concentration OF THE CHANNEL[mmol]')
    ax2.set_ylabel(' ce_na ', color='blue')
    ax3 = fig5.add_subplot(spec2[1, 0])
    ax3.plot(t, ce_k, 'b-')
    ax3.set_ylabel('ce_k ', color='blue')
    ax4 = fig5.add_subplot(spec2[2, 0])
    ax4.plot(t, ce_cl, 'b-')
    ax4.set_ylabel('ce_cl ', color='blue')
    ax5 = fig5.add_subplot(spec2[3, 0])
    ax5.plot(t, ce_urea, 'b-')
    ax5.set_ylabel('ce_urea ', color='blue')
    ax6 = fig5.add_subplot(spec2[4, 0])
    ax6.plot(t, ce_nh3, 'b-')
    ax6.set_ylabel('ce_nh3', color='blue')
    ax7 = fig5.add_subplot(spec2[5, 0])
    ax7.plot(t, ce_nh4, 'b-')
    ax7.set_ylabel('ce_nh4', color='blue')
    ax8 = fig5.add_subplot(spec2[6, 0])
    ax8.plot(t, ce_h, 'b-')
    ax8.set_ylabel('ce_h ', color='blue')
    ax8.set_xlabel('time (s) ', color='blue')
    fig3 = plt.figure(constrained_layout=True)
    spec2 = gridspec.GridSpec(ncols=1, nrows=8, figure=fig3)
    ax2 = fig3.add_subplot(spec2[0, 0])
    ax2.plot(t, ce_gluc, 'b-')
    ax2.set_ylabel('ce_gluc ', color='blue')
    ax2.set_title('Solute Concentration OF THE CELL[mmol] ')
    ax3 = fig3.add_subplot(spec2[1, 0])
    ax3.plot(t, ce_hco3, 'b-')
    ax3.set_ylabel('ce_hco3 ', color='blue')
    ax4 = fig3.add_subplot(spec2[2, 0])
    ax4.plot(t, ce_h2co3, 'b-')
    ax4.set_ylabel('ce_h2co3 ', color='blue')
    ax5 = fig3.add_subplot(spec2[3, 0])
    ax5.plot(t, ce_hpo4, 'b-')
    ax5.set_ylabel('ce_hpo4 ', color='blue')
    ax6 = fig3.add_subplot(spec2[4, 0])
    ax6.plot(t, ce_h2po4, 'b-')
    ax6.set_ylabel('ce_h2po4', color='blue')
    ax7 = fig3.add_subplot(spec2[5, 0])
    ax7.plot(t, ce_hco2, 'b-')
    ax7.set_ylabel('ce_hco2', color='blue')
    ax8 = fig3.add_subplot(spec2[6, 0])
    ax8.plot(t, ce_h2co2, 'b-')
    ax8.set_ylabel('ce_h2co2 ', color='blue')
    ax8.set_xlabel('time (s) ', color='blue')

    fig6 = plt.figure(constrained_layout=True)
    spec2 = gridspec.GridSpec(ncols=1, nrows=6, figure=fig6)
    ax5 = fig6.add_subplot(spec2[0, 0])
    ax5.plot(t, imp, 'r-', linewidth=2.0)
    ax5.set_xlim(-5, tf)
    ax5.set_ylabel('imp(mm) ', color='blue')
    ax7 = fig6.add_subplot(spec2[1, 0])
    ax7.plot(t, pe, 'b-', linewidth=2.0)
    ax7.set_ylabel('pressure(mm.Hg) ', color='blue')
    ax7.set_xlabel('time (s) ', color='blue')
    ax7.set_xlim(-5, tf)
    ax6 = fig6.add_subplot(spec2[2, 0])
    ax6.set_title('MEMBRANE POTENTIAL[mVolt]')
    ve = np.asarray(ve)

    ax6.plot(t, ve, 'b-', linewidth=2.0)
    ax6.set_ylabel('voltage_em (mv) ', color='blue')
    ax6.set_xlabel('time (s) ', color='blue')
    ax6.set_xlim(-5, tf)
    ax7 = fig6.add_subplot(spec2[3, 0])
    vie = np.asarray(vi) - np.asarray(ve)
    v_ie = np.asarray(vie)
    ax7.plot(t, v_ie, 'b-', linewidth=2.0)
    ax7.set_ylabel('voltage_ie (mv) ', color='blue')
    ax7.set_xlabel('time (s) ', color='blue')
    ax7.set_xlim(-5, tf)
    ax7 = fig6.add_subplot(spec2[4, 0])
    vm = np.asarray(vm)

    ax7.plot(t, vm, 'b-', linewidth=2.0)
    ax7.set_ylabel('voltage_m (mv) ', color='blue')
    ax7.set_xlabel('time (s) ', color='blue')
    ax7.set_xlim(-5, tf)
    ax7 = fig6.add_subplot(spec2[5, 0])
    vi = np.asarray(vi)
    ax7.plot(t, vi, 'b-', linewidth=2.0)
    ax7.set_ylabel('voltage_i (mv) ', color='blue')
    ax7.set_xlabel('time (s) ', color='blue')
    ax7.set_xlim(-5, tf)
    plt.show()
