# author: leyla noroozbabaee
# created on fri Jul  5 14:02:47 2019

ps = 9.00
vs = 0.00
rte = 2.57  # rte-gas const. times temp  (Joule/mmol)
rt = 1.93e+4 # rt-gas const. times temp  (mmhg.ml/mmol)
f = 0.965e+5  # f-faraday ()
rt_f = 0.0266362647
# z-   valence of i'th solute
z_na = 1
z_k = 1
z_cl = -1
z_hco3 = -1
z_h2co3 = 0
z_co2 = 0
z_hpo4 = -2
z_h2po4 = -1
z_urea = 0
z_nh3 = 0
z_nh4 = 1
z_h = 1
z_hco2 = -1
z_h2co2 = 0
z_gluc = 0


tau = 0
tlim = 205
chop = 80
epsi = 0.2e-9
tl = 0.1e+1
fvm = 0.012
# rm_rigid = 0.1250e-02    # rigid tubule radius (cm)
rm0 = 0.1060e-02 # compliant tubule radius
mum = 0.1000e+01
eta = 0.6400e-05
khy_5 = 0.1450e+04
kdhy_5 = 0.4960e+06
# area of different membrane
ame = 0.2000e-03
ae0 = 0.2000e-01
mua = 0.1000e+00
khy_4 = 0.1450e+04
kdhy_4 = 0.4960e+06

l0 = 0.1000e-02
chvl0 = 0.7000e-04
muv = 0.1000e+00
lpme = 0.2000e+02
lpes = 0.6000e+01

aie = 0.3600e+02
ami = 0.3600e+02
ais = 0.1000e+01

clvl0 = 0.1000e-02
imp0 = 0.6000e-01
zimp = -0.1000e+01

khy = 0.1450e+04
kdhy = 0.4960e+06
tbuf = 0.6000e-01
pkb = 0.7500e+01
qiamm = 0.700e-7
impe = 0
imps = 0.002
impm = 0
cs_na = 0.140
cs_k = 0.0049
cs_cl = 0.11322499
cs_hco3 = 0.024
cs_h2co3 = 0.000004412
cs_co2 = 0.0015
cs_hpo4 = 0.002972159
cs_h2po4 = 0.000927841
cs_urea = 0.005
cs_nh3 = 0.000002821
cs_nh4 = 0.000197179
cs_hco2 = 0.001000029
cs_h2co2 = 0.000000285
cs_gluc = 0.005


cm_k = 0.0049
cm_cl = 0.11322499
cm_hco3 = 0.02400
cm_h2co3 = 0.000004412
cm_co2 = 0.0015
cm_hpo4 = 0.002972159
cm_h2po4 = 0.000927841
cm_urea = 0.005
cm_nh3 = 0.000002821
cm_nh4 = 0.000197179
cm_hco2 = 0.001000029
cm_h2co2 = 0.000000285
cm_gluc = 0.005



sme_na = 0.750
ses_na = 0.000
hme_na = 0.1300e+01
hes_na = 0.5000e-01

sme_k = 0.600
ses_k = 0.000
hme_k = 0.1450e+01
hes_k = 0.7000e-01

sme_cl = 0.300
ses_cl = 0.000
hme_cl = 0.1000e+01
hes_cl = 0.6000e-01

sme_hco3 = 0.900
ses_hco3 = 0.000
hme_hco3 = 0.4000e+00
hes_hco3 = 0.5000e-01

sme_h2co3 = 0.900
ses_h2co3 = 0.000
hme_h2co3 = 0.4000e+00
hes_h2co3 = 0.5000e-01

sme_co2 = 0.900
ses_co2 = 0.000
hme_co2 = 0.4000e+00
hes_co2 = 0.5000e-01

sme_hpo4 = 0.900
ses_hpo4 = 0.000
hme_hpo4 = 0.2000e+00
hes_hpo4 = 0.4000e-01

sme_h2po4 = 0.900
ses_h2po4 = 0.000
hme_h2po4 = 0.2000e+00
hes_h2po4 = 0.4000e-01

sme_urea = 0.700
ses_urea = 0.000
hme_urea = 0.4000e+00
hes_urea = 0.8000e-01

sme_nh3 = 0.300
ses_nh3 = 0.000
hme_nh3 = 0.2500e+01
hes_nh3 = 0.2000e+00

sme_nh4 = 0.600
ses_nh4 = 0.000
hme_nh4 = 0.2500e+01
hes_nh4 = 0.2000e+00

sme_h = 0.200
ses_h = 0.000
hme_h = 0.3000e+02
hes_h = 0.3000e+02

sme_hco2 = 0.300
ses_hco2 = 0.000
hme_hco2 = 0.7000e+00
hes_hco2 = 0.5000e-01

sme_h2co2 = 0.700
ses_h2co2 = 0.000
hme_h2co2 = 0.1400e+01
hes_h2co2 = 0.9000e-01

sme_gluc = 1.000
ses_gluc = 0.000
hme_gluc = 0.8000e-01
hes_gluc = 0.3000e-01



n_p = 0.3000e-06
knh4 = 0.1000e+01
nphk = 0.0000e+00
# epithelial model lhp = 0.5e-7
lhp = 0.5e-7
xihp = 0.4000e+00
xhp = 0.1450e+01
# luminal model nnhe3 = 0.5500e-08
# epithelial model nnhe3 = 0.27500e-08
nae1 = 0.0000e+00
ntsc = 0.0000e+00
nnhe3 = 0.27500e-08
# luminal model  lpmi = 0.4000e-03
# epithelial model lpmi = 0.2000e-03
# luminal model  lpis = 0.4000e-03
# epithelial model lpis = 0.2000e-03
lpmi = 0.2000e-03
lpis = 0.2000e-03
# luminal model  his_na = 0.7800e-08
# epithelial model his_na = 0.3900e-08

smi_na = 1.000
sis_na = 1.000
hmi_na = 0.0000e+00
his_na = 0.3900e-08
# luminal model  hmi_k = 0.500e-06  * his_k = 0.4000e-05
# epithelial model hmi_k = 0.2500e-06  * his_k = 0.2000e-05
smi_k = 1.000
sis_k = 1.000
hmi_k = 0.2500e-06
his_k = 0.2000e-05

smi_cl = 1.000
sis_cl = 1.000
hmi_cl = 0.0000e+00
his_cl = 0.0000e+00
# luminal model  hmi_hco3 = 0.2000e-07
# epithelial model hmi_hco3 = 0.1000e-07
## epithelial model hmi_hco3 = 0.1000e-07
smi_hco3 = 1.000
sis_hco3 = 1.000
hmi_hco3 = 0.1000e-07
his_hco3 = 0.0000e+00
# luminal model   hmi_h2co3 = 0.130e-02* his_h2co3 = 0.130e-02
# epithelial model  hmi_h2co3 = 0.750e-01*  his_h2co3 = 0.750e-01
smi_h2co3 = 1.000
sis_h2co3 = 1.000
hmi_h2co3 = 0.650e-03
his_h2co3 = 0.650e-03
# luminal model   hmi_co2 = 0.2400e-01  his_co2 = 0.2400e-01
# epithelial model  hmi_co2 = 0.1200e-01 his_co2 = 0.1200e-01
smi_co2 = 1.000
sis_co2 = 1.000
hmi_co2 = 0.7500e-01
his_co2 = 0.7500e-01
# luminal model  hmi_hpo4 = 0.1900e-07 his_hpo4 = 0.450e-07
# epithelial model hmi_hpo4 = 0.950e-08 his_hpo4 = 0.2250e-07
smi_hpo4 = 1.000
sis_hpo4 = 1.000
hmi_hpo4 = 0.950e-08
his_hpo4 = 0.2250e-07
# luminal model  hmi_h2po4 = 0.0000e+00 his_h2po4 = 0.6600e-06
# epithelial model hmi_h2po4 = 0.0000e+00 his_h2po4 = 0.3300e-06
smi_h2po4 = 1.000
sis_h2po4 = 1.000
hmi_h2po4 = 0.0000e+00
his_h2po4 = 0.3300e-06
# luminal model  hmi_urea = 0.2100e-05 his_urea = 0.2000e-05
# epithelial model hmi_urea = 0.1050e-05 his_urea = 0.1000e-05
smi_urea = 0.950
sis_urea = 0.950
hmi_urea = 0.10500e-05
his_urea = 0.1000e-05
# luminal model  hmi_nh3 = 0.1700e-02  his_nh3 = 0.2000e-02
# epithelial model hmi_nh3 = 0.850e-03  his_nh3 = 0.1000e-02
smi_nh3 = 0.500
sis_nh3 = 0.500
hmi_nh3 = 0.850e-03
his_nh3 = 0.1000e-02
# luminal model  hmi_nh4 = 0.430e-06 his_nh4 = 0.12000e-05
# epithelial model hmi_nh4 = 0.215e-06 his_nh4 = 0.6000e-06
smi_nh4 = 1.000
sis_nh4 = 1.000
hmi_nh4 = 0.2150e-06
his_nh4 = 0.600e-06
# luminal model   hmi_h = 0.1700e-01 his_h = 0.1700e-01
# epithelial model hmi_h = 0.850e-02 his_h = 0.850e-02
smi_h = 1.000
sis_h = 1.000
hmi_h = 0.850e-02
his_h = 0.850e-02
# luminal model   hmi_hco2 = 0.0000e+00 his_hco2 = 0.3800e-06
# epithelial model hmi_hco2 = 0.0000e+00 his_hco2 = 0.1900e-06
smi_hco2 = 1.000
sis_hco2 = 1.000
hmi_hco2 = 0.0000e+00
his_hco2 = 0.1900e-06
# luminal model   hmi_h2co2 = 0.1000e00  his_h2co2 = 0.12000e0
# epithelial model  hmi_h2co2 = 0.0500e00  his_h2co2 = 0.6000e-01
smi_h2co2 = 0.950
sis_h2co2 = 0.950
hmi_h2co2 = 0.05000e00
his_h2co2 = 0.600e-1

# luminal model   hmi_gluc = 0.0000e+00  his_gluc = 0.1500e-04
# epithelial model  hmi_gluc = 0.0000e+00  his_gluc = 0.750e-05
smi_gluc = 1.000
sis_gluc = 1.000
hmi_gluc = 0.0000e+00
his_gluc = 0.750e-05

# boundry condition
zimpm = 1.00
zimps = -0.00
pkc = 3.57
pkf = 3.76
pkn = 9.15
pkp = 6.8
kco2 = 340.

vm = 0.0000
vs = 0.0000
vs_2 = 0.0000
pm = 15.000
ps = 9.0000

rmt0 = 1.25e-3
dgam = 1.e-4
# luminal model  net cotransporters mi
# lmi_clhco3 = 0.4000e-08  # (3,4)
# lmi_clhco2 = 0.1000e-07  # (3,13)
# lmi_nah2po4 = 0.2500e-08  # (1,8)
# lmi_nagluc = 0.1500e-07  # (1,15)

# epithelial model net cotransporters mi
lmi_clhco3 = 0.2000e-08
lmi_clhco2 = 0.5000e-08
lmi_nah2po4 = 0.1250e-08
lmi_nagluc = 0.7500e-08

# multisegment model net cotransporters mi
# lmi_nah2co3 = 0.1500e-07
# lmi_clhco3 = 0.14400e-08
# lmi_clhco2 = 0.36000e-08
# lmi_nah2po4 = 9.0e-08
# lmi_nagluc = 5.40e-08

# luminal model net cotransporters is
# lis_kcl = 0.1000e-07  # (2,3) #must be compared with the source vertically
# lis_nahco3 = 0.300e-07  # (1,4)
# lis_na_clhco3 = 0.1400e-06  # (1,3,4)
# lis_nacl = 0.700e-07  # (1,3)

# epithelial model net cotransporters is
lis_kcl = 0.50e-08
lis_nahco3 = 0.1500e-07
lis_na_clhco3 = 0.1400e-06
lis_nacl = 0.350e-07

# multisegment model net cotransporters mi
## lis_kcl = 0.3460e-08
# lis_na_clhco3 = 0.024200e-08
# lis_nahco3 = 0.1500e-07
# lis_na_clhco3 = 0.300e-07


# author: leyla noroozbabaee
# created on Monday 14/12/2020
# The soulte permeabilities for the paper Flow dependence transport for the tight junction is five times smaller than
# that in their other work of Weinstein (check the PCT_GLOB_new.py).
# The area for tight junction is  incresaed from 0.2e-3 to 0.001. The get the soulte permeabilities represented in page five ,one needs to multiply those numbers by factor of 1.e-5 and then devide them by allocated area. the soulte
# permeabilities will be the same as what is presented in PCT_GLOB_new.py file except for the tight junction memberane.
#
# num_solute = 15
# ps = 9.00
# vs = 0.00
# # z-   valence of i'th solute
# z_na = 1
# z_k = 1
# z_cl = -1
# z_hco3 = -1
# z_h2co3 = 0
# z_co2 = 0
# z_hpo4 = -2
# z_h2po4 = -1
# z_urea = 0
# z_nh3 = 0
# z_nh4 = 1
# z_h = 1
# z_hco2 = -1
# z_h2co2 = 0
# z_gluc = 0
# # Constant parameters for PCT model
# rte = 2.57  # rte-gas const. times temp  [Joule/mmol]
# rt = 1.93e+4 # rt-gas const. times temp  [mmhg.ml/mmol]
# f = 0.965e+5 # f-faraday [C/mol]
# f_current =  0.965e+5# f to be used for electrical balance [C/mol]
# pkc = 3.57
# pkf = 3.76
# pkn = 9.15
# pkp = 6.8
# khy = 0.1450e+04
# kdhy = 0.4960e+06
#
# # To get the values represented in Chloride transport in a mathematical model
# # of the rat proximal tubule, we need to multiply these numbers by related
# # membrane area and then devide by two
# # epithelial model net cotransporters mi
#
# # example : the number for na-glucose coupled transport coefficients (Flow independent transport 2007) is 0.27e-6.
# # we multiply lmi_nagluc = 0.7500e-08 by ami area, the result will be 0.27e-6. the number represented in Chloride
# # transport in a mathematical model of the rat proximal tubule, is
# lmi_nagluc = 0.7500e-08  # coupled transport coefficients [mmol2/Joule.s]
# lmi_nah2po4 = 0.1250e-08
# lmi_clhco3 = 0.2000e-08
# lmi_clhco2 = 0.5000e-08
#
# # Weintein introduced two simple exchangers in their work mentioned above and it seems to us these two
# # simple exchangers are equivalent to NHE transporter that they introduced in 1995.
# lmi_nah = 0.225e-7
# lmi_nanh4 = 0.15e-8
#
# # epithelial model net cotransporters is
# lis_kcl = 0.50e-08
# lis_nahco3 = 0.1500e-07
# lis_na_clhco3 = 0.1400e-06
#
# lie_kcl = 0.50e-08
# lie_nahco3 = 0.1500e-07
# lie_na_clhco3 = 0.1400e-06
# lie_nacl = 0.350e-07
# # H-ATPase constant parameters
# # xihp : the steepness coefficient
# # xhp : the potential difference for half maximal flux
# # lhp : a maximal rate of transport
# # luminal model lhp = 3.6e-6
# lhp = 0.5e-7
# xihp = 0.4000e+00
# xhp = 0.1450e+01
#
#
#
# tbuf = 0.6000e-01
# pkb = 0.7500e+01
#
#
#
# tau = 0
# tlim = 205
# chop = 80
# epsi = 0.2e-9
# tl = 0.1e+1
# fvm = 0.012
# # rm_rigid = 0.1250e-02    # rigid tubule radius (cm)
# # rm0 = 0.1060e-02 # compliant tubule radius
#
# eta = 0.6400e-05
#
# # area of different membrane with the unit of [cm2/cm2. epithelium]
# ame = 0.001000
# ae0 = 0.2000e-01 #aes0
# aie = 0.3600e+02
# ami = 0.3600e+02
# ais = 0.1000e+01
#
# mua = 0.1000e+00 # compliant area coefficient in LIS basement membrane [1/mmHg]
# muv = 0.1000e+00 # compliant volume coefficient in Interspace Volume [1/mmHg]
# mum = 0.1000e+01
#
# l0 = 0.1000e-02
# chvl0 = 0.7000e-04 # reference cell volume [cm3/cm2. epithelium]
# clvl0 = 0.1000e-02 # reference interspace volume [cm3/cm2. epithelium]
#
#
# imp0 = 0.6000e-01 # reference cell impermeant concentration[mmol/cm3]
# impe = 0          # interspace (e) impermeant concentration[mmol/cm3]
# imps = 0.002      # peritubular (s) impermeant concentration[mmol/cm3]
# impm = 0          # lumen (m) impermeat concentration[mmol/cm3]
# # boundry condition
# zimp = -0.100e+01 # cell impermeant valence
# zimpm = 1.00
# zimps = 0.00
# kco2 = 340.
# pm = 15.000  #[mmHg]
#
#
#
# tbuf = 0.6000e-01
# pkb = 0.7500e+01
# qiamm = 0.700e-7  # Ammonia generation[mmol/s.cm2]
#
# khy = 0.1450e+04  # hydration constant for CO2
# kdhy = 0.4960e+06 # dehydration constant for CO2
#
# cs_na = 0.140
# cs_k = 0.0049
# cs_cl = 0.11322499
# cs_hco3 = 0.024
# cs_h2co3 = 0.000004412
# cs_co2 = 0.0015
# cs_hpo4 = 0.002972159
# cs_h2po4 = 0.000927841
# cs_urea = 0.005
# cs_nh3 = 0.000002821
# cs_nh4 = 0.000197179
# cs_hco2 = 0.001000029
# cs_h2co2 = 0.000000285
# cs_gluc = 0.005
#
# cm_na = 0.140
# cm_k = 0.0049
# cm_cl = 0.11322499
# cm_hco3 = 0.02400
# cm_h2co3 = 0.000004412
# cm_co2 = 0.0015
# cm_hpo4 = 0.002972159
# cm_h2po4 = 0.000927841
# cm_urea = 0.005
# cm_nh3 = 0.000002821
# cm_nh4 = 0.000197179
# cm_hco2 = 0.001000029
# cm_h2co2 = 0.000000285
# cm_gluc = 0.005
#
#
# # luminal model  his_na = 0.7800e-08
# # epithelial model his_na = 0.3900e-08
#
# smi_na = 1.000
# sis_na = 1.000
# sie_na = 1.000
# sme_na = 0.750
# ses_na = 0.000
#
#
# def convert_permeability(x):
#     fact = 1.e-5
#     a ={"me":0.001, "es":0.02, "mi":36, "ie": 36, "is": 1}
#     dic3= dict()
#     for key in x:
#         if key in a:
#            dic3[key] =  fact*x[key] /a[key]
#     return dic3
# # hme_na = 0.1300e+01
# # hes_na = 0.5000e-01
# # hmi_na = 0.0000e+00
# # hie_na = 0.3900e-08
# # his_na = 0.0
#
# hme_na = 0.26
# hes_na = 0.5000e-01
# hmi_na = 0.0000e+00
# hie_na = 0.3900e-08
# his_na = 0.0
# h_na = {"me":26, "es":100, "mi":0.0, "ie": 0.014, "is": 0.0}
# print('h_na', convert_permeability(h_na))
#
# # luminal model  hmi_k = 0.500e-06  * his_k = 0.4000e-05
# # epithelial model hmi_k = 0.2500e-06  * his_k = 0.2000e-05
# smi_k = 1.000
# sis_k = 1.000
# sie_k = 1.000
# sme_k = 0.600
# ses_k = 0.000
# # number in PCT_GLOB_new
# # hme_k = 0.1450e+01
# # hes_k = 0.7000e-01
# # hmi_k = 0.2500e-06
# # hie_k = 0.2000e-05
# # his_k = 0.2000e-05
#
# # Flow_dependent
# hme_k = 0.29
# hes_k = 0.7000e-01
# hmi_k = 0.2500e-06
# hie_k = 0.2000e-05
# his_k = 0.2000e-05
#
# h_k = {"me":29,"es": 140,"mi":0.9, "ie":7.2,"is":0.2}
#
# print('h_k',convert_permeability(h_k))
#
# smi_cl = 1.000
# sis_cl = 1.000
# sie_cl = 1.000
# sme_cl = 0.300
# ses_cl = 0.000
#
# hme_cl = 0.2
# hes_cl = 0.6000e-01
# hmi_cl = 0.0000e+00
# hie_cl = 0.0000e+00
# his_cl = 0.0000e+00
#
# def convert_permeability(x):
#     fact = 1.e-5
#     a ={"me":0.001, "es":0.02, "mi":36, "ie": 36, "is": 1}
#     dic3= dict()
#     for key in x:
#         if key in a:
#            dic3[key] =  fact*x[key] /a[key]
#     return dic3
# h_cl = {"me":20,"es": 120,"mi":0.0, "ie":0.0,"is":0.0}
#
# print('h_cl',convert_permeability(h_cl))
# # luminal model  hmi_hco3 = 0.2000e-07
# # epithelial model hmi_hco3 = 0.1000e-07
# smi_hco3 = 1.000
# sis_hco3 = 1.000
# sie_hco3 = 1.000
# sme_hco3 = 0.900
# ses_hco3 = 0.000
#
# hme_hco3 = 0.4000e+00#/5
# hes_hco3 = 0.5000e-01
# hmi_hco3 = 0.1000e-07
# hie_hco3 = 0.0000e+00
# his_hco3 = 0.0000e+00
#
#
# # h_co3 = {"me":0.8,"es": 100,"mi":0.036, "ie":0.0,"is":0.0}
# #
# # print('h_co3',convert_permeability(h_co3))
#
# # luminal model   hmi_h2co3 = 0.130e-02 * his_h2co3 = 0.130e-02
# # epithelial model  hmi_h2co3 = 0.750e-01 *  his_h2co3 = 0.750e-01
# smi_h2co3 = 1.000
# sis_h2co3 = 1.000
# sie_h2co3 = 1.000
# sme_h2co3 = 0.900
# ses_h2co3 = 0.000
#
# hmi_h2co3 = 0.650e-03
# his_h2co3 = 0.650e-03
# hie_h2co3 = 0.650e-03
# hme_h2co3 = 0.4000e+00#/5
# hes_h2co3 = 0.5000e-01
#
# # luminal model   hmi_co2 = 0.2400e-01  his_co2 = 0.2400e-01
# # epithelial model  hmi_co2 = 0.1200e-01 his_co2 = 0.1200e-01
# smi_co2 = 1.000
# sis_co2 = 1.000
# sie_co2 = 1.000
# sme_co2 = 0.900
# ses_co2 = 0.000
#
# hmi_co2 = 0.7500e-01
# his_co2 = 0.7500e-01
# hie_co2 = 0.7500e-01
# hme_co2 = 0.4000e+00#/5
# hes_co2 = 0.5000e-01
# # luminal model  hmi_hpo4 = 0.1900e-07 his_hpo4 = 0.450e-07
# # epithelial model hmi_hpo4 = 0.950e-08 his_hpo4 = 0.2250e-07
# smi_hpo4 = 1.000
# sis_hpo4 = 1.000
# sie_hpo4 = 1.000
# sme_hpo4 = 0.900
# ses_hpo4 = 0.000
#
# hmi_hpo4 = 0.950e-08
# his_hpo4 = 0.2250e-07
# hie_hpo4 = 0.2250e-07
# hme_hpo4 = 0.2000e+00#/5
# hes_hpo4 = 0.4000e-01
# # luminal model  hmi_h2po4 = 0.0000e+00 his_h2po4 = 0.6600e-06
# # epithelial model hmi_h2po4 = 0.0000e+00 his_h2po4 = 0.3300e-06
# smi_h2po4 = 1.000
# sis_h2po4 = 1.000
# sie_h2po4 = 1.000
# sme_h2po4 = 0.900
# ses_h2po4 = 0.000
#
# hmi_h2po4 = 0.0000e+00
# his_h2po4 = 0.3300e-06
# hie_h2po4 = 0.3300e-06
# hme_h2po4 = 0.2000e+00#/5
# hes_h2po4 = 0.4000e-01
#
# # luminal model  hmi_urea = 0.2100e-05 his_urea = 0.2000e-05
# # epithelial model hmi_urea = 0.1050e-05 his_urea = 0.1000e-05
# smi_urea = 0.950
# sis_urea = 0.950
# sie_urea = 0.950
# sme_urea = 0.700
# ses_urea = 0.000
#
# hmi_urea = 0.10500e-05
# his_urea = 0.1000e-05
# hie_urea = 0.1000e-05
# hme_urea = 0.4000e+00#/5
# hes_urea = 0.8000e-01
#
# # luminal model  hmi_nh3 = 0.1700e-02  his_nh3 = 0.2000e-02
# # epithelial model hmi_nh3 = 0.850e-03  his_nh3 = 0.1000e-02
# smi_nh3 = 0.500
# sis_nh3 = 0.500
# sie_nh3 = 0.500
# sme_nh3 = 0.300
# ses_nh3 = 0.000
#
# hmi_nh3 = 0.850e-03
# his_nh3 = 0.1000e-02
# hie_nh3 = 0.1000e-02
# hme_nh3 = 0.2500e+01#/5
# hes_nh3 = 0.2000e+00
#
# # luminal model  hmi_nh4 = 0.430e-06 his_nh4 = 0.12000e-05
# # epithelial model hmi_nh4 = 0.215e-06 his_nh4 = 0.6000e-06
# smi_nh4 = 1.000
# sis_nh4 = 1.000
# sie_nh4 = 1.000
# sme_nh4 = 0.600
# ses_nh4 = 0.000
#
# hmi_nh4 = 0.2150e-06
# his_nh4 = 0.600e-06
# hie_nh4 = 0.600e-06
# hme_nh4 = 0.2500e+01#/5
# hes_nh4 = 0.2000e+00
#
# # luminal model   hmi_h = 0.1700e-01 his_h = 0.1700e-01
# # epithelial model hmi_h = 0.850e-02 his_h = 0.850e-02
# smi_h = 1.000
# sis_h = 1.000
# sie_h = 1.000
# sme_h = 0.200
# ses_h = 0.000
#
# hmi_h = 0.850e-02
# his_h = 0.850e-02
# hie_h = 0.850e-02
# hme_h = 0.3000e+02#/5
# hes_h = 0.3000e+02
#
# # luminal model   hmi_hco2 = 0.0000e+00 his_hco2 = 0.3800e-06
# # epithelial model hmi_hco2 = 0.0000e+00 his_hco2 = 0.1900e-06
# smi_hco2 = 1.000
# sis_hco2 = 1.000
# sie_hco2 = 1.000
# sme_hco2 = 0.300
# ses_hco2 = 0.000
#
# hmi_hco2 = 0.0000e+00
# his_hco2 = 0.1900e-06
# hie_hco2 = 0.1900e-06
# hme_hco2 = 0.7000e+00#/5
# hes_hco2 = 0.5000e-01
#
# # luminal model   hmi_h2co2 = 0.1000e00  his_h2co2 = 0.12000e0
# # epithelial model  hmi_h2co2 = 0.0500e00  his_h2co2 = 0.6000e-01
# smi_h2co2 = 0.950
# sis_h2co2 = 0.950
# sie_h2co2 = 0.950
# sme_h2co2 = 0.700
# ses_h2co2 = 0.000
#
# hmi_h2co2 = 0.05000e00
# his_h2co2 = 0.600e-1
# hie_h2co2 = 0.600e-1
# hme_h2co2 = 0.1400e+01#/5
# hes_h2co2 = 0.9000e-01
#
# # luminal model   hmi_gluc = 0.0000e+00  his_gluc = 0.1500e-04
# # epithelial model  hmi_gluc = 0.0000e+00  his_gluc = 0.750e-05
# smi_gluc = 1.000
# sis_gluc = 1.000
# sie_gluc = 1.000
# sme_gluc = 1.000
# ses_gluc = 0.000
#
# hmi_gluc = 0.0000e+00
# his_gluc = 0.750e-05
# hie_gluc = 0.750e-05
# hme_gluc = 0.8000e-01#/5
# hes_gluc = 0.3000e-01
#
# # Na-K pumps located on Peritubular Membrane which has the following constants
# n_p = 0.3000e-06
# knh4 = 0.1000e+01
# nphk = 0.0000e+00
# # epithelial model lhp = 0.5e-7
# # lhp = 0.5e-7
# # xihp = 0.4000e+00
# # xhp = 0.1450e+01
# # luminal model nnhe3 = 0.5500e-08
# # epithelial model nnhe3 = 0.27500e-08
# nae1 = 0.0000e+00
# ntsc = 0.0000e+00
# nnhe3 = 0.27500e-08
# # luminal model  lpmi = 0.4000e-03
# # epithelial model lpmi = 0.2000e-03
# # luminal model  lpis = 0.4000e-03
# # epithelial model lpis = 0.2000e-03
#
# # the numbers for water permeabilities are a multipication of water permeabilities and r*T.
# # lpme = 0.2000e+02/5  # lpme* rt
# # lpes = 0.6000e+01
# # lpmi = 0.2000e-03
# # lpie = 0.2000e-03
# # lpis = 0.2000e-03
# fact_lp = 1.0
# lpme = fact_lp*0.2000e+02#/5  # lpme* rt
# lpes = fact_lp*0.6000e+01
# lpmi = fact_lp*0.2000e-03
# lpie = fact_lp*0.2000e-03
# lpis = fact_lp*0.2000e-03
#
# def convert_Waterpermeability_inverse(x):
#     # converting from PCT_GLOB_new to flow-dependent
#     fact = 1
#     a ={"me":0.001, "es":0.02, "mi":36, "ie": 36, "is": 1}
#     dic3= dict()
#     for key in x:
#         if key in a:
#            dic3[key] = fact* x[key] *a[key]
#     return dic3
# lp = {"me":0.2e+2, "es":0.6e+1, "mi":0.2e-3, "ie": 0.2e-3, "is": 0.2e-3}
# print('lp',convert_Waterpermeability_inverse(lp))
# # lp {'me': 0.02, 'es': 0.12, 'mi': 0.00720, 'ie': 0.00720, 'is': 0.0002}
# def convert_Waterpermeability_PA_inverse(x):
#     # converting from flow-dependent(first definition in the table) to PCT_GLOB_new
#     fact = 1/rt
#     a = {"me": 0.001, "es": 0.02, "mi": 36, "ie": 36, "is": 1}
#     dic3 = dict()
#     for key in x:
#         if key in a:
#             dic3 [ key ] = fact * x [ key ] *a [ key ]
#     return dic3
# lp = {"me":0.2e+2, "es":0.6e+1, "mi":0.2e-3, "ie": 0.2e-3, "is": 0.2e-3}
# print('lp',convert_Waterpermeability_PA_inverse(lp))
#
#
# lp_fdt = {"me": 4.0e-3, "es": 1.2e-1, "mi": 7.2e-3, "ie": 7.2e-3, "is": 0.2e-3}
# #print('lp', convert_Waterpermeability(lp_fdt))
# def convert_Waterpermeability_PA(x):
#     # converting from flow-dependent to PCT_GLOB_new
#     fact = 1/rte
#     a = {"me": 0.001, "es": 0.02, "mi": 36, "ie": 36, "is": 1}
#     dic3 = dict()
#     for key in x:
#         if key in a:
#             dic3 [ key ] = fact * x [ key ] *a [ key ]
#     return dic3
#
#
# lp_fdt = {"me": 4.0e-3, "es": 1.2e-1, "mi": 7.2e-3, "ie": 7.2e-3, "is": 0.2e-3}
# #print('lp', convert_Waterpermeability(lp_fdt))
#
#
#
# rmt0 = 1.25e-3
# dgam = 1.e-4
# # luminal model  net cotransporters mi
# # lmi_clhco3 = 0.4000e-08  # (3,4)
# # lmi_clhco2 = 0.1000e-07  # (3,13)
# # lmi_nah2po4 = 0.2500e-08  # (1,8)
# # lmi_nagluc = 0.1500e-07  # (1,15)
#
#
#
# # multisegment model net cotransporters mi
# # lmi_nah2co3 = 0.1500e-07
# # lmi_clhco3 = 0.14400e-08
# # lmi_clhco2 = 0.36000e-08
# # lmi_nah2po4 = 9.0e-08
# # lmi_nagluc = 5.40e-08
#
# # luminal model net cotransporters is
# # lis_kcl = 0.1000e-07  # (2,3) #must be compared with the source vertically
# # lis_nahco3 = 0.300e-07  # (1,4)
# # lis_na_clhco3 = 0.1400e-06  # (1,3,4)
# # lis_nacl = 0.700e-07  # (1,3)
#
# # multisegment model net cotransporters mi
# ## lis_kcl = 0.3460e-08
# # lis_na_clhco3 = 0.024200e-08
# # lis_nahco3 = 0.1500e-07
# # lis_na_clhco3 = 0.300e-07
# #scale = 1
# # print('scale', scale) #scale =1e6
# # Nfac = 10 # input('Nfac=')
# # t0 = 0
# # tf = 560
# # T = int(560*Nfac)
# # # print('T',T)
# # dt = float(tf - t0) / float(T)
# # rtau = 1./dt
# import math
# lchm = pkc + math.log10(cm_hco3 / cm_h2co3)
# lchs = pkc + math.log10(cs_hco3 / cs_h2co3)
# cm_h = 10. ** (-lchm)
# cs_h = 10. ** (-lchs)
# C_m = [ cm_na, cm_k, cm_cl, cm_hco3, cm_h2co3, cm_co2, cm_hpo4, cm_h2po4, cm_urea, cm_nh3, cm_nh4, cs_h, cm_hco2,
#             cm_h2co2, cm_gluc ]
# C_s = [ cs_na, cs_k, cs_cl, cs_hco3, cs_h2co3, cs_co2, cs_hpo4, cs_h2po4, cs_urea, cs_nh3, cs_nh4, cs_h, cs_hco2,
#         cs_h2co2, cs_gluc ]
# # Reflection coefficients
# S_em = [ sme_na, sme_k, sme_cl, sme_hco3, sme_h2co3, sme_co2, sme_hpo4, sme_h2po4, sme_urea, sme_nh3, sme_nh4,
#          sme_h,
#          sme_hco2, sme_h2co2, sme_gluc ]
# S_es = [ ses_na, ses_k, ses_cl, ses_hco3, ses_h2co3, ses_co2, ses_hpo4, ses_h2po4, ses_urea, ses_nh3, ses_nh4,
#          ses_h,
#          ses_hco2, ses_h2co2, ses_gluc ]
# S_mi = [ smi_na, smi_k, smi_cl, smi_hco3, smi_h2co3, smi_co2, smi_hpo4, smi_h2po4, smi_urea, smi_nh3, smi_nh4,
#          smi_h, smi_hco2, smi_h2co2, smi_gluc ]
# S_is = [ sis_na, sis_k, sis_cl, sis_hco3, sis_h2co3, sis_co2, sis_hpo4, sis_h2po4, sis_urea, sis_nh3, sis_nh4,
#          sis_h, sis_hco2, sis_h2co2, sis_gluc ]
# # Permeability coefficients
# h_em = [ hme_na, hme_k, hme_cl, hme_hco3, hme_h2co3, hme_co2, hme_hpo4, hme_h2po4, hme_urea, hme_nh3, hme_nh4,
#          hme_h,
#          hme_hco2, hme_h2co2, hme_gluc ]
# h_es = [ hes_na, hes_k, hes_cl, hes_hco3, hes_h2co3, hes_co2, hes_hpo4, hes_h2po4, hes_urea, hes_nh3, hes_nh4,
#          hes_h,
#          hes_hco2, hes_h2co2, hes_gluc ]
# h_mi = [ hmi_na, hmi_k, hmi_cl, hmi_hco3, hmi_h2co3, hmi_co2, hmi_hpo4, hmi_h2po4, hmi_urea, hmi_nh3, hmi_nh4,
#          hmi_h, hmi_hco2, hmi_h2co2, hmi_gluc ]
# h_is = [ his_na, his_k, his_cl, his_hco3, his_h2co3, his_co2, his_hpo4, his_h2po4, his_urea, his_nh3, his_nh4,
#          his_h, his_hco2, his_h2co2, his_gluc ]
# # z-   valence of i'th solute
# z = [ z_na, z_k, z_cl, z_hco3, z_h2co3, z_co2, z_hpo4, z_h2po4, z_urea, z_nh3, z_nh4, z_h, z_hco2, z_h2co2, z_gluc ]
