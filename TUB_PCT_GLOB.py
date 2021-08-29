# author: leyla noroozbabaee
# created on mon April 20 11:11 2020
# globvar = 0
#   rt-   gas const. times temp  (mmhg.ml/mmol)
rte = 2.57  # rte-    gas const. times temp  (Joule/mmol)
rt = 1.93e+4
f = 0.965e+5  # f-      faraday ()
#f = 96.5
rt_f = 0.0266362647
z_na = 1
z_k = 1
z_cl = -1
z_hco3 = -1  # z-   valence of i'th solute
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

# TUBULE MODEL PARAMETERS
tau = 1
tlim = 205
chop = 80
#tubule epsi = 0.2e-8
epsi = 0.2e-8
# epithelial epsi = 0.2e-9
tl = 0.1e+1
#
rm0 = 0.1250e-02 # rigid tubule radiusin tubule (cm)
#rm0 = 0.1060e-02 # compliant tubule radius in epithelial model
mum = 0.1000e+01
eta = 0.6400e-05
khy_5 = 0.1450e+04
kdhy_5 = 0.4960e+06
#
ame = 0.2000e-03
ae0 = 0.2000e-01
# # ae0 = 0.2000e-01 that is the real value
mua = 0.1000e+00
khy_4 = 0.1450e+04
kdhy_4 = 0.4960e+06

l0 = 0.1000e-02
chvl0 = 0.7000e-04
muv = 0.1000e+00
lpme = 0.2000e+02
lpes = 0.6000e+01

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


#luminal model np = 0.6000e-06
#epithelial model np = 0.3000e-06
n_p = 0.6000e-06
knh4 = 0.1000e+01
nphk = 0.0000e+00

#luminal model lhp = 0.1e-6
#epithelial model lhp = 0.5e-7
lhp = 0.1e-6
xihp = 0.4000e+00
xhp = 0.1450e+01
#luminal model nnhe3 = 0.5500e-08
#epithelial model nnhe3 = 0.27500e-08
nae1 = 0.0000e+00
ntsc = 0.0000e+00
nnhe3 = 0.5500e-08

#luminal model  lpmi = 0.4000e-03
#epithelial model lpmi = 0.2000e-03
lpmi = 0.4000e-03
#luminal model  lpis = 0.4000e-03
#epithelial model lpis = 0.2000e-03
lpis = 0.4000e-03

smi_na = 1.000
sis_na = 1.000
hmi_na = 0.0000e+00
#luminal model  his_na = 0.7800e-08
#epithelial model his_na = 0.3900e-08
his_na = 0.7800e-08

smi_k = 1.000
sis_k = 1.000
#luminal model  hmi_k = 0.500e-06  * his_k = 0.4000e-05
#epithelial model hmi_k = 0.2500e-06  * his_k = 0.2000e-05
hmi_k = 0.500e-06
his_k = 0.4000e-05

smi_cl = 1.000
sis_cl = 1.000
hmi_cl = 0.0000e+00
his_cl = 0.0000e+00

smi_hco3 = 1.000
sis_hco3 = 1.000
#luminal model  hmi_hco3 = 0.2000e-07
#epithelial model hmi_hco3 = 0.1000e-07
hmi_hco3 = 0.2000e-07
his_hco3 = 0.0000e+00

smi_h2co3 = 1.000
sis_h2co3 = 1.000
#luminal model   hmi_h2co3 = 0.130e-02* his_h2co3 = 0.130e-02
#epithelial model  hmi_h2co3 = 0.750e-01*  his_h2co3 = 0.750e-01
hmi_h2co3 = 0.130e-02
his_h2co3 = 0.130e-02

smi_co2 = 1.000
sis_co2 = 1.000
#luminal model   hmi_co2 = 0.2400e-01  his_co2 = 0.2400e-01
#epithelial model  hmi_co2 = 0.2400e-01 his_co2 = 0.2400e-01
hmi_co2 = 0.2400e-01
his_co2 = 0.2400e-01

smi_hpo4 = 1.000
sis_hpo4 = 1.000
#luminal model  hmi_hpo4 = 0.1900e-07 his_hpo4 = 0.450e-07
#epithelial model hmi_hpo4 = 0.950e-08 his_hpo4 = 0.2250e-07
hmi_hpo4 = 0.1900e-07
his_hpo4 = 0.450e-07

smi_h2po4 = 1.000
sis_h2po4 = 1.000
#luminal model  hmi_h2po4 = 0.0000e+00 his_h2po4 = 0.6600e-06
#epithelial model hmi_h2po4 = 0.0000e+00 his_h2po4 = 0.3300e-06
hmi_h2po4 = 0.0000e+00
his_h2po4 = 0.6600e-06

smi_urea = 0.950
sis_urea = 0.950
#luminal model  hmi_urea = 0.2100e-05 his_urea = 0.2100e-05
#epithelial model hmi_urea = 0.1050e-05 his_urea = 0.1000e-05
hmi_urea = 0.2100e-05
his_urea = 0.2000e-05

smi_nh3 = 0.500
sis_nh3 = 0.500
#luminal model  hmi_nh3 = 0.1700e-02  his_nh3 = 0.2000e-02
#epithelial model hmi_nh3 = 0.850e-03  his_nh3 = 0.1000e-02
hmi_nh3 = 0.1700e-02
his_nh3 = 0.2000e-02

smi_nh4 = 1.000
sis_nh4 = 1.000
#luminal model  hmi_nh4 = 0.430e-06 his_nh4 = 0.12000e-05
#epithelial model hmi_nh4 = 0.215e-06 his_nh4 = 0.6000e-06
hmi_nh4 = 0.430e-06
his_nh4 = 0.12000e-05

smi_h = 1.000
sis_h = 1.000
#luminal model   hmi_h = 0.1700e-01 his_h = 0.1700e-01
#epithelial model hmi_h = 0.850e-02 his_h = 0.850e-02
hmi_h = 0.1700e-01
his_h = 0.1700e-01

smi_hco2 = 1.000
sis_hco2 = 1.000
#luminal model   hmi_hco2 = 0.0000e+00 his_hco2 = 0.3800e-06
#epithelial model hmi_hco2 = 0.0000e+00 his_hco2 = 0.1900e-06
hmi_hco2 = 0.0000e+00
his_hco2 = 0.3800e-06

smi_h2co2 = 0.950
sis_h2co2 = 0.950
#luminal model   hmi_h2co2 = 0.1000e00  his_h2co2 = 0.12000e00
#epithelial model  hmi_h2co2 = 0.1000e00  his_h2co2 = 0.6000e-01
hmi_h2co2 = 0.1000e00
his_h2co2 = 0.12000

#luminal model   hmi_gluc = 0.0000e+00  his_gluc = 0.1500e-04
#epithelial model  hmi_gluc = 0.0000e+00  his_gluc = 0.750e-05
smi_gluc = 1.000
sis_gluc = 1.000
hmi_gluc = 0.0000e+00
his_gluc = 0.1500e-04


# luminal model  net cotransporters mi
lmi_clhco3 = 0.4000e-08  #(3,4)
lmi_clhco2 = 0.1000e-07  #(3,13)
lmi_nah2po4 = 0.2500e-08 #(1,8)
lmi_nagluc = 0.1500e-07  #(1,15)


# epithelial model net cotransporters mi
#lmi_clhco3 = 0.2000e-08
#lmi_clhco2 = 0.5000e-08
#lmi_nah2po4 = 0.1250e-08
#lmi_nagluc = 0.7500e-08

##multisegment model net cotransporters mi
# lmi_nah2co3 = 0.1500e-07
# lmi_clhco3 = 0.14400e-08
# lmi_clhco2 = 0.36000e-08
# lmi_nah2po4 = 9.0e-08
# lmi_nagluc = 5.40e-08

# luminal model net cotransporters is
lis_kcl = 0.1000e-07#(2,3) #must be compared with the source vertically
lis_nahco3 = 0.300e-07#(1,4)
lis_na_clhco3 = 0.1400e-06#(1,3,4)
lis_nacl = 0.700e-07#(1,3)

# epithelial model net cotransporters mi
# lis_kcl = 0.50e-08
# lis_nahco3 = 0.1500e-07
#lis_na_clhco3 = 0.700e-07
#lis_nacl =  0.350e-07

# multisegment model net cotransporters mi
# lis_kcl = 0.3460e-08
# lis_na_clhco3 = 0.024200e-08
# lis_nahco3 = 0.1500e-07

# TUBULE MODEL BOUNDRY CONDITION

zimpm=1.00
zimps=-0.00
pkc=3.57
pkf=3.76
pkn=9.15
pkp=6.8
kco2=340.
#
rmt0 = 1.25e-3



dgam = 1.e-4

