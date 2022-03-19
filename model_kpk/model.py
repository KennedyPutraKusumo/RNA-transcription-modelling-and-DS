import numpy as np


def dcdt(y, t):
    c_rnatot = y[0]
    c_ppitot = y[1]
    c_ntptot = y[2]
    c_htot = y[3]
    c_t7rnaptot = y[4]
    c_mgtot = y[5]
    c_hepestot = y[6]

    # E3
    c_mgntp =

    # equation 8
    V_tr = k_app * c_t7rnaptot * (c_mg * c_mgntp) / ()


    # equation 1
    drna_dt = V_tr - V_deg
    # equation 2
    dppi_dt = (Nall - 1) * (V_tr - V_precip)
    # equation 3
    dntp_dt = -Nall * V_tr
    # equation 4
    dh_dt = (Nall - 1) * V_tr
    # equation 5
    dt7rnap_dt = -kd * c_t7rnaptot
    # equation 6
    dmg_dt = -2 * V_precip
    # equation 7
    dhepes_dt = 0
    return


