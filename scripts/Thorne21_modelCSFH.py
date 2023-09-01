import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
rng = np.random.default_rng()


dlogM = 0.01


def load_data():

    thorne21_gsmf = pd.read_table(
        "/Users/22252335/Documents/PROJ2_JWST_CSFH/data/literature/Thorne2021_SuppMaterial/TableD1_GSMFValues.txt"
    )

    thorne21_sfms = pd.read_table(
        "/Users/22252335/Documents/PROJ2_JWST_CSFH/data/literature/Thorne2021_SuppMaterial/TableD3_SFMSValues.txt"
    )

    print(thorne21_gsmf.keys())

    test_mstar = np.arange(5.0, 12.0, dlogM)
    return(
        {
            "test_mstar" : test_mstar,
            "thorne_gsmf" : thorne21_gsmf,
            "thorne_sfms" : thorne21_sfms
        }
    )


def double_schechter(x, parm):

    phi1, phi2, a1, a2, mstar = parm

    mu = pow(10, x) / pow(10, mstar)

    phi = np.log(10) * np.exp(-mu) * ( phi1*pow(mu, a1+1) + phi2*pow(mu, a2+1) )

    return(phi)

def double_powerlaw_sfms(x, parm):

    S0, M0, a, b = parm

    mu = pow(10, x) / pow(10, M0)

    logsfr = S0 - np.log10( pow(mu, -1*a) + pow(mu, -1*b) )

    return(logsfr)

def do_stuff(data):

    thorne_gsmf = data["thorne_gsmf"]
    thorne_sfms = data["thorne_sfms"]
    test_mstar = data["test_mstar"]

    CSFH_thorne = []
    CSFH_thorne_hi = []
    CSFH_thorne_lo = []

    for i in range(18):
        print(thorne_gsmf["zmed"][i])


        sample_appends = []
        for j in range(100):
            gsmf_parm = rng.normal(
                [
                    thorne_gsmf["log10(phi1)"][i],
                    thorne_gsmf["log10(phi2)"][i],
                    thorne_gsmf["alpha1"][i],
                    thorne_gsmf["alpha2"][i],
                    thorne_gsmf["Mstar"][i]
                ],
                [
                    thorne_gsmf["log10(phi1)_sd"][i],
                    thorne_gsmf["log10(phi2_sd)"][i],
                    thorne_gsmf["alpha1_sd"][i],
                    thorne_gsmf["alpha2_sd"][i],
                    thorne_gsmf["Mstar_sd"][i]
                ]
            )
            gsmf_parm[0] = pow(10, gsmf_parm[0])
            gsmf_parm[1] = pow(10, gsmf_parm[1])

            foo = double_schechter(test_mstar,
                               gsmf_parm)

            sfms_parm = rng.normal(
                [
                    thorne_sfms["S0"][i],
                    thorne_sfms["M0"][i],
                    thorne_sfms["alpha"][i],
                    thorne_sfms["beta"][i]
                ],
                [
                    thorne_sfms["S0_sd"][i],
                    thorne_sfms["M0_sd"][i],
                    thorne_sfms["alpha_sd"][i],
                    thorne_sfms["beta_sd"][i]
                ]
            )
            bar = double_powerlaw_sfms(test_mstar,
                               sfms_parm)

            csfh_sample = np.sum(
                foo * pow(10, bar)
            ) * dlogM

            sample_appends.append(csfh_sample)

        csfh=np.median(sample_appends)
        csfh_hi = np.quantile(sample_appends, 0.84)
        csfh_lo = np.quantile(sample_appends, 0.16)
        CSFH_thorne.append(csfh)
        CSFH_thorne_hi.append(csfh_hi)
        CSFH_thorne_lo.append(csfh_lo)


    thorne_csfh_fit_data = pd.DataFrame(
        {
            "z" : thorne_gsmf["zmed"],
            "csfh" : np.log10(CSFH_thorne),
            "csfh_lo" : np.log10(CSFH_thorne_lo),
            "csfh_hi" : np.log10(CSFH_thorne_hi),
        }
    )
    
    thorne_csfh_fit_data.to_csv(
        "/Users/22252335/Documents/PROJ2_JWST_CSFH/data/literature/thorne21_fitting_csfh.csv"
    )

if __name__ == "__main__":
    data = load_data()
    do_stuff(data)