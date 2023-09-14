import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.gridspec as grplt
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy import stats

plt.rcParams.update({'font.size': 20})
plt.rcParams.update({'legend.fontsize':14})
plt.rcParams.update({'xtick.direction' : 'in'})
plt.rcParams.update({'ytick.direction' : 'in'})
plt.rcParams.update({'lines.markersize': 10})

main_stub = "/Users/22252335/Documents/PROJ2_JWST_CSFH/"

import splotch as splt

def quantile50(x):
    return(
        np.nanquantile(x, 0.5)
    )

def quantile16(x):
    return (
        np.nanquantile(x, 0.16)
    )
def quantile84(x):
    return (
        np.nanquantile(x, 0.84)
    )

def load_data():
    """load the data"""

    data_drive = main_stub + "data/save/"
    catalogue_drive = main_stub + "data/catalogues/"
    literature_drive = main_stub + "data/literature/"
    sed_drive = main_stub + "data/sed_data/"

    CSFH_data = pd.read_csv(
        data_drive + "/jwst_csfh_data.csv"
    )

    gama_devils_data = pd.read_csv(
        literature_drive + "/CSFH_DSILVA+23_Ver_Final.csv"
    )

    devilsz5 = pd.read_csv(
        literature_drive + "/devils_super_z5.csv"
    )

    withAGN_catalogue = pd.read_csv(
        catalogue_drive + "/ProSpect_highz_withAGN.csv"
    )

    withoutAGN_catalogue = pd.read_csv(
        catalogue_drive + "/ProSpect_highz_withoutAGN.csv"
    )

    smf_data = pd.read_csv(
        data_drive + "/smf_data.csv"
    )
    sfs_fits = pd.read_csv(
        data_drive + "/sfs_fits.csv"
    )
    sfs_z5 = pd.read_csv(
        data_drive + "/sfs_z5.csv"
    )
    sfs_z8 = pd.read_csv(
        data_drive + "/sfs_z8.csv"
    )
    sfs_z10 = pd.read_csv(
        data_drive + "/sfs_z10.csv"
    )

    mstar_z5 = pd.read_csv(
        data_drive + "/mstar_z5.csv"
    )
    mstar_z8 = pd.read_csv(
        data_drive + "/mstar_z8.csv"
    )
    mstar_z10 = pd.read_csv(
        data_drive + "/mstar_z10.csv"
    )

    mstar_fits = pd.read_csv(
        data_drive + "/mstar_mstar_fits.csv"
    )

    redshift_info = pd.read_csv(
        data_drive + "/redshift_info.csv"
    )

    harikane2023 = pd.read_csv(
        literature_drive + "/harikane_JWST_2023.csv"
    )
    bouwens2015 = pd.read_csv(
        literature_drive + "/bouwens_2015.csv"
    )
    bouwens2023 = pd.read_csv(
        literature_drive + "/bouwens_2023.csv"
    )
    adams2023_ = np.loadtxt(
        literature_drive + "/Adams23_CSFRDCompilation.dat"
    )
    adams2023 = pd.DataFrame(
        {
            "redshift" : adams2023_[[1,5,11,17],0],
            "csfh" : adams2023_[[1,5,11,17],1],
            "lo" : adams2023_[[1,5,11,17],2],
            "hi" : adams2023_[[1,5,11,17],3]
        }
    )

    thorne2021_fitting = pd.read_csv(
        literature_drive + "/thorne21_fitting_csfh.csv"
    )

    sed_data = {
        "phot" : pd.read_csv(sed_drive + "/phot_galaxy.csv"),
        "sed_withAGN" : pd.read_csv(sed_drive + "/sed_withAGN.csv"),
        "sed_withoutAGN" : pd.read_csv(sed_drive + "/sed_withoutAGN.csv"),
        "stars_withAGN" : pd.read_csv(sed_drive + "/stars_withAGN.csv"),
        "stars_withoutAGN" : pd.read_csv(sed_drive + "/stars_withoutAGN.csv"),
        "agn_withAGN" : pd.read_csv(sed_drive + "/agnlum_withAGN.csv"),
        "transmission_curves" : pd.read_csv(sed_drive + "/transmission_curves.csv")
    }
    return(
        {
            "CSFH" : CSFH_data,
            "gama_devils" : gama_devils_data,

            "mstar_data" : [mstar_z5, mstar_z8, mstar_z10],
            "mstar_fits" : mstar_fits,

            "redshift_info" : redshift_info,

            "withAGN" : withAGN_catalogue,
            "withoutAGN" : withoutAGN_catalogue,

            "sfs_fits" : sfs_fits,
            "smfs" : smf_data,
            "sfs" : [sfs_z5,sfs_z8,sfs_z10],

            "devilsz5" : devilsz5,

            "harikane2023": harikane2023,
            "bouwens2015" : bouwens2015,
            "bouwens2023" : bouwens2023,
            "adams2023" : adams2023,

            "thorne2021_fitting" : thorne2021_fitting,

            "sed" : sed_data
        }
    )

def plot_csfh(data):
    """Plot the CSFH"""

    print(data["CSFH"].keys())

    colour_palette = {
        "main": "#4A4E69",
        "lines": "#C9ADA7",
        "points": "#9A8C98",
        "errors": "#22223B"
    }

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(8, 5), constrained_layout=True)

    ## plot md14 and harikane constant SFE
    zvec = np.linspace(0, 20)

    def harikane_csfh(z):
        csfh = (61.7 * (1 + z) ** -3.13 +
                1.0 * 10 ** (0.22 * (1 + z)) +
                2.4 * 10 ** (0.5 * (1 + z) - 3.0)) ** -1
        return (
                np.log10(csfh) + np.log10(1 / 1.53)
        )

    def md2014_fit(z):
        T0 = 0.015 * (1 + z) ** 2.7
        T1 = 1 + ((1 + z) / 2.9) ** 5.6
        return np.log10(T0 / T1) + np.log10(1 / 1.53)  # 0.63 convert to Chabrier

    ax.plot(zvec, md2014_fit(zvec), color=colour_palette["lines"], linestyle="--", linewidth=2, label="M&D+14")
    ax.plot(zvec, harikane_csfh(zvec), color=colour_palette["lines"], linewidth=2, label="Harikane+22")
    ax.fill_between(
        data["thorne2021_fitting"]["z"],
        data["thorne2021_fitting"]["csfh_lo"],
        data["thorne2021_fitting"]["csfh_hi"],
        color=colour_palette["lines"],
        alpha = 0.6,
        label = "Thorne+21"
    )

## plot my old gama+devils csfh
    ax.errorbar(
        x=data["gama_devils"]["redshift"],
        y=data["gama_devils"]["csfrd"],
        yerr=[data["gama_devils"]["err_down"], data["gama_devils"]["err_up"]],
        color="grey",
        fmt=".",
        markersize=8,
        label="D'Silva+23"
    )
## plot bouwens+15 data
    ax.errorbar(
        data["bouwens2015"]["redshift"],
        data["bouwens2015"]["csfrd"],
        yerr=[data["bouwens2015"]["err_down"], data["bouwens2015"]["err_up"]],
        color=colour_palette["points"],
        fmt="x",
        label="Bouwens+15",
        alpha = 0.7
    )
## plot new harikane data
    ax.errorbar(
        data["harikane2023"]["z"],
        data["harikane2023"]["csfh"],
        yerr = [data["harikane2023"]["lo"], data["harikane2023"]["hi"]],
        color=colour_palette["points"],
        fmt = "s",
        label = "Harikane+23",
        alpha = 0.7,
    )
## plot bouwens+23 data
    ax.errorbar(
        data["bouwens2023"]["redshift"][0:2],
        data["bouwens2023"]["csfh"][0:2],
        yerr=[data["bouwens2023"]["lo"][0:2], data["bouwens2023"]["hi"][0:2]],
        color=colour_palette["points"],
        fmt="d",
        label="Bouwens+23",
        alpha = 0.7
    )

## plot Adams+23 points
    ax.errorbar(
        data["adams2023"]["redshift"],
        data["adams2023"]["csfh"],
        yerr=[data["adams2023"]["csfh"] - data["adams2023"]["lo"], data["adams2023"]["hi"] - data["adams2023"]["csfh"]],
        color=colour_palette["points"],
        fmt="*",
        label="Adams+23",
        alpha = 0.7
    )

## plot the new csfh
    csfh = data["CSFH"]

    err_up = csfh["CSFH_without_AGN_q84"]
    ax.errorbar(
        csfh["z"], csfh["CSFH_without_AGN"],
        yerr=[csfh["CSFH_without_AGN_q16"], [err_up[0], err_up[1], 0.0]],
        xerr=[csfh["z"] - csfh["z16"], csfh["z84"] - csfh["z"]],
        markersize=12,
        markerfacecolor="tab:blue",
        fmt="^",
        capsize=3,
        label="Pro Stellar",
        alpha=1.0,
        ecolor="tab:blue",
        markeredgecolor="black",
        zorder=200,
    )

    ax.errorbar(
        csfh["z"], csfh["CSFH_without_AGN"],
        yerr=[0,0,0.2],
        xerr=[csfh["z"] - csfh["zmin"], csfh["zmax"] - csfh["z"]],
        color="tab:blue",
        fmt="none",
        alpha=0.7,
        zorder=200,
        lolims = [False, False, True]
    )

    err_up = csfh["CSFH_with_AGN_q84"]
    ax.errorbar(
        csfh["z"]+0.05, csfh["CSFH_with_AGN"],
        yerr=[csfh["CSFH_with_AGN_q16"], [err_up[0], err_up[1], 0]],
        xerr=[csfh["z"] - csfh["z16"], csfh["z84"] - csfh["z"]],
        markersize = 12,
        markerfacecolor="tab:red",
        fmt="v",
        capsize=3,
        label="Pro Stellar+AGN",
        alpha=1.0,
        ecolor="tab:red",
        markeredgecolor="black",
        zorder=200
    )
    ax.errorbar(
        csfh["z"]+0.05, csfh["CSFH_with_AGN"],
        yerr = [0,0,0.2],
        xerr=[csfh["z"] - csfh["zmin"], csfh["zmax"]-csfh["z"]],
        color = "tab:red",
        fmt = "none",
        alpha = 0.7,
        zorder=200,
        lolims = [False, False, True]
    )

    ax.set_xlim([0,16.5])
    ax.set_ylim([-4.55, -0.95])
    ax.set_yticks([-4.0, -3.0, -2.0, -1.0])
    ax.set_yticklabels([r'$10^{{{:n}}}$'.format(i) for i in ax.get_yticks()])

    ax.set_xlabel("Redshift")
    ax.set_ylabel("$\\rm{ \\rho_{SFR} \\, / \\, M_{\\odot} \\, yr^{-1} \\, Mpc^{-3} }$")

    ax.legend(ncol=1, frameon=False, loc="lower left")

    fig.savefig(
        main_stub + "plots/csfh.pdf"
    )

def plot_mstar(data):
    """Plot stellar mass against stellar mass"""

    fig, ax = plt.subplots(nrows=2, ncols=3, figsize=(10, 7), constrained_layout=True, height_ratios=(2,1), sharex="col", sharey="row")

    mstar_data = data["mstar_data"]

    mstar_fits = data["mstar_fits"]

    print(mstar_fits.keys())

    norm = plt.Normalize()
    cm = plt.colormaps["RdYlBu_r"]

    withAGN_catalogue = data["withAGN"]
    redshift_info = data["redshift_info"]

    redshift_edges =[3.5, 6.5, 9.5, 12.5]
    for i in range(3):

        redshift_idx = np.where(
            (withAGN_catalogue["z"] >= redshift_edges[i]) &
            (withAGN_catalogue["z"] < redshift_edges[i+1])
        )
        AGNlum = np.log10(
            np.array(withAGN_catalogue["AGNlum50"])[redshift_idx]
        )

        ax[0, i].set_title(
            str(redshift_edges[i]) + "$\\leq z <$" + str(redshift_edges[i+1])
        )
        ax[0, i].errorbar(
            mstar_data[i]["with_agn"],
            mstar_data[i]["without_agn"],
            xerr=mstar_data[i]["with_agn_err"],
            yerr=mstar_data[i]["without_agn_err"],
            fmt = "none",
            ecolor= "black",
            alpha = 0.5
        )
        ax0=ax[0, i].scatter(
            mstar_data[i]["with_agn"],
            mstar_data[i]["without_agn"],
            edgecolors="black",
            c=AGNlum,
            vmin=min(AGNlum),
            vmax=max(AGNlum),
            cmap=cm,
            zorder = 200
        )

        ax[0, i].plot([0,100], [0,100], ls="--", color="black")

        ax[0, i].axvspan(
            0,8, hatch = "\\", facecolor="grey", edgecolor="black", alpha = 0.3
        )

        ax[0, i].set_xlim([3.5, 12.5])
        ax[0, i].set_ylim([4.5, 12.5])
        ax[0, i].set_xticks([4, 6, 8, 10, 12])
        ax[0, i].set_yticks([6, 8, 10, 12])
        ax[0,i].set_yticklabels([r'$10^{{{:n}}}$'.format(i) for i in ax[0,i].get_yticks()])
        ax[0, i].set_xticklabels([r'$10^{{{:n}}}$'.format(i) for i in ax[0, i].get_xticks()])

        mstar_bins = np.arange(3.5, 12.5, 1.5)
        mstar_mid = mstar_bins[1:] - 0.5/2.0
        q50_delta = stats.binned_statistic(
            mstar_data[i]["with_agn"],
            mstar_data[i]["delta_mstar"],
            bins = mstar_bins,
            statistic = quantile50
        )[0]
        q16_delta = stats.binned_statistic(
            mstar_data[i]["with_agn"],
            mstar_data[i]["delta_mstar"],
            bins=mstar_bins,
            statistic=quantile16
        )[0]
        q84_delta = stats.binned_statistic(
            mstar_data[i]["with_agn"],
            mstar_data[i]["delta_mstar"],
            bins=mstar_bins,
            statistic=quantile84
        )[0]

        # ax[1, i].plot(
        #     mstar_mid,
        #     q50_delta,
        #     color = "black",
        #     linewidth = 2
        # )
        # ax[1, i].fill_between(
        #     mstar_mid,
        #     q16_delta,
        #     q84_delta,
        #     color="black",
        #     alpha = 0.4
        # )

        ax[1, i].errorbar(
            mstar_data[i]["with_agn"],
            mstar_data[i]["delta_mstar"],
            xerr=mstar_data[i]["with_agn_err"],
            yerr=mstar_data[i]["delta_mstar_err"],
            fmt="none",
            ecolor = "black",
            alpha = 0.5
        )
        ax[1, i].scatter(
            mstar_data[i]["with_agn"],
            mstar_data[i]["delta_mstar"],
            edgecolors="black",
            c=AGNlum,
            vmin=min(AGNlum),
            vmax=max(AGNlum),
            alpha = 1.0,
            cmap=cm,
            zorder=200
        )

        ax[1, i].axhline(0.0, color="black", ls = "--")
        ax[1, i].set_xlim([3.5, 12.5])
        ax[1, i].set_ylim([-5, 1.0])
        ax[1, i].set_xticks([4, 6, 8, 10, 12])
        ax[1, i].set_yticks([-4, -2, 0])
        ax[1, i].set_xticklabels([r'$10^{{{:n}}}$'.format(i) for i in ax[1, i].get_xticks()])
        ax[1,i].xaxis.set_tick_params(pad=5)

    ax[0, 0].plot(
        mstar_fits["test_mstar"], mstar_fits["z5_fit"],
        color="blueviolet",
        linewidth = 2,
    )
    ax[0, 0].fill_between(
        mstar_fits["test_mstar"], mstar_fits["z5_q16"], mstar_fits["z5_q84"],
        color = "blueviolet",
        alpha = 0.5,
    )
    ax[0, 1].plot(
        mstar_fits["test_mstar"], mstar_fits["z8_fit"],
        color="blueviolet",
        linewidth=2
    )
    ax[0, 1].fill_between(
        mstar_fits["test_mstar"], mstar_fits["z8_q16"], mstar_fits["z8_q84"],
        color="blueviolet",
        alpha=0.5
    )
    ax[0, 2].plot(
        mstar_fits["test_mstar"], mstar_fits["z10_fit"],
        color="blueviolet",
        linewidth=2
    )
    ax[0, 2].fill_between(
        mstar_fits["test_mstar"], mstar_fits["z10_q16"], mstar_fits["z10_q84"],
        color="blueviolet",
        alpha=0.5
    )

    inset_axes_font_size = 16.0
    ax0_inset = inset_axes(ax[1,2], width="80%", height="5%",
                           loc="center", bbox_to_anchor=(0.0, 0.0, 1, 1),
                           bbox_transform=ax[1,2].transAxes
                           )
    cbar = fig.colorbar(ax0, cax=ax0_inset, orientation='horizontal', ticks=[35, 37, 39, 41, 43, 45])
    cbar.ax.tick_params(labelsize=inset_axes_font_size, pad=5)
    ax0_inset.text(4.5, -4.0, "$\\rm{ \\log_{10}(L_{AGN} \\, / \\, erg \\, s^{-1}) }$",
                   transform=ax[1,2].transData, fontsize=inset_axes_font_size)

    ax[0,0].set_ylabel(
        "$\\rm{ M_{\\star}^{Pro \\, Stellar} \\, / \\, M_{\\odot} }$",

    )
    ax[1,0].set_ylabel(
        "$\\rm{\\Delta M_{\\star}}$"
    )
    fig.supxlabel(
        "$\\rm{ M_{\\star}^{Pro \\, Stellar+AGN} \\, / \\, M_{\\odot} }$"
    )

    fig.savefig(
        main_stub + "plots/mstar_mstar.pdf"
    )

def plot_sfms(data):
    """Plot stellar mass against stellar mass"""

    fig, ax = plt.subplots(nrows=2, ncols=3, figsize=(10, 7), constrained_layout=True, sharex="col", sharey="row")

    sfs_fits = data["sfs_fits"]
    smfs = data["smfs"]
    sfs_data = data["sfs"]
    devilsz5 = data["devilsz5"]


    norm = plt.Normalize()
    cm = plt.colormaps["RdYlBu_r"]

    withAGN_catalogue = data["withAGN"]
    redshift_info = data["redshift_info"]

    redshift_edges =[3.5, 6.5, 9.5, 12.5]
    for i in range(3):

        ## DEVILSD10 z>5
        ax[0, i].scatter(
            np.log10(devilsz5["Mstar"]),
            np.log10(devilsz5["SFR"]),
            marker=",",
            color="grey",
            alpha = 0.5,
            s = 5.0,
            label = "DEVILSD10 $z \\geq 5$"
        )

        redshift_idx = np.where(
            (withAGN_catalogue["z"] >= redshift_edges[i]) &
            (withAGN_catalogue["z"] < redshift_edges[i+1])
        )
        AGNlum = np.log10(
            np.array(withAGN_catalogue["AGNlum50"])[redshift_idx]
        )

        ax[0, i].set_title(
            str(redshift_edges[i]) + "$\\leq z <$" + str(redshift_edges[i+1])
        )
        ax[0, i].errorbar(
            sfs_data[i]["withAGN_mstar"],
            sfs_data[i]["withAGN_sfr"],
            xerr=sfs_data[i]["withAGN_mstar_err"],
            yerr=sfs_data[i]["withAGN_sfr_err"],
            fmt = "none",
            ecolor="tab:red",
            alpha = 0.3
        )
        ax0=ax[0, i].scatter(
            sfs_data[i]["withAGN_mstar"],
            sfs_data[i]["withAGN_sfr"],
            edgecolors="tab:red",
            c="darkred",
            alpha = 0.7,
            label = "Pro Stellar+AGN"
        )

        ax[0, i].errorbar(
            sfs_data[i]["withoutAGN_mstar"],
            sfs_data[i]["withoutAGN_sfr"],
            xerr=sfs_data[i]["withoutAGN_mstar_err"],
            yerr=sfs_data[i]["withoutAGN_sfr_err"],
            fmt="none",
            ecolor="tab:blue",
            alpha = 0.3
        )
        ax[0, i].scatter(
            sfs_data[i]["withoutAGN_mstar"],
            sfs_data[i]["withoutAGN_sfr"],
            edgecolors="tab:blue",
            c="navy",
            alpha = 0.7,
            label="Pro Stellar"
        )

        # ax[0, i].plot([-100,100], [-100,100], ls="--", color="grey")
        ax[0, i].axvspan(
            0, 8, hatch="\\", facecolor="grey", edgecolor="black", alpha=0.3
        )


        ax[0, i].set_xlim([4.5, 12.5])
        ax[0, i].set_ylim([-4.5, 2.5])
        ax[0, i].set_xticks([6, 8, 10, 12])
        ax[0, i].set_yticks([-3.0, -1.5, 0.0, 1.5])
        ax[0,i].set_yticklabels([r'$10^{{{:n}}}$'.format(i) for i in ax[0,i].get_yticks()])
        ax[0, i].set_xticklabels([r'$10^{{{:n}}}$'.format(i) for i in ax[0, i].get_xticks()])

    ax[0, 0].plot(
        sfs_fits["mstar"], sfs_fits["z5_withAGN"],
        color="tab:red",
        linewidth = 2
    )
    ax[0, 0].fill_between(
        sfs_fits["mstar"], sfs_fits["z5_withAGN_q16"], sfs_fits["z5_withAGN_q84"],
        color = "tab:red",
        alpha = 0.5
    )
    ax[0, 1].plot(
        sfs_fits["mstar"], sfs_fits["z8_withAGN"],
        color="tab:red",
        linewidth=2
    )
    ax[0, 1].fill_between(
        sfs_fits["mstar"], sfs_fits["z8_withAGN_q16"], sfs_fits["z8_withAGN_q84"],
        color="tab:red",
        alpha=0.5
    )
    ax[0, 2].plot(
        sfs_fits["mstar"], sfs_fits["z10_withAGN"],
        color="tab:red",
        linewidth=2
    )
    ax[0, 2].fill_between(
        sfs_fits["mstar"], sfs_fits["z10_withAGN_q16"], sfs_fits["z10_withAGN_q84"],
        color="tab:red",
        alpha=0.5
    )

    ax[0, 0].plot(
        sfs_fits["mstar"], sfs_fits["z5_withoutAGN"],
        color="tab:blue",
        linewidth=2
    )
    ax[0, 0].fill_between(
        sfs_fits["mstar"], sfs_fits["z5_withoutAGN_q16"], sfs_fits["z5_withoutAGN_q84"],
        color="tab:blue",
        alpha=0.5
    )
    ax[0, 1].plot(
        sfs_fits["mstar"], sfs_fits["z8_withoutAGN"],
        color="tab:blue",
        linewidth=2
    )
    ax[0, 1].fill_between(
        sfs_fits["mstar"], sfs_fits["z8_withoutAGN_q16"], sfs_fits["z8_withoutAGN_q84"],
        color="tab:blue",
        alpha=0.5
    )
    ax[0, 2].plot(
        sfs_fits["mstar"], sfs_fits["z10_withoutAGN"],
        color="tab:blue",
        linewidth=2
    )
    ax[0, 2].fill_between(
        sfs_fits["mstar"], sfs_fits["z10_withoutAGN_q16"], sfs_fits["z10_withoutAGN_q84"],
        color="tab:blue",
        alpha=0.5
    )

    ax[1, 0].plot(
        smfs["mstar"],
        np.log10(
            smfs["z5_smf"] * pow(10, sfs_fits["z5_withAGN"])
        ),
        c = "tab:red",
        linewidth = 2
    )
    ax[1, 0].fill_between(
        smfs["mstar"],
        np.log10(
            smfs["z5_smf"] * pow(10, sfs_fits["z5_withAGN_q16"])
        ),
        np.log10(
            smfs["z5_smf"] * pow(10, sfs_fits["z5_withAGN_q84"])
        ),
        alpha = 0.5,
        color="tab:red",
        linewidth=2
    )
    ax[1, 1].plot(
        smfs["mstar"],
        np.log10(
            smfs["z8_smf"] * pow(10, sfs_fits["z8_withAGN"])
        ),
        c="tab:red",
        linewidth=2
    )
    ax[1, 1].fill_between(
        smfs["mstar"],
        np.log10(
            smfs["z8_smf"] * pow(10, sfs_fits["z8_withAGN_q16"])
        ),
        np.log10(
            smfs["z8_smf"] * pow(10, sfs_fits["z8_withAGN_q84"])
        ),
        alpha=0.5,
        color="tab:red",
        linewidth=2
    )
    ax[1, 2].plot(
        smfs["mstar"],
        np.log10(
            smfs["z10_smf"] * pow(10, sfs_fits["z10_withAGN"])
        ),
        c="tab:red",
        linewidth=2
    )
    ax[1, 2].fill_between(
        smfs["mstar"],
        np.log10(
            smfs["z10_smf"] * pow(10, sfs_fits["z10_withAGN_q16"])
        ),
        np.log10(
            smfs["z10_smf"] * pow(10, sfs_fits["z10_withAGN_q84"])
        ),
        alpha=0.5,
        color="tab:red",
        linewidth=2
    )

    ax[1, 0].plot(
        smfs["mstar"],
        np.log10(
            smfs["z5_smf"] * pow(10, sfs_fits["z5_withoutAGN"])
        ),
        c="tab:blue",
        linewidth=2
    )
    ax[1, 0].fill_between(
        smfs["mstar"],
        np.log10(
            smfs["z5_smf"] * pow(10, sfs_fits["z5_withoutAGN_q16"])
        ),
        np.log10(
            smfs["z5_smf"] * pow(10, sfs_fits["z5_withoutAGN_q84"])
        ),
        alpha=0.5,
        color="tab:blue",
        linewidth=2
    )
    ax[1, 1].plot(
        smfs["mstar"],
        np.log10(
            smfs["z7_smf"] * pow(10, sfs_fits["z8_withoutAGN"])
        ),
        c="tab:blue",
        linewidth=2
    )
    ax[1, 1].fill_between(
        smfs["mstar"],
        np.log10(
            smfs["z7_smf"] * pow(10, sfs_fits["z8_withoutAGN_q16"])
        ),
        np.log10(
            smfs["z7_smf"] * pow(10, sfs_fits["z8_withoutAGN_q84"])
        ),
        alpha=0.5,
        color="tab:blue",
        linewidth=2
    )
    ax[1, 2].plot(
        smfs["mstar"],
        np.log10(
            smfs["z10_smf"] * pow(10, sfs_fits["z10_withoutAGN"])
        ),
        c="tab:blue",
        linewidth=2
    )
    ax[1, 2].fill_between(
        smfs["mstar"],
        np.log10(
            smfs["z10_smf"] * pow(10, sfs_fits["z10_withoutAGN_q16"])
        ),
        np.log10(
            smfs["z10_smf"] * pow(10, sfs_fits["z10_withoutAGN_q84"])
        ),
        alpha=0.5,
        color="tab:blue",
        linewidth=2
    )

    ax[0, 2].legend()

    for i in range(3):
        ax[1, i].set_xlim([4.5, 12.5])
        ax[1, i].set_ylim([-6.0, -1.0])
        ax[1, i].set_xticks([6, 8, 10, 12])
        ax[1, i].set_yticks([-5.0, -3.0, -1.0])
        ax[1, i].set_yticklabels([r'$10^{{{:n}}}$'.format(i) for i in ax[1, i].get_yticks()])
        ax[1, i].set_xticklabels([r'$10^{{{:n}}}$'.format(i) for i in ax[1, i].get_xticks()])
        ax[1, i].xaxis.set_tick_params(pad=5)


    ax[0,0].set_ylabel(
        "$\\rm{ SFR \\, / \\, M_{\\odot} \\, yr^{-1} }$",
        fontsize = 16
    )

    ax[1, 0].set_ylabel(
        "$\\rm{ \\phi(M_{\\star}) \\times SFMS \\, / \\, M_{\\odot} \\, yr^{-1} \\, Mpc^{-3} \\, dex^{-1} }$",
        fontsize = 16
    )

    fig.supxlabel(
        "$\\rm{ M_{\\star} \\, / \\, M_{\\odot} }$"
    )

    fig.savefig(
        main_stub + "plots/sfs.pdf"
    )

def plot_delta_sfr(data):

    withAGN = data["withAGN"]
    withoutAGN = data["withoutAGN"]

    norm = plt.Normalize()
    cm = plt.colormaps["RdYlBu_r"]

    fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(12, 5), constrained_layout=True, sharey=True)

    sfs_data = data["sfs"]
    mstar_bins = np.array([4.5, 6.5, 8.5, 10.5])

    sfr_withAGN_err = (np.array(withAGN["SFRBurst16"]) + np.array(withAGN["SFRBurst84"])) / (
                np.log(10) * np.array(withAGN["SFRBurst50"]))
    sfr_withoutAGN_err = (np.array(withoutAGN["SFRBurst16"]) + np.array(withoutAGN["SFRBurst84"])) / (
                np.log(10) * np.array(withoutAGN["SFRBurst50"]))

    for i in range(3):
        idx = np.where(
            (np.log10(withAGN["Mstar50"]) >= mstar_bins[i]) &
            (np.log10(withAGN["Mstar50"]) < mstar_bins[i+1])
        )

        redshift = np.array(withAGN["z"])[idx]

        ax[i].set_title(
            str(mstar_bins[i]) + "$ \\rm {\\leq \\log_{10} (M_{\\star}^{AGN+Stellar}/M_{\\odot})< }$" + str(mstar_bins[i+1]),
            fontsize = 14
        )
        axx=ax[i].scatter(
            np.log10(np.array(withAGN["AGNlum50"]))[idx],
            np.log10(np.array(withAGN["SFRBurst50"]))[idx] - np.log10(np.array(withoutAGN["SFRBurst50"]))[idx],
            c=redshift,
            vmin=min(redshift),
            vmax=max(redshift),
            cmap=cm,
            edgecolor = "black",
            zorder = 200
        )

        ax[i].errorbar(
            np.log10(np.array(withAGN["AGNlum50"]))[idx],
            np.log10(np.array(withAGN["SFRBurst50"]))[idx] - np.log10(np.array(withoutAGN["SFRBurst50"]))[idx],
            yerr = np.sqrt(
                pow(sfr_withAGN_err[idx],2) + pow(sfr_withoutAGN_err[idx],2)
            ),
            xerr = ( np.array(withAGN["AGNlum16"])[idx] + np.array(withAGN["AGNlum84"])[idx] ) / (np.log(10)*np.array(withAGN["AGNlum50"])[idx]),
            fmt = "none",
            color = "black",
            alpha = 0.5
        )

        ax[i].set_xlim([34.5, 47.5])
        ax[i].set_xticks([35, 39, 43, 47])
        ax[i].xaxis.set_tick_params(pad=5)

        ax[i].set_xticklabels([r'$10^{{{:n}}}$'.format(i) for i in ax[i].get_xticks()])

        ax[i].set_ylim(
            -6.0, 0.5
        )
        ax[i].set_yticks([-4, -2, 0])
        ax[i].axhline(0, color="black", linestyle="--")


    ax[1].set_ylim(
        -6, 1
    )
    ax[2].set_ylim(
        -6, 1
    )

    # cbar = fig.colorbar(axx, ax=ax[2], orientation='vertical', ticks=[3.5, 6.5, 9.5, 12.5], label="Redshift")
    # cbar.ax.tick_params(labelsize=12, pad=5)
    fig.colorbar(axx, orientation='vertical', label = "Redshift")

    fig.supxlabel(
        "$\\rm{ L_{AGN} \\, / \\, erg \\, s^{-1} }$",
        fontsize = 24
    )

    fig.supylabel(
        "$\\rm{ \\Delta SFR_{AGN+Stellar - Stellar} }$",
        fontsize = 24
    )

    # plt.tight_layout()

    # fig.subplots_adjust(right=0.75)
    # cbar_ax = fig.add_axes([0.83, 0.15, 0.05, 0.7])
    # cbar = fig.colorbar(axx, cax=cbar_ax, orientation='vertical',
    #                     ticks=[3.5, 6.5, 9.5, 12.5]
    #                     )
    # cbar.set_label(label="Redshift", fontsize = 24)

    fig.savefig(
        main_stub + "plots/delta_sfr.pdf"
    )

def plot_sed(data):

    sed_data = data["sed"]
    transmission_curves = sed_data["transmission_curves"]

    phot = sed_data["phot"]
    sed_withAGN = sed_data["sed_withAGN"]
    sed_withoutAGN = sed_data["sed_withoutAGN"]

    stars_withAGN = sed_data["stars_withAGN"]
    stars_withoutAGN = sed_data["stars_withoutAGN"]

    agn_withAGN = sed_data["agn_withAGN"]

    fig, ax = plt.subplots(nrows = 2, ncols=1, figsize=(8,6), constrained_layout=True, height_ratios=(3,1), sharex=True)
    ax[0].set_title("CEERS_1019, z=8.679")

    cmap = cm.get_cmap(name='rainbow')
    print(transmission_curves.keys())
    len_filt = len(transmission_curves.keys()) - 1
    for i in range(len_filt):
        temp = transmission_curves.iloc[:, i + 1].to_numpy()
        idx = (temp > 0.01)
        transmission = temp[idx] / np.nanmax(temp[idx])
        ax[0].plot(
            transmission_curves["lambda"][idx] / 1e4,
            (transmission) * 0.2 - 0.21,
            color=cmap(i / len_filt),
            alpha=0.5,
            ls = "-"
        )

    ax[0].errorbar(
        phot["lambda"]/1e4,
        phot["flux"]/1e-19,
        phot["flux_err"]/1e-19,
        fmt = "o",
        color = "black",
        capsize = 3,
        uplims = np.array([np.nan, np.nan,1,np.nan, 1,np.nan, np.nan, 0,0,0,0,0,0,0,0,0,0])
    )
    ax[0].plot(
        sed_withAGN["wave"]/1e4,
        sed_withAGN["flux"]/1e-19,
        linewidth = 3,
        color = "tab:red",
        linestyle = "-",
        alpha = 0.7,
        label = "Pro Stellar+AGN (LP = -6.28)"
    )
    ax[0].plot(
        sed_withoutAGN["wave"] / 1e4,
        sed_withoutAGN["flux"] / 1e-19,
        linewidth=3,
        color="tab:blue",
        linestyle="-",
        alpha = 0.7,
        label="Pro Stellar (LP = -7.58)"
    )
    ax[0].plot(
        agn_withAGN["wave"]/1e4,
        agn_withAGN["flux"]/1e-19,
        linewidth=3,
        color = "tab:purple",
        alpha = 0.7,
        label = "AGN flux"
    )

    # ax[0].text(
    #     -0.3, 0.65,
    #     "Larson+23 $\\rm{L_{AGN}}$: \n $\\rm{45.10 \\pm 0.20}$",
    #     fontsize = 16
    # )
    # ax[0].text(
    #     -0.3, 0.5,
    #     "ProSpect $\\rm{L_{AGN}}$: \n $\\rm{45.02 \\pm 2.35}$",
    #     fontsize = 16
    # )

    ax[0].legend(fontsize = 18)

    ax[0].set_xlim([0.0,5.1])
    ax[0].set_ylim([-0.2,0.9])
    ax[0].set_yticks([0,0.4,0.8])
    # ax[0].set_yscale("log")

    ax[0].set_ylabel("$\\rm{F_{\\lambda} \\, / \\, 10^{-19} erg \\, s^{-1} \\, cm^{-2} \\, \\AA^{-1}}$")


    ax[1].scatter(
        phot["lambda"]/1e4,
        (phot["flux"] - phot["sed_withAGN"]) / phot["flux_err"],
        color = "tab:red",
        alpha = 0.7
    )
    ax[1].scatter(
        phot["lambda"] / 1e4,
        (phot["flux"] - phot["sed_withoutAGN"]) / phot["flux_err"],
        color="tab:blue",
        alpha = 0.7
    )
    ax[1].axhline(0, ls = "--", color="grey")
    ax[1].set_ylim([-10,10])
    ax[1].set_yticks([-5,0,5])
    ax[1].set_ylabel("$\\rm{\\frac{Data-Model}{Err}}$")

    fig.supxlabel("Wavelength/$\\rm{\\mu m}$")

    fig.savefig(main_stub + "/plots/SED_galaxy.pdf")

def main():
    print("Main script")

    data = load_data()
    # plot_csfh(data)
    # plot_mstar(data)
    # plot_sfms(data)
    # plot_delta_sfr(data)
    plot_sed(data)


if __name__ == "__main__":
    main()