def plot_everything(data):
    """Just plot the CSFH"""


    fig, ax = plt.subplots(nrows=3, ncols=1, figsize = (7, 8), constrained_layout = True)

    ## plot the SFR-Mstar
    norm = plt.Normalize()
    cm = plt.colormaps["RdYlBu_r"]

    ax0=ax[0].scatter(
        np.log10(data["catalogue"]["Mstar"]),
        np.log10(data["catalogue"]["SFRburst"]),
        edgecolors = "black",
        s = 10 * (1.1*np.log10(data["catalogue"]["AGNlum"]) - np.min(np.log10(data["catalogue"]["AGNlum"]))),
        c = data["catalogue"]["z"],
        vmin = min(data["catalogue"]["z"]),
        vmax = max(data["catalogue"]["z"]),
        cmap = cm
    )

    inset_axes_font_size = 10.0
    ax0_inset = inset_axes(ax[0], width="30%", height="5%",
                        loc = "center", bbox_to_anchor=(-0.05, -0.314, 1, 1),
                        bbox_transform=ax[0].transAxes
                        )
    cbar=fig.colorbar(ax0, cax=ax0_inset, orientation='horizontal', ticks = [3.50,5.75,8.00,10.25,12.50])
    cbar.ax.tick_params(labelsize=inset_axes_font_size, pad = 5)
    ax0_inset.text(9.15, -0.2, "Redshift", transform = ax[0].transData, fontsize = inset_axes_font_size)

    ax0_inset1 = inset_axes(ax[0], width="30%", height="30%",
                            loc="center", bbox_to_anchor=(0.30, -0.20, 1.0, 1.0),
                            bbox_transform=ax[0].transAxes
                            )
    ax0_inset1.text(9.75, -0.2, "$\\rm{log_{10}(L_{AGN} \\, / \\, erg \\, s^{-1})}$", transform=ax[0].transData, fontsize = inset_axes_font_size)

    sizes_vector = 10 * (1.1*np.linspace(35,48,5) - np.min(np.linspace(35,48,5)))
    ax0_inset1.patch.set_alpha(0.0)
    ax0_inset1.scatter(sizes_vector, -np.zeros(5), s = sizes_vector,
                       facecolor="grey", edgecolor="black")
    ax0_inset1.set_ylim([-0.02, 0.1])
    ax0_inset1.set_xlim([min(sizes_vector)*0.8, max(sizes_vector)*1.1])
    for spine in ax0_inset1.spines.values():
        spine.set_visible(False)
    ax0_inset1.set_yticks([])
    ax0_inset1.set_xticks(sizes_vector)
    ax0_inset1.set_xticklabels(np.linspace(35,48,5), fontsize = inset_axes_font_size)

    ax[0].set_ylim([-1.0, 2.0])
    ax[0].set_yticks([-0.5, 0.5, 1.5])
    ax[0].set_yticklabels([r'$10^{{{:n}}}$'.format(i) for i in [-0.5, 0.5, 1.5]])

    ax[0].set_xticklabels([r'$10^{{{:n}}}$'.format(i) for i in ax[0].get_xticks()])

    ax[0].set_xlabel("$\\rm{M_{\\star} \\, / \\, M_{\\odot}}$")
    ax[0].set_ylabel("$\\rm{SFR \\, / \\, M_{\\odot} \\, yr^{-1}}$")

    colour_palette = ["#001219", "#ffbd00", "#0A9396", "#9B2226", "#E9D8A6", "#005F73"]
    ## now FINALLY do the CSFH plot
    ax[1].set_ylim(-3.5, -0.5)
    ax[1].set_xlim(0.0, 12.5)
    ax[1].set_yticklabels([r'$10^{{{:n}}}$'.format(i) for i in ax[1].get_yticks()])

    zvec = np.linspace(0, 20)
    def harikane_csfh(z):
        return (
            np.log10(1 / (61.7 * (1 + z) ** -3.13 + 1.0 * 10 ** (0.22 * (1 + z)) + 22.4 * 10 ** (0.5 * (1 + z) - 3.0)))
        )

    def md2014_fit(z):
        T0 = 0.015 * (1 + z) ** 2.7
        T1 = 1 + ((1 + z) / 2.9) ** 5.6
        return np.log10( T0 / T1 )

    ax[1].plot(zvec, harikane_csfh(zvec), color=colour_palette[5], label = "Harikane+22")
    ax[1].plot(zvec, md2014_fit(zvec), color=colour_palette[5], linestyle = "--", label = "M&D+14")

    ax[1].errorbar(
        x=data["gama_devils"]["z"],
        y=data["gama_devils"]["csfh"],
        yerr=[data["gama_devils"]["csfh_err_down"], data["gama_devils"]["csfh_err_up"]],
        color=colour_palette[4],
        fmt=".",
        markersize=8,
        label = "D'Silva+23"
    )

    ax[1].errorbar(
        x=data["data"]["z"],
        y=data["data"]["csfh_AGN"],
        yerr=[data["data"]["csfh_AGN_err_down"], data["data"]["csfh_AGN_err_up"]],
        color=colour_palette[3],
        fmt="v",
        alpha=0.5,
        label = "SFR [Pro AGN+Stellar]"
    )
    ax[1].errorbar(
        x=data["data"]["z"],
        y=data["data"]["csfh_noAGN"],
        yerr=[data["data"]["csfh_noAGN_err_down"], data["data"]["csfh_noAGN_err_up"]],
        color=colour_palette[2],
        fmt="^",
        alpha=0.5,
        label = "SFR [Pro Stellar]"
    )

    ax[1].errorbar(
        x=data["data"]["z"],
        y=data["data"]["csfh_replace"],
        xerr=[data["data"]["zmin"], data["data"]["zmax"]],
        yerr=0,
        color=colour_palette[0],
        fmt="none",
        capsize = 3,
        alpha = 0.3
    )
    ax[1].errorbar(
        x=data["data"]["z"],
        y=data["data"]["csfh_replace"],
        xerr=[data["data"]["z16"], data["data"]["z84"]],
        yerr=[data["data"]["csfh_replace_err_down"], data["data"]["csfh_replace_err_up"]],
        markerfacecolor=colour_palette[1],
        fmt="o",
        capsize = 3,
        label = "SFR [Hybrid]",
        alpha = 1.0,
        ecolor = colour_palette[0],
        markeredgecolor = colour_palette[0]
    )
    ax[1].legend(fontsize = 10, ncol = 1, frameon = False, loc = "lower left", bbox_to_anchor = (0.02, 0.0, 1.0, 1.0))
    ax[1].set_xlabel("Redshift")
    ax[1].set_ylabel("$\\rm{CSFH\\,/\\,M_{\\odot}\\,yr^{-1}\\,Mpc^{-3}}$")

    ## now do the CAGNH plot
    ax[2].set_ylim(40.0, 42.0)
    ax[2].set_xlim(0.0, 12.5)
    ax[2].set_yticklabels([r'$10^{{{:n}}}$'.format(i) for i in ax[2].get_yticks()])


    ax[2].errorbar(
        x=data["gama_devils"]["z"],
        y=data["gama_devils"]["cagnh"],
        yerr=[data["gama_devils"]["cagnh_err_down"], data["gama_devils"]["cagnh_err_up"]],
        color=colour_palette[4],
        fmt=".",
        markersize=8,
        label = "D'Silva+23"
    )

    ax[2].errorbar(
        x=data["data"]["z"],
        y=data["data"]["cagnh_shen"],
        yerr=[data["data"]["cagnh_shen_err_down"], data["data"]["cagnh_shen_err_up"]],
        color=colour_palette[3],
        fmt="v",
        alpha=0.5,
        label = "AGN [Pro]"
    )

    ax[2].errorbar(
        x=data["data"]["z"],
        y=data["data"]["cagnh_shen_replace"],
        xerr=[data["data"]["zmin"], data["data"]["zmax"]],
        yerr=0,
        color=colour_palette[0],
        fmt="none",
        capsize = 3,
        alpha = 0.3
    )
    ax[2].errorbar(
        x=data["data"]["z"],
        y=data["data"]["cagnh_shen_replace"],
        xerr=[data["data"]["z16"], data["data"]["z84"]],
        yerr=[data["data"]["cagnh_shen_replace_err_down"], data["data"]["cagnh_shen_replace_err_up"]],
        markerfacecolor=colour_palette[1],
        fmt="o",
        capsize = 3,
        label = "AGN [Pro LB]",
        alpha=1.0,
        ecolor=colour_palette[0],
        markeredgecolor=colour_palette[0]
    )

    ax[2].legend(fontsize = 10, ncol = 1, frameon = False, loc="lower left")
    ax[2].set_xlabel("Redshift")
    ax[2].set_ylabel("$\\rm{CAGNH \\, / \\, erg \\, s^{-1} \\, Mpc^{-3}}$")
    fig.savefig("./plots/master_plot_csfh_cagnh.pdf")

def plot_everythingV2(data):
    """Just plot only the relevant stuff. No replacements"""


    fig, ax = plt.subplots(nrows=3, ncols=1, figsize = (7, 8), constrained_layout = True)

    ## plot the SFR-Mstar
    norm = plt.Normalize()
    cm = plt.colormaps["RdYlBu_r"]

    ax0=ax[0].scatter(
        np.log10(data["catalogue"]["Mstar"]),
        np.log10(data["catalogue"]["SFRburst"]),
        edgecolors = "black",
        # s = 10 * (1.1*np.log10(data["catalogue"]["AGNlum"]) - np.min(np.log10(data["catalogue"]["AGNlum"]))),
        s=10 * (1.1 * np.log10(data["catalogue"]["AGNlum"]) - 35),
        c = data["catalogue"]["z"],
        vmin = min(data["catalogue"]["z"]),
        vmax = max(data["catalogue"]["z"]),
        cmap = cm
    )

    inset_axes_font_size = 12.0
    ax0_inset = inset_axes(ax[0], width="30%", height="5%",
                        loc = "center", bbox_to_anchor=(-0.07, -0.314, 1, 1),
                        bbox_transform=ax[0].transAxes
                        )
    cbar=fig.colorbar(ax0, cax=ax0_inset, orientation='horizontal', ticks = [4,8,12])
    cbar.ax.tick_params(labelsize=inset_axes_font_size, pad = 5)
    ax0_inset.text(9.2, -0.2, "Redshift", transform = ax[0].transData, fontsize = inset_axes_font_size)

    ax0_inset1 = inset_axes(ax[0], width="30%", height="30%",
                            loc="center", bbox_to_anchor=(0.30, -0.20, 1.0, 1.0),
                            bbox_transform=ax[0].transAxes
                            )
    ax0_inset1.text(9.85, -0.2, "$\\rm{log_{10}(L_{AGN} \\, / \\, erg \\, s^{-1})}$", transform=ax[0].transData, fontsize = inset_axes_font_size)

    sizes_vector = 10 * (1.1*np.linspace(35,48,3) - np.min(np.linspace(35,48,3)))
    ax0_inset1.patch.set_alpha(0.0)
    ax0_inset1.scatter(sizes_vector, np.zeros(3), s = sizes_vector,
                       facecolor="grey", edgecolor="black")
    ax0_inset1.set_ylim([-0.02, 0.1])
    ax0_inset1.set_xlim([min(sizes_vector)*0.8, max(sizes_vector)*1.1])
    for spine in ax0_inset1.spines.values():
        spine.set_visible(False)
    ax0_inset1.set_yticks([])
    ax0_inset1.set_xticks(sizes_vector)
    ax0_inset1.set_xticklabels(np.round(np.linspace(35,48,3),2), fontsize = inset_axes_font_size)

    ax[0].set_ylim([-1.0, 2.0])
    ax[0].set_yticks([-0.5, 0.5, 1.5])
    ax[0].set_yticklabels([r'$10^{{{:n}}}$'.format(i) for i in [-0.5, 0.5, 1.5]])

    ax[0].set_xlim([8.25,10.75])
    ax[0].set_xticklabels([r'$10^{{{:n}}}$'.format(i) for i in ax[0].get_xticks()])

    ax[0].set_xlabel("$\\rm{M_{\\star} \\, / \\, M_{\\odot}}$")
    ax[0].set_ylabel("$\\rm{SFR \\, / \\, M_{\\odot} \\, yr^{-1}}$")

    colour_palette = {
        "main" : "#4A4E69",
        "lines" : "#C9ADA7",
        "points" : "#9A8C98",
        "errors" : "#22223B"
    }
    ## now FINALLY do the CSFH plot
    ax[1].set_ylim(-4.0, -0.5)
    ax[1].set_xlim(0.0, 12.5)
    ax[1].set_yticklabels([r'$10^{{{:n}}}$'.format(i) for i in ax[1].get_yticks()])

    zvec = np.linspace(0, 20)
    def harikane_csfh(z):
        return (
            np.log10(1 / (61.7 * (1 + z) ** -3.13 + 1.0 * 10 ** (0.22 * (1 + z)) + 22.4 * 10 ** (0.5 * (1 + z) - 3.0))) + np.log10(0.63)
        )

    def md2014_fit(z):
        T0 = 0.015 * (1 + z) ** 2.7
        T1 = 1 + ((1 + z) / 2.9) ** 5.6
        return np.log10(T0 / T1 ) + np.log10(0.63) #0.63 convert to Chabrier

    ax[1].plot(zvec, harikane_csfh(zvec), color=colour_palette["lines"], linewidth = 2, label = "Harikane+22")
    ax[1].plot(zvec, md2014_fit(zvec), color=colour_palette["lines"], linestyle = "--", linewidth = 2, label = "M&D+14")

    ax[1].errorbar(
        x=data["gama_devils"]["z"],
        y=data["gama_devils"]["csfh"],
        yerr=[data["gama_devils"]["csfh_err_down"], data["gama_devils"]["csfh_err_up"]],
        color=colour_palette["points"],
        fmt=".",
        markersize=8,
        label = "D'Silva+23"
    )

    ax[1].errorbar(
        x=data["data"]["z"],
        y=data["data"]["csfh_replace"],
        xerr=[data["data"]["zmin"], data["data"]["zmax"]],
        yerr=0,
        color="grey",
        fmt="none",
        capsize = 3,
        alpha = 0.3
    )
    ax[1].errorbar(
        x=data["data"]["z"],
        y=data["data"]["csfh_replace"],
        xerr=[data["data"]["z16"], data["data"]["z84"]],
        yerr=[data["data"]["csfh_replace_err_down"], data["data"]["csfh_replace_err_up"]],
        markerfacecolor=colour_palette["main"],
        fmt="s",
        capsize = 3,
        label = "SFR",
        alpha = 1.0,
        ecolor = colour_palette["errors"],
        markeredgecolor = colour_palette["errors"]
    )
    ax[1].legend(ncol = 1, frameon = False, loc = "lower left", bbox_to_anchor = (0.02, 0.0, 1.0, 1.0))
    ax[1].set_xlabel("Redshift")
    ax[1].set_ylabel("$\\rm{CSFH\\,/\\,M_{\\odot}\\,yr^{-1}\\,Mpc^{-3}}$")

    ## now do the CAGNH plot
    ax[2].set_ylim(39.0, 42.5)
    ax[2].set_xlim(0.0, 12.5)
    ax[2].set_yticklabels([r'$10^{{{:n}}}$'.format(i) for i in ax[2].get_yticks()])


    ax[2].errorbar(
        x=data["gama_devils"]["z"],
        y=data["gama_devils"]["cagnh"],
        yerr=[data["gama_devils"]["cagnh_err_down"], data["gama_devils"]["cagnh_err_up"]],
        color=colour_palette["points"],
        fmt=".",
        markersize=8,
        label = "D'Silva+23"
    )

    ax[2].errorbar(
        x=data["data"]["z"],
        y=data["data"]["cagnh_shen_replace"],
        xerr=[data["data"]["zmin"], data["data"]["zmax"]],
        yerr=0,
        color="grey",
        fmt="none",
        capsize = 3,
        alpha = 0.3
    )
    ax[2].errorbar(
        x=data["data"]["z"],
        y=data["data"]["cagnh_shen_replace"],
        xerr=[data["data"]["z16"], data["data"]["z84"]],
        yerr=[data["data"]["cagnh_shen_replace_err_down"], data["data"]["cagnh_shen_replace_err_up"]],
        markerfacecolor=colour_palette["main"],
        fmt="s",
        capsize = 3,
        label = "AGN",
        alpha=1.0,
        ecolor=colour_palette["errors"],
        markeredgecolor=colour_palette["errors"]
    )

    ax[2].legend(ncol = 1, frameon = False, loc="lower left")
    ax[2].set_xlabel("Redshift")
    ax[2].set_ylabel("$\\rm{CAGNH \\, / \\, erg \\, s^{-1} \\, Mpc^{-3}}$")
    fig.savefig("./plots/master_plot_csfh_cagnh.pdf")