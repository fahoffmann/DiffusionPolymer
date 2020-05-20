#!/usr/bin/env python
# coding: utf-8

import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
import statsmodels.api as sm
from scipy import stats
from my_math import stretched_exponential
from my_geometry import total_number_of_segments
from my_MonteCarlo import monte_carlo
from my_constants import Constants
from my_plotting import prepare_subplots, plot_fit_err, plot_fit, plot_ax
import sys

plt.rc('text', usetex=False)

def initialization():
    """ Loads the constants, calculates the volume and number of walkers"""
    constants = Constants()
    constants.volume = total_number_of_segments(constants.dimension, constants.radius)
    constants.number_of_walkers = int(constants.concentration * constants.volume)
    b_tau = []
    radii = [10] #[5,6,7, 8]
    diffusion_constants = [0.02, 0.05, 0.1, 0.5, 1]
    return constants, b_tau, radii, diffusion_constants

def convert(mc_time_all, constants):
    """ Post-Processing of the simulations.
        Needs DataFrame mc_time_all and constants
        Calculates fitting parameter b and tau."""
    M_t_all = pd.DataFrame([constants.number_of_walkers - i for i in range(len(mc_time_all))[::-1]])
    M_t_all = M_t_all.divide(M_t_all.iloc[-1]).multiply(constants.concentration)
    M_t_norm = M_t_all #1. - M_t_all
    mc_time_avg = mc_time_all.mean(axis=1).to_frame()
    mc_time_std = mc_time_all.std(axis=1).to_frame()
    M_t_mc_time = pd.concat([M_t_norm, mc_time_avg], axis=1)
    M_t_mc_time.columns = ['M_t', 't']
    print(M_t_mc_time, len(M_t_norm), len(mc_time_avg), constants.number_of_walkers)
    popt, popv = curve_fit(stretched_exponential, M_t_mc_time['t'], M_t_mc_time['M_t'], p0=(0.6, 150*constants.radius ** 2 / constants.D_r))
    b_tau.append(popt)
    M_t_mc_time = pd.concat([M_t_norm, mc_time_avg, stretched_exponential(mc_time_avg, *popt), mc_time_std], axis=1)
    M_t_mc_time.columns = ['M_t', 't', 'exp', 'std']
    return b_tau, M_t_mc_time

def convert_MC(mc_time_all, constants):
    """ Post-Processing of the simulations.
            Needs DataFrame mc_time_all and constants
            Calculates fitting parameter b and tau."""
    M_t_all = pd.DataFrame([constants.number_of_walkers - i for i in range(len(mc_time_all))[::-1]])
    M_t_all = M_t_all.divide(M_t_all.iloc[-1]).multiply(constants.concentration)
    M_t_norm = M_t_all #1. - M_t_all
    mc_time_avg = mc_time_all.mean(axis=1).to_frame()
    mc_time_std = mc_time_all.std(axis=1).to_frame()
    M_t_mc_time = pd.concat([M_t_norm, mc_time_avg], axis=1)
    M_t_mc_time.columns = ['M_t', 't']
    print(M_t_mc_time, len(M_t_norm), len(mc_time_avg), constants.number_of_walkers)
    popt, popv = curve_fit(stretched_exponential, M_t_mc_time['t'], M_t_mc_time['M_t'], p0=(0.6, 150*constants.radius ** 2 / constants.D_r))
    b_tau.append(popt)
    M_t_mc_time = pd.concat([M_t_norm, mc_time_avg, stretched_exponential(mc_time_avg, *popt), mc_time_std], axis=1)
    M_t_mc_time.columns = ['M_t', 't', 'exp', 'std']
    return b_tau, M_t_mc_time

def postprocess(b_tau, values, description):
    """Post-Processing of the simulations.
       Needs b, tau and values.
       Produces DataFrame with values."""
    b_df = pd.DataFrame([b_tau[i][0] for i in range(len(b_tau))], [i for i in range(1,len(values) + 1)])
    tau_df = pd.DataFrame([b_tau[i][1] for i in range(len(b_tau))], [i for i in range(1, len(values) + 1)])
    values_df = pd.DataFrame([values[i] for i in range(len(values))], [i for i in range(1, len(values) + 1)])
    b_tau_df = pd.concat([b_df, tau_df, values_df], axis=1)
    b_tau_df.columns = ['b', 'tau', description]
    b_tau_df = b_tau_df.set_index(description)
    print(b_tau_df)
    return b_tau_df

if __name__ == "__main__":
    calc_radius=True
    calc_diffusion=False
    constants, b_tau, radii, diffusion_constants = initialization()
    if calc_radius:
        xmax = 0
        ystdmax = 0
        plt, ax1, ax2 = prepare_subplots(plt)
        for radius in radii:
            print(radius)
            constants.radius = radius
            if constants.diffusion:
                mc_time_all, mc_time_pa_all, mc_time_pe_all, mc_walker_all = monte_carlo(constants)
                column_headers = ["t_"+str(i) for i in range(constants.number_of_iterations)]+["l_"+str(i) for i in
                                                                                               range(constants.number_of_iterations)]
                mc_time_all.columns = column_headers
                mc_time_pa_all.columns = column_headers
                mc_time_pe_all.columns = column_headers
                method="parallel_n1d2a5time_walker"
                mc_time_all.to_csv("diffusion_"+method+"_newmolecule.csv")
                mc_time_pa_all.to_csv("diffusion_"+method+"_parallel.csv")
                mc_time_pe_all.to_csv("diffusion_"+method+"_perpendicular.csv")
                column_headers = ["d_" + str(i) for i in range(constants.number_of_iterations)] + ["t_" + str(i) for i
                                                                                                   in range(
                        constants.number_of_iterations)]
                mc_walker_all.columns = column_headers
                mc_walker_all.to_csv("diffusion_"+method+"_walker.csv")
                for i in range(constants.number_of_iterations):
                    ax1.scatter(mc_time_all["t_"+str(i)], mc_time_all["l_"+str(i)],marker='o')
                plt.savefig("diffusion_"+method+"_buildup.png")
                #sys.exit()
                #b_tau, M_t_mc_time = convert_MC(mc_time_all, constants)
            else:
                mc_time_all = monte_carlo(constants)
                b_tau, M_t_mc_time = convert(mc_time_all, constants)
                M_t_mc_time.to_csv("allsides_R"+str(radius)+"c05_iteration1.csv")
                plot_fit_err(M_t_mc_time['t'], M_t_mc_time['M_t'], M_t_mc_time['exp'], M_t_mc_time['std'], ax1, ax2, "R=" + str(constants.radius))
                xmax = max(xmax, max(M_t_mc_time['t']))
                ystdmax = max(ystdmax, max(M_t_mc_time['std']))
                ax1 = plot_ax(ax1,0.,xmax,0.,1.1,"MC time","M_t/M_inf")
                ax2 = plot_ax(ax2,0.,xmax,0.,1.1*ystdmax,"MC time","Standard deviation")
                #plt.show()
                plt.savefig("buildup_radii.png")
        plt.close()
        sys.exit()
        if not(constants.diffusion):
            b_tau_df = postprocess(b_tau, radii, 'radius')
            plt, ax1, ax2 = prepare_subplots(plt)
            model = sm.OLS(b_tau_df['tau'], b_tau_df.index ** 2)
            fit_param = model.fit()
            plot_fit(b_tau_df.index ** 2,b_tau_df['tau'],fit_param.params[0] * b_tau_df.index ** 2,ax1,"tau = " + str(round(fit_param.params[0], 2)) + " R^2", None)
            ax1 = plot_ax(ax1, 0., 1.1 * max(b_tau_df.index) ** 2, 0., 1.1 * max(b_tau_df['tau']),"R^2","tau")
            slope, intercept, r_value, p_value, std_err = stats.linregress(1.0 / b_tau_df.index, b_tau_df['b'])
            plot_fit(1.0 / b_tau_df.index, b_tau_df['b'], slope / b_tau_df.index + intercept, ax2, None,"b = " + str(round(slope, 2)) + "/R + " + str(round(intercept, 2)))
            ax2 = plot_ax(ax2, 0., 1.1 * max(1.0 / b_tau_df.index), 0.9 * min(b_tau_df['b']), 1.1 * max(b_tau_df['b']), "1/R", "b")
            tau_av = gamma(1.0 / b_tau_df['b']) / b_tau_df['b'] * b_tau_df['tau']
            model = sm.OLS(tau_av, b_tau_df.index ** 2)
            results2 = model.fit()
            ax1.plot(b_tau_df.index ** 2, tau_av, label="tau_av = " + str(round(results2.params[0], 2)) + " R^2")
            ax1.legend()
            plt.savefig("parameters_radii.png")
            plt.close()
    if calc_diffusion:
        constants.concentration = 0.5
        constants.radius = 7
        xmax = 0
        ystdmax = 0
        plt, ax1, ax2 = prepare_subplots(plt)
        for D_r in diffusion_constants:
            print(D_r)
            constants.D_r = D_r
            mc_time_all = monte_carlo(constants)
            b_tau, M_t_mc_time = convert(mc_time_all, constants)
            plot_fit_err(M_t_mc_time['t'], M_t_mc_time['M_t'], M_t_mc_time['exp'], M_t_mc_time['std'], ax1, ax2, "D_r=" + str(constants.D_r))
            xmax = max(xmax, max(M_t_mc_time['t']))
            ystdmax = max(ystdmax, max(M_t_mc_time['std']))
        ax1 = plot_ax(ax1, 0., xmax, 0., 1.1, "MC time", "M_t/M_inf")
        ax2 = plot_ax(ax2, 0., xmax, 0., 1.1 * ystdmax, "MC time", "Standard deviation")
        plt.savefig("buildup_diffusion.png")
        plt.close()
        b_tau_df = postprocess(b_tau, diffusion_constants, 'diffusion')
        plt, ax1, ax2 = prepare_subplots(plt)
        model = sm.OLS(b_tau_df['tau'], constants.radius ** 2 / b_tau_df.index)
        fit_param = model.fit()
        plot_fit(1.0 / b_tau_df.index, b_tau_df['tau'], fit_param.params[0] * constants.radius ** 2 / b_tau_df.index, ax1,
                 "tau = " + str(round(fit_param.params[0], 2)) + " R^2/D_r", None)
        ax1 = plot_ax(ax1, 0., 1.1 * max(1.0 / b_tau_df.index), 0., 1.1 * max(b_tau_df['tau']), "1/D_r", "tau")
        slope, intercept, r_value, p_value, std_err = stats.linregress(b_tau_df.index, b_tau_df['b'])
        plot_fit(b_tau_df.index, b_tau_df['b'], slope * b_tau_df.index + intercept, ax2, None,
                 "b = " + str(round(slope, 2)) + "*R + " + str(round(intercept, 2)))
        ax2 = plot_ax(ax2, 0., 1.1 * max(b_tau_df.index), 0.9 * min(b_tau_df['b']), 1.1 * max(b_tau_df['b']),
                      "D_r", "b")
        plt.savefig("parameters_diffusion.png")
        plt.close()
        plt, ax1, ax2 = prepare_subplots(plt)
        plot_fit(M_t_mc_time['t'], M_t_mc_time['M_t'], M_t_mc_time['exp'], ax1, "R=" + str(constants.radius), None)
        ax1 = plot_ax(ax1, 0., max(M_t_mc_time['t']), 0., max(M_t_mc_time['exp']),"MC time", "xxx")
        plt.show()
        plt.close()