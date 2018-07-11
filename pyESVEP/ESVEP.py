# This file is part of pyESVEP for running ESVEP model
# Copyright 2018 Radoslaw Guzinski and contributors listed in the README.md file.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

'''
Created on May 30 2018
@author: Radoslaw Guzinski (rmgu@dhigroup.com)

DESCRIPTION
===========
This package contains the main routine inherent of End-member-based Soil and Vegetation Energy
Patrtitioning `ESVEP` model.
Additional functions needed in ESVEP, such as computing of net radiation or estimating the
resistances to heat and momentum transport are imported from pyTSEB modules.

* :doc:`pyTSEB.TSEB` for ancilary model functions (e.g. estimation of G).
* :doc:`pyTSEB.net_radiation` for the estimation of net radiation and radiation partitioning.
* :doc:`pyTSEB.clumping_index` for the estimatio of canopy clumping index.
* :doc:`pyTSEB.meteo_ytils` for the estimation of meteorological variables.
* :doc:`pyTSEB.resistances` for the estimation of the resistances to heat and momentum transport.
* :doc:`pyTSEB.MO_similarity` for the estimation of the Monin-Obukhov length and MOST-related variables.

PACKAGE CONTENTS
================

ESVEP models
-----------
* :func:`ESVEP` Implementation of the ESVEP model.

'''

import numpy as np

import pyTSEB.meteo_utils as met
import pyTSEB.resistances as res
import pyTSEB.MO_similarity as MO
import pyTSEB.net_radiation as rad
import pyTSEB.clumping_index as CI
import pyTSEB.TSEB as tseb

# ==============================================================================
# List of constants used in ESVEP model and sub-routines
# ==============================================================================
# Change threshold in  Monin-Obukhov lengh to stop the iterations
L_thres = tseb.L_thres
# mimimun allowed friction velocity
u_friction_min = tseb.u_friction_min
# Maximum number of interations
ITERATIONS = tseb.ITERATIONS
# kB coefficient
kB = 2.2
# Invalid quality flag
F_INVALID = tseb.F_INVALID

def ESVEP(
    Tr_K,
    vza,
    T_A_K,
    u,
    ea,
    p,
    Sn_C,
    Sn_S,
    L_dn,
    LAI,
    emis_C,
    emis_S,
    z_0M,
    d_0,
    z_u,
    z_T,
    z0_soil=0.001,
    x_LAD=1,
    f_c=1.0,
    f_g=1.0,
    w_C=1.0,
    calcG_params=[
        [1],
        0.35],
    UseL=False):
    '''ESVEP

    Calculates soil and vegetation energy fluxes using Soil and Vegetation Energy Patrtitioning
    `ESVEP` model and a single observation of composite radiometric temperature.

    Parameters
    ----------
    Tr_K : float
        Radiometric composite temperature (Kelvin).
    vza : float
        View Zenith Angle (degrees).
    T_A_K : float
        Air temperature (Kelvin).
    u : float
        Wind speed above the canopy (m s-1).
    ea : float
        Water vapour pressure above the canopy (mb).
    p : float
        Atmospheric pressure (mb), use 1013 mb by default.
    Sn_C : float
        Canopy net shortwave radiation (W m-2).
    Sn_S : float
        Soil net shortwave radiation (W m-2).
    L_dn : float
        Downwelling longwave radiation (W m-2).
    LAI : float
        Effective Leaf Area Index (m2 m-2).
    emis_C : float
        Leaf emissivity.
    emis_S : flaot
        Soil emissivity.
    z_0M : float
        Aerodynamic surface roughness length for momentum transfer (m).
    d_0 : float
        Zero-plane displacement height (m).
    z_u : float
        Height of measurement of windspeed (m).
    z_T : float
        Height of measurement of air temperature (m).
    z0_soil : float, optional
        bare soil aerodynamic roughness length (m).
    x_LAD : float, optional
        Campbell 1990 leaf inclination distribution function chi parameter.
    f_c : float, optional
        Fractional cover.
    f_g : float, optional
        Fraction of vegetation that is green.
    w_C : float, optional
        Canopy width to height ratio.
    calcG_params : list[list,float or array], optional
        Method to calculate soil heat flux,parameters.
            * [[1],G_ratio]: default, estimate G as a ratio of Rn_S, default Gratio=0.35.
            * [[0],G_constant] : Use a constant G, usually use 0 to ignore the computation of G.
            * [[2,Amplitude,phase_shift,shape],time] : estimate G from Santanello and Friedl with G_param list of parameters (see :func:`~TSEB.calc_G_time_diff`).
    UseL : float or None, optional
        If included, its value will be used to force the Moning-Obukhov stability length.

    Returns
    -------
    flag : int
        Quality flag, see Appendix for description.
    T_S : float
        Soil temperature  (Kelvin).
    T_C : float
        Canopy temperature  (Kelvin).
    T_sd: float
        End-member temperature of dry soil (Kelvin)
    T_vd: float
        End-member temperature of dry vegetation (Kelvin)
    T_sw: float
        End-member temperature of saturated soil (Kelvin)
    T_vw: float
        End-member temperature of well-watered vegetation (Kelvin)
    T_star: float
        Critical surface temperature (Kelvin)
    L_nS : float
        Soil net longwave radiation (W m-2)
    L_nC : float
        Canopy net longwave radiation (W m-2)
    LE_C : float
        Canopy latent heat flux (W m-2).
    H_C : float
        Canopy sensible heat flux (W m-2).
    LE_S : float
        Soil latent heat flux (W m-2).
    H_S : float
        Soil sensible heat flux (W m-2).
    G : float
        Soil heat flux (W m-2).
    r_vw: float
        Canopy resistance to heat transport of well-watered vegetation (s m-1)
    r_vd: float
        Canopy resistance to heat transport of vegetation with zero soil water avaiability (s m-1)
    r_av: float
        Aerodynamic resistance to heat transport of the vegetation (s m-1)
    r_as: float
        Aerodynamic resistance to heat transport of the soil (s m-1)
    L : float
        Monin-Obuhkov length (m).
    n_iterations : int
        number of iterations until convergence of L.

    References
    ----------
    .. [Tang2017] Tang, R., and Z. L. Li. An End-Member-Based Two-Source Approach for Estimating  
        Land Surface Evapotranspiration From Remote Sensing Data. IEEE Transactions on Geoscience 
        and Remote Sensing 55, no. 10 (October 2017): 5818–32. 
        https://doi.org/10.1109/TGRS.2017.2715361.
    '''

    # Convert input float scalars to arrays and check parameters size
    Tr_K = np.asarray(Tr_K)
    (vza,
     T_A_K,
     u,
     ea,
     p,
     Sn_C,
     Sn_S,
     L_dn,
     LAI,
     emis_C,
     emis_S,
     z_0M,
     d_0,
     z_u,
     z_T,
     z0_soil,
     x_LAD,
     f_c,
     f_g,
     w_C,
     calcG_array) = map(tseb._check_default_parameter_size,
                        [vza,
                         T_A_K,
                         u,
                         ea,
                         p,
                         Sn_C,
                         Sn_S,
                         L_dn,
                         LAI,
                         emis_C,
                         emis_S,
                         z_0M,
                         d_0,
                         z_u,
                         z_T,
                         z0_soil,
                         x_LAD,
                         f_c,
                         f_g,
                         w_C,
                         calcG_params[1]],
                        [Tr_K] * 21)

    # Create the output variables
    [flag, T_S, T_C, T_sd, T_vd, T_sw, T_vw, T_star, Ln_S, Ln_C, LE_C, H_C, LE_S, H_S, G, r_vw,
     r_vd, r_av, r_as, iterations] = [np.zeros(Tr_K.shape)+np.NaN for i in range(20)]

    # iteration of the Monin-Obukhov length
    if isinstance(UseL, bool):
        # Initially assume stable atmospheric conditions and set variables for
        L = np.asarray(np.zeros(T_S.shape) + np.inf)
        max_iterations = ITERATIONS
    else:  # We force Monin-Obukhov lenght to the provided array/value
        L = np.asarray(np.ones(T_S.shape) * UseL)
        max_iterations = 1  # No iteration

    # Calculate the general parameters
    rho = met.calc_rho(p, ea, T_A_K)  # Air density
    c_p = met.calc_c_p(p, ea)  # Heat capacity of air
    z_0H = res.calc_z_0H(z_0M, kB=kB)  # Roughness length for heat transport
    z_0H_soil = res.calc_z_0H(z0_soil, kB=kB)  # Roughness length for heat transport
    s = met.calc_delta_vapor_pressure(T_A_K) * 10  # slope of the saturation pressure curve (mb C-1)
    lbd = met.calc_lambda(T_A_K)  # latent heat of vaporisation (MJ./kg)
    gama = met.calc_psicr(p, lbd)  # psychrometric constant (mb C-1)
    vpd = met.calc_vapor_pressure(T_A_K) - ea  # vapor pressure deficit (mb)

    # Calculate LAI dependent parameters for dataset where LAI > 0
    omega0 = CI.calc_omega0_Kustas(LAI, f_c, x_LAD=x_LAD, isLAIeff=True)
    F = np.asarray(LAI / f_c)  # Real LAI
    # Fraction of vegetation observed by the sensor
    f_theta = tseb.calc_F_theta_campbell(vza, F, w_C=w_C, Omega0=omega0, x_LAD=x_LAD)

    # Initially assume stable atmospheric conditions and set variables for
    # iteration of the Monin-Obukhov length
    u_friction = MO.calc_u_star(u, z_u, L, d_0, z_0M)
    u_friction = np.asarray(np.maximum(u_friction_min, u_friction))
    u_friction_s = MO.calc_u_star(u, z_u, L, np.zeros(d_0.shape), z0_soil)
    u_friction_s = np.asarray(np.maximum(u_friction_min, u_friction_s))
    L_old = np.ones(Tr_K.shape)
    L_diff = np.asarray(np.ones(Tr_K.shape) * float('inf'))

    # First assume that canopy temperature equals the minumum of Air or
    # radiometric T
    T_C = np.asarray(np.minimum(Tr_K, T_A_K))
    flag, T_S = tseb.calc_T_S(Tr_K, T_C, f_theta)

    # Loop for estimating stability.
    # Stops when difference in consecutives L is below a given threshold
    for n_iterations in range(max_iterations):
        i = flag != F_INVALID
        if np.all(L_diff[i] < L_thres):
            if L_diff[i].size == 0:
                print("Finished iterations with no valid solution")
            else:
                print("Finished interations with a max. L diff: " + str(np.max(L_diff[i])))
            break
        i = np.logical_and(L_diff >= L_thres, flag != F_INVALID)
        print("Iteration " + str(n_iterations) + ", max. L diff: " + str(np.max(L_diff[i])))
        iterations[i] = n_iterations

        # Calculate net longwave radiation with current values of T_C and T_S
        Ln_C[i], Ln_S[i] = rad.calc_L_n_Kustas(
            T_C[i], T_S[i], L_dn[i], LAI[i], emis_C[i], emis_S[i])
        Rn_C = Sn_C + Ln_C
        Rn_S = Sn_S + Ln_S

        # Compute Soil Heat Flux
        G[i] = tseb.calc_G([calcG_params[0], calcG_array], Rn_S, i)

        # Calculate aerodynamic resistances
        r_vw[i] = 100.0/LAI[i]
        r_vd[i] = np.zeros(LAI[i].shape) + 2000.0
        r_av_params = {"z_T": z_T[i], "u_friction": u_friction[i], "L": L[i],
                       "d_0": d_0[i], "z_0H": z_0H[i]}
        r_av[i] = tseb.calc_resistances(tseb.KUSTAS_NORMAN_1999, {"R_A": r_av_params})[0]
        r_as_params = {"z_T": z_T[i], "u_friction": u_friction_s[i], "L": L[i],
                       "d_0": np.zeros(d_0[i].shape), "z_0H": z_0H_soil[i]}
        r_as[i] = tseb.calc_resistances(tseb.KUSTAS_NORMAN_1999, {"R_A": r_as_params})[0]

        # Estimate the surface temperatures of the end-members
        # Eq 8a
        T_sd[i] = r_as[i] * (Rn_S[i] - G[i]) / (rho[i] * c_p[i]) + T_A_K[i]
        # Eq 8b
        T_vd[i] = r_av[i] * Rn_C[i] / (rho[i] * c_p[i]) *\
                  gama[i] * (1 + r_vd[i] / r_av[i]) / (s[i] + gama[i] * (1 + r_vd[i] / r_av[i])) -\
                  vpd[i] / (s[i] + gama[i] * (1 + r_vd[i] / r_av[i])) + T_A_K[i]
        # Eq 8c
        T_sw[i] = r_as[i] * (Rn_S[i] - G[i]) / (rho[i] + c_p[i]) *\
                  gama[i] / (s[i] + gama[i]) - vpd[i] / (s[i] + gama[i]) + T_A_K[i]
        # Eq 8d
        T_vw[i] = r_av[i] * Rn_C[i] / (rho[i] * c_p[i]) *\
                  gama[i] * (1 + r_vw[i] / r_av[i]) / (s[i] + gama[i] * (1 + r_vw[i] / r_av[i])) -\
                  vpd[i] / (s[i] + gama[i] * (1 + r_vw[i] / r_av[i])) + T_A_K[i]

        # Estimate critical surface temperature - eq 10
        T_star[i] = (T_sd[i]**4 * (1 - f_theta[i]) + T_vw[i]**4 * f_theta[i])**0.25

        # Estimate latent heat fluxes when water in the top-soil is avaiable for evaporation
        j = np.logical_and(Tr_K <= T_star, i)
        # Eq 12a
        T_C[j] = T_vw[j]
        # Eq 12b
        T_S[j] = ((Tr_K[j]**4 - f_theta[j] * T_C[j]**4) / (1 - f_theta[j]))**0.25
        # Eq 13a
        LE_C[j] = (s[j] * Rn_C[j] + rho[j] * c_p[j] * vpd[j] / r_av[j]) /\
                  (s[j] + gama[j] * (1 + r_vw[j] / r_av[j]))
        # Eq 13b
        LE_S[j] = (T_sd[j] - T_S[j]) / (T_sd[j] - T_sw[j]) *\
                  ((s[j] * (Rn_S[j] - G[j]) + rho[j] * c_p[j] * vpd[j] / r_as[j]) /
                   (s[j] + gama[j]))

        # Estimate latent heat fluxes when no water in the top-soil is avaiable for evaporation
        j = np.logical_and(Tr_K > T_star, i)
        # Eq 14a
        T_S[j] = T_sd[j]
        # Eq 14b
        T_C[j] = ((Tr_K[j]**4 - (1 - f_theta[j]) * T_S[j]**4) / f_theta[j])**0.25
        # Eq 15a
        LE_C[j] = (T_vd[j] - T_C[j]) / (T_vd[j] - T_vw[j]) *\
                  ((s[j] * Rn_C[j] + rho[j] * c_p[j] * vpd [j] / r_av[j]) /
                   (s[j] + gama[j] * (1 + r_vw[j] / r_av[j])))
        # Eq 15b
        LE_S[j] = 0

        # Estimate sensible heat fluxes as residuals
        H_C[i] = Rn_C[i] - LE_C[i]
        H_S[i] = Rn_S[i] - G[i] - LE_S[i]

        # Calculate total fluxes
        H = np.asarray(H_C + H_S)
        LE = np.asarray(LE_C + LE_S)

        # Now L can be recalculated and the difference between iterations
        # derived
        if isinstance(UseL, bool):
            L[i] = MO.calc_L(u_friction[i], T_A_K[i], rho[i], c_p[i], H[i], LE[i])
            L_diff = np.asarray(np.fabs(L - L_old) / np.fabs(L_old))
            L_diff[np.isnan(L_diff)] = float('inf')
            L_old = np.array(L)
            L_old[L_old == 0] = 1e-36
            # Calculate again the friction velocity with the new stability
            # correctios
            u_friction[i] = MO.calc_u_star(u[i], z_u[i], L[i], d_0[i], z_0M[i])
            u_friction[i] = np.asarray(np.maximum(u_friction_min, u_friction[i]))
            u_friction_s[i] = MO.calc_u_star(u[i], z_u[i], L[i], np.zeros(d_0[i].shape),
                                          np.zeros(z_0M[i].shape) + 0.005)
            u_friction_s[i] = np.asarray(np.maximum(u_friction_min, u_friction_s[i]))

    return [flag, T_S, T_C,  T_sd, T_vd, T_sw, T_vw, T_star, Ln_S, Ln_C, LE_C, H_C, LE_S, H_S, G, 
            r_vw, r_vd, r_av, r_as, u_friction, L, n_iterations]