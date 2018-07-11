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
'''

from collections import OrderedDict

import numpy as np

from pyTSEB.PyTSEB import PyTSEB, S_N, S_P, S_A

from pyESVEP import ESVEP
from pyTSEB import TSEB


class PyESVEP(PyTSEB):

    def __init__(self, parameters):
        if "resistance_form" not in parameters.keys():
            parameters["resistance_form"] = 0
        super().__init__(parameters)

    def _get_input_structure(self):
        ''' Input fields' names for ESVEP model.  Only relevant for image processing mode.

        Parameters
        ----------
        None

        Returns
        -------
        outputStructure: string ordered dict
            Names (keys) and descriptions (values) of TSEB_2T input fields.
        '''

        input_fields = super()._get_input_structure()
        del input_fields["leaf_width"]
        del input_fields["alpha_PT"]
        del input_fields["KN_b"]
        del input_fields["KN_c"]
        del input_fields["KN_C_dash"]
        return input_fields

    def _set_special_model_input(self, field, dims):
        ''' Special processing for setting certain input fields. Only relevant for image processing
        mode.

        Parameters
        ----------
        field : string
            The name of the input field for which the special processing is needed.
        dims : int list
            The dimensions of the output parameter array.

        Returns
        -------
        success : boolean
            True is the parameter was succefully set, false otherwise.
        array : float array
            The set parameter array.
        '''

        # Those fields are optional for ESVEP.
        if field in ["x_LAD", "f_c", "f_g", "w_C"]:
            success, val = self._set_param_array(field, dims)
            if not success:
                val = np.ones(dims)
                success = True
        else:
            success = False
            val = None
        return success, val

    def _get_output_structure(self):
        ''' Output fields' names for ESVEP model.

        Parameters
        ----------
        None

        Returns
        -------
        output_structure: ordered dict
            Names of the output fields as keys and instructions on whether the output
            should be saved to file as values.
        '''

        output_structure = OrderedDict([
            # Energy fluxes
            ('R_n1', S_P),   # net radiation reaching the surface
            ('R_ns1', S_A),  # net shortwave radiation reaching the surface
            ('R_nl1', S_A),  # net longwave radiation reaching the surface
            ('delta_R_n1', S_A),  # net radiation divergence in the canopy
            ('Sn_S1', S_N),  # Shortwave radiation reaching the soil
            ('Sn_C1', S_N),  # Shortwave radiation intercepted by the canopy
            ('Ln_S1', S_N),  # Longwave radiation reaching the soil
            ('Ln_C1', S_N),  # Longwave radiation intercepted by the canopy
            ('H_C1', S_A),  # canopy sensible heat flux (W/m^2)
            ('H_S1', S_N),  # soil sensible heat flux (W/m^2)
            ('H1', S_P),  # total sensible heat flux (W/m^2)
            ('LE_C1', S_A),  # canopy latent heat flux (W/m^2)
            ('LE_S1', S_N),  # soil latent heat flux (W/m^2)
            ('LE1', S_P),  # total latent heat flux (W/m^2)
            ('LE_partition', S_A),  # Latent Heat Flux Partition (LEc/LE)
            ('G1', S_P),  # ground heat flux (W/m^2)
            # temperatures (might not be accurate)
            ('T_C1', S_A),  # Canopy temperature (Kelvin)
            ('T_S1', S_A),  # Soil temperature  (Kelvin)
            ('T_sd', S_A),  # End-member temperature of dry soil (Kelvin)
            ('T_vd', S_A),  # End-member temperature of dry vegetation (Kelvin)
            ('T_sw', S_A),  # End-member temperature of saturated soil (Kelvin)
            ('T_vw', S_A),  # End-member temperature of well-watered vegetation (Kelvin)
            ('T_star', S_A),  # Critical surface temperature (Kelvin)
            # resistances
            ('r_vw', S_A),  # Canopy resistance to heat transport of well-watered vegetation (s m-1)
            ('r_vd', S_A),  # Canopy resistance to heat transport of vegetation with zero soil water avaiability (s m-1)
            ('r_av', S_A),  # Aerodynamic resistance to heat transport of the vegetation (s m-1)
            ('r_as', S_A),  # Aerodynamic resistance to heat transport of the soil (s m-1)
            # miscaleneous
            ('omega0', S_N),  # nadir view vegetation clumping factor
            ('L', S_A),  # Monin Obukhov Length
            ('u_friction', S_A),  # Friction velocity (m/s)
            ('theta_s1', S_N),  # Sun zenith angle
            ('F', S_N),  # Leaf Area Index
            ('z_0M', S_N),  # Aerodynamic roughness length for momentum trasport (m)
            ('d_0', S_N),  # Zero-plane displacement height (m)
            ('Skyl', S_N),
            ('flag', S_A),  # Quality flag
            ('n_iterations', S_N)])  # Number of iterations before model converged to stable value

        return output_structure

    def _call_flux_model_veg(self, in_data, out_data, model_params, i):
        ''' Call ESVEP model to calculate fluxes for data points containing vegetation.

        Parameters
        ----------
        in_data : dict
            The input data for the model.
        out_data : dict
            Dict containing the output data from the model which will be updated. It also contains
            previusly calculated shortwave radiation and roughness values which are used as input
            data.

        Returns
        -------
        None
        '''

        [out_data['flag'][i], out_data['T_S1'][i], out_data['T_C1'][i], out_data['T_sd'][i],
         out_data['T_vd'][i], out_data['T_sw'][i], out_data['T_vd'][i], out_data['T_star'][i],
         out_data['Ln_S1'][i], out_data['Ln_C1'][i], out_data['LE_C1'][i], out_data['H_C1'][i],
         out_data['LE_S1'][i], out_data['H_S1'][i], out_data['G1'][i], out_data['r_vw'][i],
         out_data['r_vd'][i], out_data['r_av'][i], out_data['r_as'][i], out_data['u_friction'][i],
         out_data['L'][i], out_data['n_iterations'][i]] = \
                 ESVEP.ESVEP(in_data['T_R1'][i],
                             in_data['VZA'][i],
                             in_data['T_A1'][i],
                             in_data['u'][i],
                             in_data['ea'][i],
                             in_data['p'][i],
                             out_data['Sn_C1'][i],
                             out_data['Sn_S1'][i],
                             in_data['L_dn'][i],
                             in_data['LAI'][i],
                             in_data['emis_C'][i],
                             in_data['emis_S'][i],
                             out_data['z_0M'][i],
                             out_data['d_0'][i],
                             in_data['z_u'][i],
                             in_data['z_T'][i],
                             z0_soil=in_data['z0_soil'][i],
                             x_LAD=in_data['x_LAD'][i],
                             f_c=in_data['f_c'][i],
                             f_g=in_data['f_g'][i],
                             w_C=in_data['w_C'][i],
                             calcG_params=model_params["calcG_params"])

    def _call_flux_model_soil(self, in_data, out_data, model_params, i):
        ''' Call a OSEB model to calculate soil fluxes for data points containing no vegetation.
        NOTE: For now ESVEP does not process pixels with no vegetation.

        Parameters
        ----------
        in_data : dict
            The input data for the model.
        out_data : dict
            Dict containing the output data from the model which will be updated. It also contains
            previusly calculated shortwave radiation and roughness values which are used as input
            data.

        Returns
        -------
        None
        '''

        [out_data['flag'][i],
         out_data['Ln_S1'][i],
         out_data['LE_S1'][i],
         out_data['H_S1'][i],
         out_data['G1'][i],
         out_data['r_as'][i],
         out_data['u_friction'][i],
         out_data['L'][i],
         out_data['n_iterations'][i]] = TSEB.OSEB(in_data['T_R1'][i],
                                                  in_data['T_A1'][i],
                                                  in_data['u'][i],
                                                  in_data['ea'][i],
                                                  in_data['p'][i],
                                                  out_data['Sn_S1'][i],
                                                  in_data['L_dn'][i],
                                                  in_data['emis_S'][i],
                                                  out_data['z_0M'][i],
                                                  out_data['d_0'][i],
                                                  in_data['z_u'][i],
                                                  in_data['z_T'][i],
                                                  calcG_params=model_params["calcG_params"])
