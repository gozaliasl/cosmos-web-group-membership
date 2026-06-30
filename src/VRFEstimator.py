#
# Copyright (C) 2012-2020 Euclid Science Ground Segment
#
# This library is free software; you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 3.0 of the License, or (at your option)
# any later version.
#
# This library is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this library; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
#

"""
:file: python/LE3_SIGV_CL/VRFEstimator.py

:date: 07/24/25
:author:
        Z. Ghaffari
        A. Biviano

"""
import os
import sys
import logging
import numpy as np
from astropy import units as uu
from astropy.table import Table
import matplotlib.pyplot as plt
from astropy import constants as const

# --- Configuration for logging ---
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class VRFEstimator:
    """
    Class for iterative outlier rejection and rest frame velocity dispersion
    estimation.
    """

    def __init__(self, cfg, gals, detcat):
        """
        Initializes the VRFEstimator with parameters from a config file.

        Parameters
        ----------
        '''Constructor of the class

        Parameters
        ----------
        cfg : configurator instance
            Contains the information for computing the redshift (see
            setup method for more information)
        gals : instance of galaxy_catalog
            Contains the information on the galaxies


        """
        logger.info("VRFEstimator class initialized with parameters.")
        self._cfg = cfg

        self._gals = gals
        self._detcat = detcat

        # Parameters from a config file
        self.anis = self._cfg['anis']
        self.delta = self._cfg['delta']
        self.itmax = self._cfg['itmax']
        self.evel = self._cfg['evel']
        self.nscut = self._cfg['nscut']
        self.rvguess = self._cfg['rvguess']

        # Cosmology parameters
        self.H0 = self._cfg['H0']
        self.Om0 = self._cfg['Om0']
        self.Ode0 = self._cfg['Ode0']

    def midmean(self, x):
        """
        Calculates the midmean (MIDMEAN = mean between 25% and 75% quartiles)
        """
        x_np = np.asarray(x)
        if x_np.size == 0:
            error_message = "Input array is empty. Cannot compute midmean."
            logger.error(error_message)
            raise ValueError(error_message)

        if x_np.size == 1:
            error_message = f"Input array has only one element ({x_np[0]}). Midmean requires a range of values."
            logger.error(error_message)
            raise ValueError(error_message)

        xs = np.sort(x_np)
        n = len(xs) - 1

        q25 = int(0.25 * n + 0.5)
        q75 = int(0.75 * n + 0.5)

        # Ensure indices are valid
        q25 = max(0, min(q25, n))
        q75 = max(0, min(q75, n))

        if q25 > q75:
            return np.mean(xs)

        logger.info(
            f"mean between 25% and 75% quartiles: {np.mean(xs[q25: q75 + 1])}")
        return np.mean(xs[q25: q75 + 1])

    def even_odd_median(self, vec, nd):
        """Compute the median of a vector with known length."""
        if nd % 2 == 0:
            return (vec[nd // 2 - 1] + vec[nd // 2]) / 2.0
        else:
            return vec[nd // 2]

    # --- gapclean function ---

    def gapclean(self, vz, evel, wgmax, wgap=False):
        """
        Performs robust velocity outlier rejection based on gap analysis.
        look for weighted gaps larger than wgmax; if there are no large weighted gaps, return 1 for all gals, otherwise selects galaxies in the largest among the groups separated by large gaps, with velocities within 4000 km/s of <v>=0
        Args:
            vz : Input velocities
            evel : Velocity Error
            wgmax : Weighted gaps
            wgap : If True, performs weighted gap analysis.
                                    If False, performs simple velocity clipping.

        Returns:
            An array of 0s and 1s, where 1 indicates a 'clean' velocity.
                        Returns np.nan if a critical numerical issue occurs (e.g., midmean returns nan).


        """
        # Sort velocities and get count
        vz_np = np.asarray(vz)
        vs = np.sort(vz_np)
        n = vs.size

        # --- invalid inputs ---
        if n == 0:
            error_message = "Input 'vz' array is empty. Gap analysis requires at least two velocities."
            logger.error(error_message)
            raise ValueError(error_message)
        if n == 1:
            error_message = f"Input 'vz' array has only one element ({vz_np[0]}). Gap analysis requires at least two velocities."
            logger.error(error_message)
            raise ValueError(error_message)

        # Initialize min/max velocities and clean flag array
        vmin = vs[1]
        vmax = vs[n-1]
        iclean = np.zeros(n, dtype=int)

        if wgap:
            # Calculate weighted gaps based on Beers+91, element-wise subtraction
            gap = vs[1:] - vs[:-1]

            # Filter gaps smaller than error threshold
            wg0_indices = np.where(gap < np.sqrt(2.) * evel)[0]
            if wg0_indices.size > 0:
                gap[wg0_indices] = 0.
            logger.info(
                f"The measured gaps between the ordered velocities: {gap} ")
            # Calculate weights for gaps
            wt = (np.arange(n-1) + 1) * (n - (np.arange(n-1) + 1))
            logger.info(
                f"A set of approximately Gaussian weights: {wt} ")

            wtdgap = np.power(wt * gap, 0.5)
            logger.info(
                f"A weighted gap  : {wtdgap} ")

            # Normalize weighted gaps by midmean
            y = wtdgap / self.midmean(wtdgap)

            logger.info(
                f"Normalized weighted gap : {y} ")

            # Find significant gaps
            wg_idx = np.where(y > wgmax)[0]

            nwg = wg_idx.size
            if nwg == 0:
                # All velocities are clean if no significant gaps
                iclean[:] = 1
            else:
                # Isolate largest coherent group
                jg = np.concatenate(([-1], wg_idx, [n - 1]))

                dva = np.zeros(nwg + 1, dtype=int)

                for j in range(nwg + 1):
                    dva[j] = jg[j+1] - (jg[j] + 1)

                jmin = np.argmax(dva)
                vmin = vs[jg[jmin] + 1]
                vmax = vs[jg[jmin+1]]
        else:  # Simple velocity clipping
            vmin = -wgmax
            vmax = +wgmax

        # Apply cleaning flags to original velocities
        clean_indices = np.where((vz_np >= vmin) & (vz_np <= vmax))[0]
        iclean[clean_indices] = 1

        logger.info(
            f"Clean inices are : {iclean} ")

        return iclean

    def cofmvir(self, mvir, auth):
        """
        Calculates the concentration parameter for LCDM halos based on virial mass
        and a specified author/model relation
        (Maccio, Dutton & van den Bosch 08, relaxed halos, WMAP5, Delta=200)
        ;(Maccio, Dutton & van den Bosch 08, relaxed halos, WMAP5, Delta=95)


        Args:
            mvir : The virial mass of the halo.
            auth : the author/model relation.
                        Accepted values: 'MDvdB200', 'MDvdB'.

        Returns:
            The calculated concentration parameter.
        """
        if auth == 'MDvdB200':
            slope = -0.098
            norm = 6.76
        elif auth == 'MDvdB':
            slope = -0.094
            norm = 9.35
        else:
            raise ValueError(f"Unknown 'auth' parameter: {auth}. "
                             "Accepted values are 'MDvdB200' or 'MDvdB'.")

        #
        concentration = norm * np.power(mvir, slope)
        logger.info(
            f"Concentration  : {concentration} ")
        return concentration

    def sapovervv(self, c, anis):
        """
        Calculates the ratio of aperture velocity dispersion to virial velocity
        for an NFW model at R=r_vir, based on concentration and anisotropy.
        Args:
            c : The concentration parameter of the NFW halo.
            anis : The anisotropy model. Accepted values: 'isotropic', 'ML'.

        Returns:
            The calculated ratio (v_vir / sigma_ap).
        """

        if anis == 'isotropic':
            # Coefficients for isotropic velocity dispersion
            cfs = np.array([0.344165, 0.312747, -0.336547, 0.0390675])
        elif anis == 'ML':
            # Coefficients for Mamon-Lokas type anisotropy
            cfs = np.array([0.270654, 0.362998, -0.312684, 0.0161601])
        else:
            raise ValueError(f"Unknown 'anis' parameter: {anis}. "
                             "Accepted values are 'isotropic' or 'ML'.")

        # Calculate log base 10 of concentration
        if isinstance(c, uu.Quantity):
            c = c.value
        else:
            c = c
        lc = np.log10(c)

        poly = cfs[0] + cfs[1]*lc + cfs[2] * \
            np.power(lc, 2) + cfs[3]*np.power(lc, 3)

        # Convert back from log space
        c2nfwapx = np.power(10., poly)

        # Calculate the final ratio (v_vir / sigma_ap)
        calculated_sapovervv = 1. / np.sqrt(c2nfwapx)
        logger.info(
            f"Aperture velocity dispersion of NFW model at R=r_vir  : {calculated_sapovervv} , anisotropy = {anis}")
        return calculated_sapovervv

    def siglosnfwml1(self, x):
        """
        Calculates the approximate line-of-sight velocity dispersion for an NFW model
        with ML anisotropy at a specific scaled radius.

        Args:
            x : The scaled radius, typically (concentration * r_proj / r_vir),
                equivalent to r_proj / r_s.

        Returns:
            The calculated scaled line-of-sight velocity dispersion.
            Returns NaN if input x is invalid (e.g., <= 0).
        """
        # Coefficients for the 7th-order polynomial fit
        cfs = np.array([-0.14783, -0.110877, -0.135747, 0.00194757,
                        0.0231745, 0.000631017, -0.00323355, -0.000637035])

        # Ensure x is positive to avoid log(0) or log(negative)
        if isinstance(x, np.ndarray):
            valid_indices = x.value > 0
            lx = np.full(x.shape, np.nan, dtype=float)  # Initialize with NaN
            if np.any(valid_indices):
                lx[valid_indices] = np.log10(x[valid_indices].value)
        else:
            if x <= 0:
                return np.nan
            lx = np.log10(x.value)

        # Evaluate the 7th-order polynomial
        poly = (cfs[0] + cfs[1]*lx + cfs[2]*np.power(lx, 2) + cfs[3]*np.power(lx, 3) +
                cfs[4]*np.power(lx, 4) + cfs[5]*np.power(lx, 5) + cfs[6]*np.power(lx, 6) +
                cfs[7]*np.power(lx, 7))

        # Convert back from log space to linear scale
        calculated_siglosnfwml1 = np.power(10., poly)
        logger.info(
            f"sigma_los approx for NFW with ML anisotropy : {calculated_siglosnfwml1}")
        return calculated_siglosnfwml1

    def massnfw(self, x, a, rvir):
        """
        Calculates the NFW mass function, representing the fraction of mass
        enclosed within radius 'x' relative to the mass within 'rvir'.

        Args:
            x : The radius at which to calculate the enclosed mass.
            a : The scale radius (r_s) of the NFW profile.
            rvir : The virial radius (r_vir) of the halo.

        Returns:
        The fraction of the total halo mass within radius 'x'.
        """
        if isinstance(x, uu.Quantity):
            x = x.value
        else:
            x = x
        if isinstance(a, uu.Quantity):
            a = a.value
        else:
            a = a

        # Check if rvir is a quantity and get its value
        if isinstance(rvir, uu.Quantity):
            rvir = rvir.value
        else:
            rvir = rvir
        
        numerator = np.log(1. + x / a) - x / (x + a)
        denominator = np.log(1. + rvir / a) - rvir / (rvir + a)

        if np.isclose(denominator, 0.0):
            return np.nan

        calculated_massnfw = numerator / denominator

        logger.info(
            f"The fraction of the total halo mass within radius 'x' : {calculated_massnfw}")
        return calculated_massnfw

    def interlopmbmiter(self, rorig, zgal, zcl, *args, **kwargs):
        """
        Interloper selection using the algorithm of MBM10 with iterations.

        Parameters
        ----------
        rorig :
            Radial distances of galaxies from the detection center.
        zgal :
            Galaxies' redshifts (zph or zspec).
        zcl :
            Cluster redshift.
        *args, **kwargs: Additional arguments.



        itmax : Max number of iterations.
        *args: This will capture the positional output arguments (imembers, mdelta, cleanr, cleanvz, flag).
        **kwargs: This will capture all other optional keyword arguments.

        Notes:
        if the keyword widegap is set, then check for wide gaps on 1st iteration of size widegap (typical value is widegap=4) if the keyword plot is set, then make an r,v plot of members and interlopers.

        if the keyword nsigma is set, then cut at nsigma rather than at 2.7 (default)

        if the keywords v200, conc, vmed are set, use the input v200, conc, vmed rather than computing v200 from sigmav, conc from v200 and vmed from the velocities.

        if a galaxy is located at r=0 it would be excluded, since the predicted velocity dispersion at r=0 would be 0, given than M(r=0)=0; so we offset by 1 kpc

        Returns:
            tuple: (imembers, mdelta, cleanr, cleanvz, flag)
            imembers : 1s for members, 0s for interlopers.
            mdelta : A (rough) mass estimate in Msun.
            cleanr : Cleaned radii (members).
            cleanvz : Cleaned velocities (members).
            flag : 0 if algorithm succeeded, 1 if failed

        """

        # Initialize variables from kwargs and set defaults
        logger.info("--- Starting interlopmbmiter ---")
        logger.info(
            f"Initial inputs: rorig_len={len(rorig)}, vz_len={len(zgal)}, evel={self.evel}, delta={self.delta}, zcl={zcl}, itmax={self.itmax}")

        flag = 0  # if flag=1 on return, the algorithm has failed

        # Keyword arguments
        widegap = kwargs.get('widegap', None)
        physgap = kwargs.get('physgap', None)
        plot_flag = kwargs.get('plot', False)
        nsigma = kwargs.get('nsigma', 2.7)
        c200 = kwargs.get('c200', None)
        v200 = kwargs.get('v200', None)
        conc = kwargs.get('conc', None)
        vmed = kwargs.get('vmed', None)
        noiter = kwargs.get('noiter', False)  # no iterations beyond first pass

        c = None
        if c200 is not None:
            c = c200

        # Define r
        r = np.asarray(rorig, dtype=float)

        # Offset by 1 kpc if r is less than 1, to avoid issues where r=0
        w1 = np.where(r < 1.0)[0]
        if len(w1) > 0:
            r[w1] = 1.0
        logger.info(f"Radii (r) adjusted for r < 1kpc: {r}")

        nscut = 2.7
        if nsigma is not None:
            nscut = nsigma

        # initialize value of virial radius
        rvguess = 1.e14

        # Gravitational constant, 43 km2 Mpc / (solMass s2)
        grav = const.G.to((uu.km/uu.s)**2 * uu.Mpc / uu.M_sun) * 1e10
        # logger.info(f"The gravitational constant : {grav}")

        # Speed of light, 2.99792458e5 km/s
        clight = const.c.to(uu.km/uu.s).value
        # logger.info(f"The value of clight constant : {clight}")

        # Omega_0
        Omega0 = 0.3

        # Omega_Lambda
        OmegaLambda = 0.7

        logger.info(
            f"Constants: grav={grav}, clight={clight}, Omega0={Omega0}, OmegaLambda={OmegaLambda}")

        # velocities in rest-frame
        vrf = np.asarray(zgal, dtype=float)

        logger.info('Calculate rest-frame velocities')
        vrf = (zgal - zcl) * clight / (1.0 + zcl)

        logger.info(f"Initial vrf (rest-frame velocities): {vrf}")

        # check for wide gaps on 1st step
        nd = len(vrf)
        imembers = np.zeros(nd, dtype=int)
        logger.info(f"{nd} data before gapping")

        cleanr = np.copy(r)
        cleanvz = np.copy(vrf)

        if widegap is not None:
            logger.info(f"Applying widegap cleaning with wgmax={widegap}")
            iclean = self.gapclean(vrf, self.evel, widegap, wgap=True)
            cleanr = r[np.where(iclean)[0]]
            cleanvz = vrf[np.where(iclean)[0]]

        if physgap is not None:
            logger.info(f"Applying physgap cleaning with wgmax={physgap}")
            # wgap=False for physgap
            iclean = self.gapclean(vrf, self.evel, physgap, wgap=False)
            cleanr = r[np.where(iclean)[0]]
            cleanvz = vrf[np.where(iclean)[0]]

        nd = len(cleanvz)
        logger.info(f"{nd} data after gapping")
        logger.info(f" Gapped data set is: {np.sort(cleanvz)}")

        # iterations
        ipass = 0

        while True:  # Equivalent of IDL's GOTO loop
            logger.info(f"\n--- Iteration no. {ipass} ---")

            # r,v plot if requested
            y1 = min(-6000., np.min(cleanvz)) - 500.
            y2 = max(6000., np.max(cleanvz)) + 500.

            if plot_flag:
                plt.figure()
                plt.plot(r, vrf, 'o', markerfacecolor='none',
                         markeredgecolor='gray', markersize=8, label='Original Data')

                plt.plot(cleanr, cleanvz, 'o', color='blue',
                         markersize=8, label='Cleaned Members')

                if ipass > 0:
                    # Need to sort r1 for plotting lines correctly
                    # r1 comes from the previous iteration's cleanr
                    s1 = np.argsort(r1)
                    plt.plot(r1[s1], vzmax[s1]+vzmed,
                             'r--', label='Upper Limit')
                    plt.plot(r1[s1], -vzmax[s1]+vzmed,
                             'r--', label='Lower Limit')

                plt.xlabel('R (kpc)')
                plt.ylabel('Velocity (km/s)')
                plt.title(f'Iteration {ipass} - Members and Limits')
                plt.grid(True)
                plt.legend()
                plt.show()

            # Define the median velocity
            nd = len(cleanvz)
            logger.info(f"{nd} members selected")
            if nd == 0:
                logger.warning("No members left - algorithm failure")
                flag = 1
                # Exit loop
                break

            vzs = np.sort(cleanvz)
            vzmed = self.even_odd_median(vzs, nd)
            logger.info(f"Median velocity is {vzmed}")

            if vmed is not None:
                vzmed = vmed
                logger.info(f"Adopting the input vmed: {vzmed}")

            # h(z)
            h0 = 0.7
            zmed = vzmed / clight + zcl
            hz = h0 * np.sqrt(Omega0 * (1.0 + zmed)**3 + OmegaLambda)
            logger.info(f"zmed={zmed}, H(z)={hz*100.} km/s/Mpc")

            # Define the std dispersion within rvguess
            wrv = np.where(cleanr <= rvguess)[0]

            # At least 2 elements needed for sigma
            if len(wrv) < 2:
                logger.warning(
                    "Few elements left for sigma calculation - algorithm failure")
                flag = 1
                # Exit loop
                break

            # Sigma within rvguess
            sigmav = np.std(cleanvz[wrv], ddof=0)
            logger.info(f"Sigma_v (within rvguess) = {sigmav}")

            # Define the median absolute deviation MAD
            dvzmed = np.abs(vzs - vzmed)
            dvzs = np.sort(dvzmed)
            mad = self.even_odd_median(dvzs, nd)

            # As in BFG90
            sigmavmad = mad / 0.6745

            logger.info(f"sigmav_MAD is {sigmavmad}, sigmav is {sigmav}")

            # If sigmavmad is 0 [may happen for groups with low-res redshifts]
            # take sigmavmad=sigmav
            if sigmavmad == 0:
                sigmavmad = sigmav
                logger.warning("sigmavmad was 0, set to sigmav.")

            if ipass == 0:
                if c200 is None:
                    c = 4.0
                    sigmaap = np.sqrt(np.maximum(
                        sigmavmad**2 - self.evel**2, 1.0))
                    rvguessold = 1.0
            else:
                rvguessold = rvguess
                # Uses sigmav for subsequent iterations
                sigmaap = np.sqrt(np.maximum(sigmav**2 - self.evel**2, 1.0))

            logger.info(
                f"ipass={ipass}, sigmaap={sigmaap}, rvguessold={rvguessold}")

            # guessing rvir using sigma_ap/v_vir as predicted
            # from NFW with concentration of c (=4 on 1st pass)
            anis = 'ML'
            vvguess = sigmaap / self.sapovervv(c, anis)

            if v200 is not None and (noiter or ipass == 0):
                vvguess = v200
                logger.info(f"Adopting the input v200: {vvguess}")
            # check again just in case, though checked above
            if len(wrv) <= 1:
                logger.warning(
                    "Too few elements for NFW calculations - algorithm failure")
                flag = 1
                # Exit loop
                break

            logger.info(
                f"Using {len(wrv)} galaxies, sigma_ap(<rvir) is {sigmaap}, vvguess is {vvguess}")

            # Ensure vvguess is not zero or problematic before division
            if vvguess == 0:
                logger.warning(
                    "vvguess is zero. Cannot compute rvguess and mvguess. Setting to NaN.")
                rvguess = np.nan
                mvguess = np.nan
            else:
                # in kpc
                rvguess = np.sqrt(2.0 / self.delta) / (0.001 * hz) * \
                    vvguess / 100.0
                # in Msun
                mvguess = 1.e11 * rvguess * (vvguess / 100.0)**2 / grav

            logger.info(
                f"Guess for rv is: {rvguess} kpc, mvguess is: {mvguess:.5e} Msun")

            # LCDM concentration
            auth = 'MDvdB'
            if self.delta == 200:
                auth = 'MDvdB200'

            if c200 is None:
                # Ensure positive input for cofmvir
                if mvguess / 1.e12 <= 0:
                    logger.warning(
                        f"Non-positive input ({mvguess / 1.e12}) to cofmvir. Setting c to NaN.")
                    c = np.nan
                else:
                    c = self.cofmvir(mvguess / 1.e12, auth)
            mdelta = mvguess
            logger.info(
                f"Mass (mdelta): {mdelta:.5e}, Concentration (c): {c}, Virial Radius (rvguess): {rvguess}")

            if conc is not None:
                c = conc
                logger.info(f"Adopting the input concentration: {c}")

            # check for convergence
            if ipass > 0 and not np.isnan(rvguess) and not np.isnan(rvguessold) and rvguessold != 0 and abs(rvguess / rvguessold - 1) < 0.001:
                logger.info('Converged!')
                # Exit loop
                break

            if ipass > self.itmax:
                logger.info('Max number of iterations reached')
                # Exit loop
                break

            # Define the median velocity if required
            nd1 = nd
            r1 = np.copy(cleanr)
            vz1 = np.copy(cleanvz)

            vzs = np.sort(vz1)

            if vmed is None:
                vzmed = self.even_odd_median(vzs, nd1)

            # MBM interloper rejection
            # Prepare inputs for siglosnfwml1 and massnfw
            siglos_x = c * r1 / rvguess
            massnfw_x = 1.0 / c
            massnfw_a = 1.0 / c
            massnfw_rvir = 1.0

            # Calculate components of vzmax
            siglos_term = self.siglosnfwml1(siglos_x)
            massnfw_term = self.massnfw(
                massnfw_x, massnfw_a, massnfw_rvir)

            # Check for NaN/Inf
            if np.any(np.isinf(siglos_term)) or np.any(np.isnan(siglos_term)):
                logger.warning(
                    f"siglosnfwml1 returned Inf/NaN. siglos_x_input: {siglos_x}. Result: {siglos_term}. vzmax will be problematic.")
            if np.any(np.isinf(massnfw_term)) or np.any(np.isnan(massnfw_term)):
                logger.warning(
                    f"massnfw returned Inf/NaN. massnfw_inputs: {massnfw_x}, {massnfw_a}, {massnfw_rvir}. Result: {massnfw_term}. vzmax will be problematic.")

            # Ensure argument to sqrt is >=0
            sqrt_term = np.sqrt(np.maximum(massnfw_term * c, 0))

            # Vzmax Calculation
            vzmax = nscut * siglos_term * sqrt_term * vvguess

            # Store for plotting
            vzmax = np.copy(vzmax)

            logger.info(f"vzmax first 5 = {vzmax[:5]}, vzmed = {vzmed}")

            #
            # Check if a is a quantity an
            if isinstance(vzmax, uu.Quantity):
                vzmax = vzmax.value
            else:
                vzmax = vzmax

            if isinstance(r1, uu.Quantity):
                r1 = r1.value
            else:
                r1 = r1

            if isinstance(r, uu.Quantity):
                r = r.value
            else:
                r = r

            vzmaxo = np.interp(r, r1, vzmax)
            wmbm = np.where(np.abs(vrf - vzmed) < vzmaxo)[0]

            # logger.info(f" wmbm : {wmbm}")
            cleanvz = vrf[wmbm]
            cleanr = r[wmbm]

            # logger.info(f" cleanvz : {cleanvz}")
            # logger.info(f" cleanr : {cleanr}")
            ipass += 1

            # Re-evaluate imembers for the full original dataset based on current cleanr/cleanvz
            imembers_current_iter = np.zeros(len(vrf), dtype=int)
            for i in range(len(vrf)):
                is_member = np.any((cleanvz == vrf[i]) & (cleanr == r[i]))
                if is_member:
                    imembers_current_iter[i] = 1
            # Update the imembers output array for this iteration
            imembers = imembers_current_iter

            # If no_iteration keyword is set, break after first pass
            if noiter:
                logger.info(
                    "No iteration keyword set. Exiting after first pass.")
                break

        logger.info("--- interlopmbmiter finished ---")
        logger.info(f"Final imembers: {imembers}")
        logger.info(f"Final mdelta (Msun): {mdelta:.5e}")
        logger.info(f"Final Algorithm Flag: {flag}")

        return imembers
