import numpy as np
from scipy import interpolate
import sys


def Gam_get(model: str, datW="../Data/Width_tot/"):
    '''
    Get Width Data for different models and write interpolation
    functions corresponding to the partial widths as a tuple.

    Parameters
    ----------
    model : str
        Model for which widths are needed.
        Available models:
        B_boson
        B-L_boson
        dark_photon
        dp (dark photon with different scaling)
        Le-Lmu_boson
        Le-Ltau_boson
        Lmu-Ltau_boson
    datW : str
        Optional, Location of the width data files (default is
        '../Data/Width_tot/') The file corresponding to the width data
        sould have a structure as (Mass, Full Width, Partial to el,
        Partial to mu, partial to tau, partial to hadrons) and should be
        named 'width_<model>.csv' where model is the same name as enters
        as an input to this file. For width calculation, the coupling
        should be one. If there is no activity in a particular channel,
        it should be full of zeros and not missing.

    Returns
    -------
    List (of interpolation functions)
        List of interpolation functions (in the same order as in the
        data file). The last column gives the function for the invisible
        partial width (calculated as total - sum of visible partials).

    Raises
    ------
    FileNotFoundError
        If data file speciffied does not exist

    IndexError
        If data file does not have all necessary columns
    '''

    try:
        data = np.genfromtxt(datW+'width_'+model+".csv", skip_header=1).T
    except FileNotFoundError:
        print("Decay Data file not found in the location specified. "
              "Refer to the docs to figure out correct placement of data file."
              )
        sys.exit(1)

    try:
        Ifuncs = []
        for ii in range(1, 7):
            Ifuncs.append(interpolate.interp1d(data[0], data[ii]))
        inv_width = data[1] - np.sum(data[2:7, :], axis=0)
        Ifuncs.append(interpolate.interp1d(data[0], inv_width))
    except IndexError:
        print("Decay data file does not contain all the information "
              "necessary. Check docs to find out all the columns "
              "needed."
              )
        sys.exit(1)

    return Ifuncs
