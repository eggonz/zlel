#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

.. module:: zlel_p3.py
    :synopsis: Put yours
 
.. moduleauthor:: Put yours


"""
import numpy as np
import sys

if __name__ == "__main__":
    import zlel_p1 as zl1
    import zlel_p2 as zl2
else:
    import zlel.zlel_p1 as zl1
    import zlel.zlel_p2 as zl2



def diode_NR (I0, nD, Vdj):
    """ https://documentation.help/Sphinx/math.html
        Calculates the g and the I of a diode for a NR discrete equivalent
        Given,
        
        :math:`Id = I_0(e^{(\\frac{V_d}{nV_T})}-1)`
        
        The NR discrete equivalent will be,
        
        :math:`i_{j+1} + g v_{j+1} = I`
        
        where,
        
        :math:`g = -\\frac{I_0}{nV_T}e^{(\\frac{V_d}{nV_T})}`
        
        and
        
        :math:`I = I_0(e^{(\\frac{V_{dj}}{nV_T})}-1) + gV_{dj}`
        
    Args:
        I0: Value of I0.
        nD: Value of nD.
        Vdj: Value of Vd.
        
    Return:
        gd: Conductance of the NR discrete equivalent for the diode.
        Id: Current independent source of the NR discrete equivalent.
    
    """
    
    Vt = 8.6173324e-5*300*nD

    gd = - I0 / Vdj * np.exp(Vdj / Vt)
    Id = I0 * (np.exp(Vdj / Vt) - 1) + gd * Vdj
    
    return gd, Id


def transistor_NR(Ies, Ics, Vbej, Vbcj, alphaR, alphaF, nVt):
    """
    This function returns g matrix values for the transistor equations for the next step,
    ass well as Ebers-Moll currents for inhomogeneous term.

    Args:
        Ies: Value of Ies
        Ics: Value of Ics
        Vbej: Value of Vbe
        Vbcj: Value of Vbc
        alphaR: Value of alphaR
        alphaF: Value of alphaF
        nVt: Value of n*Vt

    Return:
        g: np.array of size(2,2) containing transistors equations variation coefficients
        ebmo: current Ebers-Moll currents in a size(2) np.array
    """

    g = np.zeros((2, 2))

    '''
    g[0, 0] = 
    g[1, 1] = 
    '''
    # TODO

    g[0, 1] = - alphaR * g[1, 1]
    g[1, 0] = - alphaF * g[0, 0]

    ebmo = ebers_moll_currents(Ies, Ics, Vbej, Vbcj, alphaR, alphaF, nVt)
    return g, ebmo


def ebers_moll_currents(Ies, Ics, Vbej, Vbcj, alphaR, alphaF, nVt):
    """
    This function returns the current values of Ebers-Moll currents.

    Args:
        Ies: Value of Ies
        Ics: Value of Ics
        Vbej: Value of Vbe
        Vbcj: Value of Vbc
        alphaR: Value of alphaR
        alphaF: Value of alphaF
        nVt: Value of n*Vt

    Return:
        ebmo: current Ebers-Moll currents in a size(2) np.array
    """

    ebmo = np.zeros(2)

    '''
    ebmo[0] = 
    ebmo[1] =
    '''
    # TODO

    return ebmo


def check_non_linear(cir_info):
    """
    This function checks the presence of non-linear elements in the circuit.

    Args:
        cir_info: dict containing all circuit info.

    Returns:
        nl_elem: list of non-linear elements
    """

    # TODO

    nl_elem = list()

    return nl_elem


def calculate_NR():
    """
    """

    # check nonlinear

    while True:
        # update M, N, u
        # solve circuit
        # check circuit valid, then break

    return # solution


"""   
https://stackoverflow.com/questions/419163/what-does-if-name-main-do
"""
if __name__ == "__main__":
#    start = time.perf_counter()
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        filename = "../cirs/all/2_zlel_Q.cir"
        filename = "../cirs/all/2_zlel_1D.cir"
        
    
#    end = time.perf_counter()
#    print ("Elapsed time: ")
#    print(end - start) # Time in seconds
    
    
    
    





