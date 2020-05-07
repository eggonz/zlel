#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

.. module:: zlel_p3.py
    :synopsis:
        This module solves more complex circuits containing some non-linear elements
        such as diodes or transistors using numerical methods.
 
.. moduleauthor:: Mikel Elorza (mikelelorza0327@gmail.com), Egoitz Gonzalez (egoitz.gonz@gmail.com)


"""
import time

import numpy as np
import sys

import zlel_p2 as zl2

d_parameters = dict()
q_parameters = dict()


def diode_nr(I0, nD, Vdj):
    """ https://documentation.help/Sphinx/math.html
        Calculates the g and the I of a diode for a NR discrete equivalent.
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


def update_diode_par(info, Vdj, br_name):
    """
        This function updates the values for the NR parameters of the given diode.

    Args:
        info: dict containing all circuit info.
        Vdj: current value for Vd.
        br_name: name of one of the element's branch, whose parameters will be updated.
    """
    ind = np.where(info["br"] == br_name)[0][0]  # getting index from np.array
    nD = info["br_val"][ind][1]
    I0 = info["br_val"][ind][0]
    gd, Id = diode_nr(I0, nD, Vdj)
    d_parameters[br_name.lower()] = [gd, Id]


def get_d_par(elem_name):
    """
        This function returns the list of parameters of the given diode.

    Args:
        elem_name: diode element name.

    Returns:
        list of diode NR parameters.
    """
    return d_parameters[elem_name.lower()]


def transistor_nr(Ies, Ics, Vbej, Vbcj, alphaR, alphaF, n=1):
    """
    This function returns g matrix values for the transistor equations for the next step,
    ass well as Ebers-Moll currents for inhomogeneous term.

    Args:
        Ies: Value of Ies.
        Ics: Value of Ics.
        Vbej: Value of Vbe.
        Vbcj: Value of Vbc.
        alphaR: Value of alphaR.
        alphaF: Value of alphaF.
        nVt: Value of n*Vt.

    Return:
        g: size 4 tuple containing transistors equations variation coefficients
            g11, g12, g21, g22.
    """

    nVt = n * 8.6173324e-5 * 300

    g11 = - Ies / nVt * np.exp(Vbej / nVt)
    g22 = - Ics / nVt * np.exp(Vbcj / nVt)
    g12 = - alphaR * g22
    g21 = - alphaF * g11

    return g11, g12, g21, g22


def ebers_moll_currents(Ies, Ics, Vbej, Vbcj, alphaR, alphaF, n=1):
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
        ebmoIe: current Ebers-Moll Ie current
        ebmoIc: current Ebers-Moll Ic current
    """

    nVt = n * 8.6173324e-5 * 300

    ebmoIe = Ies * (np.exp(Vbej / nVt) - 1) - alphaR * Ics * (np.exp(Vbcj / nVt) - 1)
    ebmoIc = - alphaF * Ies * (np.exp(Vbej / nVt) - 1) + Ics * (np.exp(Vbcj / nVt) - 1)

    return ebmoIe, ebmoIc


def update_transistor_par(info, Vbej, Vbcj, br_name):
    """
        This function updates the values for the NR parameters of the given transistor.

    Args:
        info: dict containing all circuit info.
        Vdj: current value for Vd.
        br_name: name of one of the element's branch, whose parameters will be updated.
    """
    ind = np.where(info["br"] == br_name)[0][0]  # getting index from np.array
    Ies = info["br_val"][ind][0]
    Ics = info["br_val"][ind][1]
    beta = info["br_val"][ind][2]
    alphaF = beta / (beta + 1)
    alphaR = Ies / Ics * alphaF

    g11, g12, g21, g22 = transistor_nr(Ies, Ics, Vbej, Vbcj, alphaR, alphaF, n=1)
    ebmoIe, ebmoIc = ebers_moll_currents(Ies, Ics, Vbej, Vbcj, alphaR, alphaF, n=1)
    q_parameters[br_name[:-3].lower()] = [g11, g12, g21, g22, ebmoIe, ebmoIc]


def get_q_par(elem_name):
    """
        This function returns the list of parameters of the given transistor.

    Args:
        elem_name: transistor element name.

    Returns:
        list of transistor NR parameters.
    """
    return q_parameters[elem_name.lower()]


def update_all_nl_par(info, tableau_sol):
    """
        This function updates all the parameters of every diode and transistor.

    Args:
        info: dict containing all circuit info.
        tableau_sol: np.array containing the last solution of Tableau equations' system.
    """
    n = len(info["nd"])
    for ind in range(len(info["br"])):
        br = info["br"][ind]
        if br.startswith("d"):
            vd = tableau_sol[n + ind]
            update_diode_par(info, vd, br)
        elif br.startswith("q") and br.endswith("_be"): # update once per element
            vbe = tableau_sol[n + ind]
            vbc = tableau_sol[n + ind + 1]
            update_transistor_par(info, vbe, vbc, br)


def set_initial_nl_par(info):
    """
        This function sets all the parameters of every diode and transistor
        to their respective initial values.

        Initial values are:
        Diode: Vd0 = 0.6
        Transistor: Vbe0 = 0.6 Vbc0 = 0.6

    Args:
        info: dict containing all circuit info.
    """

    # Initial values for NR iterations
    vd0 = vbe0 = vbc0 = 0.6

    for ind in range(len(info["br"])):
        br = info["br"][ind]
        if br.startswith("d"):
            update_diode_par(info, vd0, br)
        elif br.startswith("q") and br.endswith("_be"):  # update once per element
            update_transistor_par(info, vbe0, vbc0, br)


def is_linear(info):
    """
        This function checks the presence of diode and transistors
        (non-linear elements) in the circuit.

    Args:
        info: dict containing all circuit info.

    Returns:
        bool: True if circuit has no diode or transistors.
    """
    for br in info["br"]:
        if br.lower().startswith("d") or br.lower().startswith("q"):
            return False
    return True


def solve_nl_circuit_in_time(info, t):
    """
        Solves non-linear circuit in given time.
        -1 means time independent, amplitudes of B and Y will be used instead.

    Args:
        info: dict containing all circuit info.
        t: time in seconds.

    Returns:
        sol: np.array of size 2b+(n-1), solution for all circuit variables (e, v, i).
    """
    return nr_method(info, t=t)


def nr_method(info, precision=1e-5, t=-1):
    """
        Applies Newton-Raphson method to solve non-linear circuit
        in given time with given convergence precision.
        -1 means time independent, amplitudes of B and Y will be used instead.
        Default precision: 1e-5

    Args:
        info: dict containing all circuit info.
        t: time in seconds.

    Returns:
        sol: np.array of size(1,2b+(n-1)), solution for all circuit variables (e, v, i).
    """
    a = zl2.get_reduced_incidence_matrix(info)

    set_initial_nl_par(info)
    m, n, u = zl2.get_element_matrices(info, t)
    tableau_t, tableau_u = zl2.build_tableau_system(a, m, n, u)
    sol = np.linalg.solve(tableau_t, tableau_u)

    pr = False
    while not pr:
        # print(sol)
        prev_sol = sol
        update_all_nl_par(info, prev_sol)
        m, n, u = zl2.get_element_matrices(info, t)
        tableau_t, tableau_u = zl2.build_tableau_system(a, m, n, u)
        sol = np.linalg.solve(tableau_t, tableau_u)
        pr = check_precision(sol, prev_sol, precision)

    return sol


def check_precision(sol, prev_sol, precision):
    """
        This function compares 2 consecutive solutions of Tableau equation system
        to check convergence precision. It return whether precision is correct.

    Args:
        sol: current solution of Tableau equation system.
        prev_sol: previous solution of Tableau equation system.
        precision: requested convergence precision.

    Returns:
        bool: True if solution is accurate and convergent.
    """
    for ind in range(len(sol)):
        # print(sol-prev_sol)
        if abs(sol[ind] - prev_sol[ind]) > precision:
            return False
    return True


"""   
https://stackoverflow.com/questions/419163/what-does-if-name-main-do
"""
if __name__ == "__main__":
    start = time.perf_counter()
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        filename = "../cirs/all/2_zlel_1D.cir"

    cir_info = zl2.process_circuit(filename)
    zl2.run_commands(cir_info)
    
    end = time.perf_counter()
    print("Elapsed time: ")
    print(end - start)  # Time in seconds
    
    
    
    





