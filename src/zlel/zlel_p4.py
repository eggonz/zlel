#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

.. module:: zlel_p4.py
    :synopsis:
        This module adds a method to solve dynamic systems, by time discretization.
 
.. moduleauthor:: Mikel Elorza (mikelelorza0327@gmail.com), Egoitz Gonzalez (egoitz.gonz@gmail.com)


"""
import time
import sys

if __name__ == "zlel.zlel_p4":
    import zlel.zlel_p2 as zl2
else:
    import zlel_p2 as zl2

v_aurreko = {}
i_aurreko = {}
tau = None


def update_state(info, sol):
    """ This function updates relevant values of the previous solution to 
    execute the Euler method.

        Args:
            sol: np array with the solution of the Tableau equations
            (e_1,..,e_n-1,v_1,..,v_b,i_1,..i_b)
            info: dict containing all circuit info

    """
    n = len(info["nd"])
    for ind in range(len(info["br"])):
        elem = info["br"][ind].lower()
        if elem.startswith("c"):
            v_aurreko[elem] = sol[n+ind-1]
        if elem.startswith("l"):
            i_aurreko[elem] = sol[n+len(info["br"])+ind-1]


def initialize(info, step):
    """ This function initializes relevant values of to execute Euler method.

        Args:
            step: time step for the discrretization
            info: dict containing all circuit info
    
    """
    global tau
    tau = step
    br = info["br"]
    br_val = info["br_val"]
    for ind in range(len(br)):
        branch = br[ind].lower()
        if branch.startswith("c"):
            v_aurreko[branch] = br_val[ind, 1]
        if branch.startswith("l"):
            i_aurreko[branch] = br_val[ind, 1]


def is_dynamic(info):
    """ This function checks whether the circuit is dynamic or not.

        Args:
            info: dict containing all circuit info
            
        Returns:
            bool: True if circuit has capacitors or inductors.
    
    """
    for br in info["br"]:
        if br.lower().startswith("c") or br.lower().startswith("l"):
            return True
    return False

    
if __name__ == "__main__":
    start = time.perf_counter()
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        filename = "../cirs/all/3_zlel_RC.cir"

    cir_info = zl2.process_circuit(filename)
    zl2.run_commands(cir_info)
    
    end = time.perf_counter()
    # print("Elapsed time: ")
    # print(end - start)  # Time in seconds
