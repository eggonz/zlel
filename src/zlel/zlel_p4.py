#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

.. module:: zlel_p3.py
    :synopsis: Put yours
 
.. moduleauthor:: Put yours


"""
import time
import sys
import zlel_p2 as zl2

v_aurreko = {}
i_aurreko = {}
step = 1e-4


def update_state(info, sol):
    n = len(info["nd"])
    for ind in range(len(info["br"])):
        elem = info["br"][ind]
        if elem.startswith("c"):
            v_aurreko[elem] = sol[n+ind-1]
            print(v_aurreko)
        if elem.startswith("l"):
            i_aurreko[elem] = sol[n+len(info["br"])+ind-1]


def initialize(info, tau):
    # global step
    step = tau
    br = info["br"]
    br_val = info["br_val"]
    for ind in range(len(br)):
        branch = br[ind].lower()
        if branch.startswith("c"):
            v_aurreko[branch] = br_val[ind][1]
            print(v_aurreko)
        if branch.startswith("l"):
            i_aurreko[branch] = br_val[ind][1]


def is_not_dynamic(info):
    for br in info["br"]:
        if br.lower().startswith("c") or br.lower().startswith("l"):
            return False
    return True

    
if __name__ == "__main__":
    start = time.perf_counter()
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        filename = "../cirs/all/3_zlel_RC.cir"

    cir_info = zl2.process_circuit(filename)
    zl2.run_commands(cir_info)
    
    end = time.perf_counter()
    print("Elapsed time: ")
    print(end - start)  # Time in seconds
