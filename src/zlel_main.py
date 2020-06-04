#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
.. module:: zlel_main.py
    :synopsis:
        Given a circuit file, this module processes it and runs the specified commands.
 
.. moduleauthor:: Mikel Elorza (mikelelorza0327@gmail.com), Egoitz Gonzalez (egoitz.gonz@gmail.com)

"""

import time
import sys

import zlel.zlel_p2 as zl2

"""    
https://stackoverflow.com/questions/419163/what-does-if-name-main-do
https://stackoverflow.com/questions/19747371/
python-exit-commands-why-so-many-and-when-should-each-be-used
"""
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
