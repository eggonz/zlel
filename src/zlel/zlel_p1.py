#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
.. module:: zlel_main.py
    :synopsis:

.. moduleauthor:: YOUR NAME AND E-MAIL


"""

import numpy as np
import sys


def cir_parser(filename):
    """
        This function takes a .cir test circuit and parse it into
        4 matices.
        If the file has not the proper dimensions it warns and exit.

    Args:
        filename: string with the name of the file

    Returns:
        cir_el: np array of strings with the elements to parse. size(1,b)
        cir_nd: np array with the nodes to the circuit. size(b,4)
        cir_val: np array with the values of the elements. size(b,3)
        cir_ctrl: np array of strings with the element which branch
        controls the controlled sources. size(1,b)

    Rises:
        SystemExit

    """
    try:
        cir = np.array(np.loadtxt(filename, dtype=str))
    except ValueError:
        sys.exit("File corrupted: .cir size is incorrect.")

    # numpy usefull exmaples
    print("================ cir ==========")
    print(cir)
    print("\n======== a = np.array (cir[:,0], dtype = int) ==========")
    a = np.array(cir[:, 1], dtype=int)
    print(a)
    print("\n======== a = np.append(a,300) ==========")
    a = np.append(a, 300)
    print(a)
    print("\n======== b = a[a > 3] ==========")
    b = a[a > 3]
    print(b)
    print("\n======== c = np.unique(a) ==========")
    c = np.unique(a)
    print(c)
    print("\n======== d = np.flatnonzero(a != 0) ==========")
    d = np.flatnonzero(a != 0)
    print(d)
    print("\n======== e = np.flatnonzero(a == 0) ==========")
    e = np.flatnonzero(a == 0)
    print(e)

#    THIS FUNCTION IS NOT COMPLETE


def print_cir_info(cir_el, cir_nd, b, n, nodes, el_num):
    """ Prints the info of the circuit:
            1.- Elements info
            2.- Node info
            3.- Branch info
            4.- Variable info
    Args:
        cir_el: reshaped cir_el
        cir_nd: reshaped cir_nd. Now it will be a (b,2) matrix
        b: # of branches
        n: # number of nodes
        nodes: an array with the circuit nodes sorted
        el_num:  the # of elements.

    """
    # Element info
    print(str(el_num) + ' Elements')
    # Node info
    print(str(n) + ' Different nodes: ' +
          str(nodes))
    # Branch info
    print("\n" + str(b) + " Branches: ")

    for i in range(1, b+1):
        print("\t" + str(i) + ". branch:\t" + cir_el[i-1] +
              ",\ti" + str(i) +
              ",\tv" + str(i) +
              "=e" + str(cir_nd[i-1, 0]) +
              "-e" + str(cir_nd[i-1, 1]))

    # Variable info
    print("\n" + str(2*b + (n-1)) + " variables: ")
    # Print all the nodes but the first (0 because is sorted)
    for i in nodes[1:]:
        print("e"+str(i)+", ", end="", flush=True)
    for i in range(b):
        print("i"+str(i+1)+", ", end="", flush=True)
    # Print all the branches but the last to close it properly
    # It works because the minuimum amount of branches in a circuit must be 2.
    for i in range(b-1):
        print("v"+str(i+1)+", ", end="", flush=True)
    print("v"+str(b))

#    YOU ARE ENCOURAGED TO MODIFY THE ARGS OF THIS FUNCTION TO FIT INTO YOUR
#    CODE.


"""
https://stackoverflow.com/questions/419163/what-does-if-name-main-do
https://stackoverflow.com/questions/19747371/
python-exit-commands-why-so-many-and-when-should-each-be-used
"""
if __name__ == "__main__":
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        filename = "../cirs/examples/0_zlel_V_R_Q.cir"
    # Parse the circuit
    # [cir_el,cir_nd,cir_val,cir_ctr]=cir_parser(filename)
    cir_parser(filename)
