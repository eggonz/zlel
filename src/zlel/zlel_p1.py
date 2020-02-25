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
    This function takes a .cir test circuit and parse it into 4 matrices.

    If the file has not the proper dimensions it warns and exits.

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

    cir_el = np.array(cir[:, 0], dtype=str)
    cir_nd = np.array(cir[:, 1:5], dtype=int)
    cir_val = np.array(cir[:, 5:8], dtype=float)
    cir_ctrl = np.array(cir[:, 8], dtype=str)

    return cir_el, cir_nd, cir_val, cir_ctrl


def element_branches(element, nodes):
    """
    Separates elements connected to two or more nodes in different branches, depending on the element type.

    Args:
        element: element name
        nodes: integer np.array with nodes of the element, size(1,4)

    Returns:
        br: branches in which the element is separated
        nd: nodes which element is connected to

    """
    elem_type = element[0].lower()

    if elem_type == "q":
        br = np.array([element + "_be", element + "_bc"], dtype=str)
        nd = np.array([nodes[1:3], nodes[1::-1]], dtype=int)

    elif elem_type == "a":
        br = np.array([element + "_in", element + "_ou"], dtype=str)
        nd = np.array([nodes[0:2], nodes[2:4]], dtype=int)

    else:
        br = np.array([element], dtype=str)
        nd = np.array([nodes[0:2]], dtype=int)

    return br, nd


def span_branches(cir_el, cir_nd):
    """
    Takes parsed cir matrices and expands the list elements to get a list of all the branches.

    Elements connected to two or more nodes are replaced by different branches linking strictly two nodes each.

    Args:
        cir_el: parsed cir_el
        cir_nd: parsed cir_nd, size(b,4)

    Returns:
        branches: reshaped cir_el
        branch_nd: reshaped cir_nd, now it is a (b,2) matrix

    Raises:
        SystemExit

    """
    branches = np.empty((1, 0), dtype=str)
    branch_nd = np.empty((0, 2), dtype=int)

    for i in range(len(cir_el)):
        span_elem = element_branches(cir_el[i], cir_nd[i])
        branches = np.append(branches, span_elem[0])
        branch_nd = np.append(branch_nd, span_elem[1], axis=0)

    if len(np.flatnonzero(branch_nd == 0)) == 0:
        sys.exit("Reference node \"0\" is not defined in the circuit.")
    elif len(np.flatnonzero(branch_nd == 0)) == 1:
        sys.exit("Node 0 is floating.")

    return branches, branch_nd


def check_parallel_v(cir_el, cir_val, branches, incidence_matrix):
    """
    This methods checks whether any two different voltage sources are connected in parallel.

    Args:
        branches: reshaped cir_el, ordered list of branch names
        incidence_matrix: calculated incidence matrix

    Raises:
        SystemExit

    """

    v_sources = np.array([i for i in range(len(incidence_matrix[0, :])) if branches[i][0].lower() == "v"], dtype=int)
    for i in v_sources:
        value1 = cir_val[np.flatnonzero(cir_el == branches[i])[0]]
        for j in v_sources[v_sources > i]:
            value2 = cir_val[np.flatnonzero(cir_el == branches[j])[0]]
            if np.all(incidence_matrix[:, i] == incidence_matrix[:, j]) and value1 != value2:
                sys.exit("Parallel V sources at branches " + str(i) + " and " + str(j) + ".")
            elif np.all(incidence_matrix[:, i] == - incidence_matrix[:, j]) and value1 != - value2:
                sys.exit("Parallel V sources at branches " + str(i) + " and " + str(j) + ".")


def check_serial_i(cir_el, cir_val, branches, incidence_matrix):
    """
    This methods checks whether any two different current sources are connected in series.

    Args:
        branches: reshaped cir_el, list of branch names
        incidence_matrix: calculated incidence matrix

    Raises:
        SystemExit

    """
    i_sources = np.array([i for i in range(len(incidence_matrix[0, :])) if branches[i][0].lower() == "i"], dtype=int)
    for n in range(len(incidence_matrix[:, 0])):
        non_zero = np.flatnonzero(incidence_matrix[n, :] != 0)
        if len(non_zero) == 2 and (non_zero in i_sources):
            i, j = non_zero
            value1 = cir_val[np.flatnonzero(cir_el == branches[i])[0]]
            value2 = cir_val[np.flatnonzero(cir_el == branches[j])[0]]

            if incidence_matrix[n, i] == - incidence_matrix[n, j] and value1 != value2:
                sys.exit("I sources in series at node " + str(i) + ".")
            elif incidence_matrix[n, i] == incidence_matrix[n, j] and value1 != - value2:
                sys.exit("I sources in series at node " + str(i) + ".")

            sys.exit("I sources in series at node " + str(i) + ".")


def node_set(cir_nd):
    """
    This method gathers all nodes in a single sorted list.

    Args:
        cir_nd: parsed or reshaped cir_nd np.array, size(b,2) or size(b,4)

    Returns:
        nodes: sorted np.array of all nodes

    """
    return np.unique(cir_nd)


def build_incidence_matrix(branches, branch_nd, nodes):
    """
    Generates incidence matrix for provided branches, and their connection nodes.

    Args:
        branches: reshaped cir_el / list of branches, size(1,b) np.array
        branch_nd: reshaped cir_nd / list of branch nodes, size(b,2) np.array
        nodes: sorted np.array set of nodes, size(1,n) np.array

    Returns:
        incidence_mat: incidence matrix, size(n,b) np.array

    """
    incidence_mat = np.zeros((len(nodes), len(branches)), dtype=int)
    for br in range(len(branches)):
        row = np.flatnonzero(nodes == branch_nd[br][0])[0]
        incidence_mat[row][br] = 1
        row = np.flatnonzero(nodes == branch_nd[br][1])[0]
        incidence_mat[row][br] = -1
    return incidence_mat


def print_cir_info(cir_el, cir_nd, b, n, nodes, el_num):
    """
        Prints the info of the circuit:
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

    for i in range(1, b + 1):
        print("\t" + str(i) + ". branch:\t" + cir_el[i - 1] +
              ",\ti" + str(i) +
              ",\tv" + str(i) +
              "=e" + str(cir_nd[i - 1, 0]) +
              "-e" + str(cir_nd[i - 1, 1]))

    # Variable info
    print("\n" + str(2 * b + (n - 1)) + " variables: ")

    # Print all the nodes but the first (0 because is sorted)
    for i in nodes[1:]:
        print("e" + str(i) + ", ", end="", flush=True)
    for i in range(b):
        print("i" + str(i + 1) + ", ", end="", flush=True)

    # Print all the branches but the last to close it properly
    # It works because the minimum amount of branches in a circuit must be 2.
    for i in range(b - 1):
        print("v" + str(i + 1) + ", ", end="", flush=True)
    print("v" + str(b))


#    YOU ARE ENCOURAGED TO MODIFY THE ARGS OF THIS FUNCTION TO FIT INTO YOUR
#    CODE.


def print_incidence_matrix(mat):
    """
    Prints the incidence matrix provided.

    Args:
        mat: incidence matrix np.array

    """
    print("\nIncidence Matrix: ")
    print(mat)


"""
https://stackoverflow.com/questions/419163/what-does-if-name-main-do
https://stackoverflow.com/questions/19747371/
python-exit-commands-why-so-many-and-when-should-each-be-used
"""
if __name__ == "__main__":
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        filename = "../cirs/all/0_zlel_parallel_V_V.cir"

    # Parse the circuit
    [cir_el, cir_nd, cir_val, cir_ctr] = cir_parser(filename)
    el_num = len(cir_el)

    # Identify branches and nodes
    [br, br_nd] = span_branches(cir_el, cir_nd)
    nd = node_set(cir_nd)

    print(br_nd)

    mat = build_incidence_matrix(br, br_nd, nd)

    # check_parallel_v(br, mat)
    # check_serial_i(br, mat)

    # Print info
    print_cir_info(br, br_nd, len(br), len(nd), nd, el_num)
    print_incidence_matrix(mat)
