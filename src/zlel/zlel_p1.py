#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
.. module:: zlel_p1.py
    :synopsis:
        This module parses circuit info, and separates it into different data arrays.
        It also calculates incidence matrix for the given input circuit.

.. moduleauthor:: Mikel Elorza (mikelelorza0327@gmail.com), Egoitz Gonzalez (egoitz.gonz@gmail.com)


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
        number_of_branches: # of branches of the element

    """
    elem_type = element[0].lower()

    if elem_type == "q":
        br = np.array([element + "_be", element + "_bc"], dtype=str)
        nd = np.array([nodes[1:3], nodes[1::-1]], dtype=int)
        number_of_branches = 2

    elif elem_type == "a":
        br = np.array([element + "_in", element + "_ou"], dtype=str)
        nd = np.array([nodes[0:2], nodes[2:4]], dtype=int)
        number_of_branches = 2

    elif elem_type == "j":
        br = np.array([element + "_in", element + "_ou"], dtype=str)
        nd = np.array([nodes[0:2], nodes[2:4]], dtype=int)
        number_of_branches = 2

    else:
        br = np.array([element], dtype=str)
        nd = np.array([nodes[0:2]], dtype=int)
        number_of_branches = 1

    return br, nd, number_of_branches


def span_branches(cir_el, cir_nd, cir_val, cir_ctr):
    """
    Takes parsed cir matrices and expands the list elements to get a list of all the branches.

    Elements connected to two or more nodes are replaced by different branches linking strictly two nodes each.

    Args:
        cir_el: parsed cir_el
        cir_nd: parsed cir_nd, size(b,4)
        cir_val: parsed cir_val, size(b,3)
        cir_ctr: parsed cir_ctr, size(b,1)

    Returns:
        branches: reshaped cir_el
        branch_nd: reshaped cir_nd, now it is a (b,2) matrix
        branch_val: reshaped cir_val, contains values repeated in branches corresponding to same element, size(b,3)
        branch_ctr: reshaped cir_ctr, contains controls repeated in branches corresponding to same element, size(b,1)

    Raises:
        SystemExit

    """
    branches = np.empty((1, 0), dtype=str)
    branch_nd = np.empty((0, 2), dtype=int)
    branch_val = np.empty((0, 3), dtype=float)
    branch_ctr = np.empty((0, 1), dtype=str)

    for i in range(len(cir_el)):
        span_elem = element_branches(cir_el[i], cir_nd[i])
        branches = np.append(branches, span_elem[0])
        branch_nd = np.append(branch_nd, span_elem[1], axis=0)
        for _ in range(span_elem[2]):
            branch_val = np.append(branch_val, cir_val[i:i+1, :], axis=0)
            branch_ctr = np.append(branch_ctr, cir_ctr[i])

    if len(np.flatnonzero(branch_nd == 0)) == 0:
        sys.exit("Reference node \"0\" is not defined in the circuit.")
    elif len(np.flatnonzero(branch_nd == 0)) == 1:
        sys.exit("Node 0 is floating.")

    return branches, branch_nd, branch_val, branch_ctr


def check_parallel_v(branches, branches_val, incidence_matrix):
    """
    This methods checks whether any two different voltage sources are connected in parallel.

    Args:
        nd: list of ordered nodes
        branches: reshaped cir_el, list of branch names
        branches_val: reshaped cir_val, list of values of each branch's element, size(b,3) np array
        incidence_matrix: calculated incidence matrix

    Raises:
        SystemExit

    """

    v_sources = np.array([i for i in range(len(incidence_matrix[0, :])) if
                          branches[i][0].lower() == "v" or branches[i][0].lower() == "b"], dtype=int)
    for v in v_sources:
        value1 = branches_val[v, 0]
        for u in v_sources[v_sources > v]:
            value2 = branches_val[u, 0]
            if np.all(incidence_matrix[:, v] == incidence_matrix[:, u]) and value1 != value2:
                sys.exit("Parallel V sources at branches " + str(v) + " and " + str(u) + ".")
            elif np.all(incidence_matrix[:, v] == - incidence_matrix[:, u]) and value1 != - value2:
                sys.exit("Parallel V sources at branches " + str(v) + " and " + str(u) + ".")


def check_serial_i(nd, branches, branches_val, incidence_matrix):
    """
    This methods checks whether any two different current sources are connected in series.

    Args:
        nd: list of ordered nodes
        branches: reshaped cir_el, list of branch names
        branches_val: reshaped cir_val, list of values of each branch's element, size(b,3) np array
        incidence_matrix: calculated incidence matrix

    Raises:
        SystemExit

    """

    i_sources = np.array([i for i in range(len(incidence_matrix[0, :])) if
                          branches[i][0].lower() == "i" or branches[i][0].lower() == "y"], dtype=int)
    for n in range(len(incidence_matrix[:, 0])):
        non_zero = np.flatnonzero(incidence_matrix[n, :] != 0)
        if all(j in i_sources for j in non_zero):
            i_sum = 0
            for elem in non_zero:
                i_sum += branches_val[elem, 0] * incidence_matrix[n][elem]
            if i_sum != 0:
                sys.exit("I sources in series at node " + str(nd[n]) + ".")


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
        filename = "../cirs/all/0_zlel_OPAMP.cir"

    # Parse the circuit
    cir_el, cir_nd, cir_val, cir_ctr = cir_parser(filename)

    # Identify branches and nodes
    br, br_nd, br_val, br_ctr = span_branches(cir_el, cir_nd, cir_val, cir_ctr)
    nd = node_set(cir_nd)

    mat = build_incidence_matrix(br, br_nd, nd)

    check_parallel_v(br, br_val, mat)
    check_serial_i(nd, br, br_val, mat)

    # Print info
    print_cir_info(br, br_nd, len(br), len(nd), nd, len(cir_el))
    print_incidence_matrix(mat)
