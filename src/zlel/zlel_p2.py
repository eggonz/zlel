#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
.. module:: zlel_main.py0
    :synopsis:
        This module contains functions to read, manage and solve circuits.
 
.. moduleauthor:: Mikel Elorza (mikelelorza0327@gmail.com), Egoitz Gonzalez (egoitz.gonz@gmail.com)


"""

import numpy as np 
import sys
import matplotlib.pyplot as plt

if __name__ == "zlel.zlel_p2":
    import zlel.zlel_p1 as zl1
    import zlel.zlel_p3 as zl3
    import zlel.zlel_p4 as zl4
else:
    import zlel_p1 as zl1
    import zlel_p3 as zl3
    import zlel_p4 as zl4


def print_solution2(sol, b, n):
    """ This function prints the solution with format.
    
        Args:
            sol: np array with the solution of the Tableau equations 
            (e_1,..,e_n-1,v_1,..,v_b,i_1,..i_b)
            b: # of branches
            n: # of nodes
    """
    np.set_printoptions(sign=' ')
    print("\n========== Nodes voltage to reference ========")
    for i in range(1, n):
        # print("e" + str(i) + " = ", '['+','.join(['% .8f' % num for num in sol[i-1]])+']')
        print("e" + str(i) + " = ", sol[i-1])
    print("\n========== Branches voltage difference ========")
    for i in range(1, b+1):
        # print("v" + str(i) + " = ", '['+','.join(['% .8f' % num for num in sol[i+n-2]])+']')
        print("v" + str(i) + " = ", sol[i+n-2])
    print("\n=============== Branches currents ==============")
    for i in range(1, b+1):
        # print("i" + str(i) + " = ", '['+','.join(['% .8f' % num for num in sol[i+b+n-2]])+']')
        print("i" + str(i) + " = ", sol[i+b+n-2])
        
    print("\n================= End solution =================\n")


def print_solution(sol, b, n):
    """ This function prints the solution with format.

        Args:
            sol: np array with the solution of the Tableau equations
            (e_1,..,e_n-1,v_1,..,v_b,i_1,..i_b)
            b: # of branches
            n: # of nodes

    """
    # The instructor solution needs to be a numpy array of numpy arrays
    # of float. If it is not, convert it to this format.
    if sol.dtype == np.float64:
        np.set_printoptions(sign=' ')  # Only from numpy 1.14
        tmp = np.zeros([np.size(sol), 1], dtype=float)
        for ind in range(np.size(sol)):
            tmp[ind] = np.array(sol[ind])
        sol = tmp
    print("\n========== Nodes voltage to reference ========")
    for i in range(1, n):
        print("e" + str(i) + " = ", sol[i - 1])
    print("\n========== Branches voltage difference ========")
    for i in range(1, b + 1):
        print("v" + str(i) + " = ", sol[i + n - 2])
    print("\n=============== Branches currents ==============")
    for i in range(1, b + 1):
        print("i" + str(i) + " = ", sol[i + b + n - 2])

    print("\n================= End solution =================\n")
    

def build_csv_header(tvi, b, n):
    """ This function build the csv header for the output files.
        First column will be v or i if .dc analysis or t if .tr and it will
        be given by argument tvi.
        The header will be this form,
        t/v/i,e_1,..,e_n-1,v_1,..,v_b,i_1,..i_b
        
    Args:
        tvi: "v" or "i" if .dc analysis or "t" if .tran
        b: # of branches
        n: # of nodes
    
    Returns:
        header: The header in csv format as string
    """
    header = tvi
    for i in range(1, n):
        header += ",e" + str(i)
    for i in range(1, b+1):
        header += ",v" + str(i)
    for i in range(1, b+1):
        header += ",i" + str(i)
    return header


def plot_from_cvs(filename, x, y, title):
    """ This function plots the values corresponding to the x string of the 
        file filename in the x-axis and the ones corresponding to the y 
        string in the y-axis.
        The x and y strings must mach with some value of the header in the
        csv file filename.
        
    Args:
        filename: string with the name of the file (including the path).
        x: string with some value of the header of the file.
        y: string with some value of the header of the file.
    
    """
    data = np.genfromtxt(filename, delimiter=',', skip_header=0,
                         skip_footer=1, names=True)
    
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot(data[x], data[y], color='r', label=title)   
    ax1.set_xlabel(x)
    ax1.set_ylabel(y)
    plt.show()


def command_pr(info):
    """
        Function called with command .PR.
        It prints the circuit info.

    Args:
        info: dict containing all circuit info
    """
    elem_num, br, br_nd, nd = info["elem_num"], info["br"], info["br_nd"], info["nd"]
    zl1.print_cir_info(br, br_nd, len(br), len(nd), nd, elem_num)
    zl1.print_incidence_matrix(info["incidence_mat"])


def command_op(info):
    """
        Function called with command .OP.
        It solves the given circuit via Tableau equations.

    Args:
        info: dict containing all circuit info
    """
    t = -1

    sol = solve_circuit_in_time(info, t)
    print_solution(sol, len(info["br"]), len(info["nd"]))


def build_tableau_system(a, m, n, u):
    """
        This function takes a, m and n matrices and u vector and uses them to
        build the linear system that corresponds to Tableau equations.

    Args:
        a: reduced incidence matrix.
        m: voltage coefficient matrix
        n: current coefficient matrix
        u: independent terms vector

    Returns:
        T: Tableau matrix.
        U: Tableau inhomogeneous vector.
    """

    a0, a1 = np.size(a, 0), np.size(a, 1)
    m0, m1 = np.size(m, 0), np.size(m, 1)
    tableau_t = np.zeros([a0 + a1 + m0, a0 + m1 + a1], dtype=float)
    tableau_u = np.zeros((a0 + a1 + m0, 1), dtype=float)

    tableau_t[:a0, -a1:] = a
    tableau_t[a0:-m0, :a0] = -np.transpose(a)
    tableau_t[a0:-m0, a0:-a1] = np.eye(m1)
    tableau_t[-m0:, a0:-a1] = m
    tableau_t[-m0:, -a1:] = n
    tableau_u[-m0:, :] = np.reshape(u, (-1, 1))

    return tableau_t, tableau_u


def solve_circuit_in_time(info, t):
    """
        Solves circuit in given time.
        -1 means time independent, amplitudes of B and Y will be used instead.

    Args:
        info: dict containing all circuit info
        t: time in seconds

    Returns:
        sol: np.array of size 2b+(n-1), solution for all circuit variables (e, v, i)
    """
    
    if not zl3.is_linear(info):
        return zl3.solve_nl_circuit_in_time(info, t)
        
    a = get_reduced_incidence_matrix(info)
    m, n, u = get_element_matrices(info, t)
    tableau_t, tableau_u = build_tableau_system(a, m, n, u)
    if np.linalg.det(tableau_t) == 0:
        sys.exit("Error solving Tableau equations, check if det(T) != 0.")

    return np.linalg.solve(tableau_t, tableau_u)


def command_dc(info, values, control):
    """
        Function called with command .DC.
        Saves dc sweep data into csv file with same filename.
        
    Args:
        info: dict containing all circuit info
        values: value info for command, np.array size(3)
        control: control element identifier, np.array size(1)
    """
    start, end, step = values
    t = -1

    file_name = info["file_name"][:-4] + "_" + control + ".dc"  # TODO lower
    header = build_csv_header("V", len(info["br"]), len(info["nd"]))

    with open(file_name, 'w') as file:
        print(header, file=file)

        v = start
        while v < end:
            change_value_of_elem(info, control, v)
            sol = solve_circuit_in_time(info, t)

            # write in csv
            sol = np.insert(sol, 0, v)
            # sol to csv
            sol_csv = ','.join(['%.5f' % num for num in sol])
            print(sol_csv, file=file)

            v += step


def change_value_of_elem(info, elem_name, new_value):
    """
        Changes the value of the selected element.
        Element must consist of a single branch and have a single value.

    Args:
         info: dict containing all circuit info
         elem_name: element name
         new_value: new value, which will override the previous
    """
    br_lower = np.array(list(map(lambda x: x.lower(), info["br"])))
    ind = np.flatnonzero(br_lower == elem_name.lower())[0]
    info["br_val"][ind, 0] = new_value


def command_tr(info, values):
    """
        Function called with command .TR.
        Saves transient analysis data into csv file.
        
    Args:
        info: dict containing all circuit info
        values: value info for command, np.array size(3)
    """
    start, end, step = values
   
    file_name = info["file_name"][:-4] + ".tr"
    header = build_csv_header("t", len(info["br"]), len(info["nd"]))

    dynamic = zl4.is_dynamic(info)

    with open(file_name, 'w') as file:
        print(header, file=file)

        if dynamic:
            zl4.initialize(info, step)

        t = start
        sol = solve_circuit_in_time(info, t)

        # write in csv
        sol_csv = np.insert(sol, 0, t)
        # sol to csv
        sol_csv = ','.join(['%.9f' % num for num in sol_csv])
        print(sol_csv, file=file)

        t += step
        while t < end:
            if dynamic:
                zl4.update_state(info, sol)

            sol = solve_circuit_in_time(info, t)
            
            # write in csv
            sol_csv = np.insert(sol, 0, t)
            # sol to csv
            sol_csv = ','.join(['%.9f' % num for num in sol_csv])
            print(sol_csv, file=file)
            t += step


def process_circuit(filename):
    """
        This function parses and processes circuit info.
        It must be called before doing anything else, as it defines all variables.

    Args:
        filename: string with file relative path

    Returns:
        info: dictionary containing all circuit info, contains following keys:
            elem_num: # of elements in circuit
            com_el: list of command names
            com_val: list of values for each command
            com_ctr: list of control elements for each command
            br: list of branch names
            br_nd: list of nodes for each branch
            br_val: list of values for the element in each branch
            br_ctr: list of control element for the element in each branch
            nd: list of node names
            incidence_mat: incidence matrix (nd x br)
    """

    # Parse the circuit
    cir_el, cir_nd, cir_val, cir_ctr = zl1.cir_parser(filename)
    # cir_el = np.array(list(map(lambda x: x.lower(), cir_el)))
    # cir_ctr = np.array(list(map(lambda x: x.lower(), cir_ctr)))

    # Differentiate commands from elements in file
    n = sum(1 for elem in cir_el if elem.startswith("."))

    com_el, com_val, com_ctr = \
        cir_el[-n:], cir_val[-n:], cir_ctr[-n:]

    cir_el, cir_nd, cir_val, cir_ctr = \
        cir_el[:-n], cir_nd[:-n], cir_val[:-n], cir_ctr[:-n]

    # Identify branches and nodes
    br, br_nd, br_val, br_ctr = zl1.span_branches(cir_el, cir_nd, cir_val, cir_ctr)
    nd = zl1.node_set(cir_nd)

    incidence_mat = zl1.build_incidence_matrix(br, br_nd, nd)

    zl1.check_parallel_v(br, br_val, incidence_mat)
    zl1.check_serial_i(nd, br, br_val, incidence_mat)

    return {
        "file_name": filename,
        "elem_num": len(cir_el),
        "com_el": com_el,
        "com_val": com_val,
        "com_ctr": com_ctr,
        "br": br,
        "br_nd": br_nd,
        "br_val": br_val,
        "br_ctr": br_ctr,
        "nd": nd,
        "incidence_mat": incidence_mat
    }


def get_element_matrices(info, t):
    """
        Given the element equations with the form " M v + N i = u ", this function builds
        and returns matrices M and N and vector u.
        t = -1 means time independent, amplitudes of B and Y will be used instead.
        Diode and transistor branches (d and q) are empty, and must be filled afterwards.

    Args:
        info: dictionary containing info of the current circuit
        t: current time, for time dependent values/elements

    Returns:
        M: voltage coefficient matrix
        N: current coefficient matrix
        u: independent terms vector
    """

    br, br_nd, br_val, br_ctr, nd =\
        info["br"], info["br_nd"], info["br_val"], info["br_ctr"], info["nd"]

    m = np.zeros(shape=(len(br), len(br)), dtype=float)
    n = np.zeros(shape=(len(br), len(br)), dtype=float)
    u = np.zeros(shape=len(br), dtype=float)

    br_lower = np.array(list(map(lambda x: x.lower(), br)))

    for ind in range(len(br)):
        branch = br_lower[ind]

        if branch.startswith("v"):
            m[ind, ind] = 1
            u[ind] = br_val[ind, 0]

        elif branch.startswith("i"):
            n[ind, ind] = 1
            u[ind] = br_val[ind, 0]

        elif branch.startswith("r"):
            m[ind, ind] = 1
            n[ind, ind] = -br_val[ind, 0]

        elif branch.startswith("d"):
            elem = branch
            n[ind, ind] = 1
            m[ind, ind] = -zl3.get_d_par(elem)[0]
            u[ind] = zl3.get_d_par(elem)[1]

        elif branch.startswith("q"):
            elem = branch[:-3]
            n[ind, ind] = 1
            if branch.endswith("_be"):
                m[ind, ind] = -zl3.get_q_par(elem)[0]
                m[ind, ind+1] = -zl3.get_q_par(elem)[1]
                u[ind] = zl3.get_q_par(elem)[4]
            elif branch.endswith("_bc"):
                m[ind, ind] = -zl3.get_q_par(elem)[3]
                m[ind, ind-1] = -zl3.get_q_par(elem)[2]
                u[ind] = zl3.get_q_par(elem)[5]

        elif branch.startswith("a"):
            # write an equation with each branch
            ind_in = np.flatnonzero(br_lower == branch[:-2] + "in")[0]  # index of the control branch in br

            if branch.endswith("_in"):
                n[ind, ind_in] = 1
            elif branch.endswith("_ou"):
                m[ind, ind_in] = 1

        elif branch.startswith("e"):
            m[ind, ind] = 1
            ind_ctr = np.flatnonzero(br_lower == br_ctr[ind].lower())[0]  # index of the control branch in br
            m[ind, ind_ctr] = -br_val[ind, 0]

        elif branch.startswith("g"):
            n[ind, ind] = 1
            ind_ctr = np.flatnonzero(br_lower == br_ctr[ind].lower())[0]  # index of the control branch in br
            m[ind, ind_ctr] = -br_val[ind, 0]

        elif branch.startswith("h"):
            m[ind, ind] = 1
            ind_ctr = np.flatnonzero(br_lower == br_ctr[ind].lower())[0]  # index of the control branch in br
            n[ind, ind_ctr] = -br_val[ind, 0]

        elif branch.startswith("f"):
            n[ind, ind] = 1
            ind_ctr = np.flatnonzero(br_lower == br_ctr[ind].lower())[0]  # index of the control branch in br
            n[ind, ind_ctr] = -br_val[ind, 0]

        elif branch.startswith("b"):
            m[ind, ind] = 1
            v, f, p = br_val[ind, :]
            if t == -1:
                u[ind] = v
            else:
                u[ind] = v * np.sin(2*np.pi*f*t + 2*np.pi/360*p)

        elif branch.startswith("y"):
            m[ind, ind] = 1
            i, f, p = br_val[ind, :]
            if t == -1:
                u[ind] = i
            else:
                u[ind] = i * np.sin(2*np.pi*f*t + 2*np.pi/360*p)
                
        elif branch.startswith("c"):
            m[ind, ind] = 1
            n[ind, ind] = -zl4.tau / br_val[ind][0]
            u[ind] = zl4.v_aurreko[branch]

        if branch.startswith("l"):
            m[ind, ind] = -zl4.tau / br_val[ind][0]
            n[ind, ind] = 1
            u[ind] = zl4.i_aurreko[branch]

    return m, n, u


def get_reduced_incidence_matrix(info):
    """
        Takes whole incidence matrix and returns reduced version.
        If present, it deletes row corresponding to node 0.
        If not, it deletes the first row instead.

    Args:
        info: dict containing all circuit info.

    Returns:
        A: reduced incidence matrix of the circuit.
    """
    nd, mat = info["nd"], info["incidence_mat"]
    a = np.flatnonzero(nd == 0)
    if len(a) == 0:
        return mat[1:, :]
    else:
        ind = a[0]
        return np.delete(mat, ind, 0)


def run_commands(info):
    """
        Runs all commands of the circuit.

    Args:
        info: dict containing all circuit data.
    """
    for ind in range(len(info["com_el"])):
        com = info["com_el"][ind][1:].lower()
        if com == "pr":
            command_pr(info)
        elif com == "op":
            command_op(info)
        elif com == "dc":
            command_dc(info,
                       np.ndarray.tolist(info["com_val"][ind, :]),
                       info["com_ctr"][ind])
        elif com == "tr":
            command_tr(info,
                       np.ndarray.tolist(info["com_val"][ind, :]))


"""    
https://stackoverflow.com/questions/419163/what-does-if-name-main-do
https://stackoverflow.com/questions/19747371/
python-exit-commands-why-so-many-and-when-should-each-be-used
"""
if __name__ == "__main__":
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        filename = "../cirs/all/1_zlel_anpli.cir"

    cir_info = process_circuit(filename)
    
    run_commands(cir_info)
