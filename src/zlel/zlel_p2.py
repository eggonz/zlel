#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
.. module:: zlel_main.py0
    :synopsis: 
 
.. moduleauthor:: YOUR NAME AND E-MAIL 


"""

import numpy as np 
import sys
import matplotlib.pyplot as plt

if __name__ == "__main__":
    import zlel_p1 as zl1
else:
    import zlel.zlel_p1 as zl1


def print_solution(sol, b, n):
    """ This function prints the solution with format.
    
        Args:
            sol: np array with the solution of the Tableau equations 
            (e_1,..,e_n-1,v_1,..,v_b,i_1,..i_b)
            b: # of branches
            n: # of nodes
    
    """
    print("\n========== Nodes voltage to reference ========")
    for i in range(1, n):
        print("e" + str(i) + " = ", sol[i-1])
    print("\n========== Branches voltage difference ========")
    for i in range(1, b+1):
        print("v" + str(i) + " = ", sol[i+n-2])
    print("\n=============== Branches currents ==============")
    for i in range(1, b+1):
        print("i" + str(i) + " = ", sol[i+b+n-2])
        
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


def save_as_csv(b, n, filename):
    """ This function generates a csv file with the name filename.
        First it will save a header and then, it loops and saves a line in 
        csv format into the file.
        
    Args:
        b: # of branches
        n: # of nodes
        filename: string with the filename (incluiding the path)
    """
    # Sup .tr
    header = build_csv_header("t", b, n)
    with open(filename, 'w') as file:
        print(header, file=file)
        # Get the indices of the elements corresponding to the sources.
        # The freq parameter cannot be 0 this is why we choose cir_tr[0].
        t = 0
        while t < 10:
            # for t in tr["start"],tr["end"],tr["step"]
            # Recalculate the Us for the sinusoidal sources

            sol = np.full(2*b+(n-1), t+1, dtype=float)
            # Insert the time
            sol = np.insert(sol, 0, t)
            # sol to csv
            sol_csv = ','.join(['%.5f' % num for num in sol])
            print(sol_csv, file=file)
            t = t + 1


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


def count_commands(cir_el):
    """
        This function counts how many commands does the given file have,
        supposing commands begin with '.'.
        
        Args:
            cir_el: np.array list of parsed elements, the first column of
            the file
            
        Returns:
            n: # of lines with command data
                
    """
    n = 0
    for elem in cir_el:
        if elem.startswith("."):
            n += 1
    return n


def command_pr():
    """
        Function called with command .PR.
        It prints the circuit info.
    """
    zl1.print_cir_info(br, br_nd, len(br), len(nd), nd, el_num)
    zl1.print_incidence_matrix(A)


def command_op():
    """
        Function called with command .OP.
        It solves the given circuit via Tableau equations.
    """
    a0, a1 = np.size(A, 0), np.size(A, 1)
    m0, m1 = np.size(M, 0), np.size(M, 1)
    T = np.zeros([a0 + a1 + m0, a0 + m1 + a1], dtype=float)
    U = np.zeros([a0 + a1 + m0, 1], dtype=float)

    T[:a0, -a1:] = A
    T[a0:-m0, :a0] = np.transpose(A)
    T[a0:-m0, a0:-a1] = np.eye(m1)
    T[-m0:, a0:-a1] = M
    T[-m0:, -a1:] = N
    U[-m0:, :] = u
    
    sol = np.linealg.solve(T, U)
    print_solution(sol, b, n)
    

def command_dc(values, control):
    """
        Function called with command .DC.
        Saves dc sweep data into csv file with same filename.
        
    Args:
        values: value info for command, np.array size(3)
        control: control element identifier, np.array size(1)
    """
    start, end, step = values
    
    header = build_csv_header("t", b, n)
    with open(filename, 'w') as file:
        print(header, file=file)
        
        t = start
        while t < end:
            # change u
            # write in csv
            t += step


def command_tr(values):
    """
        Function called with command .TR.
        Saves transient analysis data into csv file.
        
    Args:
        values: value info for command, np.array size(3)
    """


"""    
https://stackoverflow.com/questions/419163/what-does-if-name-main-do
https://stackoverflow.com/questions/19747371/
python-exit-commands-why-so-many-and-when-should-each-be-used
"""
if __name__ == "__main__":
    if len(sys.argv) > 1:
        filename = sys.argv[1]
    else:
        filename = "../cirs/all/1_zlel_V_R_dc.cir"
    
    b = 2
    n = 2
    filename = filename[:-3] + "tr"
    save_as_csv(b, n, filename)
    plot_from_cvs(filename, "t", "e1", "")

    '''
    #COPY FROM P1:

    # Parse the circuit
    cir_el, cir_nd, cir_val, cir_ctr = cir_parser(filename)
    el_num = len(cir_el)
    
    # Differentiate commands from elements in file
    n = count_commands(cir_el)
    
    com_el, com_nd, com_val, com_ctr = \
        cir_el[-n:], cir_nd[-n:], cir_val[-n:], cir_ctr[-n:]
    
    cir_el, cir_nd, cir_val, cir_ctr = \
        cir_el[:-n], cir_nd[:-n], cir_val[:-n], cir_ctr[:-n]

    # Identify branches and nodes
    br, br_nd, br_val = zl1.span_branches(cir_el, cir_nd, cir_val)
    nd = zl1.node_set(cir_nd)

    A = zl1.build_incidence_matrix(br, br_nd, nd)

    zl1.check_parallel_v(br, br_val, A)
    zl1.check_serial_i(nd, br, br_val, A)
    
    # define M, N and u'''