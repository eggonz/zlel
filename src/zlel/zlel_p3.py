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



def diode_NR (I0,nD,Vdj):
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
        Vd: Value of Vd.
        
    Return:
        gd: Conductance of the NR discrete equivalent for the diode.
        Id: Current independent source of the NR discrete equivalent.
    
    """
    
    Vt = 8.6173324e-5*300*nD
    
    return gd, Id


    
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
    
    
    
    





