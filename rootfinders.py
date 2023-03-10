
"""
MATH2019 CW1 rootfinders module

@author: Sahil Singh
"""
import numpy as np
import matplotlib . pyplot as plt
import rootfinders as rf

class InputError(Exception):
    def __init__(self):
        pass
    def __str__(self):
        return "The differential of the function f must not be 0"
    
def fixedpoint_iteration (g , p0 , Nmax ) :
    """
    In this iteration the initual value of p0 is set by the user and then
    input into the function g where we get the new p1. The iteration is
    repeated Nmax number of times.
    Parameters
    ----------
    g : Function
    The function input by the user, which they wish to
    find at which point p, g(p) = 0 
    p0 : Integer / int
    The initual starting point where the algorithem will 
    begin to start its convergence 
    Nmax : Integer / int
    Number of iterations which will be carried out in the 
    algorithem.
    Returns
    -------
    p_array : numpy.ndarray / int array
    The array will carry the squences of the aproximation
    at which g(p)=0
    """
    p_list = []
    pn = p0
    if p0 <= 0:
        raise ValueError("The input digit must be a positive Integer")
    if Nmax <= 0 or type(Nmax)!=int:
        raise ValueError("The input for the maximum iteration needs to be posive Integer")
    for n in np.arange(0, Nmax):
        pn_plusone = g(pn)
        p_list.append(pn_plusone)
        pn = pn_plusone
    p_array = np.array(p_list)
    return p_array
#4
def fixedpoint_iteration_stop(g,p0,Nmax,TOL):
    """
    fixedpoint iteration stop returning a numpy array 
    sequence approximations which are given by the fix point iteration,
    however with this function we have a contraint of tolerance which
    allows us to make a guess as to what is a reasonable good approximate 
    root of the function g .

    Parameters:
    ----------
    g : Function
    The function input by the user, which they wish to
    find at which point p, g(p) = 0 .
    p0 : Integer / int
    The initual starting point where the algorithem will 
    begin to start its convergence  
    Nmax : Integer / int
    DESCRIPTION .
    TOL : float
    Number of iterations which will be carried out in the 
    algorithem.
    Returns
    -------
    p_array : numpy.ndarray / int array
    The array will carry the squences of the aproximation
    at which g(p)=0, with the given tolerance level .
    """
    p_list = []
    pn = p0
    if p0 <= 0:
        raise ValueError("The input digit must be a positive Integer")
    if Nmax <= 0 or type(Nmax)!=int:
        raise ValueError("The input for the maximum iteration needs to be posive Integer")
    for n in np.arange(1, Nmax):
        pn_plusone = g(pn)
        if abs(pn_plusone - pn) <= TOL:
            break 
        p_list.append(pn_plusone)
        pn = pn_plusone
    p_array = np.array(p_list)
    return p_array
#5
def newton_stop (f , dfdx , p0 , Nmax , TOL ) :
    """
    DESCRIPTION:
    Newton's Method returns a numpy array 
    sequence approximations which are given by the fix point iteration,
    however with this function we have a contraint of tolerance which
    allows us to make a guess as to what is a reasonable good approximate 
    root of the function f using dfdx which is its derivative to 
    help with a faster convergence .
    Parameters:
    ----------
    f : Function
    The function input by the user, which they wish to
    find at which point p, f(p) = 0 .
    dfdx : Function
    The function input is the first derivative of 
    the function f where it does not equal 0 .
    p0 : Integer / int
    The initual starting point where the algorithem will 
    begin to start its convergence 
    Nmax : Integer / int
    DESCRIPTION .
    TOL : float
    DESCRIPTION .
    Number of iterations which will be carried out in the 
    algorithem.
    -------
    p_array : numpy.ndarray / int array
     The array will carry the squences of the aproximation
    at which f(p)=0 .
    """
    pn = p0
    p_list = []
    if dfdx == 0:
        raise InputError()
    if p0 <= 0:
        raise ValueError("The input digit must be a positive Integer")
    if Nmax <= 0 or type(Nmax)!=int:
        raise ValueError("The input for the maximum iteration needs to be posive Integer")
    for n in np.arange(0, Nmax):
        pn_plusone = pn - f(pn)/dfdx(pn)
        p_list.append(pn_plusone)
        if abs(pn_plusone - pn) <= TOL:
            break 
        pn = pn_plusone
    p_array = np.array(p_list)
    return p_array
#6
def plot_convergence (p ,f , dfdx ,g , p0 , Nmax ):
    TOL = 10**(-20)
    # Fixed - point iteration
    p_array = fixedpoint_iteration(g,p0,Nmax)
    e_array = np.abs(p-p_array)
    n_array = 1 + np.arange(np.shape(p_array)[0])
    # Newton method
    p_newton_array = newton_stop(f,dfdx,p0,Nmax,TOL) 
    e_newton_array = np.abs(p-p_newton_array)
    n_newton_array =1+np.arange(np.shape(p_newton_array)[0])
    
    # Preparing figure , using Object - Oriented (OO) style ; see:
        # https :// matplotlib .org/ stable / tutorials / introductory / quick_start . html
    fig , ax = plt . subplots ()
    ax . set_yscale ('log')
    ax . set_xlabel ("n")
    ax . set_ylabel ("|p- p_n |")
    ax . set_title (" Convergence behaviour ")
    ax . grid ( True )
    # Plot
    ax.plot( n_array , e_array , "o", label ="FP iteration ", linestyle ='--')
    ax.plot( n_newton_array ,e_newton_array , "*", label ="Newton's Method ", linestyle = ':')
    # Add legend
    ax.legend(bbox_to_anchor = (1.395, 1))
    return fig , ax
#7
def optimize_FPmethod (f , c_array , p0 , TOL ) :
    n_iteration = []
    for i in range(len(c_array)):
        c = c_array[i,]
        g = lambda x: x - c*f(x)
        b = len(fixedpoint_iteration_stop(g,p0,100,TOL))
        n_iteration.append(b)
    n_opt = min(n_iteration)
    index = n_iteration.index(n_opt) 
    c_opt = c_array[index,]
    return c_opt,n_opt


