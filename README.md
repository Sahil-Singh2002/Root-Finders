# Linear-System
Scientific computation coursework
  
The Fix-point iteration Part 1, from here on the code is of my own. This function considers the fixed-point problem g(p) = p. To numerically find a fixed point, we consider the fixed-point iteration method (p0 given). The function is mainly writen in the method of p_n = g(p_n-1) where n belongs to a set of natural numbers. The purpose of this function is to cause us to find a point of convergence after n number of iterations which in our function is Nmax, such that g(x_n) = x_n as n aproched infinity.

    def fixedpoint_iteration (g , p0 , Nmax ) :
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
      
  The output of our function an numpy.ndarray called p_array, its dimension is (Nmax, ). This is provides an approximations p_n where n belongs to a set of naterual numbers computed by the fixed-point iteration.
  
The Fix-point iteration Part 2, from here on the code is of my own. This function considers the fixed-point problem g(p) = p with the stopping criterion. To numerically find a fixed point, for which the absolute difference between  p_k and g(p_k) is less then our tolerence. We consider the fixed-point iteration method (p0 given). The function is mainly writen in the method of p_n = g(p_n-1) where n belongs to a set of natural numbers, once the error between is below the tolerance the iteration would stop before reaching Nmax. The purpose of this function is to a value p_k which is close enough without having the function conduct a lenghty iteration using up an emense amount of run time.

    def fixedpoint_iteration_stop(g,p0,Nmax,TOL):
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
    
  The output of our function an numpy.ndarray called p_array, its dimension is (k, ). This is provides an approximations p_n where n belongs to a set of naterual numbers up to k computed by the fixed-point iteration with the stopping criterion. The value of k is the smallest integer such that stopping criterion holds, unless Nmax iteration have been performed, in that case k = Nmax. 
  
The Newton's Method, we now consider the rootfinding problem f(p) = 0. To numerically find the root, we consider Newton's Method (p0 given).
Where p_(n+1) = p_n - f(p_n)/f'(p_n) where n = 0,1,2,... . The stopping criterion causes a premature stop in our iteration |p_(k+1) - p_k| is <= TOL.

    def newton_stop (f , dfdx , p0 , Nmax , TOL ) :
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
    
  The output of our function an numpy.ndarray called p_array, its dimension is (k, ). This is provides an approximations p_n where n belongs to a set of naterual numbers up to k computed by Newton's Method with the stopping criterion. The value of k is the smallest integer such that stopping criterion holds, unless Nmax iteration have been performed, in that case k = Nmax.
  
The Convergence behaviour: Plotting the error. We start the the problem with the assumption that we know what exaclty p is. We plot the absolute value of the error, 
e_n = | p - p_n |, at each iteration n = 1,2,3,... . This function will just plot the error using matplotlib.pyplot. We see on of the plots the convergence of both iterative methods such as Fixed-point iteration and Newton's Method. In this we see clearly that Newton's method has a faster convergence rate in a low number of iterations, then Fix-point iteration does. Both starting from the same p0. 

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
    #Set Axies and Title
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
    
  As we see, the function creates a figure with a plot of the errors for the fixed-point iteration method, using the function fixedpoint_iteration. Also ploting, in the same axies, the errors for Newton's method. Computing the Newton-method approximations using your newton_stop, based on the input for (f,fddx,p0,Nmax) and setting TOL as very small value such as TOL = 1.0e-16. Then we we ensure the graphs can be distinguished from one another.
  
Optimizing the Fix-point iteration Method. In this function we return the fixed-point iteration method with stopping crition and aim to optimise the parameter c within the functuion g where g(x) = x - cf(x). Our starting off peramters are (f,c_array,p0,TOL) for optimize_FPmethod. 
    
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
    
   The function should return a real number/ float c_opt and the corresponding integer n_opt, which corresponds to the optional value c in the following way: the optimal c - c_opt is the element inside of our ndarray c_array, for which the fastest convergence takes place given the Tolerance value TOL, using the initual value p0. So | p_k - g(p_k) | reaches the TOL in the least number of iterations n_opt. This is so that we can a reduced run time for the machine, for a fastest system.
    
  
  
  
  
  
  
  
  
  
  
  
  
  
