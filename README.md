# linear-system
Scientific computation coursework

The bisection method, the code is note of my own as this is code is given to us by our lecturer to start of the next codes proceding forth. However i will do justice to this code by giving a breef explanation to this method. This method in my opinion is inspired by looking the mean value theroem where if we have an intervil of [a,b] and we find that f(a) > 0  and f(b) <0 or vise versa. Then there must exist c which is in (a,b) for which f(c) = 0. 

  The purpouse of this method is to find the root of our function by repeatedly dividing the intervil. In our function we are given the permeters of (f,a,b,Nmax).
    
    The function f is the input which the user will put in to find the root for a given interval.
    The Float a and b are two distinct values, a != b, which we are testing the convergence in. 
    The integer Nmax the maximum number of iterations which our function will run in a for loop.
   
  The output of our function an numpy.ndarray called p_array, its dimension is (Nmax, ). This is provides an approximations p_n where n belongs to a 
  set of naterual numbers computed by the bisector method.
  
The Fix-point iteration Part 1, from here on the code is of my own. This function considers the fixed-point problem g(p) = p. To numerically find a fixed point, we consider the fixed-point iteration method (p0 given). The function is mainly writen in the method of p_n = g(p_n-1) where n belongs to a set of natural numbers. The purpose of this function is to cause us to find a point of convergence after n number of iterations which in our function is Nmax, such that g(x_n) = x_n as n aproched infinity.
  
  The purpouse of this method is to find a point x_n of our function by repeatedly iterating the value of x_n-1 into our function g. 
  In our function we are given the permeters of (g,p0,Nmax).
  
    The function g is the input which the user will put in to find the convergence of our fixpoint.
    The Float p0 is the point in our domain, which represents the initual approximation to start the fix-point iteration. 
    The integer Nmax the maximum number of iterations which our function will run in a for loop.
      
  The output of our function an numpy.ndarray called p_array, its dimension is (Nmax, ). This is provides an approximations p_n where n belongs to a 
  set of naterual numbers computed by the fixed-point iteration.
  
The Fix-point iteration Part 2, from here on the code is of my own. This function considers the fixed-point problem g(p) = p with the stopping criterion. To numerically find a fixed point, for which the absolute difference between  p_k and g(p_k) is less then our tolerence. We consider the fixed-point iteration method (p0 given). The function is mainly writen in the method of p_n = g(p_n-1) where n belongs to a set of natural numbers, once the error between is below the tolerance the iteration would stop before reaching Nmax. The purpose of this function is to a value p_k which is close enough without having the function conduct a lenghty iteration using up an emense amount of run time.

  The purpouse of this method is to find a point x_n of our function by repeatedly iterating the value of x_n-1 into our function g.Once |p_k - g(p_k)| <= TOL. 
  The iteration stops. In our function we are given the permeters of (g,p0,Nmax,TOL).
  
    The function g is the input which the user will put in to find the convergence of our fixpoint.
    The Float p0 is the point in our domain, which represents the initual approximation to start the fix-point iteration. 
    The integer Nmax the maximum number of iterations which our function will run in a for loop.
    The Float TOL is the tolerance which the function will allow the difference between our p_k and g(p_k).
    
  The output of our function an numpy.ndarray called p_array, its dimension is (k, ). This is provides an approximations p_n where n belongs to a 
  set of naterual numbers up to k computed by the fixed-point iteration with the stopping criterion. The value of k is the smallest integer such that stopping
  criterion holds, unless Nmax iteration have been performed, in that case k = Nmax.
  
The Newton's Method, we now consider the rootfinding problem f(p) = 0. To numerically find the root, we consider Newton's Method (p0 given). Where
p_(n+1) = p_n - f(p_n)/f'(p_n) where n = 0,1,2,... . The stopping criterion causes a premature stop in our iteration |p_(k+1) - p_k| is <= TOL.

    The function f is the input function which the user will put in to find the convergence of our root.
    The function fdiff is the input differentiated function of f which the user will put in to find the convergence of our root.
    The Float p0 is the point in our domain, which represents the initual approximation to start the Newton's Method. 
    The integer Nmax the maximum number of iterations which our function will run in a for loop.
    The Float TOL is the tolerance which the function will allow the difference between our p_(k+1) and p_k.
    
  The output of our function an numpy.ndarray called p_array, its dimension is (k, ). This is provides an approximations p_n where n belongs to a 
  set of naterual numbers up to k computed by Newton's Method with the stopping criterion. The value of k is the smallest integer such that stopping
  criterion holds, unless Nmax iteration have been performed, in that case k = Nmax.
  
The Convergence behaviour: Plotting the error. We start the the problem with the assumption that we know what exaclty p is. We plot the absolute value of the error, 
e_n = | p - p_n |, at each iteration n = 1,2,3,... . This function will just plot the error using matplotlib.pyplot. We see on of the plots the convergence of both iterative methods such as Fixed-point iteration and Newton's Method. In this we see clearly that Newton's method has a faster convergence rate in a low number of iterations, then Fix-point iteration does. Both starting from the same p0. 

    The function f is the input function which the user will put in to find the convergence of our root.
    The function dfdx is the input differentiated function of f which the user will put in to find the convergence of our root.
    The function g is the input which the user will put in to find the convergence of our fixpoint.
    The Float p0 is the point in our domain, which represents the initual approximation to start the Newton's Method. 
    The integer Nmax the maximum number of iterations which our function will run in a for loop.
    The Float p is the value that we assume is the actual root or point of convergence. Which we are testing for error with.
    
  As we see, the function creates a figure with a plot of the errors for the fixed-point iteration method, using the function fixedpoint_iteration. Also ploting,
  in the same axies, the errors for Newton's method. Computing the Newton-method approximations using your newton_stop, based on the input for (f,fddx,p0,Nmax)
  and setting TOL as very small value such as TOL = 1.0e-16. Then we we ensure the graphs can be distinguished from one another.
  
 Optimizing the Fix-point iteration Method. In this function we return the fixed-point iteration method with stopping crition and aim to optimise the parameter c
 within the functuion g where g(x) = x - cf(x). Our starting off peramters are (f,c_array,p0,TOL) for optimize_FPmethod. 
    
    The function f is the input function which the user will put in to find the convergence of our root.
    The np.ndarray c_array is a linspace of of interval with a set amount of intervals and these points are are what we test for the fastest convergence. 
    The Float p0 is the point in our domain, which represents the initual approximation to start the Newton's Method. 
    The Float TOL is the tolerance which the function will allow the difference between our p_(k+1) and p_k.
    
   The function should return a real number/ float c_opt and the corresponding integer n_opt, which corresponds to the optional value c in the following way:
   the optimal c - c_opt is the element inside of our ndarray c_array, for which the fastest convergence takes place given the Tolerance value TOL, using the 
   initual value p0. So | p_k - g(p_k) | reaches the TOL in the least number of iterations n_opt. This is so that we can a reduced run time for the machine,
   for a fastest system.
    
  
  
  
  
  
  
  
  
  
  
  
  
  
