# Linear-System

This file consists of codes designed to help find a root of a function using different methods such as fixed-point iteration and Newton's method. It also includes code to display the difference in convergence of two systems and optimize the fixed-point iteration.

## Fixed-Point Iteration

### Part 1: Fixed-Point Problem
Consider the fixed-point problem $g(p) = p$. To numerically find a fixed point, we use the fixed-point iteration method with an initial value $p_0$. The iteration is given by $p_n = g(p_{n-1})$ for $n$ belonging to a set of natural numbers. The purpose is to find a point of convergence after $N_{\text{max}}$ iterations, such that $g(x_n) = x_n$ as $n$ approaches infinity. The output is a numpy.ndarray called `p_array` with dimensions $(N_{\text{max}},)$, providing the approximations $p_n$.

### Part 2: Fixed-Point Problem with Stopping Criterion
Similar to Part 1, we consider the fixed-point problem $g(p) = p$ with an additional stopping criterion. The iteration stops when the absolute difference between $p_k$ and $g(p_k)$ is less than the tolerance. This avoids lengthy iterations. The output is a numpy.ndarray called `p_array` with dimensions $(k,)$, providing the approximations $p_n$ up to $k$. The value of $k$ is the smallest integer for which the stopping criterion holds, unless $N_{\text{max}}$ iterations have been performed, in which case $k = N_{\text{max}}$.

## Newton's Method

We now consider the rootfinding problem $f(p) = 0$. To numerically find the root, we use Newton's Method with an initial value $p_0$. The iteration is given by $p_{n+1} = p_n - \frac{f(p_n)}{f'(p_n)}$. The stopping criterion is when $|p_{k+1} - p_k| \leq \text{TOL}$. The output is a numpy.ndarray called `p_array` with dimensions $(k,)$, providing the approximations $p_n$ up to $k$. The value of $k$ is the smallest integer for which the stopping criterion holds, unless $N_{\text{max}}$ iterations have been performed, in which case $k = N_{\text{max}}$.

## Convergence Behavior: Plotting the Error

We plot the absolute value of the error, $e_n = | p - p_n |$, at each iteration. This function uses matplotlib.pyplot. The resulting plot shows the convergence of both fixed-point iteration and Newton's Method. We observe that Newton's method has a faster convergence rate in fewer iterations compared to fixed-point iteration, both starting from the same $p_0$. The function creates a figure with a plot of the errors for the fixed-point iteration method, using the function `fixedpoint_iteration`. It also plots the errors for Newton's method, computed using the `newton_stop` function with inputs $(f, f''(x), p_0, N_{\text{max}})$ and setting $\text{TOL}$ to a very small value such as $\text{TOL} = 1.0 \times 10^{-16}$. This ensures that the graphs can be distinguished from each other.

## Optimizing the Fixed-Point Iteration Method

This function returns the fixed-point iteration method with a stopping criterion and aims to optimize the parameter $c$ within the function $g(x) = x - cf(x)$. The inputs are $(f, c_{\text{array}}, p_0, \text{TOL})$ for `optimize_FPmethod`. The function returns a real number $c_{\text{opt}}$ and the corresponding integer $n_{\text{opt}}$. The optimal $c - c_{\text{opt}}$ is the element inside the ndarray $c_{\text{array}}$, for which the fastest convergence takes place given the tolerance value $\text{TOL}$ and the initial value $p_0$. Thus, $|p_k - g(p_k)|$ reaches $\text{TOL}$ in the least number of iterations $n_{\text{opt}}$, resulting in reduced run time.

In conclusion, these codes implement convergence algorithms learned in the module "Introduction to Scientific Computation." They aim to find solutions to nonlinear systems as fast as possible. As observed, Newton's method produces the fastest convergence.
