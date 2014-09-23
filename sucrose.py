"""
Code for precisely calculating the density of sucrose solutions.
"""

import numpy, scipy.optimize
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import sys

def get_constants():
    """Returns fit constants for calculating the density of a sucrose solution.
    
    The method returns constants (c, d, e) for the equation
    
        density = c*m**2 + d*m + e
        
    where density is in g/mL, and m is the mass ratio
    
        m = (mass of sucrose) / (mass of water)
    
    Args:
        currently none.
        
    Returns:
        A list of the three constants (c, d, e) defined above.
    """
    molalities = numpy.array([0.015, 0.030, 0.060, 0.090, 0.122, 0.154, 0.186, 0.220, 0.254, 0.289, 0.325, 0.398, 0.476, 0.556, 0.641, 0.730, 0.824, 0.923, 1.026, 1.136, 1.252, 1.375, 1.505, 1.643, 1.791, 1.948, 2.116, 2.295, 2.489, 2.697, 2.921, 4.382, 6.817, 11.686])
    densities = numpy.array([1.0002, 1.0021, 1.0060, 1.0099, 1.0139, 1.0178, 1.0218, 1.0259, 1.0299, 1.0340, 1.0381, 1.0465, 1.0549, 1.0635, 1.0722, 1.0810, 1.0899, 1.0990, 1.1082, 1.1175, 1.1270, 1.1366, 1.1464, 1.1562, 1.1663, 1.1765, 1.1868, 1.1972, 1.2079, 1.2186, 1.2295, 1.2864, 1.3472, 1.4117])
    mass_sucrose_over_mass_water = molalities * 342.2965 / 1000.0

    def error_function((c, d, e)):
        sum = 0.0
        for m, j in zip(mass_sucrose_over_mass_water, range(len(mass_sucrose_over_mass_water))):
        	sum += ((c*m**2 + d*m + e) - densities[j])**2
        return sum


    c, d, e = scipy.optimize.fmin(error_function, (3, 4, 1), full_output=1, xtol=1e-29, ftol=1e-29, maxiter=10000, maxfun=10000)[0]

    y = []
    ztable = []
    zfunc = []
    count = 0
    max_diff = 0.0

    for m, j in zip(mass_sucrose_over_mass_water, range(len(mass_sucrose_over_mass_water))):
        y.append(m)
        d_table = densities.ravel()[count]
        d_func = c*m**2 + d*m + e
        ztable.append(d_table)
        zfunc.append(d_func)
        if abs(d_table - d_func) > max_diff:
            max_diff = abs(d_table - d_func)
        count += 1

    print "Density = %0.4g r^2 + %0.4g r + %0.4f" % (c, d, e)
    print "Max diff = %0.4g" % max_diff

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(y, ztable, "ro")
    ax.plot(y, zfunc, "bo")
    ax.set_xlabel("Mass sucrose / Mass H$_2$O")
    ax.set_ylabel("Density (g/mL)")
    plt.show()

    return c, d, e



def main():
    """
    """
    if len(sys.argv) == 1:
        H2O_mass = float(raw_input("Enter H2O mass: "))
        sucrose_mass = float(raw_input("Enter sucrose mass: "))
    elif len(sys.argv) == 4:
        H2O_mass = float(sys.argv[1])
        sucrose_mass = float(sys.argv[2])
    else:
        print "USAGE:    python NaCl.py mass_of_H2O mass_of_sucrose "
        print "EXAMPLE:  python NaCl.py 10.416 1.013"
        exit()
    m = sucrose_mass / H2O_mass
    c, d, e = get_constants()
    print c*m**2 + d*m + e

if __name__ == "__main__":
    main()


    
    