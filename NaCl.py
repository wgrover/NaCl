"""
Code for precisely calculating the density of sodum chloride solutions.
"""

import numpy, scipy.optimize
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import sys

"""hello world"""

def get_constants():
    """Returns fit constants for calculating the density of a sodium chloride solution.
    
    The method returns constants (a, b, c, d, e) for the equation
    
        density = a*t**2 + b*t + c*m**2 + d*m + e
        
    where density is in g/mL, t is temperature in C, and m is the mass ratio
    
        m = (mass of NaCl) / (mass of water)
    
    Args:
        currently none.
        
    Returns:
        A list of the five constants (a, b, c, d, e) defined above.
    """
    temperatures = [20.0, 25.0, 30.0, 40.0]
    molalities = numpy.array([0.100, 0.250, 0.500, 0.750, 1.000, 2.000, 3.000, 4.000, 5.000])

    specific_volumes = numpy.array([#[0.995732, 0.989259, 0.978889, 0.968991, 0.959525, 0.925426, 0.896292, 0.870996, 0.848646],  # 0.0 C
                                    #[0.995998, 0.989781, 0.979804, 0.970256, 0.961101, 0.927905, 0.899262, 0.874201, 0.851958],  # 10.0 C
                                    [0.997620, 0.991564, 0.981833, 0.972505, 0.963544, 0.930909, 0.902565, 0.877643, 0.855469],  # 20.0 C
                                    [0.998834, 0.992832, 0.983185, 0.973932, 0.965038, 0.932590, 0.904339, 0.879457, 0.857301],  # 25.0 C
                                    [1.000279, 0.994319, 0.984735, 0.975539, 0.966694, 0.934382, 0.906194, 0.881334, 0.859185],  # 30.0 C
                                    [1.003796, 0.997883, 0.988374, 0.979243, 0.970455, 0.938287, 0.910145, 0.885276, 0.863108]]) # 40.0 C

    densities = 1.0/specific_volumes 
    mass_NaCl_over_mass_water = molalities * 58.443 / 1000.0

    def error_function((a, b, c, d, e)):
        sum = 0.0
        for t, i in zip(temperatures, range(len(temperatures))):
            for m, j in zip(mass_NaCl_over_mass_water, range(len(mass_NaCl_over_mass_water))):
                sum += ((a*t**2 + b*t + c*m**2 + d*m + e) - densities[i][j])**2
        return sum


    a, b, c, d, e = scipy.optimize.fmin(error_function, (0.0001, 2, 3, 4, 1), full_output=1, xtol=1e-29, ftol=1e-29, maxiter=10000, maxfun=10000)[0]

    x = []
    y = []
    ztable = []
    zfunc = []
    count = 0
    max_diff = 0.0
    for t, i in zip(temperatures, range(len(temperatures))):
        for m, j in zip(mass_NaCl_over_mass_water, range(len(mass_NaCl_over_mass_water))):
            x.append(t)
            y.append(m)
            d_table = densities.ravel()[count]
            d_func = a*t**2 + b*t + c*m**2 + d*m + e
            ztable.append(d_table)
            zfunc.append(d_func)
            if abs(d_table - d_func) > max_diff:
                max_diff = abs(d_table - d_func)
            count += 1

    print "Density = %0.4g r^2 + %0.4g r + %0.4e t^2 + %0.4e t + %0.4f" % (c, d, a, b, e)
    print "Max diff = %0.4g" % max_diff
    return a, b, c, d, e

def plot():
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(x, y, ztable, "ro")
    ax.plot(x, y, zfunc, "bo")
    ax.set_xlabel("Temperature (C)")
    ax.set_ylabel("Mass NaCl / Mass H$_2$O")
    ax.set_zlabel("Density (g/mL)")
    plt.show()

def main():
    """
    """
    if len(sys.argv) == 1:
        H2O_mass = float(raw_input("Enter H2O mass: "))
        NaCl_mass = float(raw_input("Enter NaCl mass: "))
        t = float(raw_input("Enter temperature (C): "))
    elif len(sys.argv) == 4:
        H2O_mass = float(sys.argv[1])
        NaCl_mass = float(sys.argv[2])
        t = float(sys.argv[3])
    else:
        print "USAGE:    python NaCl.py mass_of_H2O mass_of_NaCl temperature_in_Celsius"
        print "EXAMPLE:  python NaCl.py 10.416 1.013 37"
        exit()
    m = NaCl_mass / H2O_mass
    a, b, c, d, e = get_constants()
    print a*t**2 + b*t + c*m**2 + d*m + e

if __name__ == "__main__":
    main()


    
    