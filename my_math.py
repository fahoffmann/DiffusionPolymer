from numpy import sqrt, sin, cos, exp, pi

def distance(x, y, z):
    return sqrt(x ** 2 + y ** 2 + z ** 2)

def spherical_to_cartesian(r, theta, phi):
    x = r * sin(theta) * cos(phi)
    y = r * sin(theta) * sin(phi)
    z = r * cos(theta)
    return (x, y, z)

def stretched_exponential(t, b, tau):
    return 1 - exp(-(t / tau) ** b)

def stretched_exponential2(t, b, tau):
    return exp(-(t / tau) ** b)

def fick_solution_element(n, t):
    return 1.0 / (n ** 2) * exp(-n ** 2 * pi ** 2 * t)

def fick_solution(n_max, t):
    M = 0
    for n in range(1, n_max + 1):
        M += fick_solution_element(n, t)
    return 1 - 6.0 / (pi ** 2) * M