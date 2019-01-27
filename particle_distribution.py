import random
import numpy as np
import math as mt
import scipy.constants as sp
import matplotlib.pyplot as plt

def Dist(max_x, max_y, max_z, N, radius, density):
    i = 0
    particles = []
    while i < N:
        x = random.uniform(0, max_x)
        y = random.uniform(0, max_y)
        z = random.uniform(0, max_z)
        pos = np.zeros(3)
        pos[0] = x
        pos[1] = y
        pos[2] = z
        r = radius[i]
        loop_escape = False

        for particle in particles:
            diff = abs(pos - particle[1])
            sep = np.sqrt(diff[0]**2 + diff[1]**2 + diff[2]**2)
            if sep < (r + particle[0]):
                loop_escape = True
                break
                

        if loop_escape:
            continue

        mass = density * (4/3)*np.pi*r**3
        particles.append((r, pos, mass))
        i += 1

    file = open('particles.txt', 'w')
    j = 0

    for item in particles:
        print(j, -3, item[2]/2E30, item[0]/sp.astronomical_unit,
              item[1][0]/sp.astronomical_unit,
              item[1][1]/sp.astronomical_unit,
              item[1][2]/sp.astronomical_unit,
              0, 0, 0, 0, 0, 0, 3, sep=' ', file=file)
        j += 1
    file.close()

    return particles


def Size(N, target_radius, min_radius=1):
    m = target_radius/min_radius
    dist = np.random.pareto(3, N*5)
    dist = np.sort(dist)[3*N:4*N:1]
    dist_3_sum = sum(dist**3)
    alpha = target_radius/(dist_3_sum**(1/3))
    r_dist = alpha*dist
    min_r = min(r_dist)
    max_r = max(r_dist)
    print('min r: {:.3E}'.format(min_r))
    print('max r: {:.3E}'.format(max_r))
    R = max_r+min_r
    A = 4*np.pi*(R**2)
    a = np.pi*(min_r**2)
    f = 0.9
    n = int(np.ceil((A/a) * f))
    print('n: {}'.format(n))
    r = (sum(r_dist**3) / 0.6)**(1/3)
    print('collapsed Radius: {:.3E}m'.format(r))
    return r_dist
 
N = 1000
target_radius = 2000
target_density = 200
area = 20000
radius = Size(N, target_radius)
result = Dist(area, area, area, N, radius, target_density)
total_mass = sum([particle[2] for particle in result])
total_volume = (4/3)*np.pi*(area**3)
print('total_mass: {:.3E}kg'.format(total_mass))
total_density = total_mass / total_volume
print('total_density: {:.3E}kg/m3'.format(total_density))
collapsed_volume = sum([particle[0]**3 for particle in result])
collapsed_density = total_mass / (collapsed_volume/0.6)
print('collapsed_density: {}kg/m3'.format(collapsed_density))
dynamic_time = 1/np.sqrt(2*sp.gravitational_constant*total_density)
print('dynamic_time: {:.3E}s'.format(dynamic_time))
dDelta = dynamic_time * 0.1
print('dDelta: {:.3E}s'.format(dDelta))
print('dDelta: {:.3E}years/2pi'.format((dDelta/(2*np.pi)) / (sp.year)))

#plt.hist([particle[0] for particle in result], bins=1000)
#plt.show()