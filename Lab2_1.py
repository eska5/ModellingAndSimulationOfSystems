import numpy as np
import matplotlib.pyplot as plt

g = 10
m = 1
k = 1

x0, y0, z0 = -5, 0, 200
v0x, v0y, v0z = 2, 4, 0

def find_time(x0, y0, z0, v0x, v0y, v0z, g):

    a = (-1/2*g) - (v0y**2 + v0x**2)
    b = - (2*v0y*y0 + 2*v0x*x0 - v0z)
    c = - (x0**2 + y0**2 - z0)
    discriminant = b**2 - 4*a*c
    
    if discriminant >= 0:
        t1 = (-b + np.sqrt(discriminant)) / (2*a)
        t2 = (-b - np.sqrt(discriminant)) / (2*a)
        return t1, t2
    else:
        return None, None

# Function to calculate the normal vector at a given point (x, y)
def normal_vector(x, y):
    df_dx = 2*x
    df_dy = 2*y
    norm = np.sqrt(df_dx**2 + df_dy**2 + 1)
    return np.array([-df_dx/norm, -df_dy/norm, 1/norm])

# Function to calculate the reflected velocity vector
def reflected_velocity(v_before, n, k):
    return np.sqrt(k) * (v_before - 2 * np.dot(v_before, n) * n)

# Initialize arrays to store results
impact_points = []
times_between_collisions = []

# Simulation loop
for _ in range(5):

    t1,t2 = find_time(x0, y0, z0, v0x, v0y, v0z, g)
    if t1 is None:
        print("No real roots found. Exiting simulation.")
        break
    else:
        tI = t1 if t1 > t2 else t2

    # Wyznaczamy punkt uderzenia
    xI = x0 + v0x*tI
    yI = y0 + v0y*tI
    zI = z0 + v0z*tI - 0.5*g*tI**2
    
    # zapisujemy punkt uderzenia i czas pomiędzy uderzeniami
    impact_points.append((xI, yI, zI))
    times_between_collisions.append(tI)


    # Podmieniamy wartości dla kolejnego uderzenia
    v0 = np.array([v0x, v0y, v0z - g*tI])
    n = normal_vector(xI, yI)
    v_before = np.array([v0x, v0y, v0z - g*tI])
    v_after = reflected_velocity(v_before, n, k)

    x0, y0, z0 = xI, yI, zI 
    v0x, v0y, v0z = v_after
    
for i in range(5):
     print(f"{impact_points[i][0]:.4f}, {impact_points[i][1]:.4f}, {impact_points[i][2]:.4f}, {times_between_collisions[i]:.4f}")
