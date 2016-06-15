#! /usr/bin/env python3

# See:
#      http://www.geometrictools.com/Documentation/DistancePointEllipseEllipsoid.pdf
#
# Given the ellipse (x / ex) ** 2 + (y / ey) ** 2 == 1

import ClosedForm
import FullClosedForm
import math
import numpy as np
import matplotlib.pyplot as plt

def Sqr(x):
    return x * x

def RobustLength(v0, v1):
    return math.sqrt(v0 * v0 + v1 * v1)
    if abs(v0) == max(abs(v0), abs(v1)):
        return abs(v0) * math.sqrt(1 + (v1 / v0) ** 2)
    else:
        return abs(v1) * math.sqrt(1 + (v0 / v1) ** 2)

def GetRoot(r0, z0, z1, g):
    n0 = r0 * z0
    s0 = z1 - 1
    s1 = 0 if g < 0 else RobustLength(n0, z1) - 1
    s = 0

    maxIterations = 100

    i = 0
    while i < maxIterations:
        s = (s0 + s1) / 2
        if s == s0 or s == s1:
            break
        ratio0 = n0 / (s + r0)
        ratio1 = z1 / (s +  1)
        g = Sqr(ratio0) + Sqr(ratio1) - 1
        if g > 0:
            s0 = s
        elif g < 0:
            s1 = s
        else:
            break
        i += 1

    return s

def DistancePointEllipse(e0, e1, y0, y1):

    assert e0 >= e1 > 0

    if y0 < 0:
        (x0, x1, distance) = DistancePointEllipse(e0, e1, -y0, y1)
        return (-x0, x1, distance)

    if y1 < 0:
        (x0, x1, distance) = DistancePointEllipse(e0, e1, y0, -y1)
        return (x0, -x1, distance)

    if y1 > 0:
        if y0 > 0:
            z0 = y0 / e0
            z1 = y1 / e1
            g = Sqr(z0) + Sqr(z1) - 1
            if g != 0:
                r0 = Sqr(e0 / e1)
                sbar = GetRoot(r0, z0, z1, g)
                x0 = r0 * y0 / (sbar + r0)
                x1 = y1 / (sbar + 1)
                distance = math.sqrt(Sqr(x0 - y0) + Sqr(x1 - y1))
            else:
                x0 = y0
                x1 = y1
                distance = 0
        else: # y0 == 0
            x0 = 0
            x1 = e1
            distance = abs(y1 - e1)
    else: # y1 == 0
        numer0 = e0 * y0
        denom0 = Sqr(e0) - Sqr(e1)
        if numer0 < denom0:
            xde0 = numer0 / denom0
            x0 = e0 * xde0
            x1 = e1 * math.sqrt(1 - xde0 * xde0)
            distance = math.sqrt(Sqr(x0 - y0) + Sqr(x1))
        else:
            x0 = e0
            x1 = 0
            distance = abs(y0 - e0)

    return (x0, x1, distance)

#def DistancePointEllipse(e0, e1, y0, y1):
#    return e0 + e1 + y0 + y1

DistancePointEllipseVectorized = np.vectorize(DistancePointEllipse, otypes = [np.float64, np.float64, np.float64])

if False:

    print(DistancePointEllipse(ex, ey, 2, 4))
    print(ClosedForm.distance(ex, ey, 2, 4))

    x = np.linspace(-5, 5, 101)
    y = np.linspace(-5, 5, 101)

    xv, yv = np.meshgrid(x, y)

    zv1 = DistancePointEllipseVectorized(ex, ey, xv, yv)
    zv1 = zv1[2]

    zv2 = ClosedForm.distance(ex, ey, xv, yv)

    plt.subplot(311)
    plt.imshow(zv1)
    plt.colorbar()
    plt.subplot(312)
    plt.imshow(zv2)
    plt.colorbar()
    plt.subplot(313)
    plt.imshow(zv2 - zv1)
    plt.colorbar()
    plt.show()

#print(FullClosedForm.solutions(ex, ey, 2, 4))

for i in range(10000):
    ex = np.random.rand() * 10 - 5
    ey = np.random.rand() * 10 - 5
    px = np.random.rand() * 10 - 5
    py = np.random.rand() * 10 - 5

    solutions = FullClosedForm.solutions(ex, ey, px, py)

    best_i = None
    best_x = None
    best_y = None
    best_distance = None
    for i in range(len(solutions)):
        (sx, sy) = solutions[i]
        sx = np.complex128(sx)
        sy = np.complex128(sy)

        if abs(sx.imag) > 1e-10 or abs(sy.imag) > 1e-10:
            continue

        sx = sx.real
        sy = sy.real

        distance = np.sqrt((sx - px) ** 2 + (sy - py) ** 2)

        if best_distance is None or distance < best_distance:
            best_i = i
            best_x = sx
            best_y = sy
            best_distance = distance

    if best_distance is None:
        continue

    pred1 = (ex > 0) ^ (px > 0)
    pred2 = (ex > -ey) ^ (px > 0) ^ (ex > ey) ^ (ex > 0)

    print("         {:d} {:d} {:d} {:d} {:d} {:d}       {:d}       {:d} pred1 {:d}    [{:d}]".format(ex>0, ey>0, px>0, py>0, ex>ey, ex>-ey, best_i, (best_i & 1) // 1, pred1, (best_i & 1) // 1 != pred1))

