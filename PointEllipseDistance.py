#! /usr/bin/env python3

# See:
#      http://www.geometrictools.com/Documentation/DistancePointEllipseEllipsoid.pdf
#
# Given the ellipse (x / ex) ** 2 + (y / ey) ** 2 == 1

import ClosedForm
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

ex = 2
ey = 1

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
