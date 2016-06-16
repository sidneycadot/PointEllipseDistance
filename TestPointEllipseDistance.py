#! /usr/bin/env python3

# See:
#      http://www.geometrictools.com/Documentation/DistancePointEllipseEllipsoid.pdf
#
# Given the ellipse (x / ex) ** 2 + (y / ey) ** 2 == 1

import ClosedForm
import ReferenceImplementation

import math
import numpy as np
import matplotlib.pyplot as plt

def test1():

    ex = 2.0
    ey = 1.0

    x = np.linspace(-5, 5, 101)
    y = np.linspace(-5, 5, 101)

    xv, yv = np.meshgrid(x, y)

    print("Calculating reference ...")
    #zv1 = ReferenceImplementation.DistancePointEllipseVectorized(ex, ey, xv, yv)
    zv1 = ReferenceImplementation.DistancePointEllipseVectorized(ex, ey, xv, yv)
    zv1 = zv1[2]

    print("Calculating closed-form ...")
    #zv2 = ClosedForm.solution_vectorized(ex, ey, xv, yv)
    zv2 = ClosedForm.solution_vectorized(ex, ey, xv, yv)
    print("Done.")

    zv2 = np.hypot(zv2[0] - xv, zv2[1] - yv)

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

def test2():

    for i in range(100000):

        if i % 1000 == 0:
            print("checked", i)

        ex = np.random.rand() * 20 - 10
        ey = np.random.rand() * 20 - 10
        px = np.random.rand() * 20 - 10
        py = np.random.rand() * 20 - 10

        try:
            ref = ClosedForm.reference_solution(ex, ey, px, py)
            sol = ClosedForm.solution(ex, ey, px, py)

            error = math.sqrt((sol[0] - ref[0]) ** 2 + (sol[1] - ref[1]) ** 2)

            if error > 1e-4:
                raise RuntimeError("Error too large ({})".format(error))

        except RuntimeError as exception:
            print("Exception [{}]: ex = {}, ey = {}, px = {}, py = {}".format(exception, ex, ey, px, py))
            continue

def main():
    test1()

if __name__ == "__main__":
    main()
