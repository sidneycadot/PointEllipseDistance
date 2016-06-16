import numpy as np

def Sqrt(z):
    if z >= 0:
        return np.sqrt(z)
    else:
        return np.sqrt(-z) * 1j

def PosSqrt(z):
    #if z < 0:
    #    raise RuntimeError("Attempt sqrt({})".format(z))
    return np.sqrt(abs(z))

def reference_solution(ex, ey, px, py):

    T1 = (-ex**4+2*ex**2*ey**2-ey**4+ex**2*px**2+ey**2*py**2)
    T2 = (ex**4-2*ex**2*ey**2+ey**4)
    T3 = (ex**3*px-ex*ey**2*px)
    T4 = (ex**2-ey**2)
    T5 = (ex**2*ey**2*px**2*py**2)
    T6 = (ex**2*px**2)
    T7 = (2**(1/3)*(ex**4-2*ex**2*ey**2+ey**4-ex**2*px**2-ey**2*py**2)**2)
    T8 = (2*ex**2*px**2)

    U1 = (108*T2*T3**2-108*ex**2*px**2*T3**2+72*ex**2*T2*px**2*T1+36*T3**2*T1+2*T1**3)
    U2 = (-12*ex**2*T2*px**2 + 12*T3**2 + T1**2)
    U3 = ((-16*ex*px)/T4+(8*ex**3*px**3)/T4**3-(8*ex*px*T1)/T4**3)

    V1 = (108*T2*T3**2 - 108*ex**2*px**2*T3**2 + 72*ex**2*T2*px**2*T1 + 36*T3**2*T1 + 2*T1**3 + Sqrt(U1**2 - 4*U2**3)) ** (1/3)

    W1 = (3*T4**2*V1)
    W2 = T1/(3*T2)+T7/W1+V1/(3*2**(1/3)*T4**2)

    Z1 = Sqrt(T6/T4**2-T1/T4**2+W2)
    Z2 = Sqrt(T8/T4**2-T1/T4**2-W2-U3/(4*Z1))
    Z3 = Sqrt(T8/T4**2-T1/T4**2-W2+U3/(4*Z1))

    Z5 = ((ex*px)/(2*T4)-Z1/2-Z2/2)
    Z6 = ((ex*px)/(2*T4)+Z1/2-Z3/2)
    Z7 = ((ex*px)/(2*T4)-Z1/2+Z2/2)
    Z8 = ((ex*px)/(2*T4)+Z1/2+Z3/2)

    Z9 =  (-(ex**3*px)+ex*ey**2*px+ex**4*Z5-2*ex**2*ey**2*Z5+ey**4*Z5-ey**2*py**2*Z5+ex**3*px*Z5**2-ex*ey**2*px*Z5**2-ex**4*Z5**3+2*ex**2*ey**2*Z5**3-ey**4*Z5**3)
    Z10 = (-(ex**3*px)+ex*ey**2*px+ex**4*Z7-2*ex**2*ey**2*Z7+ey**4*Z7-ey**2*py**2*Z7+ex**3*px*Z7**2-ex*ey**2*px*Z7**2-ex**4*Z7**3+2*ex**2*ey**2*Z7**3-ey**4*Z7**3)
    Z11 = (-(ex**3*px)+ex*ey**2*px+ex**4*Z6-2*ex**2*ey**2*Z6+ey**4*Z6-ey**2*py**2*Z6+ex**3*px*Z6**2-ex*ey**2*px*Z6**2-ex**4*Z6**3+2*ex**2*ey**2*Z6**3-ey**4*Z6**3)
    Z12 = (-(ex**3*px)+ex*ey**2*px+ex**4*Z8-2*ex**2*ey**2*Z8+ey**4*Z8-ey**2*py**2*Z8+ex**3*px*Z8**2-ex*ey**2*px*Z8**2-ex**4*Z8**3+2*ex**2*ey**2*Z8**3-ey**4*Z8**3)

    sols = [
        (
            (ex*Z7)/Sqrt(Z7**2+Z10**2/T5),
            -(Z10/(ex*px*py*Sqrt(Z7**2+Z10**2/T5)))
        ),
        (
            (ex*Z5)/Sqrt(Z5**2+Z9**2/T5),
            -(Z9/(ex*px*py*Sqrt(Z5**2+Z9**2/T5)))
        ),
        (
            (ex*Z8)/Sqrt(Z8**2+Z12**2/T5),
            -(Z12/(ex*px*py*Sqrt(Z8**2+Z12**2/T5)))
        ),
        (
            (ex*Z6)/Sqrt(Z6**2+Z11**2/T5),
            -(Z11/(ex*px*py*Sqrt(Z6**2+Z11**2/T5)))
        )
    ]

    pred1 = (ex > 0) ^ (px > 0)
    pred2 = (ex > -ey) ^ (px > 0) ^ (ex > ey) ^ (ex > 0)

    select = int(pred2) * 2 + int(pred1)

    sol = sols[select]

    if abs(sol[0].imag) >= 1e-3:
        raise RuntimeError("Expected real, got {}".format(sol[0]))

    if abs(sol[1].imag) >= 1e-3:
        raise RuntimeError("Expected real, got {}".format(sol[1]))

    sol = (float(sol[0].real), float(sol[1].real))

    return sol

def solution(ex, ey, px, py):

    pred1 = (ex > 0) ^ (px > 0)
    pred3 = (ex > -ey) ^ (ex > ey)

    S1 = +1 if pred1 else -1
    S3 = +1 if pred3 else -1

    T4 = (ex - ey) * (ex + ey)
    T2 = T4**2
    T3 = T4 * (px * ex)
    T6 = (ex*px) ** 2
    T8 = (ey*py) ** 2

    T1 = T6 + T8 - T2
    T7 = 2**(1/3)*T1**2

    U1 = (108*T2*T3**2 - 108*ex**2*px**2*T3**2 + 72*ex**2*T2*px**2*T1 + 36*T3**2*T1 + 2*T1**3)
    U2 = (-12*ex**2*T2*px**2 + 12*T3**2 + T1**2)
    U3 = 2 * (ex**3*px**3 - ex*px*T1 - 2*ex*px*T2)

    V1 = (U1 + Sqrt(U1 ** 2 - 4 * U2 ** 3)) ** (1/3)

    W = (T1 + T7/V1 + V1/(2**(1/3))) / 3

    Z1 = S1 * PosSqrt(T6 - T1 + W)

    ZZ = PosSqrt(2 * T6 - T1 - W + U3 / Z1)

    X = (ex*px + Z1 + S1 * S3 * ZZ) / (2 * T4)

    Y = ( ex**3*px - ex*ey**2*px - ex**4*X + 2*ex**2*ey**2*X - ey**4*X + ey**2*py**2*X - ex**3*px*X**2 + ex*ey**2*px*X**2 + ex**4*X**3 - 2*ex**2*ey**2*X**3 + ey**4*X**3) / (ex*px*ey*py)

    denom = PosSqrt(X ** 2 + Y ** 2)

    solx = (ex * X) / denom
    soly = (ey * Y) / denom

    if abs(solx.imag) >= 1e-3:
        raise RuntimeError("Expected real, got {}".format(solx))

    if abs(soly.imag) >= 1e-3:
        raise RuntimeError("Expected real, got {}".format(soly))

    sol = (float(solx.real), float(soly.real))

    return sol

solution_vectorized = np.vectorize(solution, otypes = [np.float64, np.float64])
