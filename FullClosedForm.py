import numpy as np

def Sqrt(z):
    if z >= 0:
        return np.sqrt(z)
    else:
        return np.sqrt(-z) * 1j

def solutions(ex, ey, px, py):

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

    return sols
