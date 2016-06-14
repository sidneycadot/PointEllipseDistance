import numpy

def Sqrt(x):
    return numpy.sqrt(x)

def distance(ex, ey, px, py):

    T1 = (ex**4-2*ex**2*ey**2+ey**4)
    T2 = (-ex**4+2*ex**2*ey**2-ey**4+ex**2*px**2+ey**2*py**2)
    T3 = (ex**3*px-ex*ey**2*px)
    T4 = -T2
    T5 = (ex**2-ey**2)

    U1 = (108*T1*T3**2-108*ex**2*px**2*T3**2+72*ex**2*T1*px**2*T2+36*T3**2*T2+2*T2**3)
    U2 = (U1 + Sqrt(-4*(-12*ex**2*T1*px**2+12*T3**2+T2**2)**3+U1**2))**(1/3)
    U3 = ((-16*ex*px)/T5+(8*ex**3*px**3)/T5**3-(8*ex*px*T2)/T5**3)
    U4 = Sqrt(  (ex**2*px**2)/T5**2-T2/T5**2+T2/(3*T1)+(2**(1/3)*T4**2)/(3*T5**2*U2)+U2/(3*2**(1/3)*T5**2))
    U5 = Sqrt((2*ex**2*px**2)/T5**2-T2/T5**2-T2/(3*T1)-(2**(1/3)*T4**2)/(3*T5**2*U2)-U2/(3*2**(1/3)*T5**2)-U3/(4*U4))
    U6 = ((ex*px)/(2*T5)-U4/2+U5/2)
    U7 = (-(ex**3*px)+ex*ey**2*px+ex**4*U6-2*ex**2*ey**2*U6+ey**4*U6-ey**2*py**2*U6+ex**3*px*U6**2-ex*ey**2*px*U6**2-ex**4*U6**3+2*ex**2*ey**2*U6**3-ey**4*U6**3)

    U8 = Sqrt(U6**2 + U7**2/(ex**2*ey**2*px**2*py**2))

    return Sqrt((px-(ex*U6)/U8)**2+(py+U7/(ex*px*py*U8))**2)
