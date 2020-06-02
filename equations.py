import numpy as np
from sympy import symbols, Matrix, poly

pv1_0, pv1_1, pv1_2, pv1_4, pv1_5, pv1_6, \
        pv1_7, pv1_8, pv1_9, pv1_10,  \
        pv1_12, pv1_13, pv1_14, pv1_15, pv1_16, \
        pv1_17, pv1_18, pv1_19, pv1_20, pv1_21, \
        pv1_22, \
        pv1_24, pv1_25, pv1_26, pv1_27, \
        pv1_28, pv1_29, pv1_30, \
        pv1_31, pv1_32, pv1_33, pv1_34, \
        pv1_35, pv1_36, pv1_37, pv1_38 = symbols("pv1_0 pv1_1 pv1_2 pv1_4 pv1_5 pv1_6 pv1_7 pv1_8 pv1_9 \
                                                pv1_10 pv1_12 pv1_13 pv1_14 pv1_15 pv1_16 pv1_17 pv1_18 \
                                                pv1_19 pv1_20 pv1_21 pv1_22 pv1_24 pv1_25 pv1_26 pv1_27 pv1_28 \
                                                pv1_29 pv1_30 pv1_31 pv1_32 pv1_33 pv1_34 pv1_35 pv1_36 pv1_37 pv1_38")


pv2_0, pv2_1, pv2_2, pv2_4, pv2_5, pv2_6, \
        pv2_7, pv2_8, pv2_9, pv2_10,  \
        pv2_12, pv2_13, pv2_14, pv2_15, pv2_16, \
        pv2_17, pv2_18, pv2_19, pv2_20, pv2_21, \
        pv2_22, \
        pv2_24, pv2_25, pv2_26, pv2_27, \
        pv2_28, pv2_29, pv2_30, \
        pv2_31, pv2_32, pv2_33, pv2_34, \
        pv2_35, pv2_36, pv2_37, pv2_38 = symbols("pv2_0 pv2_1 pv2_2 pv2_4 pv2_5 pv2_6 pv2_7 pv2_8 pv2_9 \
                                                pv2_10 pv2_12 pv2_13 pv2_14 pv2_15 pv2_16 pv2_17 pv2_18 \
                                                pv2_19 pv2_20 pv2_21 pv2_22 pv2_24 pv2_25 pv2_26 pv2_27 pv2_28 \
                                                pv2_29 pv2_30 pv2_31 pv2_32 pv2_33 pv2_34 pv2_35 pv2_36 pv2_37 pv2_38")

u, v = symbols("u v")
CD11, CD12, CD21, CD22 = symbols("CD11 CD12 CD21 CD22 ")
cd = Matrix([[CD11,CD12],[CD21,CD22]])
x, y = cd*Matrix([u, v])
cd_inverse=cd**-1


tpvx = pv1_0 + pv1_1*x + pv1_2*y + pv1_4*x**2 + pv1_5*x*y + pv1_6*y**2 + \
    pv1_7*x**3 + pv1_8*x**2*y + pv1_9*x*y**2 + pv1_10*y**3 +  \
    pv1_12*x**4 + pv1_13*x**3*y + pv1_14*x**2*y**2 + pv1_15*x*y**3 + pv1_16*y**4
    # pv1_17*x**5 + pv1_18*x**4*y + pv1_19*x**3*y**2 + pv1_20*x**2*y**3 + pv1_21*x*y**4 + \
    # pv1_22*y**5 + \
    # pv1_24*x**6 + pv1_25*x**5*y + pv1_26*x**4*y**2 + pv1_27*x**3*y**3 + \
    # pv1_28*x**2*y**4 + pv1_29*x*y**5 + pv1_30*y**6 + \
    # pv1_31*x**7 + pv1_32*x**6*y + pv1_33*x**5*y**2 + pv1_34*x**4*y**3 + \
    # pv1_35*x**3*y**4 + pv1_36*x**2*y**5 + pv1_37*x*y**6 + pv1_38*y**7

tpvy = pv2_0 + pv2_1*y + pv2_2*x + pv2_4*y**2 + pv2_5*y*x + pv2_6*x**2 + \
    pv2_7*y**3 + pv2_8*y**2*x + pv2_9*y*x**2 + pv2_10*x**3 + \
    pv2_12*y**4 + pv2_13*y**3*x + pv2_14*y**2*x**2 + pv2_15*y*x**3 + pv2_16*x**4 
    # pv2_17*y**5 + pv2_18*y**4*x + pv2_19*y**3*x**2 + pv2_20*y**2*x**3 + pv2_21*y*x**4 + \
    # pv2_22*x**5 + \
    # pv2_24*y**6 + pv2_25*y**5*x + pv2_26*y**4*x**2 + pv2_27*y**3*x**3 + \
    # pv2_28*y**2*x**4 + pv2_29*y*x**5 + pv2_30*x**6 + \
    # pv2_31*y**7 + pv2_32*y**6*x + pv2_33*y**5*x**2 + pv2_34*y**4*x**3 + \
    # pv2_35*y**3*x**4 + pv2_36*y**2*x**5 + pv2_37*y*x**6 + pv2_38*x**7


print(poly(tpvx))
print("\n\n")
print(poly(tpvy))

tpvu, tpvv = cd_inverse*Matrix([tpvx, tpvy])

tpvu.expand()
tpvv.expand()
