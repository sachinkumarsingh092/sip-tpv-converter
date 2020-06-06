import numpy as np
from sympy import symbols, Matrix, poly


''' 
A better implementation is:
pv1 = symbols(pv1_0:39)
and then using individual elements like pv1[0]...etc.
'''



############## For PV->SIP #############
pv1 = symbols("pv1_0:39")
pv2 = symbols("pv2_0:39")


u, v = symbols("u v")
CD11, CD12, CD21, CD22 = symbols("CD11 CD12 CD21 CD22 ")
cd = Matrix([[CD11,CD12],[CD21,CD22]])
x, y = cd*Matrix([u, v])
cd_inverse=cd**-1

# Upto 4th order only is used in calculations.
tpvx = pv1[0] + pv1[1]*x + pv1[2]*y + pv1[4]*x**2 + pv1[5]*x*y + pv1[6]*y**2 + \
    pv1[7]*x**3 + pv1[8]*x**2*y + pv1[9]*x*y**2 + pv1[10]*y**3 +  \
    pv1[12]*x**4 + pv1[13]*x**3*y + pv1[14]*x**2*y**2 + pv1[15]*x*y**3 + pv1[16]*y**4
    # pv1[17]*x**5 + pv1[18]*x**4*y + pv1[19]*x**3*y**2 + pv1[20]*x**2*y**3 + pv1[21]*x*y**4 + \
    # pv1[22]*y**5 + \
    # pv1[24]*x**6 + pv1[25]*x**5*y + pv1[26]*x**4*y**2 + pv1[27]*x**3*y**3 + \
    # pv1[28]*x**2*y**4 + pv1[29]*x*y**5 + pv1[30]*y**6 + \
    # pv1[31]*x**7 + pv1[32]*x**6*y + pv1[33]*x**5*y**2 + pv1[34]*x**4*y**3 + \
    # pv1[35]*x**3*y**4 + pv1[36]*x**2*y**5 + pv1[37]*x*y**6 + pv1[38]*y**7

tpvy = pv2[0] + pv2[1]*y + pv2[2]*x + pv2[4]*y**2 + pv2[5]*y*x + pv2[6]*x**2 + \
    pv2[7]*y**3 + pv2[8]*y**2*x + pv2[9]*y*x**2 + pv2[10]*x**3 + \
    pv2[12]*y**4 + pv2[13]*y**3*x + pv2[14]*y**2*x**2 + pv2[15]*y*x**3 + pv2[16]*x**4 
    # pv2[17]*y**5 + pv2[18]*y**4*x + pv2[19]*y**3*x**2 + pv2[20]*y**2*x**3 + pv2[21]*y*x**4 + \
    # pv2[22]*x**5 + \
    # pv2[24]*y**6 + pv2[25]*y**5*x + pv2[26]*y**4*x**2 + pv2[27]*y**3*x**3 + \
    # pv2[28]*y**2*x**4 + pv2[29]*y*x**5 + pv2[30]*x**6 + \
    # pv2[31]*y**7 + pv2[32]*y**6*x + pv2[33]*y**5*x**2 + pv2[34]*y**4*x**3 + \
    # pv2[35]*y**3*x**4 + pv2[36]*y**2*x**5 + pv2[37]*y*x**6 + pv2[38]*x**7

# print(poly(tpvx))
# print("\n\n")
# print(poly(tpvy))

tpvu, tpvv = cd_inverse*Matrix([tpvx, tpvy])
tpvu.expand()
tpvv.expand()



############For SIP->PV#############
x, y = symbols('x y')
CD_INV11, CD_INV12, CD_INV21, CD_INV22 = symbols("CD_INV11 CD_INV12 CD_INV21 CD_INV22 ")
cd_inverse = Matrix([[CD_INV11,CD_INV12],[CD_INV21,CD_INV22]])

uprime, vprime = cd_inverse*Matrix([x, y])
usum = uprime
vsum = vprime

A=symbols('A_0_0 A_0_1 A_0_2 A_0_3 A_0_4  \
            A_1_0 A_1_1 A_1_2 A_1_3 A_2_0 \
            A_2_1 A_2_2 A_3_0 A_3_1 A_4_0')


B=symbols('B_0_0 B_0_1 B_0_2 B_0_3 B_0_4  \
            B_1_0 B_1_1 B_1_2 B_1_3 B_2_0 \
            B_2_1 B_2_2 B_3_0 B_3_1 B_4_0')

k=0
for i in range(5):
    for j in range (0, 5-i):
        usum+=A[k]*uprime**i*vprime**j
        vsum+=B[k]*uprime**i*vprime**j
        k+=1

# print(usum)

sipx, sipy = cd*Matrix([usum, vsum])

print(poly(sipx.expand()))

