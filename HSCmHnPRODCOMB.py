#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      Kevin Sartor
#
# Created:     20-03-2018
# Copyright:   (c) Kevin Sartor
# Licence:     CC BY-NC-SA
#-------------------------------------------------------------------------------

from __future__ import division
from scipy.optimize import fsolve
from CoolProp.CoolProp import PropsSI
import numpy as np
import math

def mm_p(m,n,f): #Molar mass of combustion complete products of CmHn combustion
    #Based on the work of P. Ngendakumana
    T=298.15
    P=1e5
    MM=[]
    x=[]
    for i in range(5):
           MM.append(0.0)
           x.append(0.0)
    "Molar mass of the different constituents of the combustion products: kg/kmol"
    MM[1]=PropsSI('M','T', T,'P',P,'CarbonDioxide')*1000
    MM[2]=PropsSI('M','T',T,'P',P,'Water')*1000
    MM[3]=PropsSI('M','T',T,'P',P,'Oxygen')*1000
    MM[4]=PropsSI('M','T',T,'P',P,'Nitrogen')*1000

    "Stoechiometric fuel-air ratio and excess air"
    f_st = (m*12+n*1)/((m+n/4.0)*(MM[3]+(79/21.0)*MM[4]))
    e = f_st/f-1

    "Mole fractions of the different constituents - wet basis"
    ntkmol = m+(n/2.0)+(m+n/4.0)*(e+(1+e)*(79/21.0))
    x[1] = m/ntkmol
    x[2] = (n/2)/ntkmol
    x[3] = (m+n/4.0)*e/ntkmol
    x[4] = (m+n/4.0)*(1+e)*(79/21.0)/ntkmol
    mm_p=0
    for i in range(len(x)):
            mm_p += x[i]*MM[i]
    return mm_p

def ideal(T,fluid):
    if (fluid=='N2'):
        if (T>1000):
            coef=np.array([0.02926640E+02,0.14879768E-02,-0.05684760E-05,0.10097038E-09,-0.06753351E-13,-0.09227977E+04 ,0.05980528E+02] , dtype=float)
        else:
            coef = np.array([0.03298677E+02,0.14082404E-02,-0.03963222E-04,0.05641515E-07,-0.02444854E-10,-0.10208999E+04, 0.03950372E+02] , dtype=float)
    elif(fluid=='O2'):
        if (T>1000):
            coef=np.array([3.28253784E+00,1.48308754E-03,-7.57966669E-07,2.09470555E-10,-2.16717794E-14 ,-1.08845772E+03,5.45323129E+00] , dtype=float)
        else:
            coef = np.array([3.78245636E+00,-2.99673416E-03,9.84730201E-06,-9.68129509E-09,3.24372837E-12,-1.06394356E+03,3.65767573E+00] , dtype=float)
    elif(fluid=='CO2'):
        if (T>1000):
            coef=np.array([3.85746029E+00,4.41437026E-03,-2.21481404E-06,5.23490188E-10,-4.72084164E-14,-4.87591660E+04,2.27163806E+00] , dtype=float)
        else:
            coef = np.array([2.35677352E+00,8.98459677E-03,-7.12356269E-06, 2.45919022E-09,-1.43699548E-13,-4.83719697E+04,9.90105222E+00    ] , dtype=float)
    elif(fluid=='H2O'):
        if (T>1000):
            coef=np.array([ 3.03399249E+00,2.17691804E-03,-1.64072518E-07,-9.70419870E-11 ,1.68200992E-14,-3.00042971E+04,4.96677010E+00] , dtype=float)
        else:
            coef = np.array([4.19864056E+00,-2.03643410E-03,6.52040211E-06,-5.48797062E-09,1.77197817E-12,-3.02937267E+04,-8.49032208E-01] , dtype=float)
    else:
        print("The fluid "+fluid+" does not exist")
    a1=coef[0]
    a2=coef[1]
    a3=coef[2]
    a4=coef[3]
    a5=coef[4]
    a6=coef[5]
    a7=coef[6]
    MM=PropsSI('M','T', T,'P',100,fluid)*1000
    R=8314/MM
    Cp = (a1 + a2*T + a3*T**2 + a4*T**3 + a5*T**4)*R
    H = (a1 + a2*T/2 + a3*T**2 /3+ a4*T**3 /4 + a5*T**4 /5 + a6/T)*R*T
    S  =( a1*math.log(T) + a2*T + a3*T**2 /2 + a4*T**3 /3 + a5*T**4 /4 + a7)*R
    props=[Cp,H,S]
    return props

def HS_PC(m,n,f,Temp,Pres):   #Retourne H et S en fonction de T et P
    T=Temp
    P=Pres
    MM=[]
    x=[]
    h=[]
    s=[]
    pp=[]
    for i in range(5):
        MM.append(0.0)
        x.append(0.0)
        h.append(0.0)
        s.append(0.0)
        pp.append(0.0)
    MM[1]=PropsSI('M','T', T,'P',P,'CarbonDioxide')*1000
    MM[2]=PropsSI('M','T',T,'P',P,'Water')*1000
    MM[3]=PropsSI('M','T',T,'P',P,'Oxygen')*1000
    MM[4]=PropsSI('M','T',T,'P',P,'Nitrogen')*1000
    "Stoechiometric fuel-air ratio and excess air"
    f_st = (m*12+n*1)/((m+n/4.0)*(MM[3]+(79/21.0)*MM[4]))
    e = f_st/f-1

    "Mole fractions of the different constituents - wet basis"
    ntkmol = m+(n/2.0)+(m+n/4.0)*(e+(1+e)*(79/21.0))
    x[1] = m/ntkmol
    x[2] = (n/2)/ntkmol
    x[3] = (m+n/4.0)*e/ntkmol
    x[4] = (m+n/4.0)*(1+e)*(79/21.0)/ntkmol
    mm_p=0
    for i in range(len(x)):
            mm_p += x[i]*MM[i]
    MM_p = mm_p

    "Enthalpy of the mixture : h_p in J/kg"
    h[1] = ideal(Temp,'CO2')[1] #PropsSI('H','T', Temp,'P',P,'CarbonDioxide')
    h[2] = ideal(Temp,'H2O')[1]#PropsSI('H','T', Temp,'P',P,'Water')
    h[3] = ideal(Temp,'O2')[1]#PropsSI('H','T', Temp,'P',P,'Oxygen')
    h[4] = ideal(Temp,'N2')[1]#PropsSI('H','T', Temp,'P',P,'Nitrogen')

    h_pc=0
    for i in range(len(x)):
            h_pc += x[i]*MM[i]*h[i]/MM_p

    "Entropy of the mixture : s_p in J/(kg.K)"
    for i in range(len(x)):
        pp[i] = x[i]*Pres

    s[1] = ideal(Temp,'CO2')[2]
    s[2] = ideal(Temp,'H2O')[2]
    s[3] = ideal(Temp,'O2')[2]
    s[4] = ideal(Temp,'N2')[2]
    s_p=0
    for i in range(len(x)):
        s_p += x[i]*MM[i]*s[i]/MM_p
    return h_pc,s_p

def TS_PC(m,n,f,H,Pres):
    P=Pres
    MM=[]
    x=[]
    h=[]
    s=[]
    pp=[]
    for i in range(5):
        MM.append(0.0)
        x.append(0.0)
        h.append(0.0)
        s.append(0.0)
        pp.append(0.0)
    MM[1]=PropsSI('M','H',H,'P',P,'CarbonDioxide')*1000
    MM[2]=PropsSI('M','H',H,'P',P,'Water')*1000
    MM[3]=PropsSI('M','H',H,'P',P,'Oxygen')*1000
    MM[4]=PropsSI('M','H',H,'P',P,'Nitrogen')*1000
    "Stoechiometric fuel-air ratio and excess air"
    f_st = (m*12+n*1)/((m+n/4.0)*(MM[3]+(79/21.0)*MM[4]))
    e = f_st/f-1

    "Mole fractions of the different constituents - wet basis"
    ntkmol = m+(n/2.0)+(m+n/4.0)*(e+(1+e)*(79/21.0))
    x[1] = m/ntkmol
    x[2] = (n/2)/ntkmol
    x[3] = (m+n/4.0)*e/ntkmol
    x[4] = (m+n/4.0)*(1+e)*(79/21.0)/ntkmol
    mm_p=0

    for i in range(len(x)):
            mm_p += x[i]*MM[i]
    MM_p = mm_p


    def f(T):
        return (H-(x[1]*MM[1]*ideal(T,'CO2')[1]+x[2]*MM[2]*ideal(T,'H2O')[1]+x[3]*MM[3]*ideal(T,'O2')[1]+x[4]*MM[4]*ideal(T,'N2')[1])/MM_p)

    T = fsolve(f, 1500)
    f(T)

    return T
