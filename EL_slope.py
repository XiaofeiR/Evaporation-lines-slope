
# -*- coding: utf-8 -*-
# @Author    : Xiaofei Ren
# @Time      : 2023/3/29 15:50
# @FileName  : EL_slope.py
# @Software  : PyCharm
# purpose    : Calculate the slope of EL

import numpy as np

#-----------------------------------------------------------
# Calculate α+ ： The calculation formula comes from Horita and Wesolowski(1994)
# α+ is the liquid-vapor equilibrium isotopic fractionation

def Tk(t):
    """
    Convert the temperature in Celsius to the Kelvin temperature 

    t： ℃
    """
    
    return t+273.15

def alpha_plus(atom,T):
    """
    Calculate α+

    atom : 'H' or 'O'
    T : Kelvin temperature 
    """

    T3 = T**3
    T2 = T**2

    if atom=='H':
        right = 1158.8 * T3 * 10 ** (-9) - 1620.1 * T2 * 10 ** (-6) + 794.84 * T * 10 ** (-3) - 161.04 + (
                    2.9992 * 10 ** 9) / T3
    elif atom=='O':

        right = -7.685 + 6.7123 * 10 ** 3 / T - 1.6664 * 10 ** 6 / T2 + 0.3504 * 10 ** 9 / T3

    else:
        right = 'None'

    i = right * 10 ** -3

    return np.e ** i
#------------------------------------------------------------


# -----------------------------------------------------------
# Calculate ε+   

def epsilon_plus(alpha):
    return (alpha-1)*10**3

#--------------------------------------------------------------


#--------------------------------------------------------------
# Calculate εk   相对湿度的函数 存在两种计算形式 单位 千分之

def epsilon_k(atom,n,method,h):
    """
    atom: 'H' or 'O'

    n: 0.5~1 ; 0.5 for open water bodies; 1 for soil water

    method : 0 or 1
        0 : Gibson(2008) GBC
        1 ：Paolo Benettin (2018) HESS

    h : relative humidity.  e.g. 0.65

    return : ‰
    """
    
    if method==0:
        if atom=='H':
            C = 25.0
        elif atom=='O':
            C = 28.6
        else:
            C = 'error'

        d = n*C*1*(1-h)

    elif method==1:
        if atom =='H':
            DiD = 0.9755
        elif atom=='O':
            DiD = 0.9723
        else:
            DiD = 'error'

        d = 1*n*(1-h)*(1-DiD)*10**3

    return d
#---------------------------------------------------------------


# ---------------------------------------------------------------
#  Calculate δA  ： δA=(δP-ε+)/α+
def delta_A(delta_P,epsilon,alpha):

    i = (delta_P-epsilon)/alpha
    return i
#------------------------------------------------------------------


#-----------------------------------------------------------------
# Calculate slope

#根据 论文 Recharge mechanisms of deep soil water revealed by water isotopes in deep loess deposits中的计算公式这个公式与 Haowen Liu (2022)中的公式相同 这里dX值仅仅为计算slope的氢或者氧

def dX(atm,Pre,h,t,n,method):
    """
    Pre : isotopic composition of precipitation.  e.g. [-38,-6]

    Other parameters are the same as above

    """
    
    if atm =='H':
        
        alpha1 = alpha_plus('H',Tk(t))
        
        epsilon1 =  epsilon_plus(alpha1)
        
        epsilon2 = epsilon_k('H',n,method,h)
        
        deltaA1 = delta_A(Pre[0],epsilon1,alpha1)
        
        d1 = h*(Pre[0]-deltaA1)-(1+Pre[0]*10**(-3))*(epsilon2+epsilon1/alpha1)
        
        d2 = 1-h+epsilon2*10**(-3)
        
        d = d1/d2
        
    elif atm=='O':
        
        alpha1 = alpha_plus('O',Tk(t))
        
        epsilon1 = epsilon_plus(alpha1)
        
        epsilon2 = epsilon_k('O',n,method,h)
        
        deltaA1 = delta_A(Pre[1],epsilon1,alpha1)
        
        d1 = h*(Pre[1]-deltaA1)-(1+Pre[1]*10**(-3))*(epsilon2+epsilon1/alpha1)
        
        
        d2 = 1-h+epsilon2*10**(-3)
        d = d1/d2
    else:
        d = 'None'
        
    return d
        
# dX('H',[-38,-6],0.75,20,0.75,1)

def slope_line(Pre,h,t,n,method):
    """
    slope = dX_H / dX_O

    Pre : e.g. [-36,-6]
    h   : e.g. 0.75
    t   : e.g. 20
    n   : e.g. 1
    method :  e.g. 1

    """
    d1 = dX('H',Pre,h,t,n,method)

    d2 = dX('O',Pre,h,t,n,method)

    return d1/d2
# slope_line([-38,-6],0.75,20,0.75,1)
