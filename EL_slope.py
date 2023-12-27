
# -*- coding: utf-8 -*-
# @Author    : Xiaofei Ren
# @FileName  : EL_slope.py
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
# Calculate εk   

def epsilon_k(atom,n,method,h):
    """
    atom: 'H' or 'O'

    n: 0.5~1 ; 0.5 for open water bodies; 1 for soil water

    method : 0 or 1
        0 : Gibson et al.(2008)
        1 ：Benettin et al.(2018)

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

# refer to Benettin et al.(2018), Xiang et al.(2020), and Liu et al.(2022)

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



# ---------Reference---------

# Horita, J., Wesolowski, D.J., 1994. Liquid-vapor fractionation of oxygen and hydrogen isotopes of water from the freezing to the critical temperature. Geochimica et Cosmochimica Acta 58, 3425–3437. https://doi.org/10.1016/0016-7037(94)90096-5

# Gibson, J.J., Birks, S.J., Edwards, T.W.D., 2008. Global prediction of δA and δ2H-δ18O evaporation slopes for lakes and soil water accounting for seasonality. Global Biogeochemical Cycles 22. https://doi.org/10.1029/2007GB002997

# Benettin, P., Volkmann, T.H.M., von Freyberg, J., Frentress, J., Penna, D., Dawson, T.E., Kirchner, J.W., 2018. Effects of climatic seasonality on the isotopic composition of evaporating soil waters. Hydrology and Earth System Sciences 22, 2881–2890. https://doi.org/10.5194/hess-22-2881-2018

# Xiang, W., Evaristo, J., Li, Z., 2020. Recharge mechanisms of deep soil water revealed by water isotopes in deep loess deposits. Geoderma 369, 114321. https://doi.org/10.1016/j.geoderma.2020.114321

# Liu, H., Tang, J., Chen, L., Zhang, X., Zhu, B., Liang, C., Liu, C., Wang, G., 2022. Threshold recognition for shallow groundwater recharge by precipitation using dual isotopes in a small subtropical hilly catchment. CATENA 213, 106186. https://doi.org/10.1016/j.catena.2022.106186
