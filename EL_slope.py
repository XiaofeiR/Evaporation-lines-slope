
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
    摄氏度换算开尔文温度
    t 温度摄氏度
    """
    # 摄氏度换算为卡尔文
    return t+273.15

def alpha_plus(atom,T):
    """
    atom 原子 H  O
    T 开尔文温度
    """
    # 根据经验公式计算平衡分馏系数

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

# T = tk(10)
# alpha_plus('H',T)

#------------------------------------------------------------



# -----------------------------------------------------------
# 计算ε+   富集系数 单位是千分之
# he equilibrium isotopic separation between liquid and vapor is then computed

def epsilon_plus(alpha):
    return (alpha-1)*10**3

# epsilon_plus(alpha_2H(20))
#--------------------------------------------------------------



#--------------------------------------------------------------
# 计算εk   相对湿度的函数 存在两种计算形式 单位 千分之
# εk quantifies the isotopic effects during net evaporation that are associated
# with the higher diffusivities of isotopically lighter molecules.

# n 0.5,1
# 0.5 : open water table
# 1.0 : groundwater，soli water

# method  0,1
# 0 : Gibson(2008) GBC 这里存在一个问题就是C值 Haowen Liu （CATENA）中C值与原版相差10
# 1 ：Paolo Benettin (2018) HESS

def epsilon_k(atom,n,method,h):

    # atom原子类型，H ， O

    # n 0.5,1  也有的研究取0.75
    # 0.5 : open water table
    # 1 : groundwater

    # method  0,1
    # 0 : Gibson(2008) GBC 这里存在一个问题就是C值 Haowen Liu （CATENA）中C值与原版相差10
    # 1 ：Paolo Benettin (2018) HESS

    # 返回 正常的无量纲单位，已经删除了方法1中的10**3


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
#  δA 的计算 deltaP 输入为 千分之几的 几
# 计算公式为 δA=(δP-ε+)/α+
def delta_A(delta_P,epsilon,alpha):
    i = (delta_P-epsilon)/alpha
    return i




#------------------------------------------------------------------
# method 选1 与2018 Benettin 一致
# 此处delta_star 也是与 Benettin中的 dstar_H
def delta_star(atom,delP,h,t,n,method):

    alpA = alpha_plus(atom,Tk(t))
    epsK = epsilon_k(atom, n, method, h)
    epsA = epsilon_plus(alpA)

    delA = delta_A(delP,epsA,alpA)

    d1 = h * delA + epsK + epsA/alpA

    d2 = h-10**(-3)*(epsK+(epsA/alpA))

    return d1/d2
# -----------------------------------------------------------------





#-----------------------------------------------------------------
# 计算slope

# Pre 降雨的同位素组成 [H,O]  单位为 in each thousand
# h 相对湿度 0~1
# t 温度 摄氏度

# n 依据参考文献 开放水体为0.5， 土壤水 1
# reference by <Effects of climatic seasonality on the isotopic composition of evaporating soil waters>
# n的取值范围为 0.5~1
# 0.5 ： fully turbulent transport that reduces kinetic fractionation, appropriate for lakes or saturated soil conditions
# 1 ： fully diffusive transport, appropriate for very dry soil conditions




# method ：计算epsilon_k的方法，参照 epsilon_k(atom,n,method,h)、

def slope_line(Pre,h,t,n,method):
    """
    计算EL斜率，参照：Paolo Benettin 2018 Matlab 代码的计算过程
    仅代表了作者代码的前半部分，至于具体的计算参照作者后面部分
    Pre： 降雨同位素 [H,O]  单位 ‰
    h ： 相对湿度
    t ： 温度 摄氏度℃
    n ： 依据参考文献 开放水体取0.5  土壤水取 1
    method ： 计算 epsilon_K方法， 主要参照epsilon_K的过程
    """


    d1 = delta_star('H',Pre[0],h,t,n,method)-Pre[0]

    d2 = delta_star('O',Pre[1],h,t,n,method)-Pre[1]

    return d1/d2

print('参数依次为Pre,h,t,n,method| 开放水体n取0.5， 土壤水取1')




def dX(atm,Pre,h,t,n,method):
    """
    根据 论文 Recharge mechanisms of deep soil water revealed by water isotopes in deep loess deposits中的计算公式
    这个公式与 Haowen Liu (2022)中的公式相同
    这里dX值仅仅为计算slope的氢或者氧

    slope = dX_H / dX_O

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


def slope_line2(Pre,h,t,n,method):
    """
    计算 蒸发线的斜率 具体的公式的表述见dX
    Pre： 降雨同位素 [H,O]  单位 ‰
    h ： 相对湿度
    t ： 温度 摄氏度℃
    n ： 依据参考文献 开放水体取0.5  土壤水取 1
    method ： 计算 epsilon_K方法， 主要参照epsilon_K的过程
    """
    d1 = dX('H',Pre,h,t,n,method)

    d2 = dX('O',Pre,h,t,n,method)

    return d1/d2
# slope_line2([-38,-6],0.75,20,0.75,1)
# -------------------------------------------------------------



if __name__ == 'main':
    print('Benettin示例Matlab代码一致均为')
    print(slope_line([-38,-6],0.75,20,0.75,1))
    # 与Benettin计算结果一致，可以用于后续计算



