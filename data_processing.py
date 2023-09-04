# Run in Python Console of QGIS
# 获取计算slope空间点的气象数据和同位素数据 （反距离加权）


import pickle
import pandas as pd
import numpy as np
import time


T1 = time.time()

# 读取网格气象数据
mean_EAR5 = pd.read_excel(r"D:\我的坚果云\同位素\Data\climate_station\meanEAR5.xlsx")

# 读取网格点稳定同位素数据
iso_data= pd.read_csv(r"D:\我的坚果云\同位素\Data\collection data\C-Iso-GuanZhong2.csv")


def disA2B(A,B):
    # 计算A-B两点之间的距离
    # A = [lon,lat] 经度 纬度 
    lon1,lat1 = A[0],A[1]
    lon2,lat2 = B[0],B[1]
    
    point1 = QgsPointXY(lon1, lat1)
    point2 = QgsPointXY(lon2, lat2)

    distance = QgsDistanceArea()
    distance.setEllipsoid('WGS84')
    m = distance.measureLine(point1, point2)
    return m
    

try:
    l0.removeSelection()
except:
    pass
    

def selectLayer(name):
    
    # 根据名称选择图层
    layer = None
    for lyr in QgsProject.instance().mapLayers().values():
        if lyr.name() == name:
            layer = lyr
            break
    return layer

l1 = selectLayer('EAR5')
l0 = selectLayer('EL计算点')
l2 = selectLayer('C_iso_GZ_site')

# 空间索引
index1 = QgsSpatialIndex(l1.getFeatures()) # EAR5
index2 = QgsSpatialIndex(l2.getFeatures()) # C_iso



# 获取距离A点最近几个点的坐标及之间的距离
def find_near(index1,index2,A):    
    # 从空间索引中寻找距离A点最近的几个点
    # index1 为 EAR5
    # index2 为 C_iso
    try:
        l1.removeSelection()
    except:
        pass    

    try:
        l2.removeSelection()
    except:
        pass
    
    num = 4 # 搜索num个点
    
    f0_x,f0_y = A[0],A[1]


    # 寻找气象点的最近点
    nearest1 = index1.nearestNeighbor(QgsPointXY(f0_x,f0_y), num)
    l1.select(nearest1)
    # 选择点

    selection1 = l1.selectedFeatures()
    EAR_loc=[]
    for i1 in selection1:
        d = disA2B([f0_x,f0_y],[i1[2],i1[1]])
        EAR_loc.append([i1[1],i1[2],d]) # i1[1] 为 lat
    EAR_loc = pd.DataFrame(EAR_loc,columns= ['lat','lon','dis'])
    # 计算距离的倒数
    EAR_loc['1/dis'] = 1/EAR_loc['dis']
    # 通过对反距离权重公式的展开，发现只需要计算权重即可
    a1 = EAR_loc['1/dis'].sum()
    EAR_loc['weight'] = 1/(EAR_loc['dis']*a1)
    
    
    
        
    # 寻找C_iso的最近点
    nearest2 = index2.nearestNeighbor(QgsPointXY(f0_x,f0_y), num)
    l2.select(nearest2)

    selection2 = l2.selectedFeatures()
    C_iso_loc = []
    for i2 in selection2:
        d = disA2B([f0_x,f0_y],[i2[1],i2[0]]) # 计算EL计算点至C_iso的距离
        C_iso_loc.append([i2[0],i2[1],d])
    C_iso_loc = pd.DataFrame(C_iso_loc,columns = ['lat','lon','dis'])    
    C_iso_loc['1/dis'] = 1/C_iso_loc['dis']
    
    # 权重
    a2 = C_iso_loc['1/dis'].sum()
    C_iso_loc['weight'] = 1/(C_iso_loc['dis']*a2)
    

    return EAR_loc,C_iso_loc
    
    
    
# e1,c1 = find_near(index1,index2,A)



# 匹配附近点的数据 EAR5
def match_EAR_values(df):
    # 将df_source中的每个月的值匹配到df中
    df_source = mean_EAR5.copy()
    # 给df创建表头-------------
    valuelist = []
    for i in range(1,13,1):
        s1 = 't2m'+'_'+str(i)
        s2 = 'rh'+'_'+str(i)
        valuelist.append(s1)
        valuelist.append(s2)
    for i in valuelist:
        df[i]=''
    # -------------------------
    df_columns = df.columns.tolist()


    for i in range(len(df)):
        # 第i行的点
        
        [lat1,lon1] = df.iloc[i][['lat','lon']].tolist()
        
        df_v = df_source.loc[(df_source['lat']==lat1)&(df_source['lon']==lon1)]
        
        for j in range(1,13,1):
            
            [t1,rh1] = df_v.loc[df_v['month']==j][['t2m','rh']].values.flatten().tolist()
            
            
            s1 = 't2m'+'_'+str(j)
            s2 = 'rh'+'_'+str(j)
            
            s1_index = df_columns.index(s1) # t2m
            s2_index = df_columns.index(s2) # rh
            
            df.iloc[i,s1_index] = t1
            df.iloc[i,s2_index] = rh1
    return df
        
      

# 匹配附近点的数据 iso
def match_ciso(df):
    # 匹配iso到每一个点
    df2 = df
    df2_source = iso_data.copy()


    lat2 = df2['lat'].tolist()   
    lon2 = df2['lon'].tolist()

    df2_v = df2_source.loc[(df2_source['lat'].isin(lat2))&(df2_source['lon'].isin(lon2))]

    df2 = pd.merge(df2,df2_v,on=['lon','lat'],how='inner')
    del df2['ele']
    
    return df2
# df2 = c1.copy()
# match_ciso(df2)
    

def get_A_values(A):
    # 根据A点位置获取A点附近气象点和C_iso点的值
    # 并且根据反距离权重进行计算，获得其对应值
    e1,c1 = find_near(index1,index2,A)

    df1 = e1.copy()
    df1 = match_EAR_values(df1)  
        
    df2 = c1.copy()
    df2 = match_ciso(df2)

    list1 = ['lat', 'lon', 'dis', '1/dis', 'weight']

    # 开始赋权计算--EAR5
    n1 = df1.columns.tolist().index('weight') # df1需要赋权的列起始索引

    cal_list1 = df1.columns.tolist()[n1+1:]

    for i in cal_list1:
        df1[i] = df1['weight']*df1[i]
    real_df1 = df1[cal_list1].sum()

    # 计算 iso
    n2 = df2.columns.tolist().index('weight') # df1需要赋权的列起始索引

    cal_list2 = df2.columns.tolist()[n2+1:]

    for i in cal_list2:
        df2[i] = df2['weight']*df2[i]
    real_df2 = df2[cal_list2].sum()


    # 记录数据保存为df
    data = []
    for i in range(1,13,1):
        
        
        r_t1 = real_df1['t2m_'+str(i)]
        r_rh1 = real_df1['rh_'+str(i)]
        
        H2 = real_df2['hyd'+str(i)]
        O18 = real_df2['oxy'+str(i)]
        
        j = [A[0],A[1],i,r_t1,r_rh1,H2,O18]
        
        data.append(j)
    return pd.DataFrame(data,columns=['lon','lat','month','t2m','rh','2H','18O'])
    
    

Data = []
for n in range(len(l0)):
    
    # 取消上一次的选择
    try:
        l0.removeSelection()
    except:
        pass
    
    i = n # 选择第i个点

    l0.select(i) # len(l0) = 372
    l0_feat = l0.selectedFeatures()[0]
    A = [l0_feat[1] ,l0_feat[2]]  

    x1 = get_A_values(A)
    Data.append(x1)
    
    end = len(l0)-n+1
    print('还剩余：'+str(end))
    
Data = pd.concat(Data,axis = 0)

print('用时：'+str(time.time()-T1)+'秒')

# 保存slope计算点气象数据和同位素数据，为下一步计算slope做准备    
Data.to_csv(r"D:\我的坚果云\同位素\Data\EL计算点的数据表.csv",index=False)

