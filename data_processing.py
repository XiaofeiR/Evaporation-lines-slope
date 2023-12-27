# -*- coding: utf-8 -*-
# @Author    : Xiaofei Ren
# @FileName  : Data_prepare.py

# Run in Python Console of QGIS
# Acquire climatic data and isotopic data (Inverse Distance Weight) for calculating slope of spatial points


import pandas as pd
import time

T1 = time.time()


def disA2B(A, B):
    """
    Calculating distance between A and B

    A : [longitude,latitude]
    B : [longitude,latitude]

    return distance
    """

    lon1, lat1 = A[0], A[1]
    lon2, lat2 = B[0], B[1]

    point1 = QgsPointXY(lon1, lat1)
    point2 = QgsPointXY(lon2, lat2)

    distance = QgsDistanceArea()
    distance.setEllipsoid('WGS84')
    m = distance.measureLine(point1, point2)
    return m


def selectLayer(name):
    """
    Select layer

    name : Name of layer
    """
    layer = None
    for lyr in QgsProject.instance().mapLayers().values():
        if lyr.name() == name:
            layer = lyr
            break
    return layer


def find_near(index1, index2, A):
    """
    Obtain the coordinates of the four nearest meteorological and isotope points to point A, along with the distances between them and point A.

    index1 : Spatial index of meteorological points.  e.g.  QgsSpatialIndex(l1.getFeatures())
    index2 : Spatial index of isotope points
    A      :  location of A  [longitude,latitude]


    return  EAR_loc, C_iso_loc

    """
    try:
        l1.removeSelection()
    except:
        pass

    try:
        l2.removeSelection()
    except:
        pass

    num = 4  # define the number of searching points
    f0_x, f0_y = A[0], A[1]

    # --------------Search climate point
    nearest1 = index1.nearestNeighbor(QgsPointXY(f0_x, f0_y), num)
    l1.select(nearest1)

    selection1 = l1.selectedFeatures()
    EAR_loc = []
    for i1 in selection1:
        d = disA2B([f0_x, f0_y], [i1[2], i1[1]])
        EAR_loc.append([i1[1], i1[2], d])  # i1[1] 为 lat
    EAR_loc = pd.DataFrame(EAR_loc, columns=['lat', 'lon', 'dis'])

    EAR_loc['1/dis'] = 1 / EAR_loc['dis']
    a1 = EAR_loc['1/dis'].sum()
    EAR_loc['weight'] = 1 / (EAR_loc['dis'] * a1)     # Calculating distance weight
    # ----------------------------------

    # --------------------Search isotope point
    nearest2 = index2.nearestNeighbor(QgsPointXY(f0_x, f0_y), num)
    l2.select(nearest2)

    selection2 = l2.selectedFeatures()
    C_iso_loc = []
    for i2 in selection2:
        d = disA2B([f0_x, f0_y], [i2[1], i2[0]])  
        C_iso_loc.append([i2[0], i2[1], d])
    C_iso_loc = pd.DataFrame(C_iso_loc, columns=['lat', 'lon', 'dis'])

    C_iso_loc['1/dis'] = 1 / C_iso_loc['dis']
    a2 = C_iso_loc['1/dis'].sum()
    C_iso_loc['weight'] = 1 / (C_iso_loc['dis'] * a2)     # Calculating distance weight
    # -------------------------------------

    return EAR_loc, C_iso_loc


def match_EAR_values(df):
    """
    Matching meteorological data to meteorological points
    df is the output of function find_near(index1, index2, A) for EAR_loc
    """

    df_source = mean_EAR5.copy()

    valuelist = []
    for i in range(1, 13, 1):
        s1 = 't2m' + '_' + str(i)
        s2 = 'rh' + '_' + str(i)
        valuelist.append(s1)
        valuelist.append(s2)
    for i in valuelist:
        df[i] = ''
    # -------------------------
    df_columns = df.columns.tolist()

    for i in range(len(df)):
        # 第i行的点

        [lat1, lon1] = df.iloc[i][['lat', 'lon']].tolist()

        df_v = df_source.loc[(df_source['lat'] == lat1) & (df_source['lon'] == lon1)]

        for j in range(1, 13, 1):
            [t1, rh1] = df_v.loc[df_v['month'] == j][['t2m', 'rh']].values.flatten().tolist()

            s1 = 't2m' + '_' + str(j)
            s2 = 'rh' + '_' + str(j)

            s1_index = df_columns.index(s1)  # t2m
            s2_index = df_columns.index(s2)  # rh

            df.iloc[i, s1_index] = t1
            df.iloc[i, s2_index] = rh1
    return df


def match_ciso(df):
    """
    Matching isotope data to isotope points
    df is the output of function find_near(index1, index2, A) for C_iso_loc
    """
    df2 = df
    df2_source = iso_data.copy()

    lat2 = df2['lat'].tolist()
    lon2 = df2['lon'].tolist()

    df2_v = df2_source.loc[(df2_source['lat'].isin(lat2)) & (df2_source['lon'].isin(lon2))]

    df2 = pd.merge(df2, df2_v, on=['lon', 'lat'], how='inner')
    del df2['ele']

    return df2


def get_A_values(A):
    """
    Based on the data from nearby meteorological and isotope points around point A's location, calculate the meteorological and isotope data for point A.

    Using Inverse Distance Weight method

    A    :  location of A  [longitude,latitude]

    return Meteorological and isotope data for point A
    """
    e1, c1 = find_near(index1, index2, A)

    df1 = e1.copy()
    df1 = match_EAR_values(df1)

    df2 = c1.copy()
    df2 = match_ciso(df2)

    # Calculate climate data
    n1 = df1.columns.tolist().index('weight')

    cal_list1 = df1.columns.tolist()[n1 + 1:]

    for i in cal_list1:
        df1[i] = df1['weight'] * df1[i]
    real_df1 = df1[cal_list1].sum()

    # Calculate isotope data
    n2 = df2.columns.tolist().index('weight')

    cal_list2 = df2.columns.tolist()[n2 + 1:]

    for i in cal_list2:
        df2[i] = df2['weight'] * df2[i]
    real_df2 = df2[cal_list2].sum()

    # save the data
    data = []
    for i in range(1, 13, 1):
        r_t1 = real_df1['t2m_' + str(i)]
        r_rh1 = real_df1['rh_' + str(i)]

        H2 = real_df2['hyd' + str(i)]
        O18 = real_df2['oxy' + str(i)]

        j = [A[0], A[1], i, r_t1, r_rh1, H2, O18]

        data.append(j)
    return pd.DataFrame(data, columns=['lon', 'lat', 'month', 't2m', 'rh', '2H', '18O'])


# Reading climate data
mean_EAR5 = pd.read_excel(r"D:\Data\climate_station\meanEAR5.xlsx")

#  Reading isotope data
iso_data = pd.read_csv(r"D:\Data\collection data\C-Iso-GuanZhong2.csv")


l0 = selectLayer('EL_slope')         # Slope layer
l1 = selectLayer('EAR5')            # Climate data layer
l2 = selectLayer('C_iso_GZ_site')   # isotope data layer



# Create spatial index
index1 = QgsSpatialIndex(l1.getFeatures())
index2 = QgsSpatialIndex(l2.getFeatures())


Data = []
for n in range(len(l0)):

    try:
        l0.removeSelection() # Cancel the last selected slope point
    except:
        pass

    i = n  # Select the n-th slope point

    l0.select(i)  # len(l0) = 372
    l0_feat = l0.selectedFeatures()[0]
    A = [l0_feat[1], l0_feat[2]]

    x1 = get_A_values(A)
    Data.append(x1)

    end = len(l0) - n + 1
    print('Remaining number：' + str(end))

Data = pd.concat(Data, axis=0)

print('Using time：' + str(time.time() - T1) + 'S')

# save the data to a table
Data.to_csv(r"D:\Data\EL_df.csv", index=False)


