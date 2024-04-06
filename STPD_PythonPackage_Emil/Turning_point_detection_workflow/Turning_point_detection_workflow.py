# Name : TurnPoint_workflow.py
# Description	: This script calls the Sequential Turning Point Detection (STPD) function
# and identifies turining paoints, creates plots with TPS, and update the TP_tablesl in the dtaababase

# by Mustafa Kemal Emil
# modified from : Ghaderpour, E.; Benedetta, A.; Bozzano, F.; Scarascia Mugnozza, G.; Mazzanti, P.
# A Fast and Robust Method for Detecting Trend Turning Points in InSAR Displacement
# Time Series, Computers & Geosciences, 2023.

# Contact : mustafakemal.emil@wmich.edu or mkemalemil@gmail.com

### TO DOs ###
## 

import numpy as np
import pandas as pd
import arcpy
np.set_printoptions(precision=4)
np.set_printoptions(suppress=True)
from datetime import datetime, timedelta

import csv
import numpy as np
from matplotlib import pyplot as plt
from STPD import STPD
from TPTR import TPTR

# Path to the geodatabase and the feature classes 
geodatabase_path = r'E:\InSARvice\WebGISauto\myProject\myProject.gdb'
#geodatabase_path = r"E:\InSARvice\WebGISauto\test\Qatar_InSAR.gdb"

#### the input could be Asc, Dsc or Vertical,  
## If all three is needed run it multiple times and select one each time as input feature class.

## Insar time series feature class names
## (1)'PSI_decomposed_vertical_MERGED_web'
## (2)'PSI_Asc_MERGED_web'
## (3)'PSI_Dsc_MERGED_web'

## Turn Point table name
## (1)'TP_table_vertical'
## (2)'TP_table_Asc'
## (3)'TP_table_Dsc'

point_FC_name = 'PSI_decomposed_vertical_MERGED_web'



########## The output table is selected based on the input feature class
if point_FC_name == 'PSI_Asc_MERGED_web':
    TP_table_name = 'TP_table_Asc'
elif point_FC_name == 'PSI_Dsc_MERGED_web':
    TP_table_name = 'TP_table_Dsc'
elif point_FC_name == 'PSI_decomposed_vertical_MERGED_web':
    TP_table_name = 'TP_table_vertical'
else:
    print("Wrong FC, check if PS output names match with the filter")
     
# Construct the full path to the feature class and the 
TP_table_path = f'{geodatabase_path}\\{TP_table_name}'
point_FC_path = f'{geodatabase_path}\\{point_FC_name}'

# Field names to be updated in the turining point table 'TP_table_name'
TP_table_fields = ['PointOID', 'TP', 'TPDate', 'slope', 'direction', 'SNR', 'NDRI', 'RunDate']

# create a list of field names from the feature class: eg.  ['OBJECTID', 'Shape', 'velocity', 'Vprecision', 'DEM', 'Hprecision.....
field_names = [field.name for field in arcpy.ListFields(point_FC_path)]
print (field_names)

# create a list of dispalcement field names from the feature class: eg.'D_20230110', 'D_20230113',
selected_point_fields = [s for s in field_names if 'D_' in s]
print ('Selected field names for cumulative displacement time seriers : ' + str(selected_point_fields))

# create another list of fields to be used to grab Object id of each point feture
essential_point_fields = ['OBJECTID','lon', 'lat', 'h'] + selected_point_fields

# create a list of field names for the date fields : e.g  '20230110', '20230113'
field_names_dates = [s.replace('D_', '') for s in [s for s in field_names if 'D_' in s]]
print (field_names_dates)

# create a list of datetime objects : e.g  datetime.datetime(2023, 1, 10, 0, 0), datetime.datetime(2023, 1, 13, 0, 0)
list_DateTime = [datetime.strptime(x,'%Y%m%d') for x in field_names_dates]
print (list_DateTime)

# find the earliest and the oldest date objects
latest = max(list_DateTime)
earliest = min(list_DateTime)
print ('Earliest date: ' + str(earliest)  + '  Oldest date: ' + str(latest))

# Calculte the duration between the earliest and the latest datest
time_delta = (latest - earliest).days
time_delta_years = time_delta/365.25 ## this is a rough estimate due to leap years
time_period = time_delta_years
print ("Time period between the first and the last Sentinel-1 image (years): " + str(time_delta_years))

#  create lists of days and years from the begining : e.g.  0, 3... // 0.0, 0.008213552361396304...
numericalDays = [(d - earliest).days for d in list_DateTime]
numericalYears = [x/365.25 for x in numericalDays] # rough estimate
print ("Days from the begining: " + str(numericalDays))
print ("Years from the begining: " + str(numericalYears))

# create a numpy array for the time epochs in years to be used in Sequential Turning Point Detection (STPD)
t = np.array(numericalYears)

# create a datetime ofject to store the time of ptocessing, added to turining point table
RunDate = datetime.now()

# delete all the rows in the TP table before appending new ones note run date is available
arcpy.DeleteRows_management(TP_table_path)

# Truining point analyses
# calculate TPs and save TP figures for each point in the feature class, update TP table for the FC
with arcpy.da.SearchCursor(point_FC_path, essential_point_fields) as cursor:
    for row in cursor:
        f = np.array(row[4:])
        t=t
        print ("printing time steps in years")
        print (t)
        print ("printing deformation values")
        print (f)        
        PointOID = int(row[0])
        print ("detecting TPs for the point (PintOID): " + str(PointOID))
        earliest = earliest
        
        # might need to optimize these parameters depending on the total time period
        #TPs = STPD(t,f, size=60, step=12, SNR=1, NDRI=0.3, dir_th=0, tp_th=1, margin=12, alpha=0.01) # defaults
        TPs = STPD(t,f, size=5, step=1, SNR=1, NDRI=0.1, dir_th=0, tp_th=1, margin=2, alpha=0.01) # for one year
        stats, y = TPTR(t, f, TPs)
        print (" TP stats : ") 
        print (stats)
        workspace = r'E:\InSARvice\WebGISauto\myProject\TurnPoints'
        output_name = str(PointOID)
        full_path_figure = f'{workspace}\\TPsOID_{output_name}.png'
        #==============================================================================================
        # Save plotf for each  time series and detected TPs
        plt.plot(t, f, '-ok', label = 'Time Series', linewidth=1, markersize=3)
        plt.plot(t, y, 'b', label = 'Linear trend')
        if len(TPs) > 0:
            plt.plot(t[np.array(TPs)], y[np.array(TPs)], 'db', label = 'Turning point')
        for k in range(len(TPs)):
            C = "SNR = "+ str(np.round(stats[k][3],2)) + "\nNDRI = "+ str(np.round(stats[k][4],2)) +"\nPointOID = " + str(PointOID) + "\nTPdate = " + list_DateTime[k].strftime("%m/%d/%Y")
            plt.text(t[TPs[k]], f[TPs[k]], C, bbox=dict(facecolor='b', alpha=0.1))
        plt.xlabel('Time (year)')    
        plt.ylabel('Displacement (mm)')
        plt.grid(True)
        # This is to re-write the x-axis value to values you like nicely:
        #plt.xticks([0, 1, 2, 3, 4, 5, 6, 7], ['2015','2016','2017','2018','2019','2020','2021','2022'])
        plt.legend(bbox_to_anchor = (0., 1.02, 1., .102), loc = 3, ncol= 4, 
                   mode = "expand", borderaxespad = 0)
        fig = plt.gcf()
        fig.set_size_inches(8, 4)    
        plt.savefig(full_path_figure,dpi=300, bbox_inches='tight')
        #plt.show()
        plt.clf()
        
        #=============================================================================================
        TP_info_list = []
        print ()
        for k in range(len(TPs)):
            print ("TP index for the point:" + str(PointOID))
            print (k)
            PointOID = PointOID
            TP = int(stats[k][0])
            TPdate = list_DateTime[TP]            
            slope = stats[k][1]
            direction = stats[k][2]
            SNR = stats[k][3]
            NDRI = stats[k][4]
            RunDate = RunDate
            TP_info = (PointOID, TP, TPdate, slope, direction, SNR, NDRI, RunDate)
            TP_info_list.append(TP_info)
            print ("list of TPs for the point (PintOID): " + str(PointOID))
            print (TP_info_list)
        
        cursor2 = arcpy.da.InsertCursor(TP_table_path, TP_table_fields)
        for TProw in TP_info_list:
            print (" TP values to be added to the table " + str(TProw))
            cursor2.insertRow(TProw)