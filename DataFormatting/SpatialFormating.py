#---------------------------------------------#
#----Formatting of spatial data for ISDM------#
#----Includes: PO data; DS data; bias cov;----#
#----intensity covs---------------------------#
#----Created by Matthew Farr------------------#
#---------------------------------------------#

import arcpy
from arcpy import env
from arcpy.sa import *

#Working enviornment within ArcGIS folder
arcpy.env.workspace = "C:\\Users\\farrm\Documents\\GitHub\\ISDM\\ArcGIS\\"

#Set working extent
arcpy.env.extent = "C:\\Users\\farrm\\Documents\\GitHub\\ISDM\\ArcGIS\\BaseData\\MMNR.shp"

#Convert extent into raster file
arcpy.PolygonToRaster_conversion("BaseData\\MMNR.shp", "NAME", "BaseData\\MMNR.tif", "", "", 500)

#Set snap raster
arcpy.env.snapRaster = "BaseData\\MMNR.tif"

#Activate extensions
arcpy.CheckOutExtension("Spatial")

#Set raster mask
arcpy.env.mask = "C:\\Users\\farrm\\Documents\\GitHub\\ISDM\\ArcGIS\\BaseData\\MMNR.shp"

#---------#
#-PO data-#
#---------#

def unique_values(table, field):
    with arcpy.da.SearchCursor(table, [field]) as cursor:
        return sorted({row[0] for row in cursor})

month = unique_values(r"BBJData\\PO_BBJ.shp", "Month")
    
year = unique_values(r"BBJData\\PO_BBJ.shp", "Year")

for yr in year:
    for mn in month:
        out = "BBJData\\PO_BBJ" + str(mn) + str(yr) + ".tif"
        expression = "Month = " + str(mn) + "AND Year = " + str(yr)
        arcpy.MakeFeatureLayer_management("BBJData\\PO_BBJ.shp", "lay")
        arcpy.SelectLayerByAttribute_management("lay", "NEW_SELECTION", expression)
        arcpy.PointToRaster_conversion("lay", "Count", out, "MAXIMUM", "", 500)
        arcpy.Delete_management("lay")

#---------#
#-DS data-#
#---------#

#-----------#
#-Bias data-#
#-----------#

out = "C:\\Users\\farrm\\Documents\\GitHub\\ISDM\\ArcGIS\\BiasData\\"
bias = "C:\\Users\\farrm\\Documents\\GitHub\\ISDM\\DataFormatting\\FormattedData\\Bias.csv"
arcpy.MakeXYEventLayer_management(bias, in_x_field = "UTME", in_y_field = "UTMN", out_layer = "lay1", spatial_reference = arcpy.SpatialReference("WGS 1984 UTM Zone 36S"))
arcpy.FeatureClassToFeatureClass_conversion("lay1", out, "Bias.shp")

for yr in year:
    for mn in month:
        fc = out + "Bias.shp"
        expression = "Month = " + str(mn) + "AND Year = " + str(yr)
        arcpy.MakeFeatureLayer_management(fc, "lay2")
        arcpy.SelectLayerByAttribute_management("lay2", "NEW_SELECTION", expression)
        kde = KernelDensity("lay2", "NONE", 500, 2000, "SQUARE_KILOMETERS")
        file = out + "Bias" + str(mn) + str(yr) + ".tif"
        kde.save(file)
        arcpy.Delete_management("lay2")


