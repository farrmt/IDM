#---------------------------------------------#
#----Formatting of spatial data for IDM------#
#----Includes: PO data; DS data; bias cov;----#
#----intensity covs---------------------------#
#----Created by Matthew Farr------------------#
#---------------------------------------------#

from __future__ import division
import arcpy, os
#from arcpy import env
#from arcpy.sa import *
import numpy as np

#Set working enviornment
arcpy.env.workspace = "~\\DataFormatting\\"

#Buffer 500m around MMNR (include pixels that overlap border)
arcpy.Buffer_analysis("Shapefiles\\MMNR.shp", "Shapefiles\\MMNR_B.shp", 500, "FULL", "ROUND", "", "")

#Set working extent
arcpy.env.extent = "Shapefiles\\MMNR_B.shp"

#Convert extent into raster file
arcpy.PolygonToRaster_conversion("Shapefiles\\MMNR_B.shp", "NAME", "Shapefiles\\MMNR.tif", "", "", 1000)

#Set snap raster
arcpy.env.snapRaster = "Shapefiles\\MMNR.tif"

#Activate extensions
arcpy.CheckOutExtension("Spatial")

#Set raster mask
arcpy.env.mask = "Shapefiles\\MMNR_B.shp"

#-----------#
#-Bias data-#
#-----------#

#Create bias shapefile
out = "Shapefiles\\"
bias = "RawData\\Bias.csv"
arcpy.MakeXYEventLayer_management(bias, in_x_field = "UTME",
                                  in_y_field = "UTMN", out_layer = "lay1",
                                  spatial_reference = arcpy.SpatialReference("WGS 1984 UTM Zone 36S"))
arcpy.FeatureClassToFeatureClass_conversion("lay1", out, "Bias.shp")

#Create kernel density for extent
arcpy.env.workspace = os.getcwd()
kd = KernelDensity("lay1", "NONE", 1000, 2000, "SQUARE_KILOMETERS")
file = "Shapefiles\\BiasKD.tif"
kd.save(file)

#Convert raster to numpy array
kdn = arcpy.RasterToNumPyArray(file)

#95% cutoff of values that are greater than zero
percent = (1 - (np.count_nonzero(kdn>0)/(kdn.shape[0]*kdn.shape[1])*0.95))*100
cutoff = np.percentile(kdn, percent)

#Reclassify raster to 0 and 1
lay2 = LessThan(file, cutoff)
lay3 = SetNull(lay2, "Shapefiles\\MMNR.tif", "Value = 1")
file1 = "Shapefiles\\Extent.tif"
lay3.save(file1)

#Convert raster extent to polygon extent
file = "Shapefiles\\Extent.shp"
arcpy.RasterToPolygon_conversion("Shapefiles\\Extent.tif", file, "NO_SIMPLIFY")

#Buffer 650m around DS transects
arcpy.Buffer_analysis("Shapefiles\\Transect.shp", "Shapefiles\\TransectBuffer.shp", 650, "FULL", "FLAT", "", "")

#Merge extent of PO with extent of DS
arcpy.Merge_management(["Shapefiles\\TransectBuffer.shp", "Shapefiles\\Extent.shp"], "Shapefiles\\ExtentRaw.shp")

#MANUALLY EDIT EXTENT RAW IN ARCMAP & ARCSCAN (save as StudyArea.tif)

#Convert StudyArea.tiff to shp for extent and mask
arcpy.RasterToPolygon_conversion("Shapefiles\\StudyArea.tif", "Shapefiles\\StudyArea.shp", "NO_SIMPLIFY")

#Set working extent
arcpy.env.extent = "Shapefiles\\StudyArea.shp"

#Set raster mask
arcpy.env.mask = "Shapefiles\\StudyArea.shp"

#Set snap raster
arcpy.env.snapRaster = "Shapefiles\\StudyArea.tif"

#Set output directory
out1 = "Shapefiles\\ObsBias\\"

#Define cursor for loop
def unique_values(table, field):
    with arcpy.da.SearchCursor(table, [field]) as cursor:
        return sorted({row[0] for row in cursor})

#Extract months and years
month = [7, 8, 9, 10, 11, 12], [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], [1, 2, 3]
year = unique_values(r"lay1", "Year")

for yr in year:
    for mn in month[yr - year[0]]:
        fc = out + "Bias.shp"
        expression = "Month = " + str(mn) + " AND Year = " + str(yr)
        arcpy.MakeFeatureLayer_management(fc, "lay4")
        arcpy.SelectLayerByAttribute_management("lay4", "NEW_SELECTION", expression)
        kde = KernelDensity("lay4", "NONE", 1000, 2000, "SQUARE_KILOMETERS")
        file2 = out1 + "Bias" + str(mn) + str(yr) + ".tif"
        kde.save(file2)
        arcpy.Delete_management("lay4")

#---------#
#-PO data-#
#---------#

#Make PO data into shapefile
arcpy.MakeXYEventLayer_management("RawData\\PO_BBJ.csv", in_x_field = "Easting",
                                  in_y_field = "Northing", out_layer = "lay5",
                                  spatial_reference = arcpy.SpatialReference("WGS 1984 UTM Zone 36S"))

#Save PO shapefile
arcpy.FeatureClassToFeatureClass_conversion("lay5", out, "PO_BBJ.shp")

#Create PO raster for each month and year
for yr in year:
    for mn in month[yr - year[0]]:
        out3 = "Shapefiles\\BBJ\\"  + str(mn) + str(yr) + "PO.tif"
        fc = out + "PO_BBJ.shp"
        expression = "Month = " + str(mn) + " AND Year = " + str(yr)
        arcpy.MakeFeatureLayer_management(fc, "lay6")
        arcpy.SelectLayerByAttribute_management("lay6", "NEW_SELECTION", expression)
        arcpy.PointToRaster_conversion("lay6", "Count", out3, "MAXIMUM", "", 1000)
        arcpy.Delete_management("lay6")

#---------#
#-DS data-#
#---------#

#Convert to raster
arcpy.PolygonToRaster_conversion("Shapefiles\\TransectBuffer.shp", "FID", "Shapefiles\\Transect.tif", "", "", 50)

#MANUALLY CLIP RASTER EXTENT DOWN TO PIXELS INCLUDED IN REGION B: VIA ARCSCAN (DSExtent.tif)

#Convert to 50m pixels
arcpy.RasterToPolygon_conversion("Shapefiles\\DSExtent.tif", "Shapefiles\\DSpixel.shp", "NO_SIMPLIFY")
arcpy.PolygonToRaster_conversion("Shapefiles\\DSpixel.shp", "FID", "Shapefiles\\DSpixel.tif", "", "", 50)

#Make DS data into shapefile
arcpy.MakeXYEventLayer_management("RawData\\DS_BBJ.csv", in_x_field = "UTME",
                                  in_y_field = "UTMN", out_layer = "lay7",
                                  spatial_reference = arcpy.SpatialReference("WGS 1984 UTM Zone 36S"))

#Save DS shapefile
arcpy.FeatureClassToFeatureClass_conversion("lay7", out, "DS_BBJ.shp")

#Create DS raster for each month and year
for yr in year:
    for mn in month[yr - year[0]]:
        out4 = "Shapefiles\\BBJ\\" + str(mn) + str(yr) + "DS.tif"
        fc = out + "DS_BBJ.shp"
        expression = "Month = " + str(mn) + " AND Year = " + str(yr)
        arcpy.MakeFeatureLayer_management(fc, "lay8")
        arcpy.SelectLayerByAttribute_management("lay8", "NEW_SELECTION", expression)
        arcpy.PointToRaster_conversion("lay8", "Count", out4, "SUM", "", 50)
        arcpy.Delete_management("lay8")

#Set working extent
arcpy.env.extent = "Shapefiles\\Transect.tif"

#Set raster mask
arcpy.env.mask = "Shapefiles\\Transect.tif"

#Set snap raster
arcpy.env.snapRaster = "Shapefiles\\Transect.tif"

#Distance from transect to each DS pixel
dst = EucDistance("Shapefiles\\Transect.shp", "", 50, "Shapefiles\\Direction.tif")
dst.save("Shapefiles\\Distance.tif")

#----------------#
#-NDVI covariate-#
#----------------#

#Set working extent
arcpy.env.extent = "Shapefiles\\StudyArea.tif"

#Set snap raster
arcpy.env.snapRaster = "Shapefiles\\StudyArea.tif"

#Working enviornment
arcpy.env.workspace = "Shapefiles\\NDVI\\"

#List of TIF files (i.e., NDVI)
rasterlist = arcpy.ListRasters("*" "TIF")

#Working enviornment
arcpy.env.workspace = os.getcwd()

#Coordinate system
coord = arcpy.SpatialReference(32736)

#Extract by mask
for raster in rasterlist:
    inp = "Shapefiles\\NDVI\\" + raster
    arcpy.ProjectRaster_management(inp, "tmp1", coord)
    arcpy.Resample_management("tmp1", "tmp2", "1000", "CUBIC")
    out = ExtractByMask("tmp2", "Shapefiles\\StudyArea.tif")
    file = "Shapefiles\\NDVI\\formatted" + raster
    out.save(file)
    arcpy.Delete_management("tmp1")
    arcpy.Delete_management("tmp2")

#----------------#
#-Lion covariate-#
#----------------#

#Set raster mask
arcpy.env.mask = "Shapefiles\\StudyArea.tif"

#Make Lion data into shapefile
arcpy.MakeXYEventLayer_management("RawData\\PO_Lion.csv", in_x_field = "Easting",
                                  in_y_field = "Northing", out_layer = "lay9",
                                  spatial_reference = arcpy.SpatialReference("WGS 1984 UTM Zone 36S"))

#Save Lion shapefile
out = "Shapefiles\\"
arcpy.FeatureClassToFeatureClass_conversion("lay9", out, "PO_Lion.shp")

#Create PO raster for each month and year
fc = "Shapefiles\\PO_Lion.shp"
for yr in year:
    for mn in month[yr - year[0]]:
        expression = "Month = " + str(mn) + " AND Year = " + str(yr)
        arcpy.MakeFeatureLayer_management(fc, "lay10")
        arcpy.SelectLayerByAttribute_management("lay10", "NEW_SELECTION", expression)
        kde = KernelDensity("lay10", "NONE", 1000, 2000, "SQUARE_KILOMETERS")
        file = "Shapefiles\\Lion\\" + "KDE" + str(mn) + str(yr) + ".tif"
        kde.save(file)
        arcpy.Delete_management("lay10")
