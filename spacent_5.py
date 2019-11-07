import ogr
import osr
import gdal
import numpy as np
import os
import glob

raster_paths = sorted(glob.glob('/media/venkanna/a08b01ef-2095-4d1d-94b4-d32b9c6f435f/ubuntu/data/AOI_8_Mumbai/PS-RGB'
                                '/*.tif'))
vector_paths = sorted(glob.glob('/media/venkanna/a08b01ef-2095-4d1d-94b4-d32b9c6f435f/ubuntu/data/AOI_8_Mumbai'
                                '/geojson_roads_speed/*.geojson'))
mask_directory = '/media/venkanna/a08b01ef-2095-4d1d-94b4-d32b9c6f435f/ubuntu/data/AOI_8_Mumbai/mask/'
shp_directory = '/media/venkanna/a08b01ef-2095-4d1d-94b4-d32b9c6f435f/ubuntu/data/AOI_8_Mumbai/shp/'
print(len(vector_paths), len(raster_paths))

shpdriver = ogr.GetDriverByName('ESRI Shapefile')

# Projections to add buffer in meters
geosr = osr.SpatialReference()
geosr.ImportFromEPSG(4326)
utmsr = osr.SpatialReference()
utmsr.ImportFromEPSG(32643)
geosr2utmsr = osr.CoordinateTransformation(geosr, utmsr)
utmsr2geosr = osr.CoordinateTransformation(utmsr, geosr)

for geojson, reference_image in zip(vector_paths, raster_paths):
    base = os.path.basename(geojson)
    file_name, ext = os.path.splitext(base)
    print('Rasterizing:  ' + file_name)
    test = ogr.Open(geojson, 0)
    lyrTest = test.GetLayer()
    v_proj = lyrTest.GetSpatialRef()

    # Buffer
    ds = shpdriver.CreateDataSource('buffer.shp')
    shp_filename = shp_directory + file_name + '.shp'
    ds_1 = shpdriver.CreateDataSource(shp_filename)
    lyrBuffer = ds.CreateLayer('buffer', geom_type=ogr.wkbPolygon, srs=utmsr)
    lyrBuffer1 = ds_1.CreateLayer(file_name, geom_type=ogr.wkbPolygon, srs=v_proj)
    featureDefn = lyrBuffer.GetLayerDefn()

    featureTest = lyrTest.GetNextFeature()
    while featureTest:
        geomTest = featureTest.GetGeometryRef()
        geomTest.Transform(geosr2utmsr)
        buffer_width = (float(featureTest.GetField('lanes')) * 1.5)
        geomBuffer = geomTest.Buffer(buffer_width)
        outFeature = ogr.Feature(featureDefn)
        outFeature.SetGeometry(geomBuffer)
        lyrBuffer.CreateFeature(outFeature)
        outFeature.Destroy()
        featureTest = lyrTest.GetNextFeature()

    buffer_featureTest = lyrBuffer.GetNextFeature()
    while buffer_featureTest:
        geomTest = buffer_featureTest.GetGeometryRef()
        geomTest.Transform(utmsr2geosr)
        outFeature = ogr.Feature(featureDefn)
        outFeature.SetGeometry(geomTest)
        lyrBuffer1.CreateFeature(outFeature)
        outFeature.Destroy()
        buffer_featureTest = lyrBuffer.GetNextFeature()

    # Rasterize

    # reference_image = 'SN5_roads_train_AOI_8_Mumbai_PS-RGB_chip0.tif'
    raster_path = mask_directory + file_name + '.tif'
    ref_image = gdal.Open(reference_image)
    gt = ref_image.GetGeoTransform()
    proj = ref_image.GetProjection()

    # 2) Creating the destination raster data source
    pixelWidth = pixelHeight = gt[1]  # depending how fine you want your raster ##COMMENT 1
    [cols, rows] = np.array(ref_image.GetRasterBand(1).ReadAsArray()).shape
    target_ds = gdal.GetDriverByName('GTiff').Create(raster_path, cols, rows, 1, gdal.GDT_Byte)
    target_ds.SetGeoTransform(gt)

    # 5) Adding a spatial reference
    target_ds.SetProjection(proj)
    band = target_ds.GetRasterBand(1)
    # gdal.RasterizeLayer(ds,bands,layer,burn_values, options = ["BURN_VALUE_FROM=Z"])
    gdal.RasterizeLayer(target_ds, [1], lyrBuffer1)
    # print(target_ds)
    target_ds = None

    # Deleting buffer shapefie with UTM projection
    os.remove("buffer.shp")
    os.remove("buffer.dbf")
    os.remove("buffer.shx")
    os.remove("buffer.prj")

