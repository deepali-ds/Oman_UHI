
// 1. REGION (HEX GRID)

var region = ee.FeatureCollection("projects/ee-deepalibidwai/assets/oman_h3_res6");
Map.centerObject(region, 6);


// 2. LANDSAT

var landsat = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
  .filterBounds(region)
  .filterDate("2022-01-01", "2025-12-31");

function prepL8(img) {
  return ee.Image.cat([
    img.normalizedDifference(["SR_B5","SR_B4"]).rename("NDVI"),
    img.normalizedDifference(["SR_B3","SR_B5"]).rename("NDWI"),
    img.normalizedDifference(["SR_B6","SR_B5"]).rename("NDBI"),
    img.select("ST_B10")
      .multiply(0.00341802).add(149.0).subtract(273.15)
      .rename("LST")
  ]);
}

var processed = landsat.map(prepL8);
var fallback = processed.mean();  // fallback if no data


// 3. ERA5 LAND 

var era5 = ee.ImageCollection("ECMWF/ERA5_LAND/MONTHLY_AGGR")
  .filterDate("2022-01-01", "2025-12-31");


// 4. TERRAIN (FOR MAR)

var dem = ee.Image("USGS/SRTMGL1_003");
var slope = ee.Terrain.slope(dem).rename("slope");


// 5. MONTH LOOP (48 months)

var months = ee.List.sequence(0,47);

var monthly = ee.ImageCollection.fromImages(
  months.map(function(i){

    var start = ee.Date("2022-01-01").advance(i,"month");
    var end = start.advance(1,"month");

    // Landsat safe composite
    var l8 = processed.filterDate(start,end);
    var sat = ee.Image(
      ee.Algorithms.If(l8.size().gt(0), l8.mean(), fallback)
    );

    // ERA5 monthly
    var clim = era5.filterDate(start,end).mean();

    var era = ee.Image.cat([
      clim.select("temperature_2m").rename("temp_air"),
      clim.select("total_precipitation_sum").rename("precip"),
      clim.select([
        "volumetric_soil_water_layer_1",
        "volumetric_soil_water_layer_2"
      ]).reduce(ee.Reducer.mean()).rename("soil_moist"),
      clim.expression(
        "sqrt(u*u+v*v)",
        {
          u: clim.select("u_component_of_wind_10m"),
          v: clim.select("v_component_of_wind_10m")
        }
      ).rename("wind_speed")
    ]).unmask();  //  prevents NaNs

    
    // UHI (THERMAL ANOMALY)
    
    var meanLST = sat.select("LST").reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: region,
      scale: 10000
    }).get("LST");

    var uhi = sat.select("LST")
      .subtract(ee.Number(meanLST))
      .rename("UHI");

    
    // MAR SUITABILITY MODEL 
    
    var mar_score =
      sat.select("NDWI").multiply(0.4)
      .add(era.select("soil_moist").multiply(0.3))
      .add(era.select("precip").multiply(0.2))
      .subtract(sat.select("NDBI").multiply(0.3))
      .subtract(slope.multiply(0.1))
      .rename("MAR_score");

    
    // FINAL IMAGE
    
    return sat
      .addBands(era)
      .addBands(uhi)
      .addBands(mar_score)
      .addBands(slope)
      .set("month", i)
      .clip(region);

  })
);

//  DEBUG CHECK
print("Bands check:", monthly.first().bandNames());


// 6. HEX AGGREGATION

var table = monthly.map(function(img){

  return img.reduceRegions({
    collection: region,
    reducer: ee.Reducer.mean(),
    scale: 10000,   // matches ERA5 resolution
    tileScale: 2
  }).map(function(f){
    return f.set("month", img.get("month"));
  });

}).flatten();


// 7. EXPORT

Export.table.toDrive({
  collection: table,
  description: "Oman_2022_2025_MAR_UHI",
  fileFormat: "CSV"
});