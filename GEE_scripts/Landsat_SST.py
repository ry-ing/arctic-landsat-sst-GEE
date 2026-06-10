/**
 * Earth Engine Module for Landsat Sea Surface Temperature (SST) Retrieval.
 *
 * Implements the multivariate regression algorithm derived from
 * 'derive_COEFFS_MERRA.py', which uses Tb, TCWV, and TCWV^2.
 *
 * This algorithm relies on RobustScaler parameters (median and IQR).
 */

var get_MERRA = require('users/ryaning/Landsat_plume_SST:ArcticLandsatSST/get_MERRA');
var mask = require('users/ryaning/Landsat_plume_SST:ArcticLandsatSST/water_mask');
var SST_algorithm = require('users/ryaning/Landsat_plume_SST:ArcticLandsatSST/SST_algorithm');
//===============================================================
// - GEE COLLECTION DEFINITIONS 
// ============================================================================

var COLLECTION_INFO = ee.Dictionary({
  'L4': {
    'TOA': ee.ImageCollection('LANDSAT/LT04/C02/T1_TOA'),
    'SR': ee.ImageCollection('LANDSAT/LT04/C02/T1_L2'),
    'TIR': ['B6'], // Primary (and only) thermal band
    'SR_BANDS': ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7', 'QA_PIXEL'],
  },
  'L5': {
    'TOA': ee.ImageCollection('LANDSAT/LT05/C02/T1_TOA'),
    'SR': ee.ImageCollection('LANDSAT/LT05/C02/T1_L2'),
    'TIR': ['B6'], // Primary (and only) thermal band
    'SR_BANDS': ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7', 'QA_PIXEL'],
  },
  'L7': {
    'TOA': ee.ImageCollection('LANDSAT/LE07/C02/T1_TOA'),
    'SR': ee.ImageCollection('LANDSAT/LE07/C02/T1_L2'),
    'TIR': ['B6_VCID_1'], // Primary thermal band (B6_VCID_2 excluded; QA_PIXEL comes from SR)
    'SR_BANDS': ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7', 'QA_PIXEL'],
  },
  'L8': {
    'TOA': ee.ImageCollection('LANDSAT/LC08/C02/T1_TOA'),
    'SR': ee.ImageCollection('LANDSAT/LC08/C02/T1_L2'),
    'TIR': ['B10'], // Primary thermal band (B11 is NOT used)
    'SR_BANDS': ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7', 'QA_PIXEL'],
  },
  'L9': {
    'TOA': ee.ImageCollection('LANDSAT/LC09/C02/T1_TOA'),
    'SR': ee.ImageCollection('LANDSAT/LC09/C02/T1_L2'),
    'TIR': ['B10'], // Primary thermal band (B11 is NOT used)
    'SR_BANDS': ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7', 'QA_PIXEL'],
  }
});

var scale_SR = function(image){
  // 1. Select only the optical bands (Band 1 through Band 7)
  var optical = image.select('SR_B.*');
  
  // 2. Scale ONLY the optical bands
  var scaled = optical.multiply(0.0000275).add(-0.2);
  
  // 3. Overwrite the original optical bands with the scaled ones
  //    This leaves 'QA_PIXEL' untouched (as an Integer).
  return image.addBands(scaled, null, true);
};

exports.getCollection = function(landsat, start_date, end_date, geometry) {
  
  // --- 2. Get sensor-specific info ---
  var info = ee.Dictionary(COLLECTION_INFO.get(landsat));

  // --- 3. Load and Join SR and TOA collections for THIS sensor ---
  var landsat_SR = ee.ImageCollection(info.get('SR'))
      .filter(ee.Filter.date(start_date, end_date))
      .filterBounds(geometry)
      .map(get_MERRA.addBand)
      .map(scale_SR)
      

  var landsat_TOA = ee.ImageCollection(info.get('TOA'))
      .filter(ee.Filter.date(start_date, end_date))
      .filterBounds(geometry);
    //  .map(mask.toa);

  var tir = ee.List(info.get('TIR'));
  var vis = ee.List(info.get('SR_BANDS')).add('TCWV');

  var Landsat_merged = (landsat_SR.select(vis).combine(landsat_TOA.select(tir), true));
  
  var Landsat_SST = Landsat_merged.map(SST_algorithm.addBand(landsat));
  
  return Landsat_SST
};
  
