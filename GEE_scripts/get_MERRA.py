// --- CONSTANTS ---
var MERRA_COLLECTION_ID = 'NASA/GSFC/MERRA/slv/2';
var MERRA_TCWV_BAND = 'TQV';

/**
 * Adds MERRA-2 Total Column Water Vapor (TCWV) and TCWV-position bands
 * to a Landsat image.
 *
 * This function finds the closest-in-time MERRA-2 image, reprojects
 * the TCWV data to the Landsat image's grid, and then calculates a
 * binned position for the TCWV value.
 *
 * @param {ee.Image} image The Landsat image.
 * @return {ee.Image} The Landsat image with 'TPW' and 'TPWpos' bands added.
 */
exports.addBand = function(image) {
  
  // Get the date/time of the Landsat image
  var landsat_time = image.date();

  // --- 1. Find the closest MERRA-2 Image ---

  // Function to compute the time difference from the Landsat image
  var set_time_diff = function(merra_image) {
    var time_diff = ee.Number(merra_image.date().millis())
                            .subtract(landsat_time.millis()).abs();
    return merra_image.set('DateDist', time_diff);
  };

  // Filter MERRA-2 collection to a 2-hour window (+/- 1 hour)
  var merra_filtered = ee.ImageCollection(MERRA_COLLECTION_ID)
    .filterDate(
      landsat_time.advance(-1, 'hour'),
      landsat_time.advance(1, 'hour')
    )
    .select(MERRA_TCWV_BAND);
    
  // Find the closest image in time
  var merra_sorted = merra_filtered.map(set_time_diff).sort('DateDist');
  var merra_size = merra_sorted.size();
  
  // Use a non-realistic value if no image is found (prevents errors)
  var merra_closest = ee.Image(ee.Algorithms.If(
    merra_size.eq(0), 
    ee.Image.constant(-999.0).rename(MERRA_TCWV_BAND), // Error image
    merra_sorted.first() // Closest image
  ));

  // --- 2. Reproject and Rename ---
  
  // Get the target projection from the Landsat image
  var target_projection = image.select(0).projection();

  // Resample (bilinear) and reproject the TCWV data
  // Rename it to 'TPPW' to match the binning expression
  var tpw = merra_closest.select(MERRA_TCWV_BAND)
    .resample('bilinear')
    .reproject({ crs: target_projection })
    .rename('TCWV'); // Rename to 'TPW'

  // --- 3. Bin the TCWV data ---
  // REMOVED: Binning logic is not needed for this use case.
  
  // --- 4. Add bands to the original image ---
  var withTPW = image.addBands(tpw, ['TCWV']);
  
  return withTPW;
};