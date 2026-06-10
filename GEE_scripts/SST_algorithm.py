/**
 * SST_algorithm.js
 * GEE module: Landsat SST retrieval using the interaction-term algorithm.
 *
 * Algorithm (from derive_SST-COEFFS-main.py):
 *   Tb_scaled   = (Tb_K   - MEDIAN[0]) / IQR[0]     [RobustScaler on Tb]
 *   TCWV_scaled = (TCWV   - MEDIAN[1]) / IQR[1]     [RobustScaler on TCWV]
 *   SST_K       = A * Tb_scaled + B * (Tb_scaled * TCWV_scaled) + Intercept + D
 *   SST_C       = SST_K - 273.15
 *
 * B is the interaction-term coefficient (Coeff_B_Interaction in CSV).
 * B is NEGATIVE for all sensors — MERRA-2 TQV is positively correlated
 * with warm humid conditions that cause slight overestimation, so the
 * interaction term applies a small downward correction when both Tb and
 * TCWV are large.
 *
 * Output bands added to image:
 *   'SST'  — sea surface temperature (°C), all pixels
 *   'bt_K' — brightness temperature (K) of the primary thermal band,
 *             retained for downstream extraction and quality control.
 */

var SST_params  = require('users/ryaning/Landsat_plume_SST:ArcticLandsatSST/SST_params');
var COEFFS      = SST_params.COEFFS;
var SCALER_PARAMS = SST_params.SCALER_PARAMS;


// ============================================================================
// Internal helpers
// ============================================================================

/**
 * Apply RobustScaler:  z = (x - median) / IQR
 * Input image must have bands ['Tb', 'TCWV'].
 */
var _applyScaling = function(inputFeatures, landsat_id) {
  var params = ee.Dictionary(SCALER_PARAMS.get(landsat_id));
  var medians = ee.List(params.get('MEDIAN'));
  var iqrs    = ee.List(params.get('IQR'));

  var median_img = ee.Image.constant(medians).rename(['Tb', 'TCWV']);
  var iqr_img    = ee.Image.constant(iqrs).rename(['Tb', 'TCWV']);

  // (x - median) / IQR  — epsilon prevents divide-by-zero
  return inputFeatures.subtract(median_img)
                      .divide(iqr_img.add(1e-9));
};


/**
 * Apply the interaction-term regression.
 * scaledFeatures must have bands ['Tb', 'TCWV'] (after RobustScaler).
 * Returns SST in Kelvin as a single-band image named 'SST_Kelvin'.
 */
var _calculateSst = function(scaledFeatures, landsat_id) {
  var coeffs = ee.Dictionary(COEFFS.get(landsat_id));

  var tb_s   = scaledFeatures.select('Tb');
  var tcwv_s = scaledFeatures.select('TCWV');

  // Interaction term: Tb_scaled * TCWV_scaled
  var interaction = tb_s.multiply(tcwv_s).rename('interaction');

  // SST_K = Intercept + A * Tb_s + B * (Tb_s * TCWV_s) + D
  var sst_kelvin = scaledFeatures.addBands(interaction).expression(
    'intercept + (A * Tb) + (B * interaction) + D', {
      'intercept':   ee.Number(coeffs.get('Intercept')),
      'A':           ee.Number(coeffs.get('A_Tb')),
      'B':           ee.Number(coeffs.get('B_Interaction')),
      'D':           ee.Number(coeffs.get('D')),
      'Tb':          tb_s,
      'interaction': interaction
    }
  );

  return sst_kelvin.rename('SST_Kelvin');
};


// ============================================================================
// Exported function
// ============================================================================

/**
 * Returns a mapping function that adds 'SST' (°C) and 'bt_K' bands to a
 * Landsat image. To be used as: collection.map(SST_algorithm.addBand(landsat)).
 *
 * The image must already have the relevant TOA thermal band and a 'TCWV' band
 * (added by get_MERRA.addBand).
 *
 * @param {string} landsat  Sensor ID: 'L4', 'L5', 'L7', 'L8', or 'L9'.
 */
exports.addBand = function(landsat) {

  return function(image) {

    // 1. Select primary thermal band (brightness temperature in K from T1_TOA)
    var tir_band = ee.String(
      ee.Algorithms.If(landsat === 'L9', 'B10',
      ee.Algorithms.If(landsat === 'L8', 'B10',
      ee.Algorithms.If(landsat === 'L7', 'B6_VCID_1',
      'B6')))  // L4, L5
    );

    var tb = image.select(tir_band).rename('Tb');   // Tb in Kelvin
    var tcwv = image.select('TCWV');

    // 2. Bundle and scale
    var inputFeatures  = tb.addBands(tcwv);
    var scaledFeatures = _applyScaling(inputFeatures, landsat);

    // 3. Compute SST
    var sst_kelvin  = _calculateSst(scaledFeatures, landsat);
    var sst_celsius = sst_kelvin.subtract(273.15).rename('SST');

    // 4. Add both SST (°C) and bt_K (K) to the original image.
    //    bt_K is the unscaled brightness temperature, kept for QC and export.
    return image
      .addBands(sst_celsius)
      .addBands(tb.rename('bt_K'));
  };

};
