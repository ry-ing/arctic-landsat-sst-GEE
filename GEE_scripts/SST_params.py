// ============================================================================
// SST_params.js
// Algorithm coefficients and RobustScaler parameters.
//
// Source: MASTER_ALGORITHM_PARAMETERS.csv from Arctic_landsat_SST_cleanrepo
// Derived by: derive_SST-COEFFS-main.py
// Algorithm: SST_K = A*Tb_scaled + B*(Tb_scaled * TCWV_scaled) + Intercept + D
// where:
//   Tb_scaled   = (Tb_K   - Scaler_Tb_Center)   / Scaler_Tb_Scale   [RobustScaler]
//   TCWV_scaled = (TQV    - Scaler_TCWV_Center)  / Scaler_TCWV_Scale
//   A, B, Intercept: from libRadtran radiative-transfer simulations
//   D: bias correction from 5-fold cross-validation on EasyCORA in-situ data
//
// Note: L4 borrows L5's D (no L4 in-situ validation data available).
// Note: TCWV scaler is the same for all sensors (derived from MERRA-2 distribution).
// ============================================================================

// ============================================================================
// 1. REGRESSION COEFFICIENTS
// ============================================================================
exports.COEFFS = ee.Dictionary({

  'L4': {
    'A_Tb':         3.4823219125277167,
    'B_Interaction': -0.011247282638209392,   // coefficient for Tb_scaled * TCWV_scaled
    'Intercept':    275.4503361096741,
    'D':            1.0489051256597346         // borrowed from L5 (no L4 validation data)
  },

  'L5': {
    'A_Tb':         3.4582955936783155,
    'B_Interaction': -0.021736111011342183,
    'Intercept':    275.45102979666893,
    'D':            1.0489051256597346
  },

  'L7': {
    'A_Tb':         3.465105453797611,
    'B_Interaction': -0.016949389372123317,
    'Intercept':    275.45154866249226,
    'D':            0.7687188792617585
  },

  'L8': {
    'A_Tb':         3.5180006879356767,
    'B_Interaction': -0.0068141085635177205,
    'Intercept':    275.45720026324545,
    'D':            1.3281548649017547
  },

  'L9': {
    'A_Tb':         3.518000675498537,
    'B_Interaction': -0.007285746076529987,
    'Intercept':    275.4565997320123,
    'D':            1.1523923728819236
  }

});


// ============================================================================
// 2. ROBUSTSCALER PARAMETERS
// ============================================================================
// MEDIAN = RobustScaler center_ (median of training data)
// IQR    = RobustScaler scale_  (interquartile range of training data)
// Order: [Tb, TCWV]
// Source: Scaler_Tb_Center, Scaler_Tb_Scale, Scaler_TCWV_Center, Scaler_TCWV_Scale
//         from MASTER_ALGORITHM_PARAMETERS.csv
//
// TCWV parameters are identical across sensors (MERRA-2 TQV distribution
// does not depend on Landsat sensor).

exports.SCALER_PARAMS = ee.Dictionary({

  'L4': {
    'MEDIAN': [274.58207742866665,  10.624123096466064],
    'IQR':    [3.4076957065901183,   6.399625420570372]
  },

  'L5': {
    'MEDIAN': [274.5828486755939,   10.624123096466064],
    'IQR':    [3.371546533425999,    6.399625420570372]
  },

  'L7': {
    'MEDIAN': [274.3076460803507,   10.624123096466064],
    'IQR':    [3.364911523674209,    6.399625420570372]
  },

  'L8': {
    'MEDIAN': [274.7343404955682,   10.624123096466064],
    'IQR':    [3.4227572696189554,   6.399625420570372]
  },

  'L9': {
    'MEDIAN': [274.6956193993194,   10.624123096466064],
    'IQR':    [3.421352186794479,    6.399625420570372]
  }

});
