// Water masking for Landsat Collection 2 QA_PIXEL.

// Universal water mask — returns a boolean mask image (1 = good water, 0 = masked).
// Applies all QA bits including cirrus. Safe for all sensors: cirrus bits (2, 14-15)
// are always 0 on L4/L5/L7, so those conditions are always true and have no effect.
exports.waterMask = function(image) {
  var qa = image.select('QA_PIXEL');
  var isWater = qa.bitwiseAnd(1 << 7).neq(0);
  var good = qa.bitwiseAnd(1 << 0).eq(0)          // no fill
    .and(qa.bitwiseAnd(1 << 1).eq(0))             // no dilated cloud
    .and(qa.bitwiseAnd(1 << 2).eq(0))             // no cirrus (0 for L4/L5/L7)
    .and(qa.bitwiseAnd(1 << 3).eq(0))             // no cloud
    .and(qa.bitwiseAnd(1 << 4).eq(0))             // no cloud shadow
    .and(qa.bitwiseAnd(1 << 5).eq(0))             // no snow
    .and(qa.rightShift(8).bitwiseAnd(3).lt(2))    // cloud confidence < 2
    .and(qa.rightShift(10).bitwiseAnd(3).lt(2))   // shadow confidence < 2
    .and(qa.rightShift(12).bitwiseAnd(3).lt(2))   // snow confidence < 2
    .and(qa.rightShift(14).bitwiseAnd(3).lt(2));  // cirrus confidence < 2 (0 for L4/L5/L7)
  return isWater.and(good);
};

// Full quality mask matching the Python training pipeline (get_INSITU-LANDSAT-PAIRS.py).
// Keeps only water pixels (bit 7) that pass all contamination checks.
// sensorKey: 'L4', 'L5', 'L7', 'L8', or 'L9' (cirrus flags only apply to L8/L9).
exports.sr = function(image, sensorKey) {
  var qa = image.select('QA_PIXEL');
  var isWater = qa.bitwiseAnd(1 << 7).neq(0);
  var goodQA = qa.bitwiseAnd(1 << 0).eq(0)          // no fill
    .and(qa.bitwiseAnd(1 << 1).eq(0))               // no dilated cloud
    .and(qa.bitwiseAnd(1 << 3).eq(0))               // no cloud
    .and(qa.bitwiseAnd(1 << 4).eq(0))               // no cloud shadow
    .and(qa.bitwiseAnd(1 << 5).eq(0))               // no snow
    .and(qa.rightShift(8).bitwiseAnd(3).lt(2))      // cloud confidence < 2
    .and(qa.rightShift(10).bitwiseAnd(3).lt(2))     // shadow confidence < 2
    .and(qa.rightShift(12).bitwiseAnd(3).lt(2));    // snow confidence < 2
  if (sensorKey === 'L8' || sensorKey === 'L9') {
    goodQA = goodQA
      .and(qa.bitwiseAnd(1 << 2).eq(0))             // no cirrus
      .and(qa.rightShift(14).bitwiseAnd(3).lt(2));  // cirrus confidence < 2
  }
  return image.updateMask(isWater.and(goodQA));
};

// Legacy bit-7-only mask (kept for backwards compatibility).
exports.toa = function(image) {
  var qa = image.select('QA_PIXEL');
  var is_water = qa.bitwiseAnd(1 << 7).neq(0);
  return image.updateMask(is_water);
};