// ============================================================================
// gee_plume_extractor.js
//
// Generic Landsat plume SST extractor app.
// Implements the updated interaction-term algorithm via the LandsatSST module.
//
// Output filename: landsat_pst_{site-name}_{YYYY-MM-DD}_stats
// Output columns (match data/Landsat-plume-surface-temps_Store-glacier/ CSVs):
//   system:index, SST_mean, SST_median, SST_stdDev,
//   date, glacier, orig_thermal_band, sensor,
//   sst_masked_mean, sst_masked_median, sst_masked_stdDev,
//   bt_mean_K, bt_median_K, bt_stddev_K, pixel_count, tqv_kg_m2,
//   .geo
// ============================================================================

var LandsatSST = require('users/ryaning/Landsat_plume_SST:ArcticLandsatSST/Landsat_SST');
var palettes   = require('users/gena/packages:palettes');

Map.setOptions('SATELLITE');
Map.style().set('cursor', 'crosshair');
Map.setCenter(-50.6, 70.4, 5);

// ============================================================================
// App state
// ============================================================================
var appState = {
  roi:          null,
  imageList:    [],
  currentIndex: 0,
  currentImage: null,
  glacierName:  'store-glacier',
  sensorKey:    null    // 'L5', 'L8', etc. for the current image
};

var drawingTools = Map.drawingTools();

// Sensor label → SPACECRAFT_ID format used in the output CSVs
var SENSOR_TO_SPACECRAFT = {
  'L4': 'LANDSAT_4',
  'L5': 'LANDSAT_5',
  'L7': 'LANDSAT_7',
  'L8': 'LANDSAT_8',
  'L9': 'LANDSAT_9'
};

// Sensor label → primary thermal band name
var SENSOR_TO_THERMAL = {
  'L4': 'B6',
  'L5': 'B6',
  'L7': 'B6_VCID_1',
  'L8': 'B10',
  'L9': 'B10'
};

// QA_PIXEL water + quality mask matching the Python training pipeline.
// Keeps only water pixels (bit 7) that pass all contamination checks.
function buildWaterMask(qa, sensorKey) {
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
  return isWater.and(goodQA);
}

// ============================================================================
// UI
// ============================================================================
var panel = ui.Panel({style: {width: '320px', padding: '10px'}});
ui.root.insert(0, panel);

panel.add(ui.Label({
  value: 'Landsat Plume SST Extractor',
  style: {fontSize: '18px', fontWeight: 'bold', color: '#1a73e8'}
}));

panel.add(ui.Label('1. Define Study Area', {fontWeight: 'bold', margin: '15px 0 5px 0'}));

var nameBox = ui.Textbox({
  placeholder: 'site name, e.g. store-glacier',
  value: 'store-glacier',
  onChange: function(text) { appState.glacierName = text.toLowerCase().replace(/\s+/g, '-'); }
});
panel.add(ui.Label('Site name (used in filename):'));
panel.add(nameBox);

var roiLabel = ui.Label('Draw a rectangle around your study area.', {color: 'gray'});
panel.add(roiLabel);

var setRoiBtn = ui.Button({label: 'Set Study Area & Find Images', onClick: findImages});
panel.add(setRoiBtn);

panel.add(ui.Label('2. Navigate Images', {fontWeight: 'bold', margin: '15px 0 5px 0'}));
var navPanel = ui.Panel({layout: ui.Panel.Layout.flow('horizontal')});
var prevBtn   = ui.Button('← Prev', function() { updateImage(-1); }, true);
var nextBtn   = ui.Button('Next →', function() { updateImage(1);  }, true);
navPanel.add(prevBtn);
navPanel.add(nextBtn);
panel.add(navPanel);

var countLabel = ui.Label('No images loaded.');
panel.add(countLabel);

panel.add(ui.Label('3. Draw Plume & Extract', {fontWeight: 'bold', margin: '15px 0 5px 0'}));
panel.add(ui.Label('Draw POLYGON(S) over the plume, then click below.'));

var downloadBtn = ui.Button({label: 'Download CSV (Plume Stats)', disabled: true, onClick: downloadStats});
panel.add(downloadBtn);

var msgLabel = ui.Label('');
panel.add(msgLabel);

// ============================================================================
// ROI / search
// ============================================================================
function resetDrawingTools() {
  drawingTools.clear();
  drawingTools.setShown(true);
  drawingTools.setDrawModes(['rectangle']);
  drawingTools.setShape('rectangle');
  drawingTools.addLayer([], 'StudyArea', 'red');
}
resetDrawingTools();

function findImages() {
  var layers = drawingTools.layers();
  if (layers.length() === 0 || layers.get(0).geometries().length() === 0) {
    roiLabel.setValue('Please draw a rectangle first.');
    return;
  }
  appState.roi = layers.get(0).geometries().get(0);
  roiLabel.setValue('Searching Landsat archive...');
  setRoiBtn.setDisabled(true);

  var merged = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2')
    .map(function(i){ return i.set('sensorKey', 'L9'); })
    .merge(ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').map(function(i){ return i.set('sensorKey', 'L8'); }))
    .merge(ee.ImageCollection('LANDSAT/LE07/C02/T1_L2').map(function(i){ return i.set('sensorKey', 'L7'); }))
    .merge(ee.ImageCollection('LANDSAT/LT05/C02/T1_L2').map(function(i){ return i.set('sensorKey', 'L5'); }))
    .merge(ee.ImageCollection('LANDSAT/LT04/C02/T1_L2').map(function(i){ return i.set('sensorKey', 'L4'); }))
    .filterBounds(appState.roi)
    .filter(ee.Filter.calendarRange(6, 9, 'month'))
    .filter(ee.Filter.lt('CLOUD_COVER', 30))
    .map(function(img) { return img.set('dateStr', img.date().format('YYYY-MM-dd')); });

  merged.reduceColumns({
    reducer: ee.Reducer.toList().repeat(3),
    selectors: ['dateStr', 'sensorKey', 'system:time_start']
  }).evaluate(function(result, err) {
    if (err) {
      roiLabel.setValue('Search error — check Console.');
      print('Archive search error:', err);
      setRoiBtn.setDisabled(false);
      return;
    }
    if (!result || !result.list || result.list[0].length === 0) {
      roiLabel.setValue('No images found (Jun–Sep, <30% cloud). Try a larger area.');
      setRoiBtn.setDisabled(false);
      return;
    }

    var dates     = result.list[0];
    var sensors   = result.list[1];
    var timestamps = result.list[2];

    var combined = [];
    for (var i = 0; i < dates.length; i++) {
      combined.push({date: dates[i], sensor: sensors[i], ts: timestamps[i]});
    }
    combined.sort(function(a, b) { return a.ts - b.ts; });

    appState.imageList    = combined;
    appState.currentIndex = 0;

    roiLabel.setValue('Found ' + combined.length + ' images.');
    countLabel.setValue('1 / ' + combined.length);
    prevBtn.setDisabled(true);
    nextBtn.setDisabled(combined.length === 1);

    drawingTools.clear();
    drawingTools.setShape('polygon');
    drawingTools.setDrawModes(['polygon']);
    drawingTools.addLayer([], 'PlumePolygons', 'purple');

    loadImage(0);
  });
}

// ============================================================================
// Image display
// ============================================================================
function updateImage(direction) {
  var newIndex = appState.currentIndex + direction;
  if (newIndex >= 0 && newIndex < appState.imageList.length) {
    appState.currentIndex = newIndex;
    loadImage(newIndex);
    prevBtn.setDisabled(newIndex === 0);
    nextBtn.setDisabled(newIndex === appState.imageList.length - 1);
    countLabel.setValue((newIndex + 1) + ' / ' + appState.imageList.length);
  }
}

function getSSTCollection(sensorKey, dateStr, geom) {
  var start = dateStr;
  var end   = ee.Date(dateStr).advance(1, 'day');
  return LandsatSST.getCollection(sensorKey, start, end, geom);
}

function loadImage(index) {
  var data = appState.imageList[index];
  msgLabel.setValue('Loading ' + data.date + ' (' + data.sensor + ')...');
  downloadBtn.setDisabled(true);

  // Clear previous plume drawings
  var polyLayer = drawingTools.layers().get(0);
  while (polyLayer.geometries().length() > 0) {
    polyLayer.geometries().remove(polyLayer.geometries().get(0));
  }
  Map.layers().reset();

  var col = getSSTCollection(data.sensor, data.date, appState.roi);

  col.size().evaluate(function(count) {
    if (count === 0) {
      msgLabel.setValue('Error processing ' + data.date + ' — no image returned.');
      return;
    }

    var img = col.first().clip(appState.roi);

    // True-colour RGB (scaled SR bands)
    var isModern = data.sensor === 'L8' || data.sensor === 'L9';
    var rgbBands = isModern
      ? ['SR_B4', 'SR_B3', 'SR_B2']
      : ['SR_B3', 'SR_B2', 'SR_B1'];
    var rgb = img.select(rgbBands, ['Red', 'Green', 'Blue']);
    var hsv = rgb.rgbToHsv();
    var enhanced = ee.Image.cat([
      hsv.select('hue'),
      hsv.select('saturation').multiply(1.6),
      hsv.select('value')
    ]).hsvToRgb();

    // Water-masked SST
    var qa = img.select('QA_PIXEL');
    var sstMasked = img.select('SST').updateMask(buildWaterMask(qa, data.sensor)).rename('sst_masked');

    // Store image and sensor key for download
    appState.currentImage = img.addBands(sstMasked);
    appState.sensorKey    = data.sensor;

    Map.addLayer(enhanced, {min: 0.0, max: 0.6, gamma: 1.3}, 'True Colour');
    Map.addLayer(sstMasked, {
      min: -2, max: 5,
      palette: palettes.colorbrewer.RdBu[9].reverse()
    }, 'SST — water pixels (°C)');

    msgLabel.setValue('Loaded: ' + data.date + ' (' + data.sensor + ')');
    downloadBtn.setDisabled(false);
  });
}

// ============================================================================
// Download / extraction
// ============================================================================

function downloadStats() {
  var layer      = drawingTools.layers().get(0);
  var geometries = layer.geometries();

  if (geometries.length() === 0) {
    msgLabel.setValue('Error: no plume polygon drawn.');
    return;
  }

  msgLabel.setValue('Calculating stats...');

  var data      = appState.imageList[appState.currentIndex];
  var dateStr   = data.date;
  var sensorKey = appState.sensorKey || data.sensor;

  // Metadata for the output feature
  var spacecraft_id     = SENSOR_TO_SPACECRAFT[sensorKey] || sensorKey;
  var orig_thermal_band = SENSOR_TO_THERMAL[sensorKey] || 'B6';

  // Output filename: landsat_pst_{site}_{date}_stats
  var siteName = appState.glacierName.replace(/[^a-zA-Z0-9\-]/g, '-');
  var filename  = 'landsat_pst_' + siteName + '_' + dateStr + '_stats';

  // Polygon geometry (same geometry the JS app drew — stored as .geo in CSV)
  var region = ee.Geometry(layer.toGeometry());

  // Water + quality mask matching the Python training pipeline
  var qa         = appState.currentImage.select('QA_PIXEL');
  var waterMask  = buildWaterMask(qa, sensorKey);

  // Build a multi-band image with descriptive band names for clean reducer output:
  //   SST           → SST_mean / SST_median / SST_stdDev
  //   sst_masked    → sst_masked_mean / sst_masked_median / sst_masked_stdDev
  //   bt_K_masked   → used for bt_mean_K / bt_median_K / bt_stddev_K / pixel_count
  //   TCWV          → tqv_kg_m2 (mean over polygon, MERRA-2 is spatially coarse)
  var sst_all    = appState.currentImage.select('SST');
  var sst_masked = sst_all.updateMask(waterMask).rename('sst_masked');
  var bt_masked  = appState.currentImage.select('bt_K').updateMask(waterMask);
  var tqv        = appState.currentImage.select('TCWV');

  var extractImg = sst_all
    .addBands(sst_masked)
    .addBands(bt_masked)    // bt_K
    .addBands(tqv);         // TCWV

  var reducer = ee.Reducer.mean()
    .combine({reducer2: ee.Reducer.median(), sharedInputs: true})
    .combine({reducer2: ee.Reducer.stdDev(), sharedInputs: true})
    .combine({reducer2: ee.Reducer.count(),  sharedInputs: true});

  var fc = ee.FeatureCollection([ee.Feature(region)]);

  var stats = extractImg.reduceRegions({
    collection: fc,
    reducer:    reducer,
    scale:      30          // Landsat native resolution
  });

  // Rename reducer outputs to final CSV column names
  var output = stats.map(function(f) {
    return ee.Feature(f.geometry(), {
      // SST over all pixels in the polygon
      'SST_mean':           f.get('SST_mean'),
      'SST_median':         f.get('SST_median'),
      'SST_stdDev':         f.get('SST_stdDev'),
      // SST over water pixels only (bit 7 masked)
      'sst_masked_mean':    f.get('sst_masked_mean'),
      'sst_masked_median':  f.get('sst_masked_median'),
      'sst_masked_stdDev':  f.get('sst_masked_stdDev'),
      // Brightness temperature of water pixels (K) — for algorithm verification
      'bt_mean_K':          f.get('bt_K_mean'),
      'bt_median_K':        f.get('bt_K_median'),
      'bt_stddev_K':        f.get('bt_K_stdDev'),
      // Number of valid (water) pixels contributing to the masked stats
      'pixel_count':        f.get('bt_K_count'),
      // MERRA-2 total precipitable water vapour (kg m⁻²) — mean over polygon
      'tqv_kg_m2':          f.get('TCWV_mean'),
      // Scene metadata
      'date':               dateStr,
      'glacier':            appState.glacierName,
      'sensor':             spacecraft_id,
      'orig_thermal_band':  orig_thermal_band
    });
  });

  output = ee.FeatureCollection(output);

  output.getDownloadURL({filename: filename}, function(url) {
    if (url) {
      msgLabel.setValue('Ready — click link below.');
      print('DOWNLOAD: ' + dateStr + ' (' + sensorKey + ')', url);
      panel.add(ui.Label({
        value:     'Download ' + dateStr,
        targetUrl: url,
        style:     {color: 'blue', fontWeight: 'bold', margin: '5px 0'}
      }));
    } else {
      msgLabel.setValue('Error: could not generate download link.');
    }
  });
}
