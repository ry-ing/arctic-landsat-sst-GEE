// Google Earth Engine Landsat-derived SST's - optimised for Greenland coastal waters
// Created by Ryan Ing, University of Edinburgh
// Email: ryan.ing@ed.ac.uk
//
// EDITOR VERSION — run this in the GEE Code Editor (not as a published App).
// Use sst_web_app_viewer.js for App deployment.
//
// ROI definition: click "Draw ROI", then click two corners on the map.
// Drawing tools layer API (layer.geometries / layer.toGeometry) is App-runtime-only
// and fails in Code Editor. This version uses Map.onClick instead.

// ============================================================================
// 1. IMPORTS AND CONFIGURATION
// ============================================================================

var LandsatSST = require('users/ryaning/arctic-landsat-sst:sst-retrieval/Landsat_SST');
var palettes = require('users/gena/packages:palettes');

var areaName = 'Landsat_SST_export';

Map.setOptions('SATELLITE');
Map.style().set('cursor', 'crosshair');

var currentImageToExport = null;

// ============================================================================
// 2. ROI STATE (two-click rectangle, no drawing tools layer API)
// ============================================================================

var capturedGeometry = null;
var roiCorner1 = null;
var roiMode = false;
var roiMapLayer = null;

function getUserGeometry() {
  return capturedGeometry;
}

// ============================================================================
// 3. UI SETUP
// ============================================================================

var panel = ui.Panel({style: {width: '420px', padding: '10px'}});
ui.root.insert(0, panel);

panel.add(ui.Label({
  value: 'Landsat-derived Arctic SST Viewer',
  style: {fontSize: '18px', fontWeight: 'bold', color: '#004d40', textAlign: 'center', margin: '1px 10px 1px 10px'}
}));
panel.add(ui.Label({
  value: 'Created by Ryan Ing, University of Edinburgh',
  style: {fontSize: '15px', fontWeight: 'italic', color: '#004d40', margin: '1px 10px 1px 10px', textAlign: 'center'}
}));
panel.add(ui.Label({
  value: 'Email: ryan.ing@ed.ac.uk',
  style: {fontSize: '15px', color: '#0096FF', textAlign: 'center', margin: '3px 10px 1px 10px'}
}));

// --- STEP 1: ROI ---
panel.add(ui.Label('1. Define Study Area', {fontWeight: 'bold', margin: '10px 0 5px 0'}));
var roiStatus = ui.Label(
  'Click "Draw ROI", then click two corners on the map to define a rectangle.',
  {fontSize: '12px', color: 'gray', margin: '1px 10px 5px 10px'}
);
panel.add(roiStatus);

var roiBtnRow = ui.Panel({layout: ui.Panel.Layout.flow('horizontal')});

var drawROIBtn = ui.Button({
  label: 'Draw ROI',
  style: {margin: '0 5px 0 0'},
  onClick: function() {
    roiCorner1 = null;
    capturedGeometry = null;
    roiMode = true;
    if (roiMapLayer) { Map.layers().remove(roiMapLayer); roiMapLayer = null; }
    roiStatus.setValue('Click the FIRST corner of your study area on the map.');
    msgLabel.setValue('Waiting for ROI...');
    dateSelect.items().reset([]);
    dateSelect.setDisabled(true);
    exportBtn.setDisabled(true);
  }
});

var clearROIBtn = ui.Button({
  label: 'Clear ROI',
  onClick: function() {
    roiCorner1 = null;
    capturedGeometry = null;
    roiMode = false;
    if (roiMapLayer) { Map.layers().remove(roiMapLayer); roiMapLayer = null; }
    roiStatus.setValue('Click "Draw ROI", then click two corners on the map.');
    msgLabel.setValue('Please draw a study area first.');
    dateSelect.items().reset([]);
    dateSelect.setDisabled(true);
    exportBtn.setDisabled(true);
  }
});

roiBtnRow.add(drawROIBtn);
roiBtnRow.add(clearROIBtn);
panel.add(roiBtnRow);

// --- STEP 2: CONFIGURATION ---
panel.add(ui.Label('2. Configuration', {fontWeight: 'bold', margin: '15px 0 5px 0'}));
var cloudPanel = ui.Panel({layout: ui.Panel.Layout.flow('horizontal')});
cloudPanel.add(ui.Label('Max Cloud Cover (%):', {padding: '5px'}));
var cloudInput = ui.Textbox({value: '30', style: {width: '60px'}});
cloudPanel.add(cloudInput);
panel.add(cloudPanel);

var findDatesBtn = ui.Button({label: 'Find Available Dates', onClick: findDates});
panel.add(findDatesBtn);

// --- STEP 3: SELECT IMAGE ---
panel.add(ui.Label('3. Select Image', {fontWeight: 'bold', margin: '15px 0 5px 0'}));
var dateSelect = ui.Select({
  placeholder: 'Waiting for dates...',
  disabled: true,
  onChange: function(date) { if (date) loadDate(date); }
});
panel.add(dateSelect);
var msgLabel = ui.Label('Please draw a study area first.');
panel.add(msgLabel);

// --- STEP 4: EXPORT ---
panel.add(ui.Label('4. Export', {fontWeight: 'bold', margin: '15px 0 5px 0'}));
panel.add(ui.Label('After clicking, go to the Tasks tab and click Run.', {margin: '1px 10px 1px 10px'}));
var exportBtn = ui.Button({
  label: 'Export Scene to Drive',
  disabled: true,
  onClick: function() { exportToDrive(dateSelect.getValue()); }
});
panel.add(exportBtn);

// --- STEP 5: VISUALIZATION ---
panel.add(ui.Label('5. Visualisation', {fontWeight: 'bold', margin: '15px 0 5px 0'}));
panel.add(ui.Label('Update to change SST visualisation range.', {margin: '1px 10px 1px 10px'}));
var visPanel = ui.Panel({layout: ui.Panel.Layout.flow('horizontal')});
var visMin = ui.Textbox({placeholder: 'Min', value: '-2', style: {width: '50px'}});
var visMax = ui.Textbox({placeholder: 'Max', value: '5', style: {width: '50px'}});
var visBtn = ui.Button('Update Range', updateVis);
visPanel.add(ui.Label('Min:')); visPanel.add(visMin);
visPanel.add(ui.Label('Max:')); visPanel.add(visMax);
visPanel.add(visBtn);
panel.add(visPanel);

// --- STEP 6: INSPECTOR ---
panel.add(ui.Label('6. Inspector', {fontWeight: 'bold', margin: '15px 0 5px 0'}));
var inspectorLabel = ui.Label('Click on the map to see SST values', {color: 'gray'});
panel.add(inspectorLabel);

// --- STEP 7 & 8: TIME-SERIES UI ---
panel.add(ui.Label('7 & 8. Time-series Analysis', {fontWeight: 'bold', margin: '15px 0 5px 0'}));
panel.add(ui.Label('With inspector (white hand) click on a point to calculate SST timeseries and seasonal range.', {margin: '1px 10px 1px 10px'}));
var yearPanel = ui.Panel({layout: ui.Panel.Layout.flow('horizontal')});
var startYearBox = ui.Textbox({value: '2015', style: {width: '60px'}});
var endYearBox = ui.Textbox({value: '2024', style: {width: '60px'}});
yearPanel.add(ui.Label('Start:')); yearPanel.add(startYearBox);
yearPanel.add(ui.Label('End:')); yearPanel.add(endYearBox);
panel.add(yearPanel);
var chartPanel = ui.Panel({style: {margin: '10px 0 0 0'}});
panel.add(chartPanel);

// ============================================================================
// 4. SHARED PROCESSING LOGIC
// ============================================================================

function processImage(img) {
  var satId = ee.String(img.get('SPACECRAFT_ID'));
  var isModern = satId.index('LANDSAT_8').gt(-1).or(satId.index('LANDSAT_9').gt(-1));

  var rgbImage = ee.Image(ee.Algorithms.If(isModern,
    img.select(['SR_B4', 'SR_B3', 'SR_B2'], ['Red', 'Green', 'Blue']),
    img.select(['SR_B3', 'SR_B2', 'SR_B1'], ['Red', 'Green', 'Blue'])
  ));

  var qa = img.select('QA_PIXEL');
  var maskWater = qa.bitwiseAnd(1 << 7).neq(0);

  var finalMask = ee.Image(ee.Algorithms.If(isModern,
    maskWater
      .and(qa.bitwiseAnd(1 << 3).neq(0).not())
      .and(qa.bitwiseAnd(1 << 4).neq(0).not()),
    maskWater
  )).rename('water_mask');

  var sstRaw = img.select('SST').rename('sst_unmasked');
  var sstMasked = sstRaw.updateMask(finalMask).rename('sst_masked');

  return img.addBands([rgbImage, sstMasked, sstRaw, finalMask], null, true)
    .copyProperties(img, ['system:time_start']);
}

// ============================================================================
// 5. FIND DATES
// ============================================================================

function findDates() {
  var roi = getUserGeometry();
  if (!roi) { msgLabel.setValue('ERROR: No study area defined. Click "Draw ROI" first.'); return; }
  msgLabel.setValue('Searching archive...');
  var cloudMax = parseFloat(cloudInput.getValue());

  var merged = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2')
    .merge(ee.ImageCollection('LANDSAT/LC08/C02/T1_L2'))
    .merge(ee.ImageCollection('LANDSAT/LE07/C02/T1_L2'))
    .merge(ee.ImageCollection('LANDSAT/LT05/C02/T1_L2'))
    .merge(ee.ImageCollection('LANDSAT/LT04/C02/T1_L2'))
    .filterBounds(roi)
    .filter(ee.Filter.lt('CLOUD_COVER', cloudMax));

  merged.aggregate_array('system:time_start')
    .map(function(time) { return ee.Date(time).format('YYYY-MM-dd'); })
    .distinct().sort().evaluate(function(dates) {
      if (!dates || dates.length === 0) {
        msgLabel.setValue('No images found under ' + cloudMax + '% cloud cover.');
        return;
      }
      dates.reverse();
      dateSelect.items().reset(dates);
      dateSelect.setDisabled(false);
      msgLabel.setValue('Found ' + dates.length + ' images. Select one.');
    });
}

// ============================================================================
// 6. LOAD AND DISPLAY
// ============================================================================

function getSSTCollection(sensor, dateStr, geom) {
  var col = LandsatSST.getCollection(sensor, dateStr, ee.Date(dateStr).advance(1, 'day'), geom);
  var correctId = ee.String('LANDSAT_').cat(sensor.substring(1));
  return col.map(function(img) { return img.set('SPACECRAFT_ID', correctId); });
}

function loadDate(dateStr) {
  var roi = getUserGeometry();
  if (!roi) return;
  msgLabel.setValue('Processing SST for ' + dateStr + '...');
  Map.layers().reset();
  roiMapLayer = null; // stale after reset

  // Re-draw ROI outline so it stays visible over the SST layer
  roiMapLayer = Map.addLayer(roi, {color: '2196F3', fillColor: '00000000'}, 'Study Area ROI');

  var collection = ee.ImageCollection([])
    .merge(getSSTCollection('L9', dateStr, roi))
    .merge(getSSTCollection('L8', dateStr, roi))
    .merge(getSSTCollection('L7', dateStr, roi))
    .merge(getSSTCollection('L5', dateStr, roi))
    .merge(getSSTCollection('L4', dateStr, roi));

  var processedCol = collection.map(processImage);

  processedCol.size().evaluate(function(count) {
    if (count === 0) { msgLabel.setValue('Error: No data returned for ' + dateStr + '.'); return; }

    currentImageToExport = processedCol.first();
    Map.addLayer(currentImageToExport.select(['Red', 'Green', 'Blue']),
      {min: 0.0, max: 0.5, gamma: 1.3}, 'True Color');
    Map.addLayer(currentImageToExport.select('sst_masked'),
      {min: -2, max: 5, palette: palettes.colorbrewer.RdBu[9].reverse()}, 'SST (Masked)');
    msgLabel.setValue('Displaying: ' + dateStr);
    exportBtn.setDisabled(false);
  });
}

// ============================================================================
// 7. EXPORT & VIS
// ============================================================================

function exportToDrive(dateStr) {
  if (!currentImageToExport) return;
  var fileName = areaName + '_' + dateStr;
  Export.image.toDrive({
    image: currentImageToExport.select(['Red', 'Green', 'Blue', 'sst_masked']).toFloat(),
    description: fileName,
    folder: 'GEE_SST_Exports',
    fileNamePrefix: fileName,
    region: getUserGeometry(),
    scale: 30,
    maxPixels: 1e9
  });
  msgLabel.setValue('Export Task Created! Check Tasks tab and click Run.');
  print('==============================================');
  print('EXPORT: ' + dateStr + '  →  Tasks tab → Run → ' + fileName);
  print('Output folder: Google Drive / GEE_SST_Exports');
  print('==============================================');
}

function updateVis() {
  var min = parseFloat(visMin.getValue());
  var max = parseFloat(visMax.getValue());
  Map.layers().forEach(function(layer) {
    if (layer.getName().indexOf('SST') > -1) {
      layer.setVisParams({min: min, max: max, palette: palettes.colorbrewer.RdBu[9].reverse()});
    }
  });
}

// ============================================================================
// 8. TIME-SERIES HELPERS
// ============================================================================

function getRawSensor(sensor, start, end, geom) {
  var col = LandsatSST.getCollection(sensor, start, end, geom);
  var correctId = ee.String('LANDSAT_').cat(sensor.substring(1));
  return col
    .filterBounds(geom)
    .filter(ee.Filter.lt('CLOUD_COVER', 60))
    .map(function(img) { return img.set('SPACECRAFT_ID', correctId); });
}

function extractTimeseries(collection, geom) {
  return collection.map(function(img) {
    var stats = img.reduceRegion({
      reducer: ee.Reducer.mean(), geometry: geom, scale: 30, maxPixels: 1e9
    });
    return ee.Feature(null, {
      'sst_masked': stats.get('sst_masked'),
      'system:time_start': img.date().millis(),
      'label': img.date().format('YYYY-MM-dd')
    });
  }).filter(ee.Filter.notNull(['sst_masked']));
}

function calculateRunningMean(featureCol, windowDays) {
  var range = windowDays * 24 * 60 * 60 * 1000;
  var withMatches = ee.Join.saveAll({matchesKey: 'matches'}).apply(
    featureCol, featureCol,
    ee.Filter.maxDifference({difference: range, leftField: 'system:time_start', rightField: 'system:time_start'})
  );
  return withMatches.map(function(f) {
    return f.set('sst_smooth', ee.FeatureCollection(ee.List(f.get('matches'))).aggregate_mean('sst_masked'));
  });
}

function aggregateToMonthly(featureCol) {
  var withProps = featureCol.map(function(f) {
    var date = ee.Date(f.get('system:time_start'));
    var year = date.get('year');
    var month = date.get('month');
    var block = year.divide(5).floor().multiply(5);
    var blockLabel = ee.String(ee.Number(block).format('%d'))
      .cat('-').cat(ee.String(ee.Number(block.add(4)).format('%d')));
    return f.set({
      'month': month,
      'year_block': blockLabel,
      'group_id': blockLabel.cat('_').cat(ee.Number(month).format('%02d'))
    });
  });

  var stats = withProps.reduceColumns({
    reducer: ee.Reducer.mean().group({groupField: 1, groupName: 'group_id'}),
    selectors: ['sst_masked', 'group_id']
  });

  return ee.FeatureCollection(ee.List(stats.get('groups')).map(function(item) {
    var dict = ee.Dictionary(item);
    var parts = ee.String(dict.get('group_id')).split('_');
    return ee.Feature(null, {
      'sst_mean': dict.get('mean'),
      'year_block': parts.get(0),
      'month': ee.Number.parse(parts.get(1))
    });
  }));
}

// ============================================================================
// 9. UNIFIED MAP CLICK HANDLER
// ============================================================================

Map.onClick(function(coords) {

  // --- Priority 1: ROI definition (two-click rectangle) ---
  if (roiMode && !roiCorner1) {
    roiCorner1 = [coords.lon, coords.lat];
    roiStatus.setValue(
      'First corner set (' + coords.lon.toFixed(3) + ', ' + coords.lat.toFixed(3) + '). ' +
      'Now click the OPPOSITE corner.'
    );
    if (roiMapLayer) Map.layers().remove(roiMapLayer);
    roiMapLayer = Map.addLayer(
      ee.Geometry.Point(roiCorner1), {color: 'red'}, 'ROI Corner 1'
    );
    return;
  }

  if (roiMode && roiCorner1) {
    var west  = Math.min(roiCorner1[0], coords.lon);
    var south = Math.min(roiCorner1[1], coords.lat);
    var east  = Math.max(roiCorner1[0], coords.lon);
    var north = Math.max(roiCorner1[1], coords.lat);
    capturedGeometry = ee.Geometry.Rectangle([west, south, east, north]);
    roiMode = false;
    roiCorner1 = null;
    if (roiMapLayer) Map.layers().remove(roiMapLayer);
    roiMapLayer = Map.addLayer(capturedGeometry, {color: '2196F3', fillColor: '2196F322'}, 'Study Area ROI');
    roiStatus.setValue(
      'ROI: (' + west.toFixed(2) + ', ' + south.toFixed(2) + ') → (' + east.toFixed(2) + ', ' + north.toFixed(2) + ')'
    );
    msgLabel.setValue('ROI defined. Click "Find Available Dates".');
    return;
  }

  // --- Priority 2: Inspector + time-series ---
  var point = ee.Geometry.Point([coords.lon, coords.lat]);
  var sYear = startYearBox.getValue();
  var eYear = endYearBox.getValue();

  inspectorLabel.setValue('Inspecting...');
  if (currentImageToExport) {
    currentImageToExport.select('sst_masked').reduceRegion({
      reducer: ee.Reducer.first(), geometry: point, scale: 30
    }).evaluate(function(val) {
      inspectorLabel.setValue(
        val && val.sst_masked !== null
          ? 'SST: ' + val.sst_masked.toFixed(2) + '°C'
          : 'No Data (Cloud/Land/Ice)'
      );
    });
  }

  chartPanel.clear();
  chartPanel.add(ui.Label('Extracting data... (' + sYear + '–' + eYear + ')', {color: 'gray'}));

  var searchGeom = point.buffer(500);
  var rawCol = ee.ImageCollection([])
    .merge(getRawSensor('L9', sYear + '-01-01', eYear + '-12-31', searchGeom))
    .merge(getRawSensor('L8', sYear + '-01-01', eYear + '-12-31', searchGeom))
    .merge(getRawSensor('L7', sYear + '-01-01', eYear + '-12-31', searchGeom))
    .merge(getRawSensor('L5', sYear + '-01-01', eYear + '-12-31', searchGeom))
    .merge(getRawSensor('L4', sYear + '-01-01', eYear + '-12-31', searchGeom));

  var validTable = extractTimeseries(rawCol.map(processImage), searchGeom)
    .sort('system:time_start');

  validTable.size().evaluate(function(count) {
    if (count < 1) {
      chartPanel.clear();
      chartPanel.add(ui.Label('Insufficient data (' + count + ' pts). Try another point.', {color: 'red'}));
      return;
    }

    var trendData = calculateRunningMean(validTable, 60);
    var seasonalData = aggregateToMonthly(validTable);
    chartPanel.clear();

    chartPanel.add(
      ui.Chart.feature.byFeature({
        features: trendData,
        xProperty: 'system:time_start',
        yProperties: ['sst_masked', 'sst_smooth']
      })
      .setSeriesNames(['Raw Data', 'Running Mean (60d)'])
      .setOptions({
        title: 'SST Trend (' + count + ' points)',
        vAxis: {title: 'SST (°C)', viewWindow: {min: 0, max: 25}},
        hAxis: {title: 'Date', format: 'yyyy', gridlines: {count: -1}},
        series: {
          0: {pointSize: 2, lineWidth: 0, color: '#A1C4FD'},
          1: {pointSize: 0, lineWidth: 2, color: '#0052cc'}
        }
      })
    );

    var monthTicks = [
      {v:1,f:'Jan'},{v:2,f:'Feb'},{v:3,f:'Mar'},{v:4,f:'Apr'},
      {v:5,f:'May'},{v:6,f:'Jun'},{v:7,f:'Jul'},{v:8,f:'Aug'},
      {v:9,f:'Sep'},{v:10,f:'Oct'},{v:11,f:'Nov'},{v:12,f:'Dec'}
    ];

    chartPanel.add(
      ui.Chart.feature.groups({
        features: seasonalData,
        xProperty: 'month',
        yProperty: 'sst_mean',
        seriesProperty: 'year_block'
      })
      .setOptions({
        title: 'Mean Seasonal Cycle (Monthly Aggregates)',
        curveType: 'function',
        vAxis: {title: 'SST (°C)'},
        hAxis: {title: 'Month', ticks: monthTicks, viewWindow: {min: 1, max: 12}},
        lineWidth: 3,
        pointSize: 4
      })
    );
  });
});
