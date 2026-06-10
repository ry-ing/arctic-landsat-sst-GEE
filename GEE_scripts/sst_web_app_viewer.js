// Google Earth Engine Landsat-derived SST's - optimised for Greenland coastal waters
// Created by Ryan Ing, University of Edinburgh
// Email: ryan.ing@ed.ac.uk

// ============================================================================
// 1. IMPORTS AND CONFIGURATION
// ============================================================================

var LandsatSST = require('users/ryaning/arctic-landsat-sst:sst-retrieval/Landsat_SST');
var palettes = require('users/gena/packages:palettes');

// Define filename export
var areaName = 'Landsat_SST_export';

// Initialize Map
// Map.setCenter(-50.6, 70.4, 6); 
Map.setOptions('SATELLITE');
Map.style().set('cursor', 'crosshair');

//  Current image for exporting
var currentImageToExport = null;

// ============================================================================
// 2. DRAWING TOOLS SETUP
// ============================================================================

var drawingTools = Map.drawingTools();
drawingTools.setShown(true);
drawingTools.setDrawModes(['polygon', 'rectangle']);
drawingTools.addLayer([], 'Geometry', 'blue');
drawingTools.setShape('rectangle');
drawingTools.draw();

function getUserGeometry() {
  var layers = drawingTools.layers();
  if (layers.length() > 0) {
    var layer = layers.get(0);
    var geometries = layer.geometries();
    if (geometries.length() > 0) {
      return ee.Geometry(layer.toGeometry()); 
    }
  }
  return null;
}

// ============================================================================
// 3. UI SETUP
// ============================================================================

var panel = ui.Panel({style: {width: '420px', padding: '10px'}});
ui.root.insert(0, panel);

// Title
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

// --- STEP 1: DRAW ROI ---
panel.add(ui.Label('1. Define Area', {fontWeight: 'bold', margin: '10px 0 5px 0'}));
panel.add(ui.Label('Draw geometry around region of interest.', {margin: '1px 10px 1px 10px'}));

var clearGeomBtn = ui.Button({
  label: 'Clear Geometry',
  onClick: function() {
    while (drawingTools.layers().length() > 0) {
      drawingTools.layers().remove(drawingTools.layers().get(0));
    }
    drawingTools.addLayer([], 'Geometry', 'red');
    drawingTools.draw();
    dateSelect.items().reset();
    dateSelect.setDisabled(true);
    msgLabel.setValue('Geometry cleared. Draw new area.');
  }
});
panel.add(clearGeomBtn);

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
var msgLabel = ui.Label('Please draw a geometry first.');
panel.add(msgLabel);

// --- STEP 4: EXPORT ---
panel.add(ui.Label('4. Export', {fontWeight: 'bold', margin: '15px 0 5px 0'}));
panel.add(ui.Label(
  'Export to Drive is not available in the published App. ' +
  'To export a scene, run sst_web_app_viewer_EDITOR.js directly in the GEE Code Editor.',
  {margin: '1px 10px 5px 10px', color: 'gray', fontSize: '12px'}
));
var exportBtn = ui.Button({label: 'Export Scene to Drive', disabled: true,
  onClick: function() {
    msgLabel.setValue('Export is only available in the Code Editor version.');
  }
});
panel.add(exportBtn);

// --- STEP 5: VISUALIZATION ---
panel.add(ui.Label('5. Visualisation', {fontWeight: 'bold', margin: '15px 0 5px 0'}));
panel.add(ui.Label('Update to change SST visulisation range.', {margin: '1px 10px 1px 10px'}));
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
panel.add(ui.Label('With inspector (white hand) click on a point to calculate SST timseries and seasonal range.', {margin: '1px 10px 1px 10px'}));
var yearPanel = ui.Panel({layout: ui.Panel.Layout.flow('horizontal')});
var startYearBox = ui.Textbox({value: '2015', style: {width: '60px'}});
var endYearBox = ui.Textbox({value: '2024', style: {width: '60px'}});
yearPanel.add(ui.Label('Start:')); yearPanel.add(startYearBox);
yearPanel.add(ui.Label('End:')); yearPanel.add(endYearBox);
panel.add(yearPanel);
var chartPanel = ui.Panel({style: {margin: '10px 0 0 0'}});
panel.add(chartPanel);

// ============================================================================
// 4. SHARED PROCESSING LOGIC (CRITICAL FIX)
// ============================================================================

// This function applies the EXACT SAME masking logic to single images and time-series
function processImage(img) {
  var satId = ee.String(img.get('SPACECRAFT_ID'));
  var isModern = satId.index('LANDSAT_8').gt(-1).or(satId.index('LANDSAT_9').gt(-1));
  
  // 1. RGB Logic
  var rgbImage = ee.Image(ee.Algorithms.If(isModern,
    img.select(['SR_B4', 'SR_B3', 'SR_B2'], ['Red', 'Green', 'Blue']),
    img.select(['SR_B3', 'SR_B2', 'SR_B1'], ['Red', 'Green', 'Blue'])
  ));

  // 2. MASKING LOGIC (Replicated from your loadDate function)
  var qa = img.select('QA_PIXEL');
  var maskWater = qa.bitwiseAnd(1 << 7).neq(0); // Standard Water Bit
  
  var finalMask = ee.Image(ee.Algorithms.If(isModern,
    // Modern: Water AND NOT Cloud(3) AND NOT Shadow(4)
    maskWater
      .and(qa.bitwiseAnd(1 << 3).neq(0).not())
      .and(qa.bitwiseAnd(1 << 4).neq(0).not()),
    // Old: Water Only
    maskWater
  )).rename('water_mask');
  
  // 3. Band Renaming
  var sstRaw = img.select('SST').rename('sst_unmasked');
  var sstMasked = sstRaw.updateMask(finalMask).rename('sst_masked');
  
  // 4. Return processed image with all bands and time properties
  return img.addBands([rgbImage, sstMasked, sstRaw, finalMask], null, true)
    .copyProperties(img, ['system:time_start']); 
}

// ============================================================================
// 5. LOGIC: FIND DATES
// ============================================================================

function findDates() {
  var roi = getUserGeometry();
  if (!roi) { msgLabel.setValue('ERROR: No geometry found.'); return; }
  msgLabel.setValue('Searching archive...');
  var cloudMax = parseFloat(cloudInput.getValue());
  
  var l9 = ee.ImageCollection("LANDSAT/LC09/C02/T1_L2");
  var l8 = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2");
  var l7 = ee.ImageCollection("LANDSAT/LE07/C02/T1_L2");
  var l5 = ee.ImageCollection("LANDSAT/LT05/C02/T1_L2");
  var l4 = ee.ImageCollection("LANDSAT/LT04/C02/T1_L2");
  
  var merged = l9.merge(l8).merge(l7).merge(l5).merge(l4)
    .filterBounds(roi)
    .filter(ee.Filter.lt('CLOUD_COVER', cloudMax));

  merged.aggregate_array('system:time_start')
    .map(function(time) { return ee.Date(time).format('YYYY-MM-dd'); })
    .distinct().sort().evaluate(function(dates) {
      if (!dates || dates.length === 0) {
        msgLabel.setValue('No images found.');
        return;
      }
      dates.reverse();
      dateSelect.items().reset(dates);
      dateSelect.setDisabled(false);
      msgLabel.setValue('Found ' + dates.length + ' images.');
    });
}

// ============================================================================
// 6. LOGIC: LOAD AND DISPLAY
// ============================================================================

function getSSTCollection(sensor, dateStr, geom) {
  var startDate = dateStr;
  var endDate = ee.Date(dateStr).advance(1, 'day');
  var col = LandsatSST.getCollection(sensor, startDate, endDate, geom);
  // Pre-tag spacecraft ID for the processor
  var correctId = ee.String('LANDSAT_').cat(sensor.substring(1));
  return col.map(function(img) {
    return img.set('SPACECRAFT_ID', correctId);
  });
}

function loadDate(dateStr) {
  var roi = getUserGeometry();
  if (!roi) return;
  msgLabel.setValue('Processing SST for ' + dateStr + '...');
  Map.layers().reset();
  
  var collection = ee.ImageCollection([])
    .merge(getSSTCollection('L9', dateStr, roi))
    .merge(getSSTCollection('L8', dateStr, roi))
    .merge(getSSTCollection('L7', dateStr, roi))
    .merge(getSSTCollection('L5', dateStr, roi))
    .merge(getSSTCollection('L4', dateStr, roi));
    
  // Apply the shared processor
  var processedCol = collection.map(processImage);

  processedCol.size().evaluate(function(count) {
    if (count === 0) { msgLabel.setValue('Error: No data.'); return; }
    
    currentImageToExport = processedCol.first();
    var rgbImage = currentImageToExport.select(['Red', 'Green', 'Blue']);
    var sstMasked = currentImageToExport.select('sst_masked');
    var sstRaw = currentImageToExport.select('sst_unmasked');

    Map.addLayer(rgbImage, {min: 0.0, max: 0.5, gamma: 1.3}, 'True Color');
    Map.addLayer(sstMasked, {min: -2, max: 5, palette: palettes.colorbrewer.RdBu[9].reverse()}, 'SST (Masked)');
    msgLabel.setValue('Displaying: ' + dateStr);
    exportBtn.setDisabled(false);
  });
}

// ============================================================================
// 7. EXPORT & VIS LOGIC
// ============================================================================

function exportToDrive(dateStr) {
  msgLabel.setValue('Export is only available in the Code Editor version (sst_web_app_viewer_EDITOR.js).');
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
// 8. OPTIMIZED FEATURE-BASED TIME-SERIES (FIXED)
// ============================================================================

// 1. Helper to get raw collection for a sensor
function getRawSensor(sensor, start, end, geom) {
  var col = LandsatSST.getCollection(sensor, start, end, geom);
  var correctId = ee.String('LANDSAT_').cat(sensor.substring(1));
  return col
    .filterBounds(geom)
    .filter(ee.Filter.lt('CLOUD_COVER', 60)) // Drop scenes >60% cloud (useless for SST)
    .map(function(img) {
    return img.set('SPACECRAFT_ID', correctId);
  });
}

// 2. FEATURE EXTRACTION (The Memory Saver)
function extractTimeseries(collection, geom) {
  return collection.map(function(img) {
    var stats = img.reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: geom,
      scale: 30,
      maxPixels: 1e9
    });
    
    // We explicitly save system:time_start so we can read it reliably later
    return ee.Feature(null, {
      'sst_masked': stats.get('sst_masked'),
      'system:time_start': img.date().millis(),
      'label': img.date().format('YYYY-MM-dd')
    });
  }).filter(ee.Filter.notNull(['sst_masked']));
}

// 3. TABLE SMOOTHING (Running Mean)
function calculateRunningMean(featureCol, windowDays) {
  var range = windowDays * 24 * 60 * 60 * 1000;
  var join = ee.Join.saveAll({matchesKey: 'matches'});
  var filter = ee.Filter.maxDifference({
    difference: range,
    leftField: 'system:time_start',
    rightField: 'system:time_start'
  });
  
  var withMatches = join.apply(featureCol, featureCol, filter);
  
  return withMatches.map(function(f) {
    var matches = ee.List(f.get('matches'));
    var mean = ee.FeatureCollection(matches).aggregate_mean('sst_masked');
    return f.set('sst_smooth', mean);
  });
}

// 4. MONTHLY AGGREGATOR (For Seasonal Plot)
function aggregateToMonthly(featureCol) {
  // Add grouping properties
  var withProps = featureCol.map(function(f) {
    // FIX: Explicitly read the timestamp property instead of calling .date()
    // This prevents the "Expected Image, Actual Feature" error
    var date = ee.Date(f.get('system:time_start'));
    
    var year = date.get('year');
    var month = date.get('month');
    
    // Calculate 5-Year Block (e.g., 2015-2019)
    var block = year.divide(5).floor().multiply(5);
    var blockLabel = ee.String(ee.Number(block).format('%d'))
        .cat('-').cat(ee.String(ee.Number(block.add(4)).format('%d')));
    
    return f.set({
      'month': month,
      'year_block': blockLabel,
      'group_id': blockLabel.cat('_').cat(ee.Number(month).format('%02d'))
    });
  });

  // Group by (YearBlock + Month) and calculate Mean SST
  var stats = withProps.reduceColumns({
    reducer: ee.Reducer.mean().group({groupField: 1, groupName: 'group_id'}),
    selectors: ['sst_masked', 'group_id']
  });

  // Convert dictionary back to FeatureCollection for plotting
  var groupList = ee.List(stats.get('groups'));
  var monthlyFeatures = groupList.map(function(item) {
    var dict = ee.Dictionary(item);
    var groupId = ee.String(dict.get('group_id'));
    var parts = groupId.split('_');
    
    return ee.Feature(null, {
      'sst_mean': dict.get('mean'),
      'year_block': parts.get(0),
      'month': ee.Number.parse(parts.get(1))
    });
  });

  return ee.FeatureCollection(monthlyFeatures);
}

// 5. MAP CLICK HANDLER
Map.onClick(function(coords) {
  var point = ee.Geometry.Point([coords.lon, coords.lat]);
  var sYear = startYearBox.getValue();
  var eYear = endYearBox.getValue();
  
  // A. Inspector Update
  inspectorLabel.setValue('Inspecting...');
  if (currentImageToExport) {
    var sample = currentImageToExport.select('sst_masked').reduceRegion({
      reducer: ee.Reducer.first(), geometry: point, scale: 30
    });
    sample.evaluate(function(val) {
      inspectorLabel.setValue(val && val.sst_masked !== null ? 'SST: ' + val.sst_masked.toFixed(2) + '°C' : 'No Data');
    });
  }

  // B. Charts
  chartPanel.clear();
  chartPanel.add(ui.Label('Extracting Data... (Processing ' + sYear + '-' + eYear + ')', {color: 'gray'}));
  
  var searchGeom = point.buffer(500); // 1km buffer
  var startD = sYear + '-01-01';
  var endD = eYear + '-12-31';

  // 1. Fetch & Process Images
  var rawCol = ee.ImageCollection([])
    .merge(getRawSensor('L9', startD, endD, searchGeom))
    .merge(getRawSensor('L8', startD, endD, searchGeom))
    .merge(getRawSensor('L7', startD, endD, searchGeom))
    .merge(getRawSensor('L5', startD, endD, searchGeom))
    .merge(getRawSensor('L4', startD, endD, searchGeom));
    
  var processedCol = rawCol.map(processImage);

  // 2. Convert to Table
  var table = extractTimeseries(processedCol, searchGeom);
  
  // 3. Filter Valid Data & Sort Chronologically
  var validTable = table
    .sort('system:time_start');

  validTable.size().evaluate(function(count) {
    if (count < 1) { 
      chartPanel.clear(); 
      chartPanel.add(ui.Label('Insufficient data ('+count+' pts). Try another point.', {color: 'red'})); 
      return; 
    }

    // 4. Generate Data for Plots
    var trendData = calculateRunningMean(validTable, 60); // 60-Day Smoothing
    var seasonalData = aggregateToMonthly(validTable);    // Monthly Means

    chartPanel.clear();

// --- CHART 1: TREND (Clean Lines) ---
    var trendChart = ui.Chart.feature.byFeature({
      features: trendData,
      xProperty: 'system:time_start',
      // SWAPPED HERE: Put 'sst_smooth' LAST so it draws ON TOP
      yProperties: ['sst_masked', 'sst_smooth'] 
    }).setSeriesNames(['Raw Data', 'Running Mean (60d)']) // Names match the order above
      .setOptions({
        title: 'SST Trend (' + count + ' points)',
        vAxis: {title: 'SST (°C)', viewWindow: {min: 0, max: 25}},
        hAxis: {title: 'Date', format: 'yyyy', gridlines: {count: -1}},
        series: {
          // Series 0 is now 'sst_masked' (Raw Data) -> Plot as Dots
          0: {pointSize: 2, lineWidth: 0, color: '#A1C4FD'}, 
          // Series 1 is now 'sst_smooth' (Running Mean) -> Plot as Line on Top
          1: {pointSize: 0, lineWidth: 2, color: '#0052cc'}  
        }
      });
    chartPanel.add(trendChart);

    // --- CHART 2: SEASONAL (Smooth Monthly Curves) ---
    var monthTicks = [
      {v: 1, f: 'Jan'}, {v: 2, f: 'Feb'}, {v: 3, f: 'Mar'},
      {v: 4, f: 'Apr'}, {v: 5, f: 'May'}, {v: 6, f: 'Jun'},
      {v: 7, f: 'Jul'}, {v: 8, f: 'Aug'}, {v: 9, f: 'Sep'},
      {v: 10, f: 'Oct'}, {v: 11, f: 'Nov'}, {v: 12, f: 'Dec'}
    ];
 
    var seasonalChart = ui.Chart.feature.groups({
      features: seasonalData,
      xProperty: 'month',
      yProperty: 'sst_mean',
      seriesProperty: 'year_block'
    }).setOptions({
        title: 'Mean Seasonal Cycle (Monthly Aggregates)',
        curveType: 'function',
        vAxis: {title: 'SST (°C)'},
        hAxis: {
          title: 'Month', 
          ticks: monthTicks,
          viewWindow: {min: 1, max: 12}
        },
        lineWidth: 3,
        pointSize: 4
    });
    chartPanel.add(seasonalChart);
  });
});