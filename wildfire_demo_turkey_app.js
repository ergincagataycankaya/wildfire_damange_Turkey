// =====================================================================
//  WILDFIRE DEMO IN TURKEY – Earth Engine App
//  Sentinel-2 burn scar mapping & burn severity (dNBR)
//  No imports needed. Publish via: Apps  ->  New App
// =====================================================================

// ------------------------- FIRE PRESETS -------------------------
var FIRES = {
  'İzmir – Karabağlar/Menderes (Aug 2019)': {
    aoi: ee.Geometry.Polygon([[
      [26.99059815747183, 38.27096885531639], [26.962789014405423, 38.261636895183976],
      [26.96484895092886, 38.24761769501632], [26.945279553956205, 38.24006777471134],
      [26.94218964917105, 38.21336712677703], [26.958325818604642, 38.202036618379246],
      [26.98785157544058, 38.2120183492158], [26.989911511964017, 38.193942318311926],
      [27.04827638012808, 38.23035958219917], [27.070592359131986, 38.27349718272393],
      [27.054112866944486, 38.3125685123244], [27.004674390381986, 38.313107280125934],
      [26.99059815747183, 38.27096885531639]]]),
    preStart: '2019-07-05', fireDate: '2019-08-17', postEnd: '2019-10-05', zoom: 12
  },
  'Antalya – Manavgat (Jul–Aug 2021)': {
    aoi: ee.Geometry.Rectangle([31.25, 36.72, 31.80, 37.12]),
    preStart: '2021-06-01', fireDate: '2021-07-27', postEnd: '2021-09-30', zoom: 10
  },
  'Muğla – Marmaris (Aug 2021)': {
    aoi: ee.Geometry.Rectangle([28.02, 36.76, 28.42, 37.02]),
    preStart: '2021-06-15', fireDate: '2021-08-01', postEnd: '2021-09-30', zoom: 10
  },
  'Custom area (draw on map)': {
    aoi: null,
    preStart: '2021-06-01', fireDate: '2021-07-27', postEnd: '2021-09-30', zoom: 10
  }
};
var FIRE_NAMES = Object.keys(FIRES);
var CUSTOM = FIRE_NAMES[3];

// ------------------------- SEVERITY CLASSES (USGS) -------------------------
var SEV_NAMES   = ['High regrowth', 'Low regrowth', 'Unburned', 'Low severity',
                   'Moderate-low', 'Moderate-high', 'High severity'];
var SEV_PALETTE = ['#7a8737', '#acbe4d', '#0ae042', '#fff70b',
                   '#ffaf38', '#ff641b', '#a41fd6'];

var rgbVis  = {bands: ['B4', 'B3', 'B2'], min: 0, max: 0.3};
var swirVis = {bands: ['B12', 'B8', 'B4'], min: 0, max: 0.4};
var dnbrVis = {min: -0.3, max: 0.9,
               palette: ['green', 'white', 'yellow', 'orange', 'red', 'purple']};

// ------------------------- PROCESSING -------------------------
// Cloud mask that works with both QA60 and the newer MSK_CLASSI_* bands
function maskL1C(img) {
  var mask = ee.Image(ee.Algorithms.If(
    img.bandNames().contains('QA60'),
    img.select('QA60').bitwiseAnd(1 << 10).eq(0)
       .and(img.select('QA60').bitwiseAnd(1 << 11).eq(0)),
    img.select('MSK_CLASSI_OPAQUE').eq(0)
       .and(img.select('MSK_CLASSI_CIRRUS').eq(0))
  ));
  return img.select('B.*').updateMask(mask).divide(10000)
            .copyProperties(img, ['system:time_start']);
}

function maskL2A(img) {
  var scl = img.select('SCL');
  var mask = scl.neq(3).and(scl.neq(8)).and(scl.neq(9)).and(scl.neq(10));
  return img.select('B.*').updateMask(mask).divide(10000)
            .copyProperties(img, ['system:time_start']);
}

function l1c(aoi, start, end) {
  return ee.ImageCollection('COPERNICUS/S2_HARMONIZED')
    .filterBounds(aoi)
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
    .filterDate(start, end)
    .map(maskL1C);
}

function analyse(aoi, preStart, fireDate, postEnd) {
  var preCol  = l1c(aoi, preStart, fireDate);
  var postCol = l1c(aoi, fireDate, postEnd);
  var preMed  = preCol.median().clip(aoi);
  var postMed = postCol.median().clip(aoi);

  var nbrPre  = preMed.normalizedDifference(['B8', 'B12']);
  var nbrPost = postMed.normalizedDifference(['B8', 'B12']);
  // mask permanent water (ESA WorldCover class 80) so the sea is not "unburned"
  var land    = ee.ImageCollection('ESA/WorldCover/v200').first().neq(80);
  var dNBR    = nbrPre.subtract(nbrPost).updateMask(land).rename('dNBR');

  var severity = ee.Image(3)
    .where(dNBR.lt(-0.25), 1)
    .where(dNBR.gte(-0.25).and(dNBR.lt(-0.1)), 2)
    .where(dNBR.gte(0.1).and(dNBR.lt(0.27)), 4)
    .where(dNBR.gte(0.27).and(dNBR.lt(0.44)), 5)
    .where(dNBR.gte(0.44).and(dNBR.lt(0.66)), 6)
    .where(dNBR.gte(0.66), 7)
    .updateMask(dNBR.mask()).clip(aoi).rename('severity');

  var burned = dNBR.gte(0.1).selfMask().rename('burned');

  return {preCol: preCol, postCol: postCol, preMed: preMed, postMed: postMed,
          dNBR: dNBR, severity: severity, burned: burned};
}

function nbrTimeSeries(aoi, start, end) {
  var col = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(aoi)
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 30))
    .filterDate(start, end)
    .map(maskL2A);
  // one image per day (AOIs can span several tiles)
  var days = col.aggregate_array('system:time_start')
    .map(function(t) { return ee.Date(t).format('YYYY-MM-dd'); }).distinct();
  return ee.ImageCollection(days.map(function(d) {
    var day = ee.Date.parse('YYYY-MM-dd', d);
    return col.filterDate(day, day.advance(1, 'day')).mosaic()
      .normalizedDifference(['B8', 'B12']).rename('NBR')
      .set('system:time_start', day.millis());
  }));
}

// ------------------------- UI: MAPS -------------------------
var leftMap  = ui.Map();
var rightMap = ui.Map();
[leftMap, rightMap].forEach(function(m) {
  m.setOptions('HYBRID');
  m.setControlVisibility({layerList: true, zoomControl: true, mapTypeControl: true,
                          scaleControl: true, drawingToolsControl: false,
                          fullscreenControl: false});
});
ui.Map.Linker([leftMap, rightMap]);

var leftTag  = ui.Label('BEFORE', {position: 'top-left', fontWeight: 'bold', fontSize: '14px'});
var rightTag = ui.Label('AFTER',  {position: 'top-right', fontWeight: 'bold', fontSize: '14px'});
leftMap.add(leftTag);
rightMap.add(rightTag);

var split = ui.SplitPanel({firstPanel: leftMap, secondPanel: rightMap,
                           wipe: true, style: {stretch: 'both'}});

// Drawing tools for custom AOI
var drawingTools = rightMap.drawingTools();
drawingTools.setShown(false);
drawingTools.setLinked(false);
while (drawingTools.layers().length() > 0) {
  drawingTools.layers().remove(drawingTools.layers().get(0));
}
var drawLayer = ui.Map.GeometryLayer({geometries: null, name: 'Custom AOI', color: 'yellow'});
drawingTools.layers().add(drawLayer);

// ------------------------- UI: CONTROL PANEL -------------------------
var S = {
  h1:   {fontSize: '22px', fontWeight: 'bold', color: '#b22222', margin: '4px 8px 0 8px'},
  h2:   {fontSize: '14px', fontWeight: 'bold', margin: '12px 8px 4px 8px'},
  txt:  {fontSize: '12px', color: '#444', margin: '2px 8px'},
  box:  {width: '110px', margin: '2px 8px'}
};

var panel = ui.Panel({style: {width: '360px', padding: '6px'}});
panel.add(ui.Label('🔥 Wildfire Demo in Turkey', S.h1));
panel.add(ui.Label('Burn scar mapping and burn severity (dNBR) from Sentinel-2 imagery. ' +
                   'Choose a fire, adjust the dates if needed and press "Run analysis".', S.txt));

panel.add(ui.Label('1. Select a fire', S.h2));
var fireSelect = ui.Select({items: FIRE_NAMES, style: {stretch: 'horizontal', margin: '2px 8px'}});
panel.add(fireSelect);

var drawPanel = ui.Panel({style: {shown: false}});
drawPanel.add(ui.Label('Draw a polygon around the burned area on the right map.', S.txt));
drawPanel.add(ui.Panel([
  ui.Button('✏️ Draw area', function() {
    drawLayer.geometries().reset();
    drawingTools.setShape('polygon');
    drawingTools.draw();
  }),
  ui.Button('Clear', function() { drawLayer.geometries().reset(); })
], ui.Panel.Layout.flow('horizontal')));
panel.add(drawPanel);

panel.add(ui.Label('2. Dates (YYYY-MM-DD)', S.h2));
var preBox  = ui.Textbox({style: S.box});
var fireBox = ui.Textbox({style: S.box});
var postBox = ui.Textbox({style: S.box});
function dateRow(label, box) {
  return ui.Panel([ui.Label(label, {width: '150px', fontSize: '12px', margin: '8px 8px'}), box],
                  ui.Panel.Layout.flow('horizontal'));
}
panel.add(dateRow('Pre-fire start', preBox));
panel.add(dateRow('Fire date (pre end / post start)', fireBox));
panel.add(dateRow('Post-fire end', postBox));

var runButton = ui.Button({label: '▶ Run analysis', style: {stretch: 'horizontal', margin: '10px 8px'}});
panel.add(runButton);

// Zoom helpers: fit the whole area with some margin, or show all of Turkey
function zoomToArea() {
  var aoi = getAOI();
  if (!aoi) return;
  var z = FIRES[fireSelect.getValue()].zoom;
  leftMap.centerObject(aoi, z);
}
var TURKEY = ee.Geometry.Rectangle([25.6, 35.8, 44.8, 42.1]);
panel.add(ui.Panel([
  ui.Button({label: '🔍 Zoom to area', onClick: zoomToArea, style: {stretch: 'horizontal'}}),
  ui.Button({label: '🌍 Zoom out (Türkiye)', style: {stretch: 'horizontal'},
             onClick: function() { leftMap.centerObject(TURKEY, 6); }})
], ui.Panel.Layout.flow('horizontal'), {margin: '0 4px'}));

// Layer switches for both maps (the left map's own layer list is hidden
// under the right map in wipe mode, so control both from here)
var LAYER_NAMES = ['True colour', 'SWIR false colour', 'dNBR', 'Burn severity', 'Area outline'];
var DEFAULT_ON = {left: [true, false, false, false, true], right: [true, false, false, true, true]};
var layerBoxes = {left: [], right: []};

panel.add(ui.Label('Map layers', S.h2));
panel.add(ui.Panel([
  ui.Label('', {width: '150px', margin: '0 8px'}),
  ui.Label('Before', {fontSize: '11px', fontWeight: 'bold', width: '60px', margin: '0'}),
  ui.Label('After', {fontSize: '11px', fontWeight: 'bold', width: '60px', margin: '0'})
], ui.Panel.Layout.flow('horizontal')));

LAYER_NAMES.forEach(function(name, i) {
  function box(side, map) {
    var cb = ui.Checkbox({label: '', value: DEFAULT_ON[side][i],
                          style: {width: '60px', margin: '0'}});
    cb.onChange(function(on) {
      var lyr = map.layers().get(i);
      if (lyr) lyr.setShown(on);
    });
    layerBoxes[side].push(cb);
    return cb;
  }
  panel.add(ui.Panel([
    ui.Label(name, {fontSize: '12px', width: '150px', margin: '4px 8px'}),
    box('left', leftMap),
    box('right', rightMap)
  ], ui.Panel.Layout.flow('horizontal')));
});

panel.add(ui.Label('3. Results', S.h2));
var results = ui.Panel();
panel.add(results);

// Legend
panel.add(ui.Label('Burn severity (dNBR)', S.h2));
SEV_NAMES.forEach(function(n, i) {
  panel.add(ui.Panel([
    ui.Label('', {backgroundColor: SEV_PALETTE[i], padding: '8px', margin: '2px 6px 2px 8px'}),
    ui.Label(n, {fontSize: '12px', margin: '4px 0'})
  ], ui.Panel.Layout.flow('horizontal')));
});

panel.add(ui.Label('Data: Copernicus Sentinel-2 (ESA). Severity thresholds: USGS (Key & Benson, 2006). ' +
                   'Results are indicative and not an official damage assessment.',
                   {fontSize: '10px', color: '#888', margin: '14px 8px 4px 8px'}));

// ------------------------- LOGIC -------------------------
fireSelect.onChange(function(name) {
  var f = FIRES[name];
  preBox.setValue(f.preStart);
  fireBox.setValue(f.fireDate);
  postBox.setValue(f.postEnd);
  drawPanel.style().set('shown', name === CUSTOM);
  if (f.aoi) {
    leftMap.centerObject(f.aoi, f.zoom);
  }
  ui.url.set('fire', FIRE_NAMES.indexOf(name));
});

function status(msg, color) {
  results.clear();
  results.add(ui.Label(msg, {fontSize: '12px', color: color || '#444', margin: '2px 8px'}));
}

function getAOI() {
  var name = fireSelect.getValue();
  if (name !== CUSTOM) return FIRES[name].aoi;
  if (drawLayer.geometries().length() === 0) return null;
  return drawLayer.toGeometry();
}

function run() {
  var name = fireSelect.getValue();
  var aoi = getAOI();
  if (!aoi) { status('Please draw an area on the right map first.', '#b22222'); return; }

  var preStart = preBox.getValue(), fireDate = fireBox.getValue(), postEnd = postBox.getValue();
  var re = /^\d{4}-\d{2}-\d{2}$/;
  if (!re.test(preStart) || !re.test(fireDate) || !re.test(postEnd) ||
      !(preStart < fireDate && fireDate < postEnd)) {
    status('Dates must be YYYY-MM-DD and in order: pre start < fire date < post end.', '#b22222');
    return;
  }

  status('⏳ Processing Sentinel-2 imagery…');
  leftTag.setValue('BEFORE  ' + preStart + ' → ' + fireDate);
  rightTag.setValue('AFTER  ' + fireDate + ' → ' + postEnd);

  var r = analyse(aoi, preStart, fireDate, postEnd);
  var outline = ee.Image().paint(aoi, 0, 2);

  // Same layer set on both sides; visibility comes from the sidebar checkboxes
  function layersFor(side, composite) {
    var images = [composite, composite, r.dNBR, r.severity, outline];
    var vis = [rgbVis, swirVis, dnbrVis, {min: 1, max: 7, palette: SEV_PALETTE}, {palette: ['yellow']}];
    return LAYER_NAMES.map(function(n, i) {
      return ui.Map.Layer(images[i], vis[i], n, layerBoxes[side][i].getValue(), i === 3 ? 0.75 : 1);
    });
  }
  leftMap.layers().reset(layersFor('left', r.preMed));
  rightMap.layers().reset(layersFor('right', r.postMed));
  zoomToArea();

  // Statistics (scale 20 m, tileScale for large areas)
  var haImg = ee.Image.pixelArea().divide(10000);
  var stats = ee.Dictionary({
    nPre:  r.preCol.size(),
    nPost: r.postCol.size(),
    aoiHa: aoi.area(10).divide(10000),
    meanDNBR: r.dNBR.reduceRegion({reducer: ee.Reducer.mean(), geometry: aoi,
                                   scale: 20, maxPixels: 1e13, tileScale: 4}).get('dNBR'),
    burnedHa: haImg.updateMask(r.burned).reduceRegion({reducer: ee.Reducer.sum(), geometry: aoi,
                                   scale: 20, maxPixels: 1e13, tileScale: 4}).get('area'),
    groups: haImg.addBands(r.severity).reduceRegion({
      reducer: ee.Reducer.sum().group({groupField: 1, groupName: 'cls'}),
      geometry: aoi, scale: 20, maxPixels: 1e13, tileScale: 4}).get('groups')
  });

  stats.evaluate(function(s, err) {
    results.clear();
    if (err) { status('Error: ' + err, '#b22222'); return; }
    if (s.nPre === 0 || s.nPost === 0) {
      status('No cloud-free images found for ' + (s.nPre === 0 ? 'the pre-fire' : 'the post-fire') +
             ' period. Try widening the dates.', '#b22222');
      return;
    }
    function row(k, v) {
      results.add(ui.Panel([
        ui.Label(k, {fontSize: '12px', width: '170px', margin: '2px 8px'}),
        ui.Label(v, {fontSize: '12px', fontWeight: 'bold', margin: '2px 0'})
      ], ui.Panel.Layout.flow('horizontal')));
    }
    var d = s.meanDNBR;
    row('Images used (pre / post)', s.nPre + ' / ' + s.nPost);
    row('Analysed area', Math.round(s.aoiHa).toLocaleString() + ' ha');
    row('Burned area (dNBR ≥ 0.1)', Math.round(s.burnedHa || 0).toLocaleString() + ' ha');
    row('Burned share', ((s.burnedHa || 0) / s.aoiHa * 100).toFixed(1) + ' %');
    row('Mean dNBR (whole area)', d === null ? 'n/a' : d.toFixed(3));

    results.add(ui.Label('Area by severity class', {fontSize: '12px', fontWeight: 'bold', margin: '8px 8px 2px 8px'}));
    (s.groups || []).forEach(function(g) {
      if (g.cls >= 4) row('  ' + SEV_NAMES[g.cls - 1], Math.round(g.sum).toLocaleString() + ' ha');
    });

    // Download link for severity GeoTIFF
    r.severity.toByte().getDownloadURL({
      name: 'burn_severity', region: aoi, scale: 20, format: 'GEO_TIFF'
    }, function(url, e) {
      if (url) results.add(ui.Label('⬇ Download severity map (GeoTIFF)',
                                    {fontSize: '12px', margin: '8px 8px'}, url));
    });

    // NBR time series chart
    var ts = nbrTimeSeries(aoi, ee.Date(preStart), ee.Date(postEnd).advance(30, 'day'));
    results.add(ui.Chart.image.series(ts, aoi, ee.Reducer.mean(), 60).setOptions({
      title: 'Mean NBR over the area (drop = fire)',
      hAxis: {title: 'Date'}, vAxis: {title: 'NBR'},
      lineWidth: 1, pointSize: 3, legend: {position: 'none'},
      series: {0: {color: 'b22222'}}
    }));
  });
}
runButton.onClick(run);

// ------------------------- START -------------------------
ui.root.clear();
ui.root.add(panel);
ui.root.add(split);

var startIdx = parseInt(ui.url.get('fire', 0), 10);
if (!(startIdx >= 0 && startIdx < FIRE_NAMES.length)) startIdx = 0;
fireSelect.setValue(FIRE_NAMES[startIdx]);   // triggers onChange (dates + zoom)
if (FIRE_NAMES[startIdx] !== CUSTOM) run();
