// =====================================================================
//  Burn analysis – Sentinel-2 (revised)
//  Requires imports (keep your existing ones): soil, tree, burned
//  FeatureCollections, each feature with property 'LC' (1/2/3).
// =====================================================================

// ================== USER INPUTS ==================
// Time series shows the fire between Aug 16 and Aug 21, 2019
var preStart  = ee.Date('2019-07-05');
var fireDate  = ee.Date('2019-08-17');   // pre window ends here, post starts here
var postEnd   = ee.Date('2019-10-05');

var geometry = ee.Geometry.Polygon(
  [[26.99059815747183,38.27096885531639],
   [26.962789014405423,38.261636895183976],
   [26.96484895092886,38.24761769501632],
   [26.945279553956205,38.24006777471134],
   [26.94218964917105,38.21336712677703],
   [26.958325818604642,38.202036618379246],
   [26.98785157544058,38.2120183492158],
   [26.989911511964017,38.193942318311926],
   [27.04827638012808,38.23035958219917],
   [27.070592359131986,38.27349718272393],
   [27.054112866944486,38.3125685123244],
   [27.004674390381986,38.313107280125934],
   [26.99059815747183,38.27096885531639]]
);

// =========== CLOUD MASK (works with QA60 or MSK_CLASSI_* bands) ===========
function maskS2clouds(img) {
  var mask = ee.Image(ee.Algorithms.If(
    img.bandNames().contains('QA60'),
    img.select('QA60').bitwiseAnd(1 << 10).eq(0)
       .and(img.select('QA60').bitwiseAnd(1 << 11).eq(0)),
    img.select('MSK_CLASSI_OPAQUE').eq(0)
       .and(img.select('MSK_CLASSI_CIRRUS').eq(0))
  ));
  return img.select('B.*')
            .updateMask(mask)
            .divide(10000)
            .copyProperties(img, ['system:time_start']);
}

// Cloud mask for L2A using the Scene Classification Layer
function maskS2sr(img) {
  var scl = img.select('SCL');
  var mask = scl.neq(3).and(scl.neq(8)).and(scl.neq(9)).and(scl.neq(10));
  return img.select('B.*')
            .updateMask(mask)
            .divide(10000)
            .copyProperties(img, ['system:time_start']);
}

// ================== IMAGE COLLECTIONS ==================
function l1c(start, end) {
  return ee.ImageCollection('COPERNICUS/S2_HARMONIZED')
    .filterBounds(geometry)
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
    .filterDate(start, end)
    .map(maskS2clouds);
}
var preCol  = l1c(preStart, fireDate);
var postCol = l1c(fireDate, postEnd);
print('Pre images:', preCol.size(), 'Post images:', postCol.size());

var srCol = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(geometry)
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20))
  .filterDate(preStart, postEnd.advance(50, 'day'))
  .map(maskS2sr);

// ============ NBR & NDVI TIME SERIES (one value per date) ============
// The AOI spans two tiles (T35SMC, T35SNC) -> mosaic images of the same day
var days = srCol.aggregate_array('system:time_start')
  .map(function(t) { return ee.Date(t).format('YYYY-MM-dd'); })
  .distinct();

var dailySR = ee.ImageCollection(days.map(function(d) {
  var day = ee.Date.parse('YYYY-MM-dd', d);
  var im = srCol.filterDate(day, day.advance(1, 'day')).mosaic();
  return im.set('system:time_start', day.millis());
}));

var indices = dailySR.map(function(im) {
  return im.normalizedDifference(['B8', 'B12']).rename('NBR')
    .addBands(im.normalizedDifference(['B8', 'B4']).rename('NDVI'))
    .copyProperties(im, ['system:time_start']);
});

print(ui.Chart.image.series(indices.select('NBR'), geometry, ee.Reducer.mean(), 10)
  .setOptions({title: 'Sentinel-2 NBR (time series)', hAxis: {title: 'Date'},
               vAxis: {title: 'NBR'}, lineWidth: 1, pointSize: 3,
               series: {0: {color: 'FF0000'}}}));
print(ui.Chart.image.series(indices.select('NDVI'), geometry, ee.Reducer.mean(), 10)
  .setOptions({title: 'Sentinel-2 NDVI (time series)', hAxis: {title: 'Date'},
               vAxis: {title: 'NDVI'}, lineWidth: 1, pointSize: 3,
               series: {0: {color: '00AA00'}}}));

// ================== COMPOSITES & NBR BEFORE / AFTER ==================
var preMed  = preCol.median().clip(geometry);
var postMed = postCol.median().clip(geometry);

// NBR = (NIR - SWIR2) / (NIR + SWIR2)  -> B8 and B12 (standard for burn severity)
var nbrPre  = preMed.normalizedDifference(['B8', 'B12']).rename('NBR_pre');
var nbrPost = postMed.normalizedDifference(['B8', 'B12']).rename('NBR_post');

// ================== dNBR & BURN SEVERITY (USGS classes) ==================
var dNBR = nbrPre.subtract(nbrPost).rename('dNBR');

// 1 high regrowth, 2 low regrowth, 3 unburned, 4 low, 5 mod-low, 6 mod-high, 7 high
var severity = ee.Image(3)
  .where(dNBR.lt(-0.25), 1)
  .where(dNBR.gte(-0.25).and(dNBR.lt(-0.1)), 2)
  .where(dNBR.gte(0.1).and(dNBR.lt(0.27)), 4)
  .where(dNBR.gte(0.27).and(dNBR.lt(0.44)), 5)
  .where(dNBR.gte(0.44).and(dNBR.lt(0.66)), 6)
  .where(dNBR.gte(0.66), 7)
  .updateMask(dNBR.mask())
  .clip(geometry)
  .rename('severity');

var sevNames   = ['High regrowth', 'Low regrowth', 'Unburned', 'Low severity',
                  'Moderate-low', 'Moderate-high', 'High severity'];
var sevPalette = ['#7a8737', '#acbe4d', '#0ae042', '#fff70b',
                  '#ffaf38', '#ff641b', '#a41fd6'];

// Mean dNBR over the AOI + textual class
var meanDNBR = ee.Number(dNBR.reduceRegion({
  reducer: ee.Reducer.mean(), geometry: geometry, scale: 20, maxPixels: 1e13
}).get('dNBR'));
print('Mean dNBR over AOI:', meanDNBR);

meanDNBR.evaluate(function(d) {
  var label =
    d > 0.66  ? 'High severity burn' :
    d > 0.44  ? 'Moderate-high severity burn' :
    d > 0.27  ? 'Moderate-low severity burn' :
    d > 0.1   ? 'Low severity burn' :
    d > -0.1  ? 'Unburned' :
    d > -0.25 ? 'Low post-fire regrowth' : 'High post-fire regrowth';
  print('AOI mean severity class: ' + label);
});

// Area per severity class (ha)
var sevArea = ee.Image.pixelArea().divide(10000).addBands(severity)
  .reduceRegion({
    reducer: ee.Reducer.sum().group({groupField: 1, groupName: 'class'}),
    geometry: geometry, scale: 20, maxPixels: 1e13
  });
print('Area (ha) by dNBR severity class (1-7):', sevArea.get('groups'));

// ================== SUPERVISED CLASSIFICATION ===========
function forceNumericLC(fc) {
  return ee.FeatureCollection(fc).map(function(f) {
    return f.set('LC', ee.Number.parse(ee.String(f.get('LC'))));
  });
}
var trainingFC = forceNumericLC(soil)            // LC = 1
  .merge(forceNumericLC(tree))                   // LC = 2
  .merge(forceNumericLC(burned));                // LC = 3

// RGB + NIR + SWIR + NBR separates burned areas far better than RGB alone
var classImage = postMed.select(['B2', 'B3', 'B4', 'B8', 'B11', 'B12'])
  .addBands(nbrPost.rename('NBR'))
  .addBands(dNBR);
var bands = classImage.bandNames();

var samples = classImage.sampleRegions({
  collection: trainingFC, properties: ['LC'], scale: 10, geometries: false
}).randomColumn('rand', 42);

var training   = samples.filter(ee.Filter.lt('rand', 0.7));
var validation = samples.filter(ee.Filter.gte('rand', 0.7));
print('Training size:', training.size(), 'Validation size:', validation.size());
print('LC histogram:', samples.aggregate_histogram('LC'));

var classifier = ee.Classifier.smileRandomForest({numberOfTrees: 200, seed: 42})
  .train({features: training, classProperty: 'LC', inputProperties: bands});

var classified = classImage.classify(classifier).clip(geometry);

// Accuracy on held-out samples
var errorMatrix = validation.classify(classifier)
  .errorMatrix('LC', 'classification');
print('Confusion matrix:', errorMatrix);
print('Overall accuracy:', errorMatrix.accuracy());
print('Kappa:', errorMatrix.kappa());

// 1 soil -> orange, 2 tree -> dark green, 3 burned -> dark red
var classPalette = ['#FFA500', '#0B6623', '#8B0000'];
var burnedOnly = classified.eq(3).selfMask();

// ================== BURNED AREA ==================
var burnedHa = ee.Number(ee.Image.pixelArea().divide(10000)
  .updateMask(burnedOnly)
  .reduceRegion({reducer: ee.Reducer.sum(), geometry: geometry,
                 scale: 10, maxPixels: 1e13})
  .get('area'));
print('Burned area – classification (ha):', burnedHa);

// ================== VIS PARAMS ==================
var rgbVis  = {bands: ['B4', 'B3', 'B2'], min: 0.0, max: 0.3};
var swirVis = {bands: ['B12', 'B8', 'B4'], min: 0.0, max: 0.4}; // burn scar pops out
var nbrVis  = {min: -0.5, max: 0.8, palette: ['black', 'white', 'green']};
var dnbrVis = {min: -0.3, max: 0.9, palette: ['green', 'white', 'yellow', 'orange', 'red', 'purple']};

// ================== SIDE-BY-SIDE UI ==================
var leftMap  = ui.Map();
var rightMap = ui.Map();

leftMap.add(ui.Label('BEFORE (' + preStart.format('YYYY-MM-dd').getInfo() + ' → ' +
                     fireDate.format('YYYY-MM-dd').getInfo() + ')',
                     {position: 'top-left', fontWeight: 'bold'}));
rightMap.add(ui.Label('AFTER (' + fireDate.format('YYYY-MM-dd').getInfo() + ' → ' +
                      postEnd.format('YYYY-MM-dd').getInfo() + ')',
                      {position: 'top-right', fontWeight: 'bold'}));

leftMap.addLayer(preMed,  rgbVis,  'RGB Before');
leftMap.addLayer(preMed,  swirVis, 'SWIR false colour Before', false);
leftMap.addLayer(nbrPre,  nbrVis,  'NBR Before', false);

rightMap.addLayer(postMed,  rgbVis,  'RGB After');
rightMap.addLayer(postMed,  swirVis, 'SWIR false colour After', false);
rightMap.addLayer(nbrPost,  nbrVis,  'NBR After', false);
rightMap.addLayer(dNBR,     dnbrVis, 'dNBR', false);
rightMap.addLayer(severity, {min: 1, max: 7, palette: sevPalette}, 'Burn severity (dNBR)', false);
rightMap.addLayer(classified, {min: 1, max: 3, palette: classPalette}, 'Classified After', false);
rightMap.addLayer(burnedOnly, {palette: ['#8B0000']}, 'Only Burned');
rightMap.addLayer(ee.Image().paint(geometry, 0, 2), {palette: ['yellow']}, 'AOI');
leftMap.addLayer(ee.Image().paint(geometry, 0, 2), {palette: ['yellow']}, 'AOI');

ui.Map.Linker([leftMap, rightMap]);
leftMap.centerObject(geometry, 13);

ui.root.widgets().reset([ui.SplitPanel({
  firstPanel: leftMap, secondPanel: rightMap, wipe: true
})]);

// ================== LEGENDS ==================
function makeLegend(title, names, colors, position) {
  var panel = ui.Panel({style: {position: position, padding: '8px'}});
  panel.add(ui.Label(title, {fontWeight: 'bold', margin: '0 0 6px 0'}));
  names.forEach(function(n, i) {
    panel.add(ui.Panel([
      ui.Label('', {backgroundColor: colors[i], padding: '8px', margin: '0 6px 4px 0'}),
      ui.Label(n, {margin: '0 0 4px 0'})
    ], ui.Panel.Layout.Flow('horizontal')));
  });
  return panel;
}
rightMap.add(makeLegend('Land cover', ['Soil', 'Tree', 'Burned'], classPalette, 'bottom-left'));
rightMap.add(makeLegend('Burn severity (dNBR)', sevNames, sevPalette, 'bottom-right'));

// ================== OPTIONAL EXPORTS ==================
// Export.image.toDrive({image: severity.toByte(), description: 'burn_severity',
//   region: geometry, scale: 20, maxPixels: 1e13});
// Export.image.toDrive({image: classified.toByte(), description: 'landcover_post',
//   region: geometry, scale: 10, maxPixels: 1e13});
