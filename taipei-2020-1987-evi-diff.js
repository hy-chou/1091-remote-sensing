var TPx = 121.55330608757961;
var TPy = 25.082922283113202;
var roi = ee
  // add .shp to asset from https://data.gov.tw/dataset/7442/
  .FeatureCollection("users/<username>/TaiwanCounty")
  .filter("COUNTYID == 'A'")
  .geometry();

// var year = 2020;
// var year = 1987;

function l5evi(year) {
  var mask = function (image) {
    var qaBand = image.select("pixel_qa");
    var water = qaBand.bitwiseAnd(1 << 2);
    var cloud = qaBand.bitwiseAnd(1 << 3).or(qaBand.bitwiseAnd(1 << 5));
    var lim = image.lt(10000).and(image.gte(0));

    return image.updateMask(cloud.not().and(water.not())).updateMask(lim);
  };
  var evimask = function (image) {
    var lim = image.lt(1).and(image.gte(-1));

    return image.updateMask(lim);
  };

  // var vpl8 = {
  //   bands: ['B4', 'B3', 'B2'],
  //   min: 0,
  //   max: 0.3,
  //   gamma: 1.8,
  // };
  // var vpl7 = {
  //   bands: ['B3', 'B2', 'B1'],
  //   min: 0,
  //   max: 0.3,
  //   gamma: 1.4,
  // };
  var vpl5 = {
    bands: ["B3", "B2", "B1"],
    min: 0,
    max: 0.3,
    gamma: 1.4,
  };

  // var imcl8 = ee.ImageCollection("LANDSAT/LC08/C01/T1_SR")
  //             .filterDate(year.toString(), (year + 1).toString())
  //             .filterBounds(ee.Geometry.Point(TPx, TPy))
  //             .map(mask)
  //             .select("B2", "B3", "B4", "B5");
  // var imcl7 = ee.ImageCollection("LANDSAT/LE07/C01/T1_SR")
  //             .filterDate(year.toString(), (year + 1).toString())
  //             .filterBounds(ee.Geometry.Point(TPx, TPy))
  //             .map(mask)
  //             .select("B1", "B2", "B3", "B4");
  var imcl5 = ee
    .ImageCollection("LANDSAT/LT05/C01/T1_SR")
    .filterDate(year.toString(), (year + 1).toString())
    .filterBounds(ee.Geometry.Point(TPx, TPy))
    .map(mask)
    .select("B1", "B2", "B3", "B4");

  // var imgl8 = imcl8.median().divide(10000);
  // var imgl7 = imcl7.median().divide(10000);
  var imgl5 = imcl5.median().divide(10000);

  // var evil8 = imgl8.expression(
  //     '2.5 * (NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1)', {
  //       'NIR': imgl8.select('B5'),
  //       'RED': imgl8.select('B4'),
  //       'BLUE': imgl8.select('B2'),
  //   });
  // var evil7 = imgl7.expression(
  //     '2.5 * (NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1)', {
  //       'NIR': imgl7.select('B4'),
  //       'RED': imgl7.select('B3'),
  //       'BLUE': imgl7.select('B1'),
  //   });
  var evil5 = imgl5.expression(
    "2.5 * (NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1)",
    {
      NIR: imgl5.select("B4"),
      RED: imgl5.select("B3"),
      BLUE: imgl5.select("B1"),
    }
  );

  Map.setCenter(TPx, TPy, 11.75);
  Map.addLayer(imgl5.clip(roi), vpl5, "imgmed" + year.toString());
  Map.addLayer(
    evil5.clip(roi),
    {
      min: 0,
      max: 1.0,
      // palette: [
      //   'FFFFFF', 'CE7E45', 'DF923D', 'F1B555',
      //   'FCD163', '99B718', '74A901', '66A000',
      //   '529400', '3E8601', '207401', '056201',
      //   '004C00', '023B01', '012E01', '011D01', '011301'
      // ],
      palette: ["white", "green", "black"],
      // palette: ["white", "black"],
    },
    "evi" + year.toString()
  );

  var rdc1 = ee.Reducer.mean();
  var rdc = rdc1
    .combine({
      reducer2: ee.Reducer.percentile([1, 25, 50, 75, 99]),
      sharedInputs: true,
    })
    .combine({ reducer2: ee.Reducer.count(), sharedInputs: true })
    .combine({ reducer2: ee.Reducer.stdDev(), sharedInputs: true });
  var stats = evimask(evil5).clip(roi).reduceRegion({
    reducer: rdc,
    geometry: roi,
    scale: 30,
    bestEffort: true,
  });

  print(year, stats);
  return evil5;
}

function l8evi(year) {
  var mask = function (image) {
    var qaBand = image.select("pixel_qa");
    var water = qaBand.bitwiseAnd(1 << 2);
    var cloud = qaBand.bitwiseAnd(1 << 3).or(qaBand.bitwiseAnd(1 << 5));
    var lim = image.lt(10000).and(image.gte(0));

    return image.updateMask(cloud.not().and(water.not())).updateMask(lim);
  };
  var evimask = function (image) {
    var lim = image.lt(1).and(image.gte(-1));

    return image.updateMask(lim);
  };

  var vpl8 = {
    bands: ["B4", "B3", "B2"],
    min: 0,
    max: 0.3,
    gamma: 1.8,
  };

  var imcl8 = ee
    .ImageCollection("LANDSAT/LC08/C01/T1_SR")
    .filterDate(year.toString(), (year + 1).toString())
    .filterBounds(ee.Geometry.Point(TPx, TPy))
    .map(mask)
    .select("B2", "B3", "B4", "B5");

  var imgl8 = imcl8.median().divide(10000);

  var evil8 = imgl8.expression(
    "2.5 * (NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1)",
    {
      NIR: imgl8.select("B5"),
      RED: imgl8.select("B4"),
      BLUE: imgl8.select("B2"),
    }
  );

  Map.setCenter(TPx, TPy, 11.75);
  Map.addLayer(imgl8.clip(roi), vpl8, "imgmed" + year.toString());
  Map.addLayer(
    evil8.clip(roi),
    {
      min: 0,
      max: 1.0,
      // palette: [
      //   'FFFFFF', 'CE7E45', 'DF923D', 'F1B555',
      //   'FCD163', '99B718', '74A901', '66A000',
      //   '529400', '3E8601', '207401', '056201',
      //   '004C00', '023B01', '012E01', '011D01', '011301'
      // ],
      palette: ["white", "green", "black"],
      // palette: ["white", "black"],
    },
    "evi" + year.toString()
  );

  var rdc1 = ee.Reducer.mean();
  var rdc = rdc1
    .combine({
      reducer2: ee.Reducer.percentile([1, 25, 50, 75, 99]),
      sharedInputs: true,
    })
    .combine({ reducer2: ee.Reducer.count(), sharedInputs: true })
    .combine({ reducer2: ee.Reducer.stdDev(), sharedInputs: true });
  var stats = evimask(evil8).clip(roi).reduceRegion({
    reducer: rdc,
    geometry: roi,
    scale: 30,
    bestEffort: true,
  });

  print(year, stats);
  return evil8;
}

var diff = l8evi(2020).subtract(l5evi(1987)).divide(2);
Map.addLayer(
  diff.clip(roi),
  {
    min: -0.1,
    max: 0.1,
    // palette: [
    //   'FFFFFF', 'CE7E45', 'DF923D', 'F1B555',
    //   'FCD163', '99B718', '74A901', '66A000',
    //   '529400', '3E8601', '207401', '056201',
    //   '004C00', '023B01', '012E01', '011D01', '011301'
    // ],
    // palette: ["white", "green", "black"],
    palette: ["white", "black"],
  },
  "diff"
);
