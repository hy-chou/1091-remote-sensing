var lonROI = 121.55330608757961;
var latROI = 25.082922283113202;
var zoomROI = 10;
var roi = ee
  // add .shp to asset from https://data.gov.tw/dataset/7442/
  .FeatureCollection("users/<username>/TaiwanCounty")
  .filter("COUNTYID == 'A'")
  .geometry();

function getL5img(year) {
  function imgMask(image) {
    var qaBand = image.select("pixel_qa");
    var water = qaBand.bitwiseAnd(1 << 2);
    var cloud = qaBand.bitwiseAnd(1 << 3).or(qaBand.bitwiseAnd(1 << 5));
    var inRange = image.lt(10000).and(image.gte(0));

    return image.updateMask(cloud.not().and(water.not())).updateMask(inRange);
  }

  var imcL5 = ee
    .ImageCollection("LANDSAT/LT05/C01/T1_SR")
    .filterDate(year.toString(), (year + 1).toString())
    .filterBounds(ee.Geometry.Point(lonROI, latROI))
    .map(imgMask)
    .select("B1", "B2", "B3", "B4");

  var imgL5 = imcL5.median().divide(10000);

  return imgL5;
}

function getL8img(year) {
  function imgMask(image) {
    var qaBand = image.select("pixel_qa");
    var water = qaBand.bitwiseAnd(1 << 2);
    var cloud = qaBand.bitwiseAnd(1 << 3).or(qaBand.bitwiseAnd(1 << 5));
    var inRange = image.lt(10000).and(image.gte(0));

    return image.updateMask(cloud.not().and(water.not())).updateMask(inRange);
  }

  var imcL8 = ee
    .ImageCollection("LANDSAT/LC08/C01/T1_SR")
    .filterDate(year.toString(), (year + 1).toString())
    .filterBounds(ee.Geometry.Point(lonROI, latROI))
    .map(imgMask)
    .select("B2", "B3", "B4", "B5");

  var imgL8 = imcL8.median().divide(10000);

  return imgL8;
}

function getL5evi(year) {
  function evimask(image) {
    var lim = image.lt(1).and(image.gte(-1));

    return image.updateMask(lim);
  }

  var imgL5 = getL5img(year);

  var eviL5 = imgL5.expression(
    "2.5 * (NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1)",
    {
      BLUE: imgL5.select("B1"),
      RED: imgL5.select("B3"),
      NIR: imgL5.select("B4"),
    }
  );
  eviL5 = eviL5.updateMask(eviL5.lt(1).and(eviL5.gte(-1)));

  return eviL5;
}

function getL8evi(year) {
  var evimask = function (image) {
    var lim = image.lt(1).and(image.gte(-1));

    return image.updateMask(lim);
  };

  var imgL8 = getL8img(year);

  var eviL8 = imgL8.expression(
    "2.5 * (NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1)",
    {
      BLUE: imgL8.select("B2"),
      RED: imgL8.select("B4"),
      NIR: imgL8.select("B5"),
    }
  );
  eviL8 = eviL8.updateMask(eviL8.lt(1).and(eviL8.gte(-1)));

  return eviL8;
}

function drawL5img(year) {
  var vpL5img = {
    bands: ["B3", "B2", "B1"],
    min: 0,
    max: 0.3,
    gamma: 1.4,
  };

  Map.setCenter(lonROI, latROI, zoomROI);
  Map.addLayer(getL5img(year).clip(roi), vpL5img, "imgL5" + year.toString());
}

function drawL8img(year) {
  var vpL8img = {
    bands: ["B4", "B3", "B2"],
    min: 0,
    max: 0.3,
    gamma: 1.4,
  };

  Map.setCenter(lonROI, latROI, zoomROI);
  Map.addLayer(getL8img(year).clip(roi), vpL8img, "imgL8" + year.toString());
}

function drawL5evi(year) {
  var vpL5evi = {
    min: 0,
    max: 1.0,
    palette: ["white", "green", "black"],
  };

  Map.setCenter(lonROI, latROI, zoomROI);
  Map.addLayer(getL5evi(year).clip(roi), vpL5evi, "eviL5" + year.toString());
}

function drawL8evi(year) {
  var vpL8evi = {
    min: 0,
    max: 1.0,
    palette: ["white", "green", "black"],
  };

  Map.setCenter(lonROI, latROI, zoomROI);
  Map.addLayer(getL8evi(year).clip(roi), vpL8evi, "eviL8" + year.toString());
}

function getL5stats(year) {
  function evimask(image) {
    var lim = image.lt(1).and(image.gte(-1));

    return image.updateMask(lim);
  }

  var rdc = ee.Reducer.mean()
    .combine({
      reducer2: ee.Reducer.percentile([0, 1, 25, 50, 75, 99, 100]),
      sharedInputs: true,
    })
    .combine({ reducer2: ee.Reducer.count(), sharedInputs: true })
    .combine({ reducer2: ee.Reducer.stdDev(), sharedInputs: true });

  var stats = evimask(getL5evi(year)).clip(roi).reduceRegion({
    reducer: rdc,
    geometry: roi,
    scale: 30,
    bestEffort: true,
  });

  return year, stats;
}

function getL8stats(year) {
  function evimask(image) {
    var lim = image.lt(1).and(image.gte(-1));

    return image.updateMask(lim);
  }

  var rdc = ee.Reducer.mean()
    .combine({
      reducer2: ee.Reducer.percentile([0, 1, 25, 50, 75, 99, 100]),
      sharedInputs: true,
    })
    .combine({ reducer2: ee.Reducer.count(), sharedInputs: true })
    .combine({ reducer2: ee.Reducer.stdDev(), sharedInputs: true });

  var stats = evimask(getL8evi(year)).clip(roi).reduceRegion({
    reducer: rdc,
    geometry: roi,
    scale: 30,
    bestEffort: true,
  });

  return year, stats;
}

// TODO

var evis = getL5evi(1987).addBands(getL8evi(2020));

var his = ui.Chart.image
  .histogram(evis, roi, 30)
  .setSeriesNames(["evi1987", "evi2020"])
  .setOptions({
    vAxis: { title: "count" },
    hAxis: { title: "EVI" },
    series: {
      0: { color: "orange" },
      1: { color: "blue" },
    },
  });
print(his);
