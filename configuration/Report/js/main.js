function parseVector(str) {
  vector = [];
  str.replace(/[-+]?([0-9]+\.?[0-9]*|\.[0-9]+)([eE][-+]?[0-9]+)?/g,
      function (value) {
        vector[vector.length] = parseFloat(value);
      });
  return vector;
}

function createChart(elementId) {
  // Set the dimensions of the canvas / graph
  var margin = {top: 30, right: 20, bottom: 30, left: 50},
      width = 600 - margin.left - margin.right,
      height = 270 - margin.top - margin.bottom;

  // Set the ranges
  var x = d3.scale.linear().range([0, width]);
  var y = d3.scale.linear().range([height, 0]);

  // Define the axes
  var xAxis = d3.svg.axis().scale(x)
      .orient("bottom").ticks(5);

  var yAxis = d3.svg.axis().scale(y)
      .orient("left").ticks(5);
      
  // Adds the svg canvas
  var svg = d3.select(elementId)
      .append("svg")
          .attr("width", width + margin.left + margin.right)
          .attr("height", height + margin.top + margin.bottom)
      .append("g")
          .attr("transform", 
                "translate(" + margin.left + "," + margin.top + ")");
  return {x: x,
          y: y,
          height: height,
          xAxis: xAxis,
          yAxis: yAxis,
          svg: svg};
}


function buildCharts(data) {

  var chart = createChart("#charts");

  // Define the line
  var valueline = d3.svg.line()
      .x(function(d) { return chart.x(d.iterationNumber); })
      .y(function(d) { return chart.y(d.timeStepSize); });

  // Scale the range of the data
  chart.x.domain(d3.extent(data, function(d) { return d.iterationNumber; }));
  chart.y.domain(d3.extent(data, function(d) { return d.timeStepSize; }));

  // Add the valueline path.
  chart.svg.append("path")
      .attr("class", "line")
      .attr("d", valueline(data));

  // Add the X Axis
  chart.svg.append("g")
      .attr("class", "x axis")
      .attr("transform", "translate(0," + chart.height + ")")
      .call(chart.xAxis);

  // Add the Y Axis
  chart.svg.append("g")
      .attr("class", "y axis")
      .call(chart.yAxis);
}

function main() {
  var url = window.location.pathname;
  var filename = url.substring(url.lastIndexOf('/') + 1);
  d3.csv(filename.substring(0, filename.indexOf('.')) + ".csv").row(function(d) {
    return {
      iterationNumber: parseInt(d.IterationNumber),
      time: parseFloat(d.Time),
      timeStepSize:  parseFloat(d.TimeStepSize),
      force1: typeof something === "undefined" ? [ 0, 0, 0 ] : parseVector(d.Force1),
      force2: typeof something === "undefined" ? [ 0, 0, 0 ] : parseVector(d.Force2),
      force3: typeof something === "undefined" ? [ 0, 0, 0 ] : parseVector(d.Force3),
    };
  }).get(function(error, rows) {
    if (error != null) {
      console.log(error);
      d3.select("#charts").append("div").append("Error occured");
    } else {
      buildCharts(rows);
    }
  });
}
