var fpfmt = function () {
  var d3fpfrm = d3.format(".5g");
  return function(value) {
    return d3fpfrm(value).replace(/0+$/gm,'');
  }
}();

function vectorToString(vector) {
  var string = ""
  for (var i = 0; i < vector.length; ++i) {
    if (i != 0) {
      string += " | ";
    }
    string += fpfmt(vector[i]);
  }
  return string;
}

function Handler(name, time, iterationsNumber) {
  this.name = name;
  this.time = time;
  this.iterationsNumber = iterationsNumber;

  this.Turek2D1Handler = function(d) {
    for (var i = 0; i < d.ibForce.length; ++i) {
      d.ibForce[i] = 500.0 * d.ibForce[i];
      d.force1[i] = 500.0 * d.force1[i];
      d.force2[i] = 500.0 * d.force2[i];
      d.force3[i] = 500.0 * d.force3[i];
    }
    return true;
  };
  this.Turek2D2Handler = function(d) {
    for (var i = 0; i < d.ibForce.length; ++i) {
      d.ibForce[i] = 20.0 * d.ibForce[i];
      d.force1[i] = 20.0 * d.force1[i];
      d.force2[i] = 20.0 * d.force2[i];
      d.force3[i] = 20.0 * d.force3[i];
    }
    return true;
  };
  this.handle = function(d) { return true; }

  this.iterationNumberSlice = Math.ceil(this.iterationsNumber / 11000);

  this.globalFilter = function(obj) {
    var doProcess = false;
    if (obj.iterationNumber == this.iterationsNumber) {
      doProcess = true;
    } else if (obj.iterationNumber % this.iterationNumberSlice == 0) {
      doProcess = true;
    }
    if (doProcess) {
      return this.handle(obj);
    }
    return false;
  }

  this.timeSlice = Math.min(0.1 * this.time, 0.1);

  this.preChartFilter = function(obj) {
    if (obj.time > this.timeSlice) {
      return true;
    }
    return false;
  }

  if (this.name == "Turek2D1"
      || this.name == "Template") {
    this.handle = this.Turek2D1Handler;
    this.timeSlice = Math.min(0.6 * this.time, 0.6);
    this.iterationNumberSlice = Math.ceil(this.iterationsNumber / 5000);
  } else if (this.name == "Turek2D2") {
    this.handle = this.Turek2D2Handler;
    this.timeSlice = Math.min(0.6 * this.time, 0.6);
    this.iterationNumberSlice = Math.ceil(this.iterationsNumber / 5000);
  }
}

function computeVectorN2(vector) {
  var norm = 0.0;
  vector.forEach(function(d) {
    norm += d * d;
  });
  norm = Math.sqrt(norm);
  return norm;
}

function parseVector(str) {
  vector = [];
  str.replace(/[-+]?([0-9]+\.?[0-9]*|\.[0-9]+)([eE][-+]?[0-9]+)?/g,
      function (value) {
        vector[vector.length] = parseFloat(value);
      });
  return vector;
}

function createInfoTable(data) {
  var table = d3.select("#info-table").append("table");
  var tbody = table.append("tbody");

  function createRow(name, value) {
    var tr = tbody.append("tr");
    tr.append("td")
      .attr("style", "text-align: left;")
      .html(name);
    tr.append("td")
      .attr("style", "text-align: left;")
      .html(value);
  }

  function wallTypeToString(type) {
    return type == 0 ? "Input" :
           type == 1 ? "ParabolicInput" :
           type == 2 ? "Output" : "Input";
  }

  createRow("Name", data.Name);
  createRow("Solver",
      data.SolverId == 1 ? "Improved Fractional-Step Finite-Difference"
      : "Simple Fractional-Step Finite-Difference");
  tr = tbody.append("tr").append("td")
    .attr("colspan", 2)
    .attr("style", "text-align: center; font-weight: bold;")
    .html("Parameters");
  createRow("Reynolds Number", data.Re);
  createRow("Gamma", data.Gamma);
  createRow("Tau", data.Tau);
  createRow("Environment force", data.G);
  createRow("IB outer layer size", data.OuterLayerSize);
  createRow("IB inner layer size", data.InnerLayerSize);
  createRow("Number of Processors", data.ProcessorSize.join(" "));
  createRow("Number of Cells in the Domain", data.GlobalCellSize);
  createRow("Number of Cells in a Subdomain", data.UniformLocalCellSize);
  createRow("Number of Cells in the Last Subdomain", data.LastLocalCellSize);
  createRow("Geometric Size of the Domain", data.Width);
  createRow("Geometric Size of a Cell", data.CellWidth);
  tr = tbody.append("tr").append("td")
    .attr("colspan", 2)
    .attr("style", "text-align: center; font-weight: bold;")
    .html("Domain Boundaries");
  createRow("Left Wall Type", wallTypeToString(data["00WallType"]));
  if (data["00WallType"] != 2) {
    createRow("Left Wall Velocity", data["00WallVelocity"]);
  }
  createRow("Right Wall Type", wallTypeToString(data["01WallType"]));
  if (data["01WallType"] != 2) {
    createRow("Right Wall Velocity", data["01WallVelocity"]);
  }
  createRow("Bottom Wall Type", wallTypeToString(data["10WallType"]));
  if (data["10WallType"] != 2) {
    createRow("Bottom Wall Velocity", data["10WallVelocity"]);
  }
  createRow("Top Wall Type", wallTypeToString(data["11WallType"]));
  if (data["11WallType"] != 2) {
    createRow("Top Wall Velocity", data["11WallVelocity"]);
  }
  if (data.ProcessorSize.length == 3) {
    createRow("Back Wall Type", wallTypeToString(data["20WallType"]));
    if (data["20WallType"] != 2) {
      createRow("Back Wall Velocity", data["20WallVelocity"]);
    }
    createRow("Fron Wall Type", wallTypeToString(data["21WallType"]));
    if (data["21WallType"] != 2) {
      createRow("Fron Wall Velocity", data["21WallVelocity"]);
    }
  }
  tr = tbody.append("tr").append("td")
    .attr("colspan", 2)
    .attr("style", "text-align: center; font-weight: bold;")
    .html("Statistics");
  createRow("Number of Iterations", data.iterationSize);
  createRow("Time", data.time);
  createRow("Time Step Size",
      "[" + fpfmt(data.minTimeStepSize) + ", " + fpfmt(data.maxTimeStepSize) + "]" +
       ", avr. " + fpfmt(data.avrTimeStepSize) );
}

function createIterationsTable(data, hasForces) {
  var table = d3.select("#iterations-table").append("table");
  var thead = table.append("thead");
  var tbody = table.append("tbody");

  var theadTr = thead.append("tr");
  theadTr.append("td").html("It. Num.");
  theadTr.append("td").html("t");
  theadTr.append("td").html("dt");
  theadTr.append("td").html("Min Velocity");
  theadTr.append("td").html("Max Velocity");
  theadTr.append("td").html("Min Pressure");
  theadTr.append("td").html("Max Pressure");
  theadTr.append("td").html("IB Force");
  theadTr.append("td").html("Force1");
  theadTr.append("td").html("Force2");
  theadTr.append("td").html("Force3");

  data.forEach(function(d) {
    var tr = tbody.append("tr");
    tr.append("td").html(d.iterationNumber);
    tr.append("td").html(d.time);
    tr.append("td").html(d.timeStepSize);
    tr.append("td").html(vectorToString(d.minVelocity));
    tr.append("td").html(vectorToString(d.maxVelocity));
    tr.append("td").html(d.minPressure);
    tr.append("td").html(d.maxPressure);
    if (hasForces) {
      tr.append("td").html(vectorToString(d.ibForce));
      tr.append("td").html(vectorToString(d.force1));
      tr.append("td").html(vectorToString(d.force2));
      tr.append("td").html(vectorToString(d.force3));
    }
  });
}

function buildChart(data, accessor) {
  // Set the dimensions of the canvas / graph
  var margin = {top: 30, right: 20, bottom: 70, left: 50},
      width = 700 - margin.left - margin.right,
      height = 300 - margin.top - margin.bottom;
      
  // Adds the svg canvas
  var svg = d3.select("#charts")
      .append("svg")
          .attr("class", "svg-graph")
          .attr("width", width + margin.left + margin.right)
          .attr("height", height + margin.top + margin.bottom)
      .append("g")
          .attr("transform", 
                "translate(" + margin.left + "," + margin.top + ")");

  // Set the ranges
  var x = d3.scale.linear().range([0, width]);
  var y = d3.scale.linear().range([height, 0]);

  // Scale the range of the data
  x.domain(d3.extent(data, function(d) { return accessor.x(d); }));
  y.domain([
      d3.min(data, function(d) { return Math.min.apply(null, accessor.y(d)); }),
      d3.max(data, function(d) { return Math.max.apply(null, accessor.y(d)); })]);

  // Define the axes
  var xAxis = d3.svg.axis().scale(x).orient("bottom").ticks(10);

  var yAxis = d3.svg.axis().scale(y).orient("left").ticks(10);

  var color = d3.scale.category10();
  var legendSpace = width / accessor.names.length;

  accessor.names.forEach(function(name, index) {
    // Add the curve
    tempColor = color(index);
    svg.append("path")
      .attr("class", "line")
      .style("stroke", tempColor)
      .attr("d",  d3.svg.line()
        .x(function(d2) { return x(accessor.x(d2)); })
        .y(function(d2) { return y(accessor.y(d2)[index]); })(data));
    var coeff = Math.ceil(10.0 / (width / data.length));
     if (index == 0) {
      svg.selectAll("dot").data(data) .enter().append("circle")
          .filter(function(d, i) { return (i % coeff) == 0; })
          .attr("r", 3.5)
          .attr("cx", function(d2) { return x(accessor.x(d2)); })
          .attr("cy", function(d2) { return y(accessor.y(d2)[index]); })
          .attr("fill-opacity", 0.0)
          .attr("stroke-width", 1.0)
          .attr("stroke", tempColor);
     } else if (index == 1) {
      svg.selectAll("dot").data(data).enter().append("rect")
          .filter(function(d, i) { return (i % coeff) == 0; })
          .attr("width", 4).attr("height", 4)
          .attr("x", function(d2) { return x(accessor.x(d2)) - 2.0; })
          .attr("y", function(d2) { return y(accessor.y(d2)[index]) - 2.0; })
          .attr("fill-opacity", 0.0)
          .attr("stroke-width", 1.0)
          .attr("stroke", tempColor);
     } else if (index == 2) {
      svg.selectAll("dot").data(data).enter().append("rect")
          .filter(function(d, i) { return (i % coeff) == 0; })
          .attr("width", 3.5).attr("height", 3.5)
          .attr("x", function(d2) { return x(accessor.x(d2)) - 1.75; })
          .attr("y", function(d2) { return y(accessor.y(d2)[index]) - 1.75; })
          .attr("fill", tempColor);
     } else if (index == 3) {
      svg.selectAll("dot").data(data).enter().append("circle")
          .filter(function(d, i) { return (i % coeff) == 0; })
          .attr("r", 2)
          .attr("cx", function(d2) { return x(accessor.x(d2)); })
          .attr("cy", function(d2) { return y(accessor.y(d2)[index]); })
          .attr("fill", tempColor);
     }
    // Add the Legend
    svg.append("text")
        .attr("x", (legendSpace/2) + index*legendSpace)
        .attr("y", height + (margin.bottom/2)+ 5)
        .attr("class", "legend")
        .style("fill", function() { return tempColor; })
        .text(name);
  });

  svg.append("g")
      .attr("class", "x axis")
      .attr("transform", "translate(0," + height + ")")
      .call(xAxis);

  svg.append("g")
      .attr("class", "y axis")
      .call(yAxis);
}

function main() {
  var url = window.location.pathname;
  var filename = url.substring(url.lastIndexOf('/') + 1);
  d3.csv(filename.substring(0, filename.indexOf('.')) + "-info.csv",
         function(error, data) {
    if (error != null) {
      console.log(error);
      d3.select("#info-table").append("div")
        .append("Error occured reading info CSV file");
      return;
    }
    var info = data[0];
    info.ProcessorSize = parseVector(info.ProcessorSize);
    proccedWithIterationResults(info);
  });
}

function proccedWithIterationResults(info) {
  var url = window.location.pathname;
  var filename = url.substring(url.lastIndexOf('/') + 1);
  filename = filename.substring(0, filename.indexOf('.'));

  info.iterationSize = 0;

  d3.csv(filename + ".csv",
         function(error, data) {
    if (error != null) {
      console.log(error);
      d3.select("#charts").append("div").append("Error occured");
      return;
    }

    var skipSize = 0;
    var hasForces = true;
    var weightCoeff = 1.0 / (data.length - skipSize);

    var rows = [];

    info.maxTimeStepSize = -Number.MAX_VALUE;
    info.minTimeStepSize = +Number.MAX_VALUE;

    info.avrTimeStepSize = 0.0;

    data.forEach(function(d, index) {
      if (index < 2) {
        if (typeof d.IbForceSum === "undefined"
            || typeof d.Force1 === "undefined"
            || typeof d.Force2 === "undefined"
            || typeof d.Force3 === "undefined") {
          hasForces = false;
        }
      } 
      if (index == (data.length - 1)) {
        if (typeof d.Time === "undefined"
            || typeof d.TimeStepSize === "undefined") {
          return;
        }
      }

      var iterationNubmer = parseInt(d.IterationNumber);
      if (iterationNubmer < skipSize) {
        return;
      }

      var it = new Object();
      it.iterationNumber = iterationNubmer;
      it.time = parseFloat(d.Time);
      it.timeStepSize = parseFloat(d.TimeStepSize);
      it.minVelocity = parseVector(d.MinVelocity);
      it.maxVelocity = parseVector(d.MaxVelocity);
      it.minPressure = parseFloat(d.MinPressure);
      it.maxPressure = parseFloat(d.MaxPressure);
      if (hasForces) {
        it.ibForce = parseVector(d.IbForceSum);
        it.force1 = parseVector(d.Force1);
        it.force2 = parseVector(d.Force2);
        it.force3 = parseVector(d.Force3);
      }

      info.maxTimeStepSize
        = Math.max(it.timeStepSize, info.maxTimeStepSize);
      info.minTimeStepSize
        = Math.min(it.timeStepSize, info.minTimeStepSize);
      info.avrTimeStepSize += weightCoeff * it.timeStepSize;

      rows.push(it);
    });

    if (info != null) {
      info.iterationSize = rows[rows.length - 1].iterationNumber;
      info.time = rows[rows.length - 1].time;
      createInfoTable(info);
    }

    var handler = new Handler(filename,
                              rows[rows.length - 1].time,
                              rows[rows.length - 1].iterationNumber);

    rows = rows.filter(function(obj){return handler.globalFilter(obj);});

    buildChart(rows, new function () {
      this.names = ["Time Step Size"]
      this.x = function(element) {
        return element.iterationNumber;
      }
      this.y = function(element) {
        return [element.timeStepSize];
      }
    });

    createIterationsTable(rows, hasForces);

    rows = rows.filter(function(obj) {return handler.preChartFilter(obj);});

    if (hasForces) {
      buildChart(rows, new function () {
        this.names = ["Max Velocity(X)", "Max Velocity(Y)"];
        this.x = function(element) {
          return element.iterationNumber;
        }
        this.y = function(element) {
          return [element.maxVelocity[0],
                  element.maxVelocity[1]];
        }
      });
      buildChart(rows, new function () {
        this.names = ["Min Velocity(X)", "Min Velocity(Y)"];
        this.x = function(element) {
          return element.iterationNumber;
        }
        this.y = function(element) {
          return [element.minVelocity[0],
                  element.minVelocity[1]];
        }
      });
      buildChart(rows, new function () {
        this.names = ["Max Pressure"];
        this.x = function(element) {
          return element.iterationNumber;
        }
        this.y = function(element) {
          return [element.maxPressure];
        }
      });
      buildChart(rows, new function () {
        this.names = ["Min Pressure"];
        this.x = function(element) {
          return element.iterationNumber;
        }
        this.y = function(element) {
          return [element.minPressure];
        }
      });
      buildChart(rows, new function () {
        this.names = ["Direct Forcing(X)", "Direct Forcing(Y)"];
        this.x = function(element) {
          return element.iterationNumber;
        }
        this.y = function(element) {
          return [element.ibForce[0],
                  element.ibForce[1]];
        }
      });
      buildChart(rows, new function () {
        this.names = ["Force1(X)", "Force1(Y)"];
        this.x = function(element) {
          return element.iterationNumber;
        }
        this.y = function(element) {
          return [element.force1[0],
                  element.force1[1]];
        }
      });
      buildChart(rows, new function () {
        this.names = ["Force2(X)", "Force2(Y)"];
        this.x = function(element) {
          return element.iterationNumber;
        }
        this.y = function(element) {
          return [element.force2[0],
                  element.force2[1]];
        }
      });
      buildChart(rows, new function () {
        this.names = ["Force3(X)", "Force3(Y)"];
        this.x = function(element) {
          return element.iterationNumber;
        }
        this.y = function(element) {
          return [element.force3[0],
                  element.force3[1]];
        }
      });
    }
  });
}
