//! class RTSModel
//! <<module>> RTSModel
//! RTSModel o-- RITSModel

var {
  robust_interrupted_time_series,
  robust_interrupted_time_series_approximated,
  existence_change_point_hypothesis
} = require("./model/model")

//@Unused
exports.robust_interrupted_time_series = robust_interrupted_time_series

//@Unused
exports.robust_interrupted_time_series_approximated = robust_interrupted_time_series_approximated

exports.existence_change_point_hypothesis = existence_change_point_hypothesis

console.log(":: LOADING RITS Models")

//! RTSModel: +fit_model(dates, values, theoretical_change_point, candidates_before, candidates_after)
const fit_model = (dates, values, theoretical_change_point, candidates_before, candidates_after, covariance_structure_type) => {
  let model = {};
  model.dates = dates;
  model.x = [...Array(model.dates.length).keys()];
  model.change_point = {};
  //Just for backcompatibility
  model.change_point.theoretical = theoretical_change_point;
  model.change_point.before = candidates_before;
  model.change_point.after = candidates_after;
  //End backcompatibility
  model.change_point.min = theoretical_change_point - candidates_before;
  model.change_point.max = theoretical_change_point + candidates_after;
  model.y = values;
  model.estimations = robust_interrupted_time_series(
    model.x, model.y,
    model.change_point.min, model.change_point.max,
    covariance_structure_type,
    );
  model.covariance_structure_type = covariance_structure_type;
  return model;
};
exports.fit_model = fit_model;

var {
  ModelDataSource,
  ModelDataSourceJson,
  ModelDataSourceTest,
  ModelDataSourceLargeTest
} = require("./model-datasources.js")
var {
  Stats,
  Model
} = require("./model-core.js")
var {
  CallbackEvent,
  ModelParameters
} = require("./model-parameters.js")
var {
  PlotTimeSeries,
  PlotLikelihood,
  PlotTimeSeriesEstimation,
  PlotPreResidualsHistogram,
  PlotPreResidualsHistogram,
  PlotPreResidualsACF,
  PlotPostResidualsHistogram,
  PlotPostResidualsACF,
  PlotCollector,
  PlotThumbnailsCollector
} = require("./model-plots.js")


//@Unused
exports.ModelDataSource = ModelDataSource;

exports.ModelDataSourceJson = ModelDataSourceJson;

//@Unused
exports.ModelDataSourceTest = ModelDataSourceTest;

exports.ModelDataSourceLargeTest = ModelDataSourceLargeTest;

//@Unused
exports.CallbackEvent = CallbackEvent;

//@Unused
exports.ModelParameters = ModelParameters;

//@Deprecated
exports.Stats = Stats;

//@Deprecated
exports.RTSModel = Model;

//exports.UnitModelPlotter = UnitModelPlotter;
exports.PlotTimeSeries = PlotTimeSeries;
exports.PlotLikelihood = PlotLikelihood;
exports.PlotTimeSeriesEstimation = PlotTimeSeriesEstimation;
exports.PlotPreResidualsHistogram = PlotPreResidualsHistogram;
exports.PlotPreResidualsACF = PlotPreResidualsACF;
exports.PlotPostResidualsHistogram = PlotPostResidualsHistogram;
exports.PlotPostResidualsACF = PlotPostResidualsACF;
exports.PlotCollector = PlotCollector;
exports.PlotThumbnailsCollector = PlotThumbnailsCollector;

