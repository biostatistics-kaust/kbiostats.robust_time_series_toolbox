import math
import sequtils
import nimstats/nimstats
import nimstats/rits_cov_struct
###########################################################
##! class RobustInterruptedCandidateModel
##! <<Model>> RobustInterruptedCandidateModel
type RobustInterruptedCandidateModel* = ref object
  ##! RobustInterruptedCandidateModel: +after_model: OLSModel
  after_model*: OLSModel
  ##! RobustInterruptedCandidateModel: +before_model: OLSModel
  before_model*: OLSModel
  ##! RobustInterruptedCandidateModel: +after_residual_model: OLSModel
  after_residual_model*: OLSModel
  ##! RobustInterruptedCandidateModel: +before_residual_model: OLSModel
  before_residual_model*: OLSModel
  ##! RobustInterruptedCandidateModel: +loglikelihood: float
  loglikelihood*: float
  ##! RobustInterruptedCandidateModel: +change_point_candidate: int
  change_point_candidate*: int
  ##! RobustInterruptedCandidateModel: +parameters_mean: vector
  parameters_mean*: vector
  #
  ##! RobustInterruptedCandidateModel: -before_model_covariance_parameters: matrix
  before_model_covariance_parameters: matrix
  ##! RobustInterruptedCandidateModel: -existence_change_point_hypothesis: HypothesisScore[ChiSquare]
  existence_change_point_hypothesis: HypothesisScore[ChiSquare]

##! class FlattenEstimator
##! <<Model>> FlattenEstimator
type FlattenEstimator* = ref object
  ##! FlattenEstimator: +mean: float
  mean*: float
  ##! FlattenEstimator: +standard_deviation: float
  standard_deviation*: float
  ##! FlattenEstimator: +confidence_interval: seq[float]
  confidence_interval*: seq[float]
  ##! FlattenEstimator: +p_value: float
  p_value*: float
  ##! FlattenEstimator: +t_value: float
  t_value*: float
  ##! FlattenEstimator: +score: float
  score*: float

##! FlattenEstimator: flatten()
proc flatten[T](estimator: SimpleEstimator[T]): FlattenEstimator =
  new result
  result.mean = estimator.estimate
  result.standard_deviation = estimator.standard_deviation
  result.confidence_interval = estimator.confidence_interval
  result.p_value = estimator.p_value
  result.score = estimator.test_score
  result.t_value = estimator.test_score
  
##! class FlattenHypothesisTest
##! <<Model>> FlattenHypothesisTest
type FlattenHypothesisTest* = ref object
  ##! FlattenHypothesisTest: -score: float
  score: float
  ##! FlattenHypothesisTest: -p_value: float
  p_value: float
  ##! FlattenHypothesisTest: -null_confidence_interval: seq[float]
  null_confidence_interval: seq[float]
  ##! FlattenHypothesisTest: -distribution: string
  distribution: string

#@Internal@IgnoreInClassDiagram
#! HypothesisScore: flatten()
proc flatten[T](test: HypothesisScore[T]): FlattenHypothesisTest =
  new result
  result.score = test.test_score
  result.p_value = test.p_value
  result.null_confidence_interval = test.null_confidence_interval
  result.distribution = $(test.distribution)

##! class RITSResiduals
##! <<Model>> RITSResiduals
type RITSResiduals* = ref object
  ##! RITSResiduals:-values: seq[float]
  values: seq[float]
  ##! RITSResiduals: -autocorrelation: seq[float]
  autocorrelation: seq[float]
  ##! RITSResiduals: -autocorrelation_null_confidence_interval: seq[float]
  autocorrelation_null_confidence_interval: seq[float]

##! class RITSPeriod
##! <<Model>> RITSPeriod
##! RITSPeriod <-- RITSResiduals
##! RITSPeriod <-- FlattenEstimator
type RITSPeriod* = ref object
  ##! RITSPeriod: +slope: FlattenEstimator
  slope*: FlattenEstimator
  ##! RITSPeriod: +intercept: FlattenEstimator
  intercept*: FlattenEstimator
  ##! RITSPeriod: +level: FlattenEstimator
  level*: FlattenEstimator
  ##! RITSPeriod: +autocorrelation: FlattenEstimator
  autocorrelation*: FlattenEstimator
  ##! RITSPeriod: +noise: FlattenEstimator
  noise*: FlattenEstimator
  ##! RITSPeriod: +residuals: RITSResiduals
  residuals*: RITSResiduals

proc newRITSPeriod(): RITSPeriod =
  new result
  new result.slope
  new result.intercept
  new result.level
  new result.autocorrelation
  new result.noise
  new result.residuals

##! class RobustInterruptedModel
##! <<Model>> RobustInterruptedModel
##! RobustInterruptedModelSet <-- RITSPeriod
##! RobustInterruptedModelSet <-- FlattenHypothesisTest
type RobustInterruptedModel* = ref object
  ##! RobustInterruptedModel: +change_point_candidates: seq[int]
  change_point_candidates*: seq[int]
  ##! RobustInterruptedModel: +loglikelihood_candidates: seq[float]
  loglikelihood_candidates*: seq[float]

  ##! RobustInterruptedModel: +change_point_index: int
  change_point_index*: int
  ##! RobustInterruptedModel: +change_point_x: float
  change_point_x*: float

  ##! RobustInterruptedModel: +initial: RITSPeriod
  initial*: RITSPeriod
  ##! RobustInterruptedModel: +increment_change: RITSPeriod
  increment_change*: RITSPeriod

  existence_change_point_hypothesis*: FlattenHypothesisTest

##! class RobustInterruptedModelSet
##! <<Model>> RobustInterruptedModelSet
##! RobustInterruptedModelSet <-- FlattenHypothesisTest
##! RobustInterruptedModelSet <-- RobustInterruptedModel
type RobustInterruptedModelSet* = ref object
  models*: seq[RobustInterruptedModel]
  existence_change_point_hypothesis*: FlattenHypothesisTest

##! RITSModel: -flatten(v: vector, d: float=NegInf)
proc clean_nan(v: vector, d: float=NegInf): vector = v.mapIt(if it.classify == fcNan: d else: it)

##! RITSModel: -flatten()
proc flatten(model: RITSModel): RobustInterruptedModel =
  new result
  result.initial = newRITSPeriod()
  result.increment_change = newRITSPeriod()

  result.change_point_candidates = model.change_point_candidates
  result.loglikelihood_candidates = model.loglikelihood_candidates.clean_nan
  #
  result.change_point_index = model.change_point_index
  result.change_point_x = model.change_point_x
  #
  let
    acf1 = model.model.before_residual_model.residuals.clean_nan.autocorrelation_function(max_lags=100)
    acf2 = model.model.after_residual_model.residuals.clean_nan.autocorrelation_function(max_lags=100)
  result.initial.slope = model.model.before_model.coefficients[0].flatten
  result.initial.intercept = model.model.before_model.coefficients[1].flatten
  result.initial.autocorrelation = model.model.before_residual_model.coefficients[0].flatten
  result.initial.noise = model.model.before_residual_model.noise_variance.flatten
  result.initial.residuals.values = model.model.before_residual_model.residuals.clean_nan
  result.initial.residuals.autocorrelation = acf1.estimators_mean.clean_nan(0.0)
  result.initial.residuals.autocorrelation_null_confidence_interval = acf1.null_distribution.confidence_interval
  result.increment_change.slope = model.model.after_model.coefficients[0].flatten
  result.increment_change.intercept = model.model.after_model.coefficients[1].flatten
  result.increment_change.autocorrelation = model.model.after_residual_model.coefficients[0].flatten
  result.increment_change.noise = model.model.after_residual_model.noise_variance.flatten
  result.increment_change.residuals.values = model.model.after_residual_model.residuals.clean_nan
  result.increment_change.residuals.autocorrelation = acf2.estimators_mean.clean_nan(0.0)
  result.increment_change.residuals.autocorrelation_null_confidence_interval = acf2.null_distribution.confidence_interval
  #
  result.increment_change.level.mean = (
    result.increment_change.slope.mean * result.change_point_x + result.increment_change.intercept.mean
  )
  result.increment_change.level.standard_deviation = (
    result.increment_change.slope.standard_deviation * result.change_point_x + result.increment_change.intercept.standard_deviation
  )
  result.increment_change.level.confidence_interval = (
    result.increment_change.slope.confidence_interval * result.change_point_x + result.increment_change.intercept.confidence_interval
  )
  #
  result.existence_change_point_hypothesis = model.model.existence_change_point_hypothesis.flatten

##! RITSModel: +existence_change_point_hypothesis(test_scores: seq[float])
proc existence_change_point_hypothesis(test_scores: seq[float]): FlattenHypothesisTest {.exportc: "existence_change_point_hypothesis".} =
  echo "test_scores:", $test_scores
  let total_score = test_scores.len.toFloat * test_scores.min
  chisquare(dof=2.0*test_scores.len.toFloat).htest_score(total_score, test_type=oneTailed).flatten

when false:
  proc existence_change_point_hypothesis(test_scores: seq[float]): FlattenHypothesisTest {.exportc: "existence_change_point_hypothesis".} =
    let total_score = test_scores.len.toFloat * test_scores.min
    chisquare(dof=2.0*test_scores.len.toFloat).htest_score(total_score, test_type=oneTailed).flatten

##! RITSModel: +existence_change_point_hypothesis(x: vector, y: vector, change_point_candidates_start: int, change_point_candidates_end: int)
proc robust_interrupted_time_series(x: vector, y: vector, change_point_candidates_start: int, change_point_candidates_end: int, covariance_structure_type: cstring): RobustInterruptedModel  {.exportc: "robust_interrupted_time_series".} =
  echo "covariance_structure_type:::", covariance_structure_type
  let model = rits_model(x, y, change_point_candidates_start, change_point_candidates_end, $covariance_structure_type)
  echo model
  model.flatten
  
##! RITSModel: +robust_interrupted_time_series_approximated(sampling: int, x: vector, y: vector, change_point_candidates_start: int, change_point_candidates_end: int)
proc robust_interrupted_time_series_approximated(sampling: int, x: vector, y: vector, change_point_candidates_start: int, change_point_candidates_end: int, covariance_structure_type: cstring): RobustInterruptedModel  {.exportc: "robust_interrupted_time_series_approximated".} =
  let model = rits_model(x.subsample(sampling), y, (change_point_candidates_start.toFloat/sampling.toFloat).toInt, (change_point_candidates_end.toFloat/sampling.toFloat).toInt, $covariance_structure_type)
  echo model
  model.flatten
  