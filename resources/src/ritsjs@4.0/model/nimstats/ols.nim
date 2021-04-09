import math
import stats
import strutils
import sequtils
import strformat
import sugar
import core
import distributions


when isMainModule:
  echo "[Initializing]"
  let X1 = @[
    @[1.0, 2.0, 3.0],
    @[4.0, 5.0, 6.0],
    @[7.0, 8.0, 9.0],
    @[10.0, 11.0, 12.0],
  ]

###########################################################
# UTILS
###########################################################
type data_transformer = proc (y: vector): vector
type data_transformer_with_index = proc (y: vector, idx: int): vector

proc data_transform*(X: matrix, dt: data_transformer): matrix =
  for r in 1 .. X.rows:
    result.add dt(X[r - 1])

proc data_transform*(X: matrix, dt: data_transformer_with_index): matrix =
  for r in 0 .. (X.rows - 1):
    result.add dt(X[r], r)

when isMainModule:
  echo "Testing ols function: data_transform"
  assert(data_transform(X1, (x) => x) == @[@[1.0, 2.0, 3.0], @[4.0, 5.0, 6.0], @[7.0, 8.0, 9.0], @[10.0, 11.0, 12.0]])
  assert(data_transform(X1, (x) => x .^ 2) == @[@[1.0, 4.0, 9.0], @[16.0, 25.0, 36.0], @[49.0, 64.0, 81.0], @[100.0, 121.0, 144.0]])
  assert(data_transform(X1, (x) => hstack(x, 1)) == @[@[1.0, 2.0, 3.0, 1.0], @[4.0, 5.0, 6.0, 1.0], @[7.0, 8.0, 9.0, 1.0], @[10.0, 11.0, 12.0, 1.0]])
  echo "  [OK]"

proc add_intercept*(X: matrix): matrix =
  for r in 1 .. X.rows:
    result.add hstack(X[r - 1], 1)

proc add_intercept*(x: vector): matrix =
  for r in x.low .. x.high:
    result.add @[x[r], 1]

when isMainModule:
  echo "Testing ols function: add_intercept"
  assert(add_intercept(X1) == @[@[1.0, 2.0, 3.0, 1.0], @[4.0, 5.0, 6.0, 1.0], @[7.0, 8.0, 9.0, 1.0], @[10.0, 11.0, 12.0, 1.0]])
  assert(add_intercept(@[1.0, 2.0, 9.0, 12.0]) == @[@[1.0, 1.0], @[2.0, 1.0], @[9.0, 1.0], @[12.0, 1.0]])
  echo "  [OK]"
    
###########################################################
# Estimators
###########################################################
type BaseEstimator*[T] = ref object
  distribution: T

proc estimator*[T](distribution: T): BaseEstimator[T] =
  new result
  result.distribution = distribution

proc estimate*[T](sd: BaseEstimator[T]): float = sd.distribution.mean

#proc standard_deviation*[T](sd: Estimator[T]): float = sd.distribution.variance.sqrt
proc standard_deviation*[T](sd: BaseEstimator[ScaledShiftedDistribution[T]]): float = sd.distribution.scale

proc confidence_interval*[T](sd: BaseEstimator[T], alpha=0.05): seq[float] = @[sd.distribution.inv(0.5 * alpha), sd.distribution.inv(1.0 - 0.5 * alpha)]

#proc score*[V](sd: BaseEstimator[V]): float = NaN

#proc score*[V](sd: BaseEstimator[ScaledShiftedDistribution[V]]): float = sd.distribution.location / sd.distribution.scale

#proc p_value*[V](sd: BaseEstimator[V]): float = 1 - sd.distribution.distribution.cdf(sd.score.abs) + sd.distribution.distribution.cdf(-sd.score.abs)

#proc t_score*(sd: BaseEstimator[ScaledShiftedDistribution[StudentT]]): float = sd.score

###########################################################
# Hypothesis score
###########################################################
type HypothesisTestType* = enum
  #twoTailed, oneTailedLeft, oneTailedRight
  twoTailed, oneTailed, rightTailed, leftTailed
  #oneTailedLeft = -1, oneTailedRight=1, twoTailed = 90

type HypothesisScore*[T] = ref object
  distribution*: T
  test_score*: float
  test_type*: HypothesisTestType

proc htest_score*[T](distribution: T, score: float, test_type: HypothesisTestType = HypothesisTestType.twoTailed): HypothesisScore[T] =
  new result
  result.distribution = distribution
  result.test_score = score
  result.test_type = test_type

proc null_confidence_interval*[T](sd: HypothesisScore[T], alpha=0.05): seq[float] =
  let
    default_ci = @[sd.distribution.inv(0.5 * alpha), sd.distribution.inv(1.0 - 0.5 * alpha)]
    right_ci = @[sd.distribution.inv(0), sd.distribution.inv(1.0 - alpha)]
    left_ci = @[sd.distribution.inv(1 - alpha), sd.distribution.inv(1.0)]
    test_score = sd.distribution.mean
  if test_score < default_ci[0]:
    return right_ci
  elif test_score > default_ci[1]:
    return left_ci
  return default_ci

proc p_value*[V](sd: HypothesisScore[V]): float =
  if sd.test_type == HypothesisTestType.twoTailed:
    return 1 - sd.distribution.cdf(sd.test_score.abs) + sd.distribution.cdf(-sd.test_score.abs)
  if sd.test_type == HypothesisTestType.leftTailed:
    return sd.distribution.cdf(sd.test_score)
  if sd.test_type == HypothesisTestType.rightTailed:
    return 1 - sd.distribution.cdf(sd.test_score)
  let p1 = 1 - sd.distribution.cdf(sd.test_score)
  let p2 = sd.distribution.cdf(sd.test_score)
  return (if p1 < p2 and p1 >= 0: p1 else: p2)

import typetraits
proc `$`*(d: HypothesisScore): string =
  result &= &"[Hypothesis test]\n* Null distribution: {d.distribution} \n"
  result &= &"      Null C.I. 95%| Test value|    p-value| \n"
  result &= &"{d.null_confidence_interval[0]:>9.5f} {d.null_confidence_interval[1]:>9.5f}|  {d.test_score:>9.5f}| {d.p_value:>9.5e}|"
  

###########################################################
type SimpleEstimator*[T] = ref object
  distribution*: ScaledShiftedDistribution[T]
  hypothesis_test*: HypothesisScore[ScaledShiftedDistribution[T]]

proc shifted_estimator*[T](distribution: T, location=0.0, scale=1.0): SimpleEstimator[T] =
  let est_distribution = distribution.scaled_shifted(location=location, scale=scale)
  let est_null_distribution = distribution.scaled_shifted(location=0.0, scale=scale)
  new result
  result.distribution = est_distribution
  result.hypothesis_test = est_null_distribution.htest_score(score=location, test_type=HypothesisTestType.twoTailed)

proc null_distribution*[T](d: SimpleEstimator[T]): float = d.hypothesis_test.distribution

proc estimate*[T](d: SimpleEstimator[T]): float = d.distribution.mean

proc standard_deviation*[T](sd: SimpleEstimator[T]): float = sd.distribution.scale

proc total_standard_deviation*[T](sd: SimpleEstimator[T]): float = sd.distribution.variance ^ 0.5

proc confidence_interval*[T](sd: SimpleEstimator[T], alpha=0.05): seq[float] = @[sd.distribution.inv(0.5 * alpha), sd.distribution.inv(1.0 - 0.5 * alpha)]

proc p_value*[T](d: SimpleEstimator[T]): float = d.hypothesis_test.p_value

proc test_score*[T](d: SimpleEstimator[T]): float = d.distribution.location / d.distribution.scale
  
proc `$`*(d: SimpleEstimator): string =
  result &= &" Estimate| Std. dev.| Confidence int. 95%|   t-value|   p-value| \n"
  result &= &"{d.estimate:>9.5f}| {d.standard_deviation:>9.5f}| {d.confidence_interval[0]:>9.5f} {d.confidence_interval[1]:>9.5f}| {d.test_score:>9.5f}|{d.p_value:>9.5e}"
  

###########################################################
# AUTOCORRELATION FUNCTION
###########################################################
#https://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/autoband.htm
type AutoCorrelationFunctionData = ref object
  estimators*: seq[SimpleEstimator[Normal]] #Under the assumption that the data are from a moving average proces
  null_distribution*: SimpleEstimator[Normal] #to assess independence (i.e., randomness) of a set of observations.
  estimators_mean*: vector
  lags*: vector

proc autocorrelation_function*(x: vector, max_lags: int=40): AutoCorrelationFunctionData =
  new result
  # Ineficient but easy to debug
  let
    x_mean = x.mean
    N = x.high
  for k in 0..min(max_lags, N):
    let acf = ((x[0..(N-k)] - x_mean) * (x[k..N] - x_mean)).sum / ((x - x.mean) .^ 2).sum
    result.lags.add k.toFloat
    result.estimators_mean.add acf
    result.estimators.add normal().shifted_estimator(location=acf, scale=(1.0 + (1.0 + 2.0 * (x[0..k] .^ 2).sum) / N.toFloat).sqrt)
  result.null_distribution = normal().shifted_estimator(scale=1.0/N.toFloat.sqrt)

###########################################################
# OLS
###########################################################
type OLSModel* = ref object
  residuals*: vector
  sum_squared_errors*: float
  degrees_of_freedom*: float
  variance_matrix_coefficients*: matrix
  #
  include_intercept*: bool #Because it could be informative. Maybe.
  #
  loglikelihood*: float
  #
  R2*: float#Coefficient of determination
  adjustedR2*: float#Coefficient of determination
  beta_hat*: vector
  #
  coefficients*: seq[SimpleEstimator[StudentT]]
  noise_variance*: SimpleEstimator[InvChiSquare]
  #https://stats.stackexchange.com/questions/256726/linear-regression-what-does-the-f-statistic-r-squared-and-residual-standard-err
  model_significance*: HypothesisScore[CentralF]
  #
  feature_names*: seq[string]

#http://www3.grips.ac.jp/~yamanota/Lecture%20Note%204%20to%207%20OLS.pdf
#https://newonlinecourses.science.psu.edu/stat414/node/280/
#http://global.oup.com/booksites/content/0199268010/samplesec3
proc ols_model*(X: matrix, y: vector, names: seq[string] = @[]): OLSModel =
  new result
  if names.len > 0 and names.len != X.cols:
    raise newException(ValueError, "incorrect number of feature names")
  let
    Y = y.as_column_vector
    XpX = (X.T * X).inverse
    beta_hat = (XpX * X.T) * Y
    Ypred = X * beta_hat
    residuals = (Y - Ypred).T[0]
    sse = norm(residuals) ^ 2
    #var total_sample_variation = norm(Y - mean(y)) ^ 2 #SST
    variance_normalization_factor = ((X.rows - X.cols - 1).toFloat / (X.rows - X.cols).toFloat)
    s2 = sse / (X.rows - X.cols - 1).toFloat * variance_normalization_factor
    var_beta_hat = s2 * XpX
    estimate_std = var_beta_hat.diag .^ 0.5
    coefficients = beta_hat.T[0]
    include_intercept = X.wise_standard_deviation(axis=0).any_val(true_val=0.0)

  var total_model_variation = norm(Ypred.ravel - Ypred.ravel.mean) ^ 2 #SSE, ESS
  if include_intercept:
    total_model_variation = norm(Ypred.ravel - y.mean) ^ 2 #RSS
  else:
    total_model_variation = norm(y - residuals) ^ 2
  let f_score = (total_model_variation / (X.cols + (if include_intercept: -1 else: 0)).toFloat) / (sse / (X.rows - X.cols).toFloat)
  
  result.include_intercept = include_intercept
  result.residuals = residuals
  result.sum_squared_errors = sse
  result.variance_matrix_coefficients = var_beta_hat
  result.R2 = total_model_variation/(total_model_variation + sse)
  result.adjustedR2 = 1.0 - (X.rows + (if include_intercept: -1 else: 0)).toFloat / (X.rows - X.cols).toFloat * (1.0 - result.R2)
  result.beta_hat = beta_hat.ravel
  #
  result.loglikelihood = normal(mean=0.0, std=(sse * (result.degrees_of_freedom + 1e-10 - 2.0) / (X.rows - X.cols).toFloat).abs.sqrt).loglikelihood(residuals)
  #Fixed according https://stats.stackexchange.com/questions/277009/why-are-the-degrees-of-freedom-for-multiple-regression-n-k-1-for-linear-reg
  #let dof = (X.rows - X.cols - 1).toFloat #DOF without intercept
  result.degrees_of_freedom = (X.rows - X.cols).toFloat#Because if includes intercept, it is in the design matrix
  
  result.noise_variance = shifted_estimator(
    distribution=invchisquare(dof=result.degrees_of_freedom + 1e-10),
    location=0.0,#result.degrees_of_freedom,
    scale=sse * (result.degrees_of_freedom + 1e-10 - 2.0) / (X.rows - X.cols).toFloat
  )
  result.coefficients = (0..estimate_std.high).mapIt(
    shifted_estimator(
      distribution=studentt(dof=result.degrees_of_freedom),
      location=coefficients[it], scale=estimate_std[it]
    )
  )
  #Note s2 and noise_variance.estimate are the same, but s2 avoids to calculate it many times
  
  let ms_model_dof = (X.cols + (if include_intercept: -1 else: 0)).toFloat
  let ms_residual_dof = (X.rows - X.cols).toFloat
  result.model_significance = central_f(df1=ms_model_dof, df2=ms_residual_dof).htest_score( 
    score=(total_model_variation/ms_model_dof)/(sse/ms_residual_dof),
    test_type=oneTailed
  )
  if names.len > 0:
    result.feature_names = names
  else:
    result.feature_names = (1..X.cols).mapIt(fmt"x{it}")

proc estimate_as_string*[T](estimator: SimpleEstimator[T], title: bool=true, name: string="x"): string =
  if title:
    result &= " ".repeat(name.len + 1) & &"| Estimate| Std. dev.| Confidence int. 95%|   t-value|   p-value| \n"
  result &= &" {name}|{estimator.estimate:>9.5f}| {estimator.standard_deviation:>9.5f}| {estimator.confidence_interval[0]:>9.5f} {estimator.confidence_interval[1]:>9.5f}| {estimator.test_score:>9.5f}| {estimator.p_value:>9.5e}|"

proc coefficients_as_string*(model: OLSModel, title: bool=true): string =
  let min_len = max(model.feature_names.mapIt(it.len))
  if title:
    result &= " ".repeat(min_len + 1) & &"| Estimate| Std. dev.| Confidence int. 95%|   t-value|   p-value| \n"
  for k in 0..model.coefficients.high:
    result &= model.coefficients[k].estimate_as_string(title=false, name=model.feature_names[k].align(min_len))
    #result &= model.feature_names[k].align(min_len) & &"|{model.coefficients[k].estimate:>9.5f}| {model.coefficients[k].standard_deviation:>9.5f}| {model.coefficients[k].confidence_interval[0]:>9.5f} {model.coefficients[k].confidence_interval[1]:>9.5f}| {model.coefficients[k].test_score:>9.5f}| {model.coefficients[k].p_value:>9.5e}|"
    if k != model.coefficients.high: result &= "\n"

proc `$`*(model: OLSModel): string =
  result = "[Ordinary Least Squares Model]\n"
  result &= fmt"* Design matrix: {model.coefficients.len}x{model.residuals.len}" & (if model.include_intercept: " (include intercept)\n" else: "\n")
  result &= "* Coefficients:\n"
  result &= model.coefficients_as_string & "\n"
  result &= "* Noise:\n"
  result &= $model.noise_variance & "\n"
  result &= &"* Residual standard error: {(model.noise_variance.estimate ^ 0.5):.3f} on {model.degrees_of_freedom.toInt} degrees of freedom\n"
  result &= &"* Multiple R-squared: {model.R2:.5f},	Adjusted R-squared: {model.adjustedR2:.5f}\n"
  #result &= &"* F-statistic: {model.model_significance.test_score:.5f} on {model.model_significance.distribution.df1.toInt} and {model.model_significance.distribution.df2.toInt} DF, p-value: {model.model_significance.p_value:.5e} \n"
  #
  result &= &"* Analysis of variance:\n"
  let anova_model_dof = model.model_significance.distribution.df1
  let anova_residual_dof = model.model_significance.distribution.df2
  let anova_residual_mse = model.sum_squared_errors / anova_residual_dof
  let anova_model_mse = model.model_significance.test_score * anova_residual_mse
  result &= "           |  D.f.|    Sum Sq.|   Mean Sq.|   F value |   p-value\n"
  result &= &" Predictors|  {anova_model_dof:>4}| {anova_model_mse*anova_model_dof:>9.5e}| {anova_model_mse:>9.5e}| {model.model_significance.test_score:>9.5e}| {model.model_significance.p_value:.5e}\n"
  result &= &"  Residuals|  {anova_residual_dof:>4}| {anova_residual_mse*anova_residual_dof:>9.5e}| {anova_residual_mse:>9.5e}|"
    
when isMainModule:
  echo "Testing ols function: ols_model"
  if true:
    const datax1 = @[0.0, 1.0, 1.9, 2.9, 4.01, 4.9]
    const datay1 = @[5.0, 6.0, 7.0, 8.0, 9.0, 10.1]
    
    if true:
      let ols_model = ols_model(data_transform(datax1.as_column_vector, (x) => hstack(x, 1)), datay1)
      echo "\n   y ~ x + 1"
      echo ols_model
      assert (ols_model.model_significance.test_score - 3250.5).abs < 1e-1
      assert (ols_model.model_significance.p_value - 5.667e-7).abs < 1e-7
      assert (ols_model.R2 - 0.9988).abs < 1e-4
      assert (ols_model.adjustedR2 - 0.9985).abs < 1e-4
      assert (ols_model.degrees_of_freedom - 4.0).abs < 1e-4
      assert ((ols_model.residuals - @[0.001699, -0.025506, 0.050008, 0.022803, -0.117395, 0.068391]).abs < 1e-5).all_val
      assert ((ols_model.coefficients.mapIt(it.estimate) - @[1.02721, 4.99830]).abs < 1e-5).all_val
      assert ((ols_model.coefficients.mapIt(it.standard_deviation) - @[0.01802, 0.05360]).abs < 1e-5).all_val
      assert ((ols_model.coefficients.mapIt(it.test_score) - @[57.01, 93.25]).abs < 1e-2).all_val
      assert ((ols_model.coefficients.mapIt(it.p_value) - @[5.67e-7, 7.93e-8]).abs < 1e-9).all_val
      assert (((ols_model.noise_variance.estimate ^ 0.5) - 0.07439).abs < 1e-5)
    
    if true:
      let ols_model = ols_model(data_transform(datax1.as_column_vector, (x) => hstack(x)), datay1)
      echo "\n   y ~ x + 0"
      echo ols_model
      assert (ols_model.model_significance.test_score - 32.09).abs < 1e-1
      assert (ols_model.model_significance.p_value - 0.002385).abs < 1e-4
      assert (ols_model.R2 - 0.8652).abs < 1e-3
      assert (ols_model.adjustedR2 - 0.8382).abs < 1e-3
      assert (ols_model.degrees_of_freedom - 5.0).abs < 1e-4
      assert ((ols_model.residuals - @[5.0000, 3.5884, 2.4180, 1.0064, -0.6705, -1.7168]).abs < 1e-3).all_val
      assert ((ols_model.coefficients.mapIt(it.estimate) - @[2.4116]).abs < 1e-3).all_val
      assert ((ols_model.coefficients.mapIt(it.standard_deviation) - @[0.4257]).abs < 1e-3).all_val
      assert ((ols_model.coefficients.mapIt(it.test_score) - @[5.664]).abs < 1e-3).all_val
      assert ((ols_model.coefficients.mapIt(it.p_value) - @[0.00238]).abs < 1e-5).all_val
      assert (((ols_model.noise_variance.estimate ^ 0.5) - 3.103).abs < 1e-3)
    
    if true:
      let ols_model = ols_model(datax1.as_column_vector, datay1, @["x"])
      echo "\n   y ~ x + 0"
      echo ols_model
      assert (ols_model.model_significance.test_score - 32.09).abs < 1e-1
      assert (ols_model.model_significance.p_value - 0.002385).abs < 1e-4
      assert (ols_model.R2 - 0.8652).abs < 1e-3
      assert (ols_model.adjustedR2 - 0.8382).abs < 1e-3
      assert (ols_model.degrees_of_freedom - 5.0).abs < 1e-4
      assert ((ols_model.residuals - @[5.0000, 3.5884, 2.4180, 1.0064, -0.6705, -1.7168]).abs < 1e-3).all_val
      assert ((ols_model.coefficients.mapIt(it.estimate) - @[2.4116]).abs < 1e-3).all_val
      assert ((ols_model.coefficients.mapIt(it.standard_deviation) - @[0.4257]).abs < 1e-3).all_val
      assert ((ols_model.coefficients.mapIt(it.test_score) - @[5.664]).abs < 1e-3).all_val
      assert ((ols_model.coefficients.mapIt(it.p_value) - @[0.00238]).abs < 1e-5).all_val
      assert (((ols_model.noise_variance.estimate ^ 0.5) - 3.103).abs < 1e-3)
  
    const datax2 = @[1.0,  2.0,  4.0,  6.0,  9.0, 10.1, 11.2]
    const datay2 = @[-7.00, -6.40, -11.40, -19.40, -54, -72.0, -90]

    if true:
      let ols_model = ols_model(data_transform(datax2.as_column_vector, (x) => hstack(@[1.0], x)), datay2, @["x", "c"])
      echo "\n   y ~ x + 1"
      echo ols_model
      assert (ols_model.model_significance.test_score - 47.85).abs < 1e-1
      assert (ols_model.model_significance.p_value - 0.000968).abs < 1e-4
      assert (ols_model.R2 - 0.9054).abs < 1e-3
      assert (ols_model.adjustedR2 - 0.8865).abs < 1e-2
      assert (ols_model.degrees_of_freedom - 5.0).abs < 1e-4
      assert ((ols_model.residuals - @[-11.984,  -3.255,   8.003,  16.262,   6.049,  -3.009, -12.067]).abs < 1e-3).all_val
      assert ((ols_model.coefficients.mapIt(it.estimate) - @[13.113, -8.129]).abs < 1e-3).all_val
      assert ((ols_model.coefficients.mapIt(it.standard_deviation) - @[8.491, 1.175]).abs < 1e-3).all_val
      assert ((ols_model.coefficients.mapIt(it.test_score) - @[1.544, -6.918]).abs < 1e-3).all_val
      assert ((ols_model.coefficients.mapIt(it.p_value) - @[0.183147, 0.000968]).abs < 1e-5).all_val
      assert (((ols_model.noise_variance.estimate ^ 0.5) - 11.61).abs < 1e-3)
  
    if true:
      let ols_model = ols_model(data_transform(datax2.as_column_vector, (x) => hstack(x .^ 2.0, x, @[1.0])), datay2, @["x^2", "x", "c"])
      echo "\n   y ~ x + x^2 + 1"
      echo ols_model
      assert (ols_model.model_significance.test_score - 1732.2).abs < 1e-1
      assert (ols_model.model_significance.p_value - 1.33e-6).abs < 1e-6
      assert (ols_model.R2 - 0.9988).abs < 1e-3
      assert (ols_model.adjustedR2 - 0.9983).abs < 1e-3
      assert (ols_model.degrees_of_freedom - 4.0).abs < 1e-4
      assert ((ols_model.residuals - @[0.6199, -0.2095, -1.8471,  1.8102,  0.2488, -0.9609,  0.3386 ]).abs < 1e-3).all_val
      assert ((ols_model.coefficients.mapIt(it.estimate) - @[-1.03685, 4.53990, -11.12294]).abs < 1e-3).all_val
      assert ((ols_model.coefficients.mapIt(it.standard_deviation) - @[0.05759, 0.71852, 1.70616]).abs < 1e-3).all_val
      assert ((ols_model.coefficients.mapIt(it.test_score) - @[-18.003, 6.318, -6.519]).abs < 1e-3).all_val
      assert ((ols_model.coefficients.mapIt(it.p_value) - @[5.6e-05, 0.00321, 0.00286]).abs < 1e-5).all_val
      assert (((ols_model.noise_variance.estimate ^ 0.5) - 1.433).abs < 1e-3)
  
    if true:
      let ols_model = ols_model(data_transform(datax2.as_column_vector, (x) => hstack(x .^ 2.0, x)), datay2)
      echo "\n   y ~ x + x^2 + 0"
      echo ols_model
      assert (ols_model.model_significance.test_score - 437.1).abs < 1e-1
      assert (ols_model.model_significance.p_value - 2.439e-06).abs < 1e-4
      assert (ols_model.R2 - 0.9943).abs < 1e-3
      assert (ols_model.adjustedR2 - 0.992).abs < 1e-3
      assert (ols_model.degrees_of_freedom - 5.0).abs < 1e-4
      assert ((ols_model.residuals - @[-6.6820, -4.2828, -1.2409, 4.7258,  2.1847, -0.7203, -1.8331 ]).abs < 1e-3).all_val
      assert ((ols_model.coefficients.mapIt(it.estimate) - @[-0.7406, 0.4226]).abs < 1e-3).all_val
      assert ((ols_model.coefficients.mapIt(it.standard_deviation) - @[0.1079, 1.0449]).abs < 1e-3).all_val
      assert ((ols_model.coefficients.mapIt(it.test_score) - @[-6.864, 0.404]).abs < 1e-3).all_val
      assert ((ols_model.coefficients.mapIt(it.p_value) - @[0.001, 0.703]).abs < 1e-3).all_val
      assert (((ols_model.noise_variance.estimate ^ 0.5) - 4.371).abs < 1e-3)
      
    if true:
      let ols_model = ols_model(data_transform(datax2.as_column_vector, (x) => hstack(x .^ 3, x .^ 2.0, x)), datay2)
      echo "\n   y ~ x + x^2 + x^3 + 0"
      echo ols_model
      assert (ols_model.model_significance.test_score - 647.8).abs < 1e-1
      assert (ols_model.model_significance.p_value - 7.904e-06).abs < 1e-4
      assert (ols_model.R2 - 0.9979).abs < 1e-3
      assert (ols_model.adjustedR2 - 0.9964).abs < 1e-3
      assert (ols_model.degrees_of_freedom - 4.0).abs < 1e-4
      assert ((ols_model.residuals - @[-3.8759, -0.5162, 0.6864, 2.7604, -1.6983, -2.0912 , 1.9648 ]).abs < 1e-3).all_val
      assert ((ols_model.coefficients.mapIt(it.estimate) - @[-0.07401, 0.40423, -3.45429]).abs < 1e-3).all_val
      assert ((ols_model.coefficients.mapIt(it.standard_deviation) - @[0.02783, 0.43646, 1.61781]).abs < 1e-3).all_val
      assert ((ols_model.coefficients.mapIt(it.test_score) - @[-2.660, 0.926, -2.135]).abs < 1e-3).all_val
      assert ((ols_model.coefficients.mapIt(it.p_value) - @[0.0564, 0.4068, 0.0996]).abs < 1e-4).all_val
      assert (((ols_model.noise_variance.estimate ^ 0.5) - 2.937).abs < 1e-3)
  
  echo "[OK]"

###########################################################
# OLS
###########################################################
type OLSPrediction = ref object
  predictions*: seq[SimpleEstimator[StudentT]]
  degrees_of_freedom*: float
  prediction_means*: vector

proc predict*(model: OLSModel, X: matrix): OLSPrediction =
  new result
  let
    Yhat = (X * model.beta_hat.as_column_vector).ravel
    #Yvar = X * model.variance_matrix_coefficients * X.T
    #estimate_std = Yvar.diag .^ 0.5
    #https://stackoverflow.com/questions/23903021/compute-only-diagonal-of-product-of-matrices-in-matlab
    Yvar = ((X * model.variance_matrix_coefficients) .* X).sum(axis=1)
    estimate_std = Yvar .^ 0.5
  result.prediction_means = Yhat
  result.degrees_of_freedom = (estimate_std.len - 1).toFloat
  result.predictions = (0..estimate_std.high).mapIt(
    shifted_estimator(
      distribution=studentt(dof=result.degrees_of_freedom),
      location=Yhat[it], scale=estimate_std[it]
    )
  )

when isMainModule:
  echo "Testing ols function: predict"
  if true:
    const datax1 = @[0.0, 1.0, 1.9, 2.9, 4.01, 4.9]
    const dataxp1 = @[6.0, 7.0, 8.0, 9.0, 10.0]
    const datay1 = @[5.0, 6.0, 7.0, 8.0, 9.0, 10.1]
    
    if true:
      let ols_model = ols_model(datax1.add_intercept, datay1)
      echo "   y ~ x + 1"
      let ols_prediction = ols_model.predict(dataxp1.add_intercept)
      assert ((ols_prediction.prediction_means - @[11.16153, 12.18874, 13.21595, 14.24315, 15.27036]).abs < 1e-3).all_val
      assert ((ols_prediction.predictions.mapIt(it.estimate) - @[11.16153, 12.18874, 13.21595, 14.24315, 15.27036]).abs < 1e-3).all_val
      assert ((ols_prediction.predictions.mapIt(it.standard_deviation) - @[0.07077635, 0.08739312, 0.10447509, 0.12182674, 0.13934735]).abs < 1e-4).all_val
      
    const datax2 = @[1.0,  2.0,  4.0,  6.0,  9.0, 10.1, 11.2]
    const dataxp2 = @[12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0, 21.0, 22.0, 23.0, 24.0]
    const datay2 = @[-7.00, -6.40, -11.40, -19.40, -54, -72.0, -90]

    if true:
      let ols_model = ols_model(data_transform(datax2.as_column_vector, (x) => hstack(x .^ 3, x .^ 2.0, x)), datay2)
      echo "   y ~ x + x^2 + x^3 + 0"
      let ols_prediction = ols_model.predict(data_transform(dataxp2.as_column_vector, (x) => hstack(x .^ 3, x .^ 2.0, x)))
      assert ((ols_prediction.prediction_means - @[-111.1373, -139.1980, -172.2233, -210.6572, -254.9439, -305.5273, -362.8516, -427.3608, -499.4991, -579.7105, -668.4390, -766.1289, -873.2240]).abs < 1e-3).all_val
      assert ((ols_prediction.predictions.mapIt(it.estimate) - @[-111.1373, -139.1980, -172.2233, -210.6572, -254.9439, -305.5273, -362.8516, -427.3608, -499.4991, -579.7105, -668.4390, -766.1289, -873.2240]).abs < 1e-3).all_val
      assert ((ols_prediction.predictions.mapIt(it.standard_deviation) - @[4.609635, 8.313027, 13.406165, 20.006098, 28.263133, 38.337475, 50.392947, 64.594908, 81.109439, 100.102991, 121.742216, 146.193880, 173.624819]).abs < 1e-4).all_val
      
  echo "[OK]"
    
