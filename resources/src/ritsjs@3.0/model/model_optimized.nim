import math
import stats
import sequtils
import strformat


###########################################################
# OPERATIONS
#const M = 120000
const M = 1012
const EPSILON = 1e-10

type Vector = object
  length: int
  data*: ref array[M, float32]

type ArrayCheckType = (ref array[M, float32]) | seq[float]

proc check_value(data: var ArrayCheckType, length: int, idx: int, w: openArray[float], w_idx: int) {.inline.} =
  let v = w[w_idx]
  #echo v.classify, "  ", v, "  ", (v * 10 / 10 != v)
  if v.classify == fcNaN or v.classify == fcInf or v.classify == fcNegInf or (v * 10 / 10 != v): #This works in JS
    if idx == 0:
      data[idx] = w[w_idx + 1]
    elif idx == length - 1:
      data[idx] = w[w_idx - 1]
    else:
      data[idx] = 0.5 * (w[w_idx - 1] + w[w_idx + 1])
  else:
    data[idx] = v

proc fix_vector(v: seq[float]): seq[float] {.exportc.} =
  let length = v.len
  echo length
  newSeq(result, length)
  var i: int = 0
  for k in v.low..v.high:
    check_value(result, length, i, v, k)
    inc(i)

proc vector(v: openArray[float]): Vector =
  new result.data
  result.length = v.len
  var i: int = 0
  for k in v.low..v.high:
    check_value(result.data, result.length, i, v, k)
    #result.data[i] = v[k];
    inc(i)
  result.length = i
  result

proc vector(v: openArray[int]): Vector = vector(v.mapIt(it.float))

proc vector_reduced(v: openArray[float], sampling: int): Vector =
  new result.data
  result.length = v.len
  var i: int = 0
  for k in v.low..v.high:
    if k mod sampling != 0:
      continue
    check_value(result.data, result.length, i, v, k)
    #result.data[i] = v[k];
    inc(i)
  result.length = i
  result

# SIZE 
proc len(v: Vector): int {.inline.} = v.length
proc high(v: Vector): int {.inline.} = v.length - 1
proc low(v: Vector): int {.inline.} = 0

# TO STRING
proc `$`(v: Vector): string {.inline.} =
  result = "["
  for k in v.low..v.high:
    result &= fmt"{v.data[k]:.04f}"
    if k != v.high:
      result &= ", "
  result &= "]"

proc vector_to_seq(v: Vector): seq[float] {.inline.} =
  let length = v.len
  newSeq(result, length)
  var i: int = 0
  for k in v.low..v.high:
    result[i] = v.data[k]
    inc(i)

# ACCESSORS
proc `[]`*(v: Vector, i: int): float {.inline.} =
  assert(i >= 0, "Negative index")
  assert(i < v.len, "Index out of bounds")
  v.data[i]

proc `[]=`*(v: var Vector, i: int, val:float) {.inline.} =
  assert(i >= 0, "Negative index")
  assert(i < M, "Index out of bounds")
  if i >= v.length:
    v.length = i + 1
  v.data[i] = val

# ACCESSORS
proc `[]`*(v: Vector, idx: HSlice[int, int]): Vector {.inline.} =
  new result.data
  var k = 0
  for i in idx:
    assert(i >= 0, "Negative index")
    assert(i < v.len, "Index out of bounds")
    result.data[k] = v.data[i]
    inc(k)
  result.length = k

proc `[]=`*(v: var Vector, idx: HSlice[int, int], val:float) {.inline.} =
  for i in idx:
    assert(i >= 0, "Negative index")
    assert(i < M, "Index out of bounds")
    if i >= v.length:
      v.length = i + 1
    v.data[i] = val

proc `[]=`*(v: var Vector, idx: HSlice[int, int], w: var Vector) {.inline.} =
  assert(idx.b - idx.a + 1 == w.len, "Assign vector in a wrong size")
  var k = 0
  for i in idx:
    assert(i >= 0, "Negative index")
    assert(i < M, "Index out of bounds")
    if i >= v.length:
      v.length = i + 1
    v.data[i] = w[k]
    inc(k)

# ELEMENT-WISE SUM VECTOR
proc `+`*(v: Vector, w: Vector): Vector {.inline.}=
  new result.data
  for i in v.low..v.high:
    result[i] = v[i] + w[i]

proc `+`*(v: Vector, k: float): Vector {.inline.} =
  new result.data
  for i in v.low..v.high:
    result[i] = v[i] + k

proc `+`*(k: float, v: Vector): Vector {.inline.} = v + k

proc `+`*(v: Vector, k: int): Vector {.inline.} = v + k.float

proc `+`*(k: int, v: Vector): Vector {.inline.} = v + k.float

# ELEMENT-WISE DIFFERENCE VECTOR
proc `-`*(v: Vector, w: Vector): Vector {.inline.}=
  new result.data
  for i in v.low..v.high:
    result[i] = v[i] - w[i]

proc `-`*(v: Vector, k: float): Vector {.inline.} =
  new result.data
  for i in v.low..v.high:
    result[i] = v[i] - k

proc `-`*(k: float, v: Vector): Vector {.inline.} =
  new result.data
  for i in v.low..v.high:
    result[i] = k - v[i]

proc `-`*(v: Vector, k: int): Vector {.inline.} = v - k.float

proc `-`*(k: int, v: Vector): Vector {.inline.} = k.float - v

# ELEMENT-WISE PRODUCT VECTOR
proc `.*`*(v: Vector, w: Vector): Vector {.inline.}=
  new result.data
  for i in v.low..v.high:
    result[i] = v[i] * w[i]

proc `.*`*(v: Vector, k: float): Vector {.inline.} =
  new result.data
  for i in v.low..v.high:
    result[i] = v[i] * k

proc `.*`*(k: float, v: Vector): Vector {.inline.} = v .* k

proc `.*`*(v: Vector, k: int): Vector {.inline.} = v .* k.float

proc `.*`*(k: int, v: Vector): Vector {.inline.} = k.float .* v

proc `*`*(v: Vector, k: float): Vector {.inline.} = v .* k
proc `*`*(k: float, v: Vector): Vector {.inline.} = k .* v
proc `*`*(v: Vector, k: int): Vector {.inline.} = v .* k
proc `*`*(k: int, v: Vector): Vector {.inline.} = k .* v

# ELEMENT-WISE DIVISION VECTOR
proc `./`*(v: Vector, w: Vector): Vector {.inline.}=
  new result.data
  for i in v.low..v.high:
    result[i] = v[i] / w[i]

proc `./`*(v: Vector, k: float): Vector {.inline.} =
  new result.data
  for i in v.low..v.high:
    result[i] = v[i] / k

proc `./`*(k: float, v: Vector): Vector {.inline.} =
  new result.data
  for i in v.low..v.high:
    result[i] = k / v[i]

proc `./`*(v: Vector, k: int): Vector {.inline.} = v ./ k.float

proc `./`*(k: int, v: Vector): Vector {.inline.} = k.float ./ v

proc `/`*(v: Vector, k: float): Vector {.inline.} = v ./ k
proc `/`*(k: float, v: Vector): Vector {.inline.} = k ./ v
proc `/`*(v: Vector, k: int): Vector {.inline.} = v ./ k
proc `/`*(k: int, v: Vector): Vector {.inline.} = k ./ v

# ELEMENT-WISE POWER VECTOR
proc `^` (v: float32, w: float32): float32 {.inline.} = pow(v, w)

proc `.^`*(v: Vector, w: Vector): Vector {.inline.}=
  new result.data
  for i in v.low..v.high:
    result[i] = v[i] ^ w[i]

proc `.^`*(v: Vector, k: float): Vector {.inline.} =
  new result.data
  for i in v.low..v.high:
    result[i] = v[i] ^ k

proc `.^`*(k: float, v: Vector): Vector {.inline.} =
  new result.data
  for i in v.low..v.high:
    result[i] = k ^ v[i]

proc `.^`*(v: Vector, k: int): Vector {.inline.} = v .^ k.float

proc `.^`*(k: int, v: Vector): Vector {.inline.} = k.float .^ v

proc `^`*(v: Vector, k: float): Vector {.inline.} = v .^ k
proc `^`*(k: float, v: Vector): Vector {.inline.} = k .^ v
proc `^`*(v: Vector, k: int): Vector {.inline.} = v .^ k
proc `^`*(k: int, v: Vector): Vector {.inline.} =  k .^ v

# DOT PRODUCT VECTOR
proc dot(v: Vector, w: Vector): float {.inline.} =
  for i in v.low..v.high:
    result += v[i] * w[i]

proc `*`*(v: Vector, w: Vector): float {.inline.} = dot(v, w)

###########################################################
# BASIC STATISTICS
proc sum*(v: Vector): float  {.inline.} = v.data[].sum
proc mean*(v: Vector): float  {.inline.} = v.sum / v.len.float
proc avg*(v: Vector): float  {.inline.} = mean(v)
proc variance*(v: Vector): float  {.inline.}= (v .* v).sum / v.len.float - v.mean ^ 2
proc sample_variance*(v: Vector): float  {.inline.} = ((v .* v).sum / v.len.float - v.mean ^ 2) * (v.len.float) / (v.len.float - 1)
proc covariance*(v: Vector, w: Vector): float  {.inline.} = (dot(v, w) / v.len.float - v.mean * w.mean)
proc sample_covariance*(v: Vector, w: Vector): float  {.inline.} = (dot(v, w) / v.len.float - v.mean * w.mean) * (v.len.float) / (v.len.float - 1)

proc student_t_cdf(t: float, df: int): float {.inline.} =
  let y = (df.float + 2.0/19.0) / (df.float + 22.0/71.0) * sqrt((df.float + 2.0/19.0) * ln(1.0 + t * t / (df.float + 13.0/50.0)))
  return 1.0/(1.0 + exp(-1.6 * y - 0.07 * y * y * y))

proc student_t_ppf_95p(df: int): float {.inline.} =
  const table_ppf = [1.0e50, 12.71, 4.303, 3.182, 2.776, 2.571, 2.447, 2.365,
  2.306, 2.262, 2.228, 2.201, 2.179, 2.160, 2.145, 2.131,
  2.120, 2.110, 2.101, 2.093, 2.086, 2.080, 2.074, 2.069,
  2.064, 2.060, 2.056, 2.052, 2.048, 2.045, 2.042, 2.021,
  2.009, 2.000, 1.990, 1.984, 1.980, 1.960]
  if df > table_ppf.high:
    return table_ppf[table_ppf.high]
  elif df <= 1:
    return table_ppf[1]
  return table_ppf[df]

proc student_t_ppf_95p(df: float): float {.inline.} =
  student_t_ppf_95p(df.int)

###########################################################
type AutocorrelationFunction = object
  autocorrelation: seq[float]
  confidence_interval: float

proc acf_confidence_interval(x: Vector): float = 
  1.96 * sqrt(1 / x.len)

proc acf(x: Vector, max_lag: int = 20): Vector = 
  new result.data
  let
    n = x.high
    u = x.mean
    s2 = x.variance
  result[0] = 1
  for h in 1..min(max_lag, n):
    let
      a = x[0..(n - h)]
      b = x[h..n]
    result[h] = (((a - u) .* (b - u)).mean + EPSILON) / (n.float * s2 + EPSILON )

proc acf_information(x: Vector, max_lag: int = 20): AutocorrelationFunction =
  result.autocorrelation = acf(x, max_lag).vector_to_seq
  result.confidence_interval = acf_confidence_interval(x)

when isMainModule:
  let
    x = vector([3,2,1,9,4,2,5,3,3,6,8,9,10])
    y = acf(x, 10)
  echo "AutocorrelationFunction"
  echo fmt"  X: {x}"
  echo fmt"  ACF: {y}"
  echo fmt"  ACF-info: {acf_information(x, 10)}"
  echo ""

###########################################################

type MeanEstimates = object
  beta_0: float
  beta_1: float
  delta_0: float
  delta_1: float
  variance_pre: float
  variance_post: float

type LinearRegressionParameters = object
  intercept: float
  intercept_variance: float 
  intercept_confidence_interval: array[2, float]
  intercept_width_confidence_interval: float 
  intercept_p_value: float
  #
  slope: float
  slope_variance: float
  slope_confidence_interval: array[2, float]
  slope_width_confidence_interval: float
  slope_p_value: float
  #
  R2: float
  #
  residual: seq[float]
  residual_sum_squares: float
  residual_variance: float
  #
  autocorrelation_function: AutocorrelationFunction
  #
  n: int

type MeanStructureDifferenceParameters = object
  intercept: float
  intercept_variance: float 
  intercept_confidence_interval: array[2, float]
  intercept_width_confidence_interval: float 
  intercept_p_value: float 
  #
  slope: float
  slope_variance: float
  slope_confidence_interval: array[2, float]
  slope_width_confidence_interval: float
  slope_p_value: float
  #
  autoregressive_slope: float
  autoregressive_slope_variance: float
  autoregressive_slope_confidence_interval: array[2, float]
  autoregressive_slope_width_confidence_interval: float
  autoregressive_slope_p_value: float
  #

type AR1Ma1ModelParameters = object
  mean_structure: LinearRegressionParameters
  autoregressive_structure: LinearRegressionParameters

type LikelihoodResult = object
  loglikelihood: float
  before_change: AR1Ma1ModelParameters
  after_change: AR1Ma1ModelParameters

type LikelihoodModel = object
  change_points: seq[int]
  loglikelihood: seq[float]
  best_loglikelihood: float
  best_likelihood: float
  best_time: float
  best_index: int

type RobustInterruptedModel = object
  before_change: AR1Ma1ModelParameters
  after_change: AR1Ma1ModelParameters
  likelihood: LikelihoodModel
  parameter_differences: MeanStructureDifferenceParameters

proc p_value_of_t_distr(mean: float, variance: float, n: int): float {.inline.} =
  2.0 - 2.0 * student_t_cdf(mean / sqrt(variance), n.int - 2)

proc model_differences(before_change: AR1Ma1ModelParameters, after_change: AR1Ma1ModelParameters): MeanStructureDifferenceParameters =
  result.intercept = after_change.mean_structure.intercept - before_change.mean_structure.intercept
  result.slope = after_change.mean_structure.slope - before_change.mean_structure.slope
  result.autoregressive_slope = after_change.autoregressive_structure.slope - before_change.autoregressive_structure.slope
  
  let
    intercept_dof = after_change.mean_structure.n + before_change.mean_structure.n - 4
    slope_dof = after_change.mean_structure.n + before_change.mean_structure.n - 4
    autoregressive_slope_dof = after_change.autoregressive_structure.n + before_change.autoregressive_structure.n - 4
  
  result.intercept_variance = after_change.mean_structure.intercept_variance + before_change.mean_structure.intercept_variance
  result.slope_variance = after_change.mean_structure.slope_variance + before_change.mean_structure.slope_variance
  result.autoregressive_slope_variance = after_change.autoregressive_structure.slope_variance + before_change.autoregressive_structure.slope_variance
  
  result.intercept_p_value = p_value_of_t_distr(result.intercept, result.intercept_variance, intercept_dof + 2);
  result.slope_p_value = p_value_of_t_distr(result.slope, result.slope_variance, slope_dof + 2);
  result.autoregressive_slope_p_value = p_value_of_t_distr(result.autoregressive_slope, result.autoregressive_slope_variance, autoregressive_slope_dof + 2);
  
  result.intercept_width_confidence_interval = 2.0 * sqrt(result.intercept_variance) * student_t_ppf_95p(intercept_dof)
  result.slope_width_confidence_interval = 2.0 * sqrt(result.slope_variance) * student_t_ppf_95p(slope_dof)
  result.autoregressive_slope_width_confidence_interval = 2.0 * sqrt(result.autoregressive_slope_variance) * student_t_ppf_95p(autoregressive_slope_dof)

  result.intercept_confidence_interval = [result.intercept - 0.5 * result.intercept_width_confidence_interval, result.intercept + 0.5 * result.intercept_width_confidence_interval]
  result.slope_confidence_interval =  [result.slope - 0.5 * result.slope_width_confidence_interval, result.slope + 0.5 * result.slope_width_confidence_interval]
  result.autoregressive_slope_confidence_interval =  [result.autoregressive_slope - 0.5 * result.autoregressive_slope_width_confidence_interval, result.autoregressive_slope + 0.5 * result.autoregressive_slope_width_confidence_interval]
  

proc `$`(params: LinearRegressionParameters): string =
  fmt"""
    slope: {params.slope:.3f} (var: {params.slope_variance:.3f}, CI width: {params.slope_width_confidence_interval:.3f})
    intercept: {params.intercept:.3f} (var: {params.intercept_variance:.3f}, CI width: {params.intercept_width_confidence_interval:.3f})
    R2: {params.R2:.3f}
    SSE: {params.residual_sum_squares:.3f}"""

# y = ax + b
proc simple_linear_regression(X: Vector, Y: Vector): LinearRegressionParameters {.exportc.} =
  let
    n = X.len.float
    Sx = X.sum
    Sxx = (X .* X).sum
    Sy = Y.sum
    Sxy = (X .* Y).sum
    Syy = (Y .* Y).sum
    t_value = student_t_ppf_95p(n - 2)

  result.n = n.int
  result.slope = (n * Sxy - Sx * Sy) / (n * Sxx - Sx ^ 2 + EPSILON)
  result.intercept = (Sy / n - result.slope * Sx / n)
  #
  result.R2 = (n * Sxy - Sx * Sy) ^ 2 / (n * Sxx - Sx ^ 2) / (n * Syy - Sy ^ 2 + EPSILON)
  #
  let residual = Y - (result.intercept + result.slope * X)
  result.residual = residual.vector_to_seq
  result.residual_sum_squares = (residual .^ 2.0).sum
  result.residual_variance = (n * Syy - Sy ^ 2 - result.slope ^ 2 * (n * Sxx - Sx ^ 2))/(n * n - 2 * n + EPSILON)
  #result.residual_variance = result.residual_sum_squares / (n - 2 + EPSILON)
  #
  result.autocorrelation_function = acf_information(residual, 100)
  #
  result.slope_variance = n * result.residual_variance / (n * Sxx - Sx ^ 2 + EPSILON)
  result.intercept_variance = result.slope_variance * Sxx / n
  #
  result.slope_p_value = p_value_of_t_distr(result.slope, result.slope_variance, n.int);
  result.intercept_p_value = p_value_of_t_distr(result.intercept, result.intercept_variance, n.int);
  #
  result.slope_width_confidence_interval = 2.0 * sqrt(result.slope_variance) * t_value
  result.slope_confidence_interval[0] = result.slope - 0.5 * result.slope_width_confidence_interval
  result.slope_confidence_interval[1] = result.slope + 0.5 * result.slope_width_confidence_interval
  result.intercept_width_confidence_interval = 2.0 * sqrt(result.intercept_variance) * t_value
  result.intercept_confidence_interval[0] = result.intercept - 0.5 * result.intercept_width_confidence_interval
  result.intercept_confidence_interval[1] = result.intercept + 0.5 * result.intercept_width_confidence_interval

# y = ax + b
proc simple_linear_regression_wo_intercept(X: Vector, Y: Vector): LinearRegressionParameters {.exportc.} =
  let
    n = X.len.float
    Sx = X.sum
    Sxx = (X .* X).sum
    Sy = Y.sum
    Sxy = (X .* Y).sum
    Sxxyy = ((X .* Y) .* (X .* Y)).sum
    Syy = (Y .* Y).sum
  result.n = n.int
  result.slope = (n * Sxy - Sx * Sy + Sx * Sy) / (n * Sxx - Sx ^ 2 + Sx * Sx + EPSILON)
  result.intercept = 0
  #
  result.R2 = (Sxy ^ 2 + EPSILON) / (Sxx * Syy + EPSILON) # Sample R2, by the model it is assumed that E[Y] = 0
  #
  let residual = Y - (result.intercept + result.slope * X)
  result.residual = residual.vector_to_seq
  result.residual_sum_squares = (residual .^ 2.0).sum
  result.residual_variance = result.residual_sum_squares / (n - 1 + EPSILON)
  #
  result.autocorrelation_function = acf_information(residual, 100)
  #
  result.slope_variance = result.residual_variance / (Sxx + EPSILON)
  result.intercept_variance = 0
  #
  result.slope_p_value = p_value_of_t_distr(result.slope, result.slope_variance, n.int);
  result.intercept_p_value = 0#p_value_of_t_distr(result.intercept, result.intercept_variance, n.int);
  #
  result.slope_width_confidence_interval = 2.0 * sqrt(result.slope_variance) * student_t_ppf_95p(n - 1)
  result.slope_confidence_interval[0] = result.slope - 0.5 * result.slope_width_confidence_interval
  result.slope_confidence_interval[1] = result.slope + 0.5 * result.slope_width_confidence_interval
  result.intercept_width_confidence_interval = 0
  result.intercept_confidence_interval[0] = 0
  result.intercept_confidence_interval[1] = 0

#######################################################################################

proc estimated_var_matrix(autocorrelation: float, variance: float, Y: Vector): array[2, array[2, float]] =
  var
    d1 = 1.0
    d2 = 1.0
    tau = - variance
    abs_autocorrelation = abs(autocorrelation)
    ln_autocorrelation = ln(autocorrelation)
    T = Y.high
    X = Y[0..(T-2)]
  var
    sumX = X.sum()
    sumX2 = (X .* X).sum()
  if abs_autocorrelation > 2.0:
    echo "Range not allowed!"
  if abs_autocorrelation > 0.1:
    d1 = 1 - 1/(10 * ln_autocorrelation) - 2/(1 + exp(-10 * ln_autocorrelation))
    d2 = 1 - 1/(10 * ln_autocorrelation) - 1/(1 + exp(-10 * ln_autocorrelation))
    tau = 1/(5 * ln_autocorrelation)
  result[0][0] = d1 * 2.0 + tau * (2.0 * T.float - 6.0) + (T.float - 4.0) * d2
  result[1][0] = d1 * (X[0] + X[T-2]) + d2 * (sumX - X[0] - X[T-2]) + tau * (2.0 * sumX - X[T-2] - X[0])
  result[0][1] = result[1][0]
  result[1][1] = d1 * (X[0] ^ X[0] + X[T-2] * X[T-2]) + d2 * (sumX2 - X[0] * X[0] - X[T-2] * X[T-2]) + tau * (2.0 * sumX2 - X[0] * X[0] - X[T-2] * X[T-2])

#######################################################################################

# OLS Method
# y[t] = a x[t] + b + r[t]
# r[t] = p r[t-1] + e
#   =>  r[t:1] = p r[t-1:0] + e[t: 1]
proc arma(X: Vector, Y: Vector): AR1Ma1ModelParameters {.exportc.} =
  result.mean_structure = simple_linear_regression(X, Y)
  let residuals = Y - (result.mean_structure.intercept + result.mean_structure.slope * X)
  result.autoregressive_structure = simple_linear_regression_wo_intercept(residuals[0..(residuals.high - 1)], residuals[1..residuals.high])

proc normal_loglikelihood(X: Vector, Y: Vector, slope: float, intercept: float, sigma2: float): float =
  let
    n = X.len.float
  result = -0.5 * n * ln(2 * PI * sigma2 + EPSILON)
  result -= ((Y - (intercept + slope * X)) .^ 2).sum / (2 * sigma2 + EPSILON)
  
proc analysis_at_point_t(X: Vector, Y: Vector, change_point: int): LikelihoodResult =
  let
    N = X.len
    #
    X_before = X[0..(change_point-1)]
    Y_before = Y[0..(change_point-1)]
    #
    X_after = X[change_point..X.high]
    Y_after = Y[change_point..Y.high]
  result.before_change = arma(X_before, Y_before)
  result.loglikelihood += normal_loglikelihood(X_before, Y_before, result.before_change.mean_structure.slope, result.before_change.mean_structure.intercept, result.before_change.autoregressive_structure.residual_variance)
  result.after_change = arma(X_after, Y_after)
  result.loglikelihood += normal_loglikelihood(X_after, Y_after, result.after_change.mean_structure.slope, result.after_change.mean_structure.intercept, result.after_change.autoregressive_structure.residual_variance)


proc robust_interrupted_time_series_model(X: Vector, Y: Vector, change_point: int, candidates_before: int, candidates_after: int): RobustInterruptedModel =
  var
    ll_result: LikelihoodResult
  result.likelihood.best_loglikelihood = -1e100
  for t in (change_point - candidates_before) .. (change_point + candidates_after):
    if t < 3 or t > X.high:
      continue
    ll_result = analysis_at_point_t(X, Y, t)
    result.likelihood.change_points.add t
    result.likelihood.loglikelihood.add ll_result.loglikelihood
    if ll_result.loglikelihood > result.likelihood.best_loglikelihood:
      result.likelihood.best_loglikelihood = ll_result.loglikelihood
      result.likelihood.best_likelihood = exp(ll_result.loglikelihood)
      result.likelihood.best_index = t
      result.likelihood.best_time = X[t]
      result.before_change = ll_result.before_change
      result.after_change = ll_result.after_change
      result.parameter_differences = model_differences(result.before_change, result.after_change)

proc robust_interrupted_time_series(X: openArray[float], Y: openArray[float], change_point: int, candidates_before: int, candidates_after: int): RobustInterruptedModel  {.exportc: "robust_interrupted_time_series".} =
  robust_interrupted_time_series_model(vector(X), vector(Y), change_point, candidates_before, candidates_after)
  
proc robust_interrupted_time_series_approximated(sampling: int, X: openArray[float], Y: openArray[float], change_point: int, candidates_before: int, candidates_after: int): RobustInterruptedModel  {.exportc: "robust_interrupted_time_series_approximated".} =
  robust_interrupted_time_series_model(vector_reduced(X, sampling), vector_reduced(Y, sampling), (change_point.float / sampling.float).int, (candidates_before.float / sampling.float).int, (candidates_after.float / sampling.float).int)
  


import random
import times
import sequtils

proc test_model_change_point(change_point: int=7, candidates: int=5, verbose: bool=true) =
  echo "Robust Interrupted Time Series Model v3.0\n"
  echo "Boot sample test"
  var
    x = @[1.0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20].mapIt(it.float)
    y = @[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200].mapIt(it.float)
  for t in 0..x.high:
    y[t] = y[t].float + 0 * rand(1000).float / 1000.0
  if verbose:
    echo fmt"X: {x}"
    echo fmt"Y: {y}"
  var
    candidates_after = candidates
    candidates_before = candidates
    model = robust_interrupted_time_series(x, y, change_point, candidates_before, candidates_after)
  if verbose:
    echo model

test_model_change_point()

when isMainModule and not defined(in_production):
  import random
  import times

  proc test_linear_regression() =
    let
      x0 = vector([1.47, 1.50, 1.52, 1.55, 1.57, 1.60, 1.63, 1.65, 1.68, 1.70, 1.73, 1.75, 1.78, 1.80, 1.83])
      y0 = vector([52.21, 53.12, 54.48, 55.84, 57.20, 58.57, 59.93, 61.29, 63.11, 64.47, 66.28, 68.10, 69.92, 72.19, 74.46])
      slr0 = simple_linear_regression(x0, y0)
    echo fmt"X: {x0}"
    echo fmt"Y: {y0}"
    echo $slr0
    assert(abs(slr0.slope - 61.272) < 0.01)
    assert(abs(slr0.slope_variance - 1.776 ^ 2) < 0.01)
    assert(abs(slr0.slope_width_confidence_interval - 7.67328) < 0.01)
    assert(abs(slr0.intercept - -39.062) < 0.01)
    assert(abs(slr0.intercept_variance - 2.938 ^ 2) < 0.01)
    assert(abs(slr0.intercept_width_confidence_interval - 12.6943) < 0.01)
    assert(abs(slr0.R2 - 0.9892) < 0.01)

  proc test_linear_regression_wo_intercept() =
    let
      x0 = vector([1.47, 1.50, 1.52, 1.55, 1.57, 1.60, 1.63, 1.65, 1.68, 1.70, 1.73, 1.75, 1.78, 1.80, 1.83])
      y0 = vector([52.21, 53.12, 54.48, 55.84, 57.20, 58.57, 59.93, 61.29, 63.11, 64.47, 66.28, 68.10, 69.92, 72.19, 74.46])
      slr0 = simple_linear_regression_wo_intercept(x0, y0)
    echo fmt"X: {x0}"
    echo fmt"Y: {y0}"
    echo $slr0
    assert(abs(slr0.slope - 37.7131) < 0.01)
    assert(abs(slr0.slope_variance - 0.4362 ^ 2) < 0.01)
    echo abs(slr0.slope_width_confidence_interval - 1.871)
    assert(abs(slr0.slope_width_confidence_interval - 1.871) < 0.01)
    assert(abs(slr0.intercept - 0) < 0.01)
    assert(abs(slr0.intercept_variance - 0) < 0.01)
    assert(abs(slr0.intercept_width_confidence_interval - 0) < 0.01)
    assert(abs(slr0.R2 - 0.9981) < 0.01)

  proc test_model[T: static[int]](change_point: int, candidates: int, verbose: bool=true) =
    var
      x: array[T, float]
      y: array[T, float]
      x0: Vector
      y0: Vector
    for t in 0..(T-1):
      x[t] = t.float
      y[t] = t.float * 10.0 + rand(1000).float / 1000.0    
    x0 = vector(x)
    y0 = vector(y)
    if verbose:
      echo fmt"X: {x0[0..10]}"
      echo fmt"Y: {y0[0..10]}"
    var
      candidates_after = candidates
      candidates_before = candidates
      model = robust_interrupted_time_series_model(x0, y0, change_point, candidates_before, candidates_after)
    if verbose:
      echo model

  proc test_model_sampled[T: static[int]](change_point: int, candidates: int, sampling:int=1, verbose: bool=true) =
    var
      x: array[T, float]
      y: array[T, float]
      x0: Vector
      y0: Vector
    for t in 0..(T-1):
      x[t] = t.float
      y[t] = t.float * 10.0 + rand(1000).float / 1000.0    
    x0 = vector(x)
    y0 = vector(y)
    if verbose:
      echo fmt"X: {x0[0..10]}"
      echo fmt"Y: {y0[0..10]}"
    var
      candidates_after = candidates
      candidates_before = candidates
      model = robust_interrupted_time_series_approximated(sampling, x, y, change_point, candidates_before, candidates_after)
    if verbose:
      echo model

  proc test_model_50_50(verbose: bool=false) =
    test_model[50](25, 12, verbose)
  
  proc test_model_100(verbose: bool=false) =
    test_model[100](50, 5, verbose)

  proc test_model_500(verbose: bool=false) =
    test_model[500](250, 5, verbose)
  
  proc test_model_approx_500_500(verbose: bool=false) =
    test_model_sampled[500](250, 249, 10, verbose)

  proc test_model_200(verbose: bool=false) =
    test_model[200](100, 99, verbose)

  proc test_model_500_100(verbose: bool=false) =
    test_model[500](250, 100, verbose)

  proc test_model_5000(verbose: bool=false) =
    test_model[5000](250, 5, verbose)

  proc test_model_50000_5(verbose: bool=false) =
    test_model[50000](1000, 5, verbose)

  proc test_model_50000_50(verbose: bool=false) =
    test_model[50000](1000, 50, verbose)

  proc test_model_50000_500(verbose: bool=false) =
    test_model[50000](1000, 500, verbose)


  test_linear_regression()
  test_linear_regression_wo_intercept()
  test_model_100(verbose=true)
  
  echo "# 50 points (50 cands.) (x10)"
  var time = getTime()
  for k in 1..10:
    test_model_50_50()
  echo "## Time taken: ", getTime() - time

  echo "# 100 points (x10)"
  time = getTime()
  for k in 1..10:
    test_model_100()
  echo "## Time taken: ", getTime() - time
  
  echo "# 200 points (200 cands.) (x10)"
  time = getTime()
  for k in 1..10:
    test_model_200()
  echo "## Time taken: ", getTime() - time
  
  echo "# 500 points (x10)"
  time = getTime()
  for k in 1..10:
    test_model_500()
  echo "## Time taken: ", getTime() - time
  
  echo "# 500 points (500 cands.) - approx (x10)"
  time = getTime()
  for k in 1..10:
    test_model_approx_500_500()
  echo "## Time taken: ", getTime() - time
  
  echo "# 500 points (100 cands.) (x10)"
  time = getTime()
  for k in 1..10:
    test_model_500_100()
  echo "## Time taken: ", getTime() - time
  
  echo "# 5000 points (x10)"
  time = getTime()
  for k in 1..10:
    test_model_5000()
  echo "## Time taken: ", getTime() - time
  
  echo "# 50000 (5 cands.) points (x10)"
  time = getTime()
  for k in 1..10:
    test_model_50000_5()
  echo "## Time taken: ", getTime() - time
  
  echo "# 50000 (50 cands.) points (x10)"
  time = getTime()
  for k in 1..10:
    test_model_50000_50()
  echo "## Time taken: ", getTime() - time
  
  echo "# 50000 (500 cands.) points (x10)"
  time = getTime()
  for k in 1..10:
    test_model_50000_500()
  echo "## Time taken: ", getTime() - time

  ##### JS
  # 100 points (x10)
  ## Time taken: 792 milliseconds
  # 200 points (200 cands.) (x10)
  ## Time taken: 8 seconds and 272 milliseconds
  # 500 points (x10)
  ## Time taken: 537 milliseconds
  # 500 points (100 cands.) (x10)
  ## Time taken: 7 seconds and 78 milliseconds
  # 5000 points (x10)
  ## Time taken: 366 milliseconds
  # 50000 (5 cands.) points (x10)
  ## Time taken: 386 milliseconds
  # 50000 (50 cands.) points (x10)
  ## Time taken: 3 seconds and 940 milliseconds
  # 50000 (500 cands.) points (x10)
  ## Time taken: 35 seconds and 328 milliseconds

  ##### C
  # 100 points (x10)
  ## Time taken: 234 milliseconds, 372 microseconds, and 900 nanoseconds
  # 200 points (200 cands.) (x10)
  ## Time taken: 4 seconds, 613 milliseconds, 653 microseconds, and 900 nanoseconds
  # 500 points (x10)
  ## Time taken: 243 milliseconds, 349 microseconds, and 100 nanoseconds
  # 500 points (100 cands.) (x10)
  ## Time taken: 4 seconds, 460 milliseconds, 95 microseconds, and 500 nanoseconds
  # 5000 points (x10)
  ## Time taken: 258 milliseconds, 361 microseconds, and 900 nanoseconds
  # 50000 (5 cands.) points (x10)
  ## Time taken: 411 milliseconds, 895 microseconds, and 600 nanoseconds
  # 50000 (50 cands.) points (x10)
  ## Time taken: 3 seconds, 821 milliseconds, 772 microseconds, and 500 nanoseconds
  # 50000 (500 cands.) points (x10)
  ## Time taken: 37 seconds, 3 milliseconds, and 978 microseconds

