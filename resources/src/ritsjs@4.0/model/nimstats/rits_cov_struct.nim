import sequtils
import strutils
import strformat
import sugar
import algorithm
import math
import core
import distributions
import ols
import ar
import ols_cov

###########################################################
##! class RITSCandidateModel
##! <<Model>> RITSCandidateModel
type RITSCandidateModel* = ref object
  ##! RITSCandidateModel: +after_model: OLSModel
  after_model*: OLSModel
  ##! RITSCandidateModel: +before_model: OLSModel
  before_model*: OLSModel
  ##! RITSCandidateModel: +after_residual_model: OLSModel
  after_residual_model*: OLSModel
  ##! RITSCandidateModel: +before_residual_model: OLSModel
  before_residual_model*: OLSModel
  ##! RITSCandidateModel: +loglikelihood: float
  loglikelihood*: float
  ##! RITSCandidateModel: +change_point_candidate: int
  change_point_candidate*: int
  ##! RITSCandidateModel: +parameters_mean: vector
  parameters_mean*: vector
  ##! RITSCandidateModel: +covariance_structure_type: string
  covariance_structure_type*: string
  #
  ##! RITSCandidateModel: -before_model_covariance_parameters: matrix
  before_model_covariance_parameters: matrix
  ##! RITSCandidateModel: -existence_change_point_hypothesis: HypothesisScore[ChiSquare]
  existence_change_point_hypothesis*: HypothesisScore[ChiSquare]

##! class RITSModel
##! <<Model>> RITSModel
##! RITSModel <-- RITSCandidateModel
type RITSModel* = ref object
  ##! RITSModel: +candidates: seq[RITSCandidateModel]
  candidates*: seq[RITSCandidateModel]
  ##! RITSModel: +change_point_candidates: seq[int]
  change_point_candidates*: seq[int]
  ##! RITSModel: +loglikelihood_candidates: seq[float]
  loglikelihood_candidates*: seq[float]
  ##! RITSModel: +loglikelihood: float
  loglikelihood*: float
  ##! RITSModel: +change_point_index: int
  change_point_index*: int
  ##! RITSModel: +change_point_x: float
  change_point_x*: float
  ##! RITSModel: +model*: RITSCandidateModel
  model*: RITSCandidateModel
  ##! RITSModel: +parameters_mean: vector
  parameters_mean*: vector

##! RITSModel: +toString()
proc `$`*(model: RITSModel): string =
  result = "[Robust Interrupted Time Series Model]\n"
  result &= "* Coefficients:\n"
  result &= model.model.before_model.coefficients[0].estimate_as_string(title=true, name="       slope (pre)") & "\n"
  result &= model.model.before_model.coefficients[1].estimate_as_string(title=false, name="   intercept (pre)") & "\n"
  result &= model.model.before_residual_model.coefficients[0].estimate_as_string(title=false, name="    autocorr (pre)") & "\n"
  result &= model.model.after_model.coefficients[0].estimate_as_string(title=false, name="      slope (post)") & "\n"
  result &= model.model.after_model.coefficients[1].estimate_as_string(title=false, name="  intercept (post)") & "\n"
  result &= model.model.after_residual_model.coefficients[0].estimate_as_string(title=false, name="   autocorr (post)") & "\n"
  result &= "* Noise estimates:\n"
  result &= ($model.model.before_residual_model.noise_variance.estimate_as_string(title=false, name="    variance (pre)")) & "\n"
  result &= ($model.model.after_residual_model.noise_variance.estimate_as_string(title=false, name="   variance (post)")) & "\n"
  result &= "* Log-likelihood evaluations:\n"
  result &= " Change-point| Log-likelihood\n"
  for k in 0..model.change_point_candidates.high:
    result &= &"       {model.change_point_candidates[k]:>5} | {model.loglikelihood_candidates[k]:.3f} " & (if model.change_point_candidates[k] == model.change_point_index: "***\n" else: "\n")
  result &= "* Wald test:\n"
  result &= &"  H0: There is not change point at {model.change_point_x}:\n"
  result &= "  H1: At least one of the after C.P. parameters is not 0:\n  "
  result &= ($model.model.existence_change_point_hypothesis).replace("\n", "\n  ")

const USE_APROXIMATION = true
##! RITSModel: ~covariance_matrix_for_before_model(model: RITSCandidateModel, x: vector): matrix 
proc covariance_matrix_for_before_model(model: RITSCandidateModel, x: vector): matrix =
  let
    T = x.len
    #phi = model.before_residual_model.beta_hat[0]
    phi = min(model.before_residual_model.beta_hat[0].abs, 1 - 1e-5)
    covariance_factor = (model.before_model.noise_variance.estimate) / (1 - phi ^ 2)
    #X = x[1..x.high].as_column_vector.data_transform((x, i) => @[1.0, x[0], (if i > model.change_point_candidate: 1.0 else: 0.0), (if i > model.change_point_candidate: x[0] else: 0.0)])
    # In our code the order is slope, intercept, slope, intercept
    X = x[1..x.high].as_column_vector.data_transform((x, i) => @[
      x[0], 
      1.0, 
      (if i > model.change_point_candidate: x[0] else: 0.0),
      (if i > model.change_point_candidate: 1.0 else: 0.0)
    ])
  #Numerical approximation
  when USE_APROXIMATION:
    var
      covariance_approx: matrix = @[]
    for r in 0..2:
      covariance_approx.add @[]
      for c in 0..2:
        covariance_approx[r].add covariance_factor * (phi ^ (abs((c - r) mod T).toFloat))
    let inverse_approximation = covariance_approx.inverse
    var inverted_matrix: matrix = @[]
    for r in 0..(x.high - 1):
      inverted_matrix.add @[]
      for c in 0..(x.high - 1):
        inverted_matrix[r].add 0.0
        if r == c and (r == 0 or r==(x.high - 1)):
          inverted_matrix[r][c] = inverse_approximation[0][0]
        elif r == c:
          inverted_matrix[r][c] = inverse_approximation[1][1]
        elif (r - c).abs <= 1:
          inverted_matrix[r][c] = inverse_approximation[1][0]
    #echo "Error (C) => ", (covariance.inverse - inverted_matrix).abs.mean
    return (X.T * inverted_matrix * X).inverse
  when false:
    if phi.abs <= 0.1:
      var inverted_matrix = zeros(T-1, T-1)
      for r in 0..(x.high - 1):
        for c in 0..(x.high - 1):
          if r == c:
            inverted_matrix[r][c] = 1.0 / covariance_factor
          elif (r - c).abs <= 1:
            inverted_matrix[r][c] = -(phi.abs) / covariance_factor
      #echo "Error (A) => ", (covariance.inverse - inverted_matrix).abs.mean
      return (X.T * inverted_matrix * X).inverse
    elif phi.abs <= 2.0:
      var inverted_matrix = zeros(T-1, T-1)
      for r in 0..(x.high - 1):
        for c in 0..(x.high - 1):
          if r == c and (r == 0 or r==(x.high - 1)):
            inverted_matrix[r][c] = (1.0 - 1.0 / (10.0 * ln(phi.abs)) - pow(phi.abs, 10.0) / (1.0 + pow(phi.abs, 10.0))) / covariance_factor
          elif r == c:
            inverted_matrix[r][c] = (1.0 - 1.0 / (10.0 * ln(phi.abs)) - 2.0 * pow(phi.abs, 10.0) / (1.0 + pow(phi.abs, 10.0))) / covariance_factor
          elif (r - c).abs <= 1:
            inverted_matrix[r][c] = 1.0 / (5.0 * ln(phi.abs))
      #echo "Error (B) => ", (covariance.inverse - inverted_matrix).abs.mean, "  ", phi.abs
      return (X.T * inverted_matrix * X).inverse
  var
    covariance = ones(T-1, T-1)
  for r in 0..(x.high - 1):
    for c in 0..(x.high - 1):
      covariance[r][c] = covariance_factor * (phi ^ (abs((c - r) mod T).toFloat))
  result = (X.T * (covariance.inverse) * X).inverse

##! RITSModel: ~model_candidate(x: vector, y: vector, change_point: int, covariance_structure_type: string): RITSCandidateModel
proc model_candidate(x: vector, y: vector, change_point: int, covariance_structure_type: string): RITSCandidateModel = 
  new result
  let
    pre_x = x[0..change_point].add_intercept
    pre_y = y[0..change_point]
    post_x = x[(change_point + 1)..x.high].add_intercept
    post_y = y[(change_point + 1)..y.high]
    pre_ols_model = ols_model(pre_x, pre_y, @["beta1", "beta0"])
    post_ols_model = ols_model(post_x, post_y - pre_ols_model.predict(post_x).prediction_means, @["delta1", "delta0"])
  var
    pre_residual_ar_model: OLSModel
    post_residual_ar_model: OLSModel
    likelihood_model: float
    
  if covariance_structure_type == "autoregressive":
    pre_residual_ar_model = pre_ols_model.residuals.AR1(include_mean=false, feature_name="phi0")
    post_residual_ar_model = post_ols_model.residuals.AR1(include_mean=false, feature_name="phi1")
    #likelihood_model = pre_residual_ar_model.loglikelihood + pre_residual_ar_model.loglikelihood
    likelihood_model = (
      normal(mean=0, std=pre_residual_ar_model.noise_variance.estimate.sqrt).loglikelihood(pre_ols_model.residuals) +
      normal(mean=0, std=post_residual_ar_model.noise_variance.estimate.sqrt).loglikelihood(post_ols_model.residuals)
    )
  elif covariance_structure_type == "independent":
    pre_residual_ar_model = pre_ols_model.residuals.simple_ls_model_cov_independient(@["phi0"])
    post_residual_ar_model = post_ols_model.residuals.simple_ls_model_cov_independient(@["phi1"])
    #likelihood_model = pre_residual_ar_model.loglikelihood + pre_residual_ar_model.loglikelihood
    likelihood_model = (
      normal(mean=0, std=pre_residual_ar_model.noise_variance.estimate.sqrt).loglikelihood(pre_ols_model.residuals) +
      normal(mean=0, std=post_residual_ar_model.noise_variance.estimate.sqrt).loglikelihood(post_ols_model.residuals)
    )
  elif covariance_structure_type == "exchangeable":
    pre_residual_ar_model = pre_ols_model.residuals.simple_ls_model_cov_exchangeable(@["phi0"])
    post_residual_ar_model = post_ols_model.residuals.simple_ls_model_cov_exchangeable(@["phi1"])
    #likelihood_model = pre_residual_ar_model.loglikelihood + pre_residual_ar_model.loglikelihood
    likelihood_model = (
      normal(mean=0, std=pre_residual_ar_model.noise_variance.estimate.sqrt).loglikelihood(pre_ols_model.residuals) +
      normal(mean=0, std=post_residual_ar_model.noise_variance.estimate.sqrt).loglikelihood(post_ols_model.residuals)
    )
  else:
    echo "Covariance structure unrecognized: ", covariance_structure_type
  when false:
    # In the original paper, the calculation of loglikelihood is made with this formula (original AR(1) model):
    normal(mean=0, std=pre_residual_ar_model.noise_variance.estimate.sqrt).loglikelihood(@[pre_ols_model.residuals[0]])
    # However in our implementation, that line is equivalent to
    pre_residual_ar_model.loglikelihood
  #echo "LL ", normal(mean=0, std=pre_residual_ar_model.noise_variance.estimate.sqrt).loglikelihood(pre_ols_model.residuals), "  ", pre_residual_ar_model.loglikelihood
  RITSCandidateModel(
    before_model: pre_ols_model,
    after_model: post_ols_model,
    before_residual_model: pre_residual_ar_model,
    after_residual_model: post_residual_ar_model,
    loglikelihood: likelihood_model,
    change_point_candidate: change_point,
    covariance_structure_type: covariance_structure_type,
    parameters_mean: hstack(pre_ols_model.beta_hat, pre_residual_ar_model.beta_hat, post_ols_model.beta_hat, post_residual_ar_model.beta_hat)
  )

when false:
  proc existence_change_point_wald_test*(models: seq[RITSCandidateModel]): HypothesisScore[ChiSquare] =
    let
      M = models.len
    var
      joint_mean_parameters = zeros(2 * M, 1)
      joint_cov_parameters = zeros(2 * M, 2 * M)
    for j in 0..models.high:
      joint_mean_parameters[2 * j][0] = models[j].after_model.beta_hat[0] # slope
      joint_mean_parameters[2 * j + 1][0] = models[j].after_model.beta_hat[1] # intercept
    for j in 0..models.high:
      let cov = models[j].before_model_covariance_parameters
      joint_cov_parameters[2 * j][2 * j + 0] = cov[2][2]
      joint_cov_parameters[2 * j][2 * j + 1] = cov[2][3]
      joint_cov_parameters[2 * j + 1][2 * j + 0] = cov[3][2]
      joint_cov_parameters[2 * j + 1][2 * j + 1] = cov[3][3]
    let
      score = (joint_mean_parameters.T * (joint_cov_parameters.inverse) * joint_mean_parameters)[0][0]
      hypothesis = ChiSquare(dof: 2.0 * M.toFloat).htest_score(score)
    echo "SCORE"
    echo joint_cov_parameters
    echo score
    echo hypothesis.p_value()
    echo hypothesis.null_confidence_interval()

##! RITSModel: ~existence_change_point_wald_test*(models: seq[RITSCandidateModel]): HypothesisScore[ChiSquare]
# Optimized and adjusted
# With unbalanced p-values, the previous method cannot work
proc existence_change_point_wald_test*(models: seq[RITSCandidateModel]): HypothesisScore[ChiSquare] =
  let
    M = models.len
    base_distribution = ChiSquare(dof: 2.0)
  var
    scores: seq[float] = @[]
  for j in 0..models.high:
    let
      cov = models[j].before_model_covariance_parameters
      beta = models[j].after_model.beta_hat.as_column_vector
      reduced_cov = @[
        @[cov[2, 2], cov[2, 3]],
        @[cov[3, 2], cov[3, 3]],
      ]
      score = (beta.T * reduced_cov.inverse * beta)[0][0]
    scores.add score
  #let
  #  total_score = M.toFloat * scores.min
  #  hypothesis = ChiSquare(dof: 2.0 * M.toFloat).htest_score(total_score, test_type=rightTailed)
  let
    hypothesis_set = scores.mapIt ChiSquare(dof: 2.0 * M.toFloat).htest_score(it, test_type=rightTailed)
    p_values = hypothesis_set.mapIt it.p_value
    number_rejected = p_values.Benjamini_Hochberg_FDR(alpha=0.05)
    #total_score = M.toFloat * scores.min
    hypothesis = hypothesis_set.sorted(cmp=(h1, h2)=>h1.p_value < h2.p_value)[0]
    #hypothesis = ChiSquare(dof: 2.0 * M.toFloat).htest_score(total_score, test_type=rightTailed)
  #
  when false:
    echo "hypothesis-set:", $hypothesis_set
    echo "hypothesis:", $hypothesis
    echo "number_rejected:", $number_rejected
  return hypothesis

##! RITSModel: ~existence_change_point_wald_test*(models: seq[RITSModel]): HypothesisScore[ChiSquare]
proc existence_change_point_wald_test*(models: seq[RITSModel]): HypothesisScore[ChiSquare] = existence_change_point_wald_test models.mapIt(it.model)

##! RITSModel: +rits_model(x: vector, y: vector, candidate_change_point_start: int, candidate_change_point_end: int, covariance_structure_type: string)
proc rits_model*(x: vector, y: vector, candidate_change_point_start: int, candidate_change_point_end: int, covariance_structure_type: string): RITSModel =
  new result
  if y.len != x.len:
    raise newException(ValueError, "x and y must have the same size")
  var likely_change_point = -1
  result.loglikelihood = NegInf
  for k in candidate_change_point_start..candidate_change_point_end:
    if k < 1 or k > y.high - 2: continue
    let candidate = model_candidate(x, y, k, covariance_structure_type)
    result.change_point_candidates.add k
    result.loglikelihood_candidates.add candidate.loglikelihood
    ####
    #echo "  ", candidate.loglikelihood, "  ", candidate.before_model.beta_hat, candidate.after_model.beta_hat
    #echo "CP:", k
    #echo " pre  ", candidate.before_residual_model.residuals.abs.sum
    #echo " post ", candidate.after_residual_model.residuals.abs.sum
    if candidate.loglikelihood.classify != fcNan and candidate.loglikelihood > result.loglikelihood:
      result.loglikelihood = candidate.loglikelihood
      result.model = candidate
      likely_change_point = k
  #echo likely_change_point, "!!", candidate_change_point_start, "  ", candidate_change_point_end
  result.parameters_mean = result.model.parameters_mean
  #Because estimate the covariance matrix is an expensive operation, we just make it for the best model:
  result.model.before_model_covariance_parameters = result.model.covariance_matrix_for_before_model(x)
  result.model.existence_change_point_hypothesis = existence_change_point_wald_test(@[result.model])
  result.change_point_index = likely_change_point
  result.change_point_x = x[likely_change_point]

when isMainModule:
  if true:
    let vector_x = arange(42)
    var vector_y = @[49.95, 53.01, 55.65, 59.17, 61.92, 65.2, 68.22, 71.16, 73.8, 76.99, 79.97, 82.41, 86.07, 89.05, 91.8, 95, 97.94, 100.96, 103.71, 106.49, 110.03, 109.13, 108.09, 108.76, 106.72, 106.17, 104.18, 104, 103.59, 101.86, 100.93, 100.53, 98.83, 98.47, 97.23, 96.07, 94.99, 93.56, 93.96, 91.36, 91.12, 89.23]
    echo "Model 1"
    let model = rits_model(vector_x, vector_y, 0, 42)
    echo model
    echo "(Self) Joint Wald Test"
    let joint_existence_changepoint_test = existence_change_point_wald_test(@[model, model])
    echo joint_existence_changepoint_test
    echo "(Self) Joint Wald Test"
    let joint_existence_changepoint_test1 = existence_change_point_wald_test(@[model, model, model, model, model, model, model])
    echo joint_existence_changepoint_test1
    echo "[OK]"
    echo "\n\n"

  if not true:
    let vector_x = @[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
    var vector_y = @[0.5, 1.0, 1.5, 2.0, 2.5, 24.0, 28.0, 32.0, 36.0, 40.0]
    for k in vector_y.low..vector_y.high:
      vector_y[k] += randn() * 0.1
    echo "Model 1"
    let model = rits_model(vector_x, vector_y, 3, 9)
    echo model
    echo "(Self) Joint Wald Test"
    let joint_existence_changepoint_test = existence_change_point_wald_test(@[model, model])
    echo joint_existence_changepoint_test
    echo "(Self) Joint Wald Test"
    let joint_existence_changepoint_test1 = existence_change_point_wald_test(@[model, model, model, model, model, model, model])
    echo joint_existence_changepoint_test1
    echo "[OK]"
    echo "\n\n"

  if not true:
    let vector_x1 = @[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
    let vector_x2 = @[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
    var vector_y1 = @[0.5, 1.0, 1.5, 2.0, 2.5, 24.0, 28.0, 32.0, 36.0, 40.0]
    var vector_y2 = vector_x2 .* 5.0
    for k in vector_y1.low..vector_y1.high:
      vector_y1[k] += randn() * 0.1
      vector_y2[k] += randn() * 0.1

    echo "Model 1"
    let model1 = rits_model(vector_x1, vector_y1, 3, 8)
    echo model1
    echo "Model 2"
    let model2 = rits_model(vector_x2, vector_y2, 3, 8)
    echo model2
    echo "Joint Wald Test"
    let joint_existence_changepoint_test = existence_change_point_wald_test(@[model1, model2])
    echo joint_existence_changepoint_test
    echo "[OK]"
    echo "\n\n"

  if not true:
    var vector_x1: vector
    var vector_y1: vector
    var vector_x2: vector
    var vector_y2: vector
    let change_point = 40
    vector_x1 = arange(0, 100)
    vector_x2 = arange(0, 100)
    vector_y1 = hstack(
      arange(0, change_point) * 10 + simulateAR1(change_point, autocorrelation=0.2, noise_variance=0.01), 
      arange(change_point, vector_x1.len) * 10 + simulateAR1(vector_x1.len - change_point, autocorrelation=0.2, noise_variance=0.01),
    )
    vector_y2 = hstack(
      arange(0, change_point) * 10 + simulateAR1(change_point, autocorrelation=0.3, noise_variance=0.01), 
      arange(change_point, vector_x2.len) * 11 + simulateAR1(vector_x1.len - change_point, autocorrelation=0.2, noise_variance=0.01),
    )
    ##
    echo "Model 1"
    let model1 = rits_model(vector_x1, vector_y1, change_point-5, change_point+5)
    echo model1
    echo "Model 2"
    let model2 = rits_model(vector_x2, vector_y2, change_point-5, change_point+5)
    echo model2
    echo "Joint Wald Test"
    let joint_existence_changepoint_test = existence_change_point_wald_test(@[model1, model2])
    echo joint_existence_changepoint_test
    echo "[OK]"
    echo "\n\n"

  if not true:
    var vector_x1: vector
    var vector_y1: vector
    var vector_x2: vector
    var vector_y2: vector
    let change_point = 40
    vector_x1 = arange(0, 100)
    vector_x2 = arange(0, 100)
    vector_y1 = hstack(
      arange(0, change_point) * 10 + simulateAR1(change_point, autocorrelation=0.2, noise_variance=0.01), 
      arange(change_point, vector_x1.len) * 20 + simulateAR1(vector_x1.len - change_point, autocorrelation=0.2, noise_variance=0.01),
    )
    vector_y2 = hstack(
      arange(0, change_point) * 20 + simulateAR1(change_point, autocorrelation=0.3, noise_variance=0.01), 
      arange(change_point, vector_x2.len) * 31 + simulateAR1(vector_x1.len - change_point, autocorrelation=0.2, noise_variance=0.01),
    )
    ##
    echo "Model 1"
    let model1 = rits_model(vector_x1, vector_y1, change_point-5, change_point+5)
    echo model1
    echo "Model 2"
    let model2 = rits_model(vector_x2, vector_y2, change_point-5, change_point+5)
    echo model2
    echo "Joint Wald Test"
    let joint_existence_changepoint_test = existence_change_point_wald_test(@[model1, model2])
    echo joint_existence_changepoint_test
    echo "[OK]"
