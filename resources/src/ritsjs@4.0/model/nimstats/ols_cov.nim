import math
import stats
import strutils
import sequtils
import strformat
import sugar
import core
import distributions
import ols

###########################################################
# OLS
###########################################################
type OLSModel2 = ref object
  residuals*: vector
  sum_squared_errors*: float
  degrees_of_freedom*: float
  variance_matrix_coefficients*: matrix
  #
  include_intercept: bool #Because it could be informative. Maybe.
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


proc simple_ls_model_cov_independient*(y: vector, names: seq[string] = @[]): OLSModel =
  new result
  #if names.len > 0 and names.len != X.cols:
  #  raise newException(ValueError, "incorrect number of feature names")
  #var total_model_variation = norm(Ypred.ravel - Ypred.ravel.mean) ^ 2 #SSE, ESS
  #total_model_variation = norm(y - residuals) ^ 2
  
  let
    residuals = y
    total_model_variation = norm(residuals) ^ 2
    degrees_of_freedom = (residuals.len - 1)
    sse = norm(residuals) ^ 2
    s2: float = sum(residuals * residuals) / (degrees_of_freedom - 2).toFloat
    var_beta_hat = eye(residuals.cols, residuals.cols) * s2
    
  result.include_intercept = true
  result.residuals = y
  result.sum_squared_errors = sse
  result.variance_matrix_coefficients = var_beta_hat
  result.R2 = -1#total_model_variation/(total_model_variation + sse)
  result.adjustedR2 = -1#1.0 - (X.rows + (if include_intercept: -1 else: 0)).toFloat / (X.rows - X.cols).toFloat * (1.0 - result.R2)
  result.beta_hat = @[]
  #
  result.degrees_of_freedom = degrees_of_freedom.toFloat
  result.loglikelihood = normal(mean=0.0, std=s2.sqrt).loglikelihood(residuals)
  
  result.noise_variance = shifted_estimator(
    distribution=invchisquare(dof=result.degrees_of_freedom + 1e-10),
    location=0.0,#result.degrees_of_freedom,
    scale=s2.sqrt, #sse * (degrees_of_freedom.toFloat + 1e-10 - 1.0) / degrees_of_freedom.toFloat
  )
  
  result.coefficients = @[
    shifted_estimator(
      distribution=studentt(dof=result.degrees_of_freedom),
      location=0.0,
      scale=s2
    )
  ]
  #Note s2 and noise_variance.estimate are the same, but s2 avoids to calculate it many times
  
  let ms_model_dof = 1.toFloat#(X.cols + (if include_intercept: -1 else: 0)).toFloat
  let ms_residual_dof = (y.len-1).toFloat#(X.rows - X.cols).toFloat
  result.model_significance = central_f(df1=ms_model_dof, df2=ms_residual_dof).htest_score( 
    score=(total_model_variation/ms_model_dof)/(sse/ms_residual_dof),
    test_type=oneTailed
  )
  if names.len > 0:
    result.feature_names = names
  else:
    result.feature_names = @["X1"]

proc simple_ls_model_cov_exchangeable*(y: vector, names: seq[string] = @[]): OLSModel =
  new result
  #if names.len > 0 and names.len != X.cols:
  #  raise newException(ValueError, "incorrect number of feature names")
  #var total_model_variation = norm(Ypred.ravel - Ypred.ravel.mean) ^ 2 #SSE, ESS
  #total_model_variation = norm(y - residuals) ^ 2
  let
    residuals = y
    T = y.len
    total_model_variation = norm(residuals) ^ 2
    degrees_of_freedom = (residuals.len - 1)
    sse = norm(residuals) ^ 2
    s2 = sum(residuals * residuals) / (degrees_of_freedom - 3).toFloat
    var_beta_hat = eye(residuals.cols, residuals.cols) * s2
  var
    rho_ij = 0.0
  for i in 1..(y.len - 1):
    for j in 1..(y.len - 1):
      if i != j:
        continue
      rho_ij += (residuals[i] * residuals[j])
  rho_ij = rho_ij / ((T - 1) * (T - 2) - 3).toFloat
  rho_ij = rho_ij.sqrt
  
  result.include_intercept = true
  result.residuals = y
  result.sum_squared_errors = sse
  result.variance_matrix_coefficients = var_beta_hat
  result.R2 = -1#total_model_variation/(total_model_variation + sse)
  result.adjustedR2 = -1#1.0 - (X.rows + (if include_intercept: -1 else: 0)).toFloat / (X.rows - X.cols).toFloat * (1.0 - result.R2)
  result.beta_hat = @[rho_ij]
  #
  result.degrees_of_freedom = degrees_of_freedom.toFloat
  result.loglikelihood = normal(mean=0.0, std=s2.sqrt).loglikelihood(residuals)
  
  result.noise_variance = shifted_estimator(
    distribution=invchisquare(dof=result.degrees_of_freedom + 1e-10),
    location=0.0,#result.degrees_of_freedom,
    scale=s2, #sse * (degrees_of_freedom.toFloat + 1e-10 - 1.0) / degrees_of_freedom.toFloat
  )
  result.noise_variance = shifted_estimator(
    distribution=invchisquare(dof=result.degrees_of_freedom + 1e-10),
    location=rho_ij,#result.degrees_of_freedom,
    scale=s2, #sse * (degrees_of_freedom.toFloat + 1e-10 - 1.0) / degrees_of_freedom.toFloat
  )
  
  let ms_model_dof = 1.toFloat#(X.cols + (if include_intercept: -1 else: 0)).toFloat
  let ms_residual_dof = (y.len-1).toFloat#(X.rows - X.cols).toFloat
  
  result.coefficients = @[
    shifted_estimator(
      distribution=studentt(dof=result.degrees_of_freedom),
      #distribution=central_f(df1=ms_model_dof, df2=ms_residual_dof),
      location=rho_ij,
      scale=1.0 # @TODO: Not used, need further revision
    )
  ]
  #Note s2 and noise_variance.estimate are the same, but s2 avoids to calculate it many times
  
  result.model_significance = central_f(df1=ms_model_dof, df2=ms_residual_dof).htest_score( 
    score=(total_model_variation/ms_model_dof)/(sse/ms_residual_dof),
    test_type=oneTailed
  )
  result.feature_names = @["x1"]


proc simple_ls_model_cov_exchangeable2*(y: vector, names: seq[string] = @[]): OLSModel =
  new result
  let
    X = constant_matrix(y.len, 1, 1.0)
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



proc `$`*(model: OLSModel2): string =
  result = "[Ordinary Least Squares Model]\n"
  result &= fmt"* Design matrix: {model.coefficients.len}x{model.residuals.len}" & (if model.include_intercept: " (include intercept)\n" else: "\n")
  result &= "* Coefficients:\n"
  #result &= model.coefficients_as_string & "\n"
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
    
