
riot.tag2('rts-data-block-summary', '<div class="model-title" style="" colspan="2"> <span class="mif-dashboard icon"></span> <span>{opts.model.unit_name}</span> </div> <div class="model-description"> <div class="model-mini-plot"> <rts-model-plot type="estimation" ref="plain" id="{\'executive-plot-\' + opts.index}" model="{opts.model}" toolbar="{false}" header="{false}" legend="{false}"></rts-model-plot> </div> </div> <div class="model-mini-info"> <ul class="v-menu model-information" style="width:100%;"> <li class="menu-title">Change-point estimate</li> <li class="change-point">{change_point(opts.model)}</li> <li class="menu-title">Change in trend estimate</li> <li class="info-before-TET"> {mean_estimate[\'slope\'](opts.model, 0)} </li> <li class="info-after-TET"> {mean_estimate[\'slope\'](opts.model, 1)} </li> <li class="info-diff"> {estimate_difference_on[\'slope\'](opts.model)} </li> <li class="info-95CI"> {confidence_interval[\'slope\'](opts.model, 1)} </li> <li class="menu-title">Change in level estimate</li> <li class="info-diff"> {estimate_difference_on[\'level\'](opts.model)} </li> <li class="info-95CI"> {confidence_interval[\'level\'](opts.model, 1)} </li> <li class="menu-title">Change in {correlation_change_label} estimate</li> <li class="info-before-TET"> {mean_estimate[\'autocorrelation\'](opts.model, 0)} </li> <li class="info-after-TET"> {mean_estimate[\'autocorrelation\'](opts.model, 1)} </li> <li class="info-diff"> {estimate_difference_on[\'autocorrelation\'](opts.model)} </li> <li class="menu-title">Supremum Wald test</li> <li class="info-estimate-TET"> {estimate_difference_on[\'existence_change_point_hypothesis\'](opts.model, 0)} </li> <li class="info-pvalue-TET"> {pvalues_difference_on[\'existence_change_point_hypothesis\'](opts.model, 0)} </li> </ul> </div> </div>', 'rts-data-block-summary *,[data-is="rts-data-block-summary"] *{ cursor: default !important; } rts-data-block-summary .model-mini-plot,[data-is="rts-data-block-summary"] .model-mini-plot{ width: 50%; position: absolute; top: 0; left: 0; } rts-data-block-summary .model-mini-info,[data-is="rts-data-block-summary"] .model-mini-info{ width: 50%; position: relative; right: 0; left: 50%; } rts-data-block-summary .model-title,[data-is="rts-data-block-summary"] .model-title{ width: 100%; background-color: #efefef; font-size: 15pt; font-weight: 300; } rts-data-block-summary .model-information li,[data-is="rts-data-block-summary"] .model-information li{ font-weight: 500; font-size: 11pt; } rts-data-block-summary .model-information .menu-title,[data-is="rts-data-block-summary"] .model-information .menu-title{ font-size: 12pt; font-weight: 200; } rts-data-block-summary .model-information .change-point-TET,[data-is="rts-data-block-summary"] .model-information .change-point-TET{ color: rgb(75, 22, 22); } rts-data-block-summary .model-information .info-before-TET,[data-is="rts-data-block-summary"] .model-information .info-before-TET{ font-weight: normal; } rts-data-block-summary .model-information .info-after-TET,[data-is="rts-data-block-summary"] .model-information .info-after-TET{ font-weight: normal; } rts-data-block-summary .model-information .info-diff,[data-is="rts-data-block-summary"] .model-information .info-diff{ color: rgb(19, 38, 66); color: rgb(15, 56, 56); } rts-data-block-summary .model-information .info-95CI,[data-is="rts-data-block-summary"] .model-information .info-95CI{ color: rgb(15, 56, 56); } rts-data-block-summary .model-information li::before,[data-is="rts-data-block-summary"] .model-information li::before{ font-size: 10pt; font-weight: normal; width: 160px; display: inline-block; } rts-data-block-summary .model-information .change-point::before,[data-is="rts-data-block-summary"] .model-information .change-point::before{ content: "Estimated change-point:"; } rts-data-block-summary .model-information .info-before-TET::before,[data-is="rts-data-block-summary"] .model-information .info-before-TET::before{ content: "Before change-point:"; } rts-data-block-summary .model-information .info-after-TET::before,[data-is="rts-data-block-summary"] .model-information .info-after-TET::before{ content: "After change-point:"; } rts-data-block-summary .model-information .info-diff::before,[data-is="rts-data-block-summary"] .model-information .info-diff::before{ content: "Difference:"; } rts-data-block-summary .model-information .info-95CI::before,[data-is="rts-data-block-summary"] .model-information .info-95CI::before{ content: "95% C.I.:"; } rts-data-block-summary .model-information .info-estimate-TET::before,[data-is="rts-data-block-summary"] .model-information .info-estimate-TET::before{ content: "Estimate:"; } rts-data-block-summary .model-information .info-pvalue-TET::before,[data-is="rts-data-block-summary"] .model-information .info-pvalue-TET::before{ content: "p-value:"; } rts-data-block-summary .model-information-table .model-title,[data-is="rts-data-block-summary"] .model-information-table .model-title{ width: 100%; background-color: #efefef; font-family: -apple-system, system-ui, BlinkMacSystemFont, "Segoe UI", "Roboto", "Ubuntu", "Helvetica Neue", sans-serif; font-size: 15pt; font-weight: 300; } rts-data-block-summary .model-mini-plot rts-model-plot>div,[data-is="rts-data-block-summary"] .model-mini-plot rts-model-plot>div{ border: none !important; box-shadow: none !important; } rts-data-block-summary .model-mini-plot rts-model-plot,[data-is="rts-data-block-summary"] .model-mini-plot rts-model-plot{ margin-to3p: 5px; margin-bot3tom: 5px; max-width: 110%; min-width: 110%; width: 110% !important; min-height: 450px; max-height: 450px; height: 450px !important; margin-left: -10%; } rts-data-block-summary .model-mini-plot rts-model-plot [data-role="appbar"],[data-is="rts-data-block-summary"] .model-mini-plot rts-model-plot [data-role="appbar"],rts-data-block-summary .model-mini-plot rts-model-plot .image-container h4,[data-is="rts-data-block-summary"] .model-mini-plot rts-model-plot .image-container h4{ display: none; }', '', function(opts) {


    const self = this;
    const config = opts;
    config.model = !!config.model ? config.model : null;
    self.change_point = (model) => moment(model.dates[model.estimations.change_point_index]).format(
      "MMM DD, YYYY");
    console.log(config)

    const covariance_structure_type = config.model == null? "autoregressive": config.model.covariance_structure_type;
    self.correlation_change_label = (
      covariance_structure_type == "autoregressive"? "adjacent correlation": (
        covariance_structure_type == "independent"? "variance": (
          covariance_structure_type == "exchangeable"? "correlation": "<<UNKNOWN TYPE>>"
        )
      )
    );

    const rnd = (f) => f.fmt("7.4g")
    const mean_estimate_difference = (model, estimator) => model.estimations.increment_change[estimator].mean
    const mean_estimate_difference_acf = (model, estimator) => model.estimations.increment_change[estimator].mean - model.estimations.initial[estimator].mean
    self.estimate_difference_on = {
      slope: (model) => rnd(mean_estimate_difference(model, "slope")),
      level: (model) => rnd(mean_estimate_difference(model, "level")),
      intercept: (model) => rnd(mean_estimate_difference(model, "intercept")),
      white_noise: (model) => rnd(mean_estimate_difference(model, "noise")),
      autocorrelation: (model) => rnd(mean_estimate_difference_acf(model, "autocorrelation")),
      existence_change_point_hypothesis: (model) => rnd(model.estimations.existence_change_point_hypothesis.score),
    };

    const pval_rnd = (f) => `${(f < 1e-20 ? "< 1e-20": f.fmt(".4g"))}`
    self.pvalues_difference_on = {
      slope: (model) => pval_rnd(model.estimations.increment_change.slope.p_value),
      level: (model) => pval_rnd(model.estimations.increment_change.level.p_value),
      intercept: (model) => pval_rnd(model.estimations.increment_change.intercept.p_value),
      white_noise: (model) => pval_rnd(model.estimations.increment_change.noise.p_value),
      autocorrelation: (model) => pval_rnd(model.estimations.increment_change.autocorrelation.p_value),
      existence_change_point_hypothesis: (model) => pval_rnd(model.estimations.existence_change_point_hypothesis.p_value),
    };

    const mean_estimates = (model, estimator) => [model.estimations.initial[estimator].mean, model.estimations.initial[estimator].mean + model.estimations.increment_change[estimator].mean]
    self.mean_estimate = {
      slope: (model, idx) => rnd(mean_estimates(model, "slope")[idx]),

      intercept: (model, idx) => rnd(mean_estimates(model, "intercept")[idx]),
      white_noise: (model, idx) => rnd(mean_estimates(model, "noise")[idx]),
      autocorrelation: (model, idx) => rnd(mean_estimates(model, "autocorrelation")[idx]),
    };

    const format_confidence_interval = (ci) => `(${rnd(ci[0])}, ${rnd(ci[1])})`
    const confidence_intervals = (model, estimator) => [model.estimations.initial[estimator].confidence_interval, model.estimations.increment_change[estimator].confidence_interval]
    self.confidence_interval = {
      slope: (model, idx) => format_confidence_interval(confidence_intervals(model, "slope")[idx]),
      level: (model, idx) => format_confidence_interval(confidence_intervals(model, "level")[idx]),
      intercept: (model, idx) => format_confidence_interval(confidence_intervals(model, "intercept")[idx]),
      white_noise: (model, idx) => format_confidence_interval(confidence_intervals(model, "noise")[idx]),
      autocorrelation: (model, idx) => format_confidence_interval(confidence_intervals(model, "autocorrelation")[idx]),

    };

    self.on("update", () => {
      self.refs.plain.update();
    });

    self.on("mount", () => {});
});



riot.tag2('rts-data-table-summary', '<div class="par-container"> <h4 class="text-light">Data summary</h4> </div> <div class="paper-style"> <ul class="v-menu" style="width:100%;"> <li class="menu-title"> <div class="row"> <div class="cell-3"> <div></div> </div> <div class="cell-3"> <div>Estimate of slope change</div> </div> <div class="cell-3"> <div>Estimate of level change</div> </div> <div class="cell-3"> <div>Estimated change in {correlation_change_label}</div> </div> </div> </li> <virtual each="{model, index in opts.models}"> <li class=""> <div class="{⁗row ⁗ + (changepoint_is_significant(model) ? ⁗unit-significant⁗ : ⁗⁗ )}"> <div class="cell-3"> <div class="centered-cell unit-name-cell"> <span class="icon"><svg class="{⁗unit-⁗ + (index + 1) + ⁗ unit-miniplot⁗}"></svg></span> <span>{model.unit_name} {changepoint_significance_stars(model)}</span> </div> </div> <virtual each="{parameter in [\'slope\', \'level\', \'autocorrelation\' ]}"> <div class="cell"> <div class="centered-cell"> <span class="mif-move-up mif-2x parameter-indication" if="{estimate_difference_on[parameter](model)>=                     0}"></span> <span class="mif-move-down mif-2x parameter-indication" if="{estimate_difference_on[parameter](model)                     < 0}"></span> <span class="parameter-estimation"> {estimate_difference_on[parameter](model)} </span> </div> </div> <div class="cell-2 parameter-estimates"> <div> <div class="parameter-changes" if="{mean_estimate[parameter]}"> <div> <span class="label">Changes:</span> <span class="values"> <span class="value-before">{mean_estimate[parameter](model, 0)}</span> <span class="mif-forward" style="position: relative;top: 2px;"></span> <span class="value-after">{mean_estimate[parameter](model, 1)}</span> </span> </div> </div> <div class="parameter-confidence-interval" if="{confidence_interval[parameter]}"> <span class="label">95% C.I.:</span> <span class="values"> {confidence_interval[parameter](model, 1)} </span> </div> </div> </div> </virtual> </div> </li> </virtual> </ul> </div>', 'rts-data-table-summary *,[data-is="rts-data-table-summary"] *{ cursor: default; } rts-data-table-summary .unit-miniplot,[data-is="rts-data-table-summary"] .unit-miniplot{ zoom: 0.35; margin-top: 23px; float: left; } rts-data-table-summary .icon+span,[data-is="rts-data-table-summary"] .icon+span{ display: block; width: calc(100% - 50px); height: 75%; position: absolute; left: 45px; top: 11%; } rts-data-table-summary .centered-cell,[data-is="rts-data-table-summary"] .centered-cell{ height: 100%; } rts-data-table-summary .parameter-indication,[data-is="rts-data-table-summary"] .parameter-indication{ float: left; } rts-data-table-summary .parameter-estimation,[data-is="rts-data-table-summary"] .parameter-estimation{ font-size: 115%; text-align: center; } rts-data-table-summary .parameter-changes,[data-is="rts-data-table-summary"] .parameter-changes{ min-height: 1em; } rts-data-table-summary .parameter-changes,[data-is="rts-data-table-summary"] .parameter-changes,rts-data-table-summary .parameter-confidence-interval,[data-is="rts-data-table-summary"] .parameter-confidence-interval{ font-size: 80%; text-align: center; } rts-data-table-summary .parameter-changes .values,[data-is="rts-data-table-summary"] .parameter-changes .values,rts-data-table-summary .parameter-confidence-interval .values,[data-is="rts-data-table-summary"] .parameter-confidence-interval .values{ font-weight: bold; } rts-data-table-summary .centered-cell,[data-is="rts-data-table-summary"] .centered-cell{ display: table; width: 100%; text-align: left; } rts-data-table-summary .centered-cell *:nth-child(2),[data-is="rts-data-table-summary"] .centered-cell *:nth-child(2){ margin-left: 5px; margin-top: 5px; display: inline-block; } rts-data-table-summary .centered-cell .parameter-estimation,[data-is="rts-data-table-summary"] .centered-cell .parameter-estimation{ margin-right: -5px; width: 5px; } rts-data-table-summary .paper-style .menu-title,[data-is="rts-data-table-summary"] .paper-style .menu-title{ height: 2.7em; line-height: 1.2; padding-bottom: 0.1em; padding-top: 0.1em; } rts-data-table-summary .unit-significant.row,[data-is="rts-data-table-summary"] .unit-significant.row{ background-color: rgba(255, 0, 0, 0.05); background-color: rgba(15, 87, 219, 0.1); } rts-data-table-summary .unit-significant .unit-name-cell,[data-is="rts-data-table-summary"] .unit-significant .unit-name-cell{ color: rgb(171, 21, 21); color: rgb(34, 93, 208); } rts-data-table-summary .cell-2.parameter-estimates > div > div,[data-is="rts-data-table-summary"] .cell-2.parameter-estimates > div > div{ min-height: 1.5em; } rts-data-table-summary .cell-2.parameter-estimates > div,[data-is="rts-data-table-summary"] .cell-2.parameter-estimates > div{ position: absolute; top: 0; bottom: 0; margin: auto; }', '', function(opts) {


    const self = this;
    const config = opts;
    config.models = !!config.models ? config.models : [];

    const covariance_structure_type = config.models.length == 0 || config.models[0] == null? "autoregressive": config.models[0].covariance_structure_type;
    console.warn("@covariance_structure_type:", covariance_structure_type)
    self.correlation_change_label = covariance_structure_type == "autoregressive"? "adjacent correlation": (
      covariance_structure_type == "independent"? "variance": (
        covariance_structure_type == "exchangeable"? "correlation": "<<UNKNOWN TYPE>>"
      )
    );

    self.changepoint_is_significant = (model) => model.estimations.existence_change_point_hypothesis.p_value < 0.05
    self.changepoint_significance_stars = (model) => (
      model.estimations.existence_change_point_hypothesis.p_value >= 0.1 ? "" : (
        model.estimations.existence_change_point_hypothesis.p_value >= 0.05 ? "+" : (
          model.estimations.existence_change_point_hypothesis.p_value >= 0.01 ? "*" : (
            model.estimations.existence_change_point_hypothesis.p_value >= 0.001 ? "**" : (
              model.estimations.existence_change_point_hypothesis.p_value >= 0.001 ? "***" :
              "****"
            )
          )
        )
      )
    )

    const rnd = (f) => Math.round(f * 100) / 100;
    const mean_estimate_difference_acf = (model, estimator) => model.estimations.increment_change[estimator].mean - model.estimations.initial[estimator].mean

    self.estimate_difference_on = {
      slope: (model) => rnd(model.estimations.increment_change.slope.mean),
      level: (model) => rnd(model.estimations.increment_change.level.mean),
      intercept: (model) => rnd(model.estimations.increment_change.intercept.mean),
      white_noise: (model) => rnd(model.estimations.increment_change.noise.mean),
      autocorrelation: (model) => rnd(mean_estimate_difference_acf(model, "autocorrelation")),
    };
    const mean_estimates = (model, estimator) => [model.estimations.initial[estimator].mean, model.estimations.initial[estimator].mean + model.estimations.increment_change[estimator].mean]
    const mean_estimates_acf = (model, estimator) => [model.estimations.initial[estimator].mean, model.estimations.increment_change[estimator].mean]
    self.mean_estimate = {
      slope: (model, idx) => rnd(mean_estimates(model, "slope")[idx]),

      intercept: (model, idx) => rnd(mean_estimates(model, "intercept")[idx]),
      white_noise: (model, idx) => rnd(mean_estimates(model, "noise")[idx]),
      autocorrelation: (model, idx) => rnd(mean_estimates_acf(model, "autocorrelation")[idx]),
    };
    const format_confidence_interval = (ci) => `(${rnd(ci[0])}, ${rnd(ci[1])})`
    const confidence_intervals = (model, estimator) => [model.estimations.initial[estimator].confidence_interval, model.estimations.increment_change[estimator].confidence_interval]
    self.confidence_interval = {
      slope: (model, idx) => format_confidence_interval(confidence_intervals(model, "slope")[idx]),
      level: (model, idx) => format_confidence_interval(confidence_intervals(model, "level")[idx]),
      intercept: (model, idx) => format_confidence_interval(confidence_intervals(model, "intercept")[idx]),
      white_noise: (model, idx) => format_confidence_interval(confidence_intervals(model, "noise")[idx]),

    };

    self.on("mount", () => {

    });
});


riot.tag2('rts-data-summary', '<div class="executive-summary hidden" style="" if="{model_exists()}"> <rts-data-table-summary models="{opts.models}"></rts-data-table-summary> <div class="par-container"> <h4 class="text-light">Relevant results</h4> </div> <div class="paper-style"> <div each="{model, index in opts.models}"> <rts-data-block-summary model="{model}" index="{index}"></rts-data-block-summary> </div> </div> </div>', 'rts-data-summary .hidden,[data-is="rts-data-summary"] .hidden{ display: none; } rts-data-summary .par-container,[data-is="rts-data-summary"] .par-container{ width: 80%; margin: auto; } rts-data-summary h4,[data-is="rts-data-summary"] h4{ margin-top: 40px; }', '', function(opts) {


    const self = this;
    const config = opts;
    config.models = !!config.models ? config.models : [];

    self._models = () => config.models;
    self._models_indexes = () => (!!config.models ? Array.from(Array(config.models.length).keys()) : []);

    self.model_exists = () => config.models != null && config.models.length > 0;
    self._trigger_event = (name) => {
      return () => self.trigger('app:request:current:plot:allunits', {
        plottype: name.toLowerCase()
      });
    };

    self.on("mount", () => {

    });
});


