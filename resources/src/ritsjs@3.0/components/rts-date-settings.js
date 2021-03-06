riot.tag2('rts-date-settings', '<div class="date-settings-container"> <div class="par-container"> <h4 ref="settings_range" class="text-light">Date range</h4> </div> <div class="paper-style"> <p>Define the date interval to analyze.</p> <div class="date-range-graph"></div> <div class="dynamic-text"> <p> <span>Start date:</span> <span class="start-date val" ref="tmin"></span> </p> <p> <span>End date:</span> <span class="end-date val" ref="tmax"></span> </p> </div> </div> <div class="par-container"> <h4 ref="settings_theoretical" class="text-light">Theoretical change-point information</h4> </div> <div class="paper-style"> <p>Define the date of the theoretical change-point (or interventation date) making click in the graph.</p> <p>Choose the range of date to be consider as candidates to change-point.</p> <p class="observation">(Wider candidate intervals can take more time to be processed.)</p> <div class="change-point-graph"></div> <div class="dynamic-text"> <p> <span>Theoretical change-point:</span> <span class="theoretical-date val" ref="theoretical"></span> </p> <p> <span>Candidate start date:</span> <span class="candidate-start-date val" ref="theoretical_min"></span> (<span class="candidate-before val" ref="candidates_before"></span> candidates before change-point) </p> <p> <span>Candidate end date:</span> <span class="candidate-end-date val" ref="theoretical_max"></span> (<span class="candidate-after val" ref="candidates_after"></span> candidates after change-point) </p> </div> </div> <div class="par-container"> <button class="button start-model" ref="start_model"> <i class="ms-Icon ms-Icon--TestAutoSolid icon"></i> <span class="caption">Model the dataset</span> </button> </div> </div>', 'rts-date-settings *,[data-is="rts-date-settings"] *{ user-select: none; } rts-date-settings .dynamic-text *,[data-is="rts-date-settings"] .dynamic-text *,rts-date-settings .dynamic-text,[data-is="rts-date-settings"] .dynamic-text{ user-select: text; } rts-date-settings .dynamic-text,[data-is="rts-date-settings"] .dynamic-text{ margin-top: 20px; } rts-date-settings .dynamic-text p span:nth-child(2),[data-is="rts-date-settings"] .dynamic-text p span:nth-child(2){ font-weight: 400; color: rgba(55, 126, 183, 1); } rts-date-settings .par-container,[data-is="rts-date-settings"] .par-container{ width: 80%; margin: auto; } rts-date-settings .date-settings-container,[data-is="rts-date-settings"] .date-settings-container{ width: 100% !important; } rts-date-settings .date-range-graph,[data-is="rts-date-settings"] .date-range-graph{ width: 100% !important; } rts-date-settings .change-point-graph,[data-is="rts-date-settings"] .change-point-graph{ width: 100% !important; } rts-date-settings h4,[data-is="rts-date-settings"] h4{ margin-top: 40px; } rts-date-settings p,[data-is="rts-date-settings"] p{ font-size: 12pt; font-weight: 200; margin-top: 0px; } rts-date-settings .observation,[data-is="rts-date-settings"] .observation{ color: rgb(160, 38, 38); font-size: 95%; } rts-date-settings .start-model:target,[data-is="rts-date-settings"] .start-model:target{ background: rgba(0, 0, 0, 0.45); border: 1px solid rgba(0, 0, 0, 0.45); transition: linear 200ms background-color; } rts-date-settings .start-model:hover,[data-is="rts-date-settings"] .start-model:hover{ background: rgba(55, 176, 238, 0.4); background: rgba(0, 0, 0, 0.25); border: 1px solid rgba(0, 0, 0, 0.3); transition: linear 200ms background-color; } rts-date-settings .start-model,[data-is="rts-date-settings"] .start-model{ background: rgba(55, 176, 238, 0.2); background: rgba(0, 0, 0, 0.1); border: 1px solid rgba(0, 0, 0, 0.2); transition: linear 200ms background-color; margin-bottom: 50px; position: absolute; right: 0px; cursor: default; }', '', function(opts) {


    const self = this;
    const config = opts;

    config.data_source = new RTSModel.ModelDataSourceLargeTest();

    self.change_point = {
      dates: {
        theoretical: null,
        candidates_before: null,
        candidates_after: null,
      },
      index: {
        theoretical: null,
        candidates_before: null,
        candidates_after: null,
      },
    }
    self.data_range = {
      dates: {
        start: null,
        end: null,
      },
      index: {
        start: null,
        end: null,
      },
    };

    const closest_date = (date) => {
      const datatable = self.datatable;
      for (let i = 0; i < datatable.length; i++) {
        if (date <= datatable[i][0]) {
          return [i, date];
        }
      }
      return [datatable.length - 1, datatable[datatable.length - 1][0]];
    };

    self.selected_datasource = () => {
      return {
        change_point: self.change_point,
        data_range: self.data_range,
        data_source: (new RTSModel.ModelDataSourceJson()).fill_from_data_table(config.data_source.units, self
          .datatable.filter((val, idx) => idx >= self.data_range.index.start && idx <= self.data_range.index.end))
      }
    };

    const closest_date_config = (date, type, subtype) => {
      const [index, _date] = closest_date(date);
      self[type].dates[subtype] = _date;
      self[type].index[subtype] = index;
    };

    const update_data_range = (tmin, tmax) => {
      closest_date_config(tmin, "data_range", "start");
      closest_date_config(tmax, "data_range", "end");
      self.refs.tmin.textContent = moment(self.data_range.dates.start).format("MMM D YYYY");
      self.refs.tmax.textContent = moment(self.data_range.dates.end).format("MMM D YYYY");
    };

    const update_change_point = (tmin, tmax, tval) => {
      closest_date_config(tval, "change_point", "theoretical");
      closest_date_config(tmin, "change_point", "candidates_before");
      closest_date_config(tmax, "change_point", "candidates_after");
      self.change_point.index.candidates_before = -self.change_point.index.candidates_before + self.change_point.index.theoretical;
      self.change_point.index.candidates_after = self.change_point.index.candidates_after - self.change_point.index.theoretical;
      self.refs.theoretical_min.textContent = moment(self.change_point.dates.candidates_before).format("MMM D YYYY");
      self.refs.theoretical_max.textContent = moment(self.change_point.dates.candidates_after).format("MMM D YYYY");
      self.refs.theoretical.textContent = moment(self.change_point.dates.theoretical).format("MMM D YYYY");
      self.refs.candidates_before.textContent = self.change_point.index.candidates_before
      self.refs.candidates_after.textContent = self.change_point.index.candidates_after
    };

    const default_graph_options = () => ({
      title: '',
      ylabel: '',
      colors: ["rgb(66, 131, 185)", "rgb(77, 175, 74)", "rgb(228, 59, 45)", "rgb(152, 78, 163)",
        "rgb(242, 126, 51)",
        "rgb(166, 86, 41)", "rgb(153, 153, 153)"
      ],
      showRangeSelector: true,
      rangeSelectorHeight: 100,
      rangeSelectorAlpha: 0.8,
      rangeSelectorForegroundLineWidth: '2px',
      rangeSelectorForegroundStrokeColor: 'rgba(36, 75, 107, 1)',
      rangeSelectorBackgroundLineWidth: '2px',
      rangeSelectorBackgroundStrokeColor: 'rgba(55, 126, 183, 0.1)',
      rangeSelectorPlotLineWidth: '2px',
      rangeSelectorPlotStrokeColor: 'rgba(55, 126, 183, 1)',
      rangeSelectorPlotFillColor: 'rgba(55, 126, 183, 0.2)',
      rangeSelectorPlotFillGradientColor: 'rgba(55, 126, 183, 0.2)',
      legend: "never",
    })

    const synchronizerMasterSlave = (master, slave, callback) => {

      master.updateOptions({
        drawCallback: (me, initial) => {
          let [xmin, xmax] = master.xAxisRange();
          update_data_range(xmin, xmax);
          let new_data = master.rawData_.map((c) => [(new Date(c[0])), ...c.slice(1)]).filter((c) => c[0] >=
            xmin && c[0] <= xmax)
          slave.updateOptions({
            file: new_data
          });
          callback();
        },
      })
    };

    const Crosshair = (function () {

      var crosshair = function (start_height = 0) {
        this.canvas_ = document.createElement("canvas");
        this.marked = {
          x: 1,
          xval: 0
        };
        this.start_height = start_height;
        this.last_height = 0;
      }

      crosshair.prototype.activate = function (g) {
        g.graphDiv.appendChild(this.canvas_);

        return {
          select: this.select,
          deselect: this.deselect
        };
      };

      crosshair.prototype.select = function (e) {
        if (this.direction_ === null) {
          return;
        }

        var width = e.dygraph.width_;
        var height = e.dygraph.height_;
        this.last_height = height;
        this.canvas_.width = width;
        this.canvas_.height = height;
        this.canvas_.style.width = width + "px";
        this.canvas_.style.height = height + "px";

        var ctx = this.canvas_.getContext("2d");
        ctx.clearRect(0, 0, width, height);
        ctx.strokeStyle = "rgba(0, 0, 0, 0.3)";
        ctx.lineWidth = 1;
        ctx.beginPath();

        var canvasx = Math.floor(e.dygraph.selPoints_[0].canvasx) + 0.5;

        ctx.moveTo(canvasx, 0);
        ctx.lineTo(canvasx, height - this.start_height);
        ctx.stroke();
        ctx.closePath();

        this.draw_mark();
      };

      crosshair.prototype.draw_mark = function () {
        var ctx = this.canvas_.getContext("2d");

        if (!!this.marked.x) {
          ctx.beginPath();
          ctx.lineWidth = 2;
          ctx.strokeStyle = "#0074ce";
          ctx.moveTo(this.marked.x, 0);
          ctx.lineTo(this.marked.x, this.last_height - this.start_height);
          ctx.stroke();
          ctx.closePath();
        }
      };

      crosshair.prototype.deselect = function (e) {

      };

      crosshair.prototype.reset = function () {
        var ctx = this.canvas_.getContext("2d");
        ctx.clearRect(0, 0, this.canvas_.width, this.canvas_.height);
        this.marked = {
          x: 1,
          xval: 0
        };
      };

      crosshair.prototype.checkBounds = function (master, slave) {
        if (!this.marked.x) return;
        var ctx = this.canvas_.getContext("2d");
        let [xmin, xmax] = master.xAxisRange();
        this.marked.xval = Math.max(xmin + 2, this.marked.xval);
        this.marked.xval = Math.min(xmax - 2, this.marked.xval);
        this.marked.x = slave.toDomXCoord(this.marked.xval);
        ctx.clearRect(0, 0, this.canvas_.width, this.canvas_.height);
        this.draw_mark();
        update_change_point(xmin, xmax, this.marked.xval);
      };

      crosshair.prototype.destroy = function () {
        this.canvas_ = null;
      };

      return crosshair;
    })();

    let permanent_info = {
      change_point: {
        dates: {
          theoretical: null,
          candidates_before: null,
          candidates_after: null,
        },
        index: {
          theoretical: null,
          candidates_before: null,
          candidates_after: null,
        },
      },
      data_range: {
        dates: {
          start: null,
          end: null,
        },
        index: {
          start: null,
          end: null,
        },
      }
    };
    self.on("update", () => {
      self.datatable = config.data_source.datatable;
      const labels = ["Time", ...config.data_source.units];

      const graph_options_master = default_graph_options();
      graph_options_master.labels = labels;
      if (!!permanent_info.data_range.dates.start) {
        graph_options_master.dateWindow = [permanent_info.data_range.dates.start, permanent_info.data_range.dates
          .end
        ];
      }

      const line_annotation_plugin = new Crosshair(start_height = 125);
      const graph_options_slave = default_graph_options();
      graph_options_slave.labels = labels;
      graph_options_slave.plugins = [line_annotation_plugin];
      if (!!permanent_info.change_point.dates.candidates_before) {
        graph_options_slave.dateWindow = [permanent_info.change_point.dates.candidates_before, permanent_info
          .change_point.dates.candidates_after
        ];
        line_annotation_plugin.marked = {
          x: 1,
          xval: permanent_info.change_point.dates.theoretical
        }
      }

      const check_bounds = () => {
        line_annotation_plugin.checkBounds(self.change_point_graph, self.change_point_graph);
      };
      graph_options_slave.interactionModel = {
        mousedown: function (event, g, context) {
          event.preventDefault();
          line_annotation_plugin.marked.x = Math.floor(g.selPoints_[0].canvasx) + 0.5;
          line_annotation_plugin.marked.xval = g.selPoints_[0].xval;
          line_annotation_plugin.marked.yval = g.selPoints_[0].yval;
          line_annotation_plugin.select({
            dygraph: g
          });
          check_bounds();
        }
      }

      self.date_range_graph = new Dygraph(
        self.root.querySelector(".date-range-graph"),
        self.datatable,
        graph_options_master
      )

      self.change_point_graph = new Dygraph(
        self.root.querySelector(".change-point-graph"),
        self.datatable,
        graph_options_slave
      )

      self.change_point_graph.updateOptions({
        drawCallback: check_bounds
      });

      setTimeout(check_bounds, 100);

      synchronizerMasterSlave(self.date_range_graph, self.change_point_graph, check_bounds);
    });

    self.on("mount", () => {
      self.update();
      self.refs.start_model.addEventListener("click", () => {
        self.trigger("app:request:update:dataset", self.selected_datasource())
      });
    });

    self.scroll_into_data_range_settings = () => {
      self.refs.settings_range.scrollIntoView();
    }

    self.scroll_into_theoretical_settings = () => {
      self.refs.settings_theoretical.scrollIntoView();
    };

    self.update_configuration = (change_point, data_range) => {
      permanent_info.change_point = change_point
      permanent_info.data_range = data_range
      self.update();
    };
});