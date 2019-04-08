riot.tag2('rts-main-window', '<window-decorator ref="window"> <yield to="title"> <span class="btn-custom" role="button"><img src="./images/app/app.png" style="width: 23px;height: 23px;position: relative;top: 7px;display: block;"></span> <span class="label">Robust Time Series Toolbox - KAUST Biostatistics Group</span> </yield> <yield to="left-title"> <span class="btn-custom mif-eyedropper load-test-data"></span> </yield> <yield to="right-title"> <span class="btn-custom mif-home open-home-page"></span> </yield> <yield to="menu-bar"> <rts-ribbon-menu-bar ref="menubar"></rts-ribbon-menu-bar> </yield> <yield to="content"> <rts-three-column-panel ref="panels"></rts-three-column-panel> </yield> </window-decorator> </div>', 'rts-main-window .unit-miniplot,[data-is="rts-main-window"] .unit-miniplot{ width: 100px !important; height: 100px !important; } rts-main-window .window-content,[data-is="rts-main-window"] .window-content{ height: 100% !important; overflow: hidden; } rts-main-window .plot .apexcharts-title-text,[data-is="rts-main-window"] .plot .apexcharts-title-text{ color: #000001 !important; font-family: -apple-system, system-ui, BlinkMacSystemFont, "Segoe UI", "Roboto", "Ubuntu", "Helvetica Neue", sans-serif !important; font-family: \'Segoe UI Web (West European)\' !important; font-weight: 200 !important; font-size: 16pt !important; }', '', function(opts) {


    const require_ = require;

    const _global = window;
    const app_state = () => _global.__RTSdata;
    const global_state = () => _global.__gRTSdata;

    const update_unit_names_in_outline = (unit_names) => {
      self.refs.window.refs.panels.change_unit_names(unit_names);
    };
    const update_models_in_summary_reports = (models) => {
      self.refs.window.refs.panels.change_summary_models(models);
      self.refs.window.refs.panels.change_report_models(models);
    };

    const raise_resize_event = () => {
      setTimeout(() => {
        window.dispatchEvent(new Event('resize'));
      }, 30);
    };

    const connect_event_with_component_action = (event_name, component_receiver, action_name, condition = null) => {
      condition = !!condition ? condition : ((d) => true);
      self.on(event_name, (data) => {
        console.log("NEW  calling => ", component_receiver, "=>", action_name, " [", data, "]");
        if (condition(data)) {
          component_receiver[action_name]();
        }
      });
    };

    const event_registry = (event_name, component_trigger, event_name_in_component) => {
      event_name_in_component = !!event_name_in_component ? event_name_in_component : event_name;
      event_name_in_component = event_name_in_component.toLowerCase();
      event_name = event_name.toLowerCase();
      component_trigger.on(event_name_in_component, (data) => {
        console.log("Default calling => ", event_name_in_component, "=>", event_name, " [", data,
          "]");
        self.trigger(event_name, data);
      });
    };

    const show_loading_infobox = (title, content, operation, ms_before_start = 100, ms_before_close = 100) => {
      Metro.infobox.create(`
                    <h3>${title}</h3>
                    <p>${content}</p>
                `, "", {
        closeButton: false,
        onOpen: function () {
          setTimeout(() => {
            const ib = jQuery(this).data("infobox");
            try {
              operation();
            } catch (e) {
              Metro.infobox.create(`
                                        <h3>Processing error</h3>
                                        <p>${e.message}</p>
                                    `, "warning");
            } finally {
              setTimeout(() => {
                ib.close();
              }, ms_before_close);
            }
          }, ms_before_start);
        }
      });
    };

    const register_events_in_toolbar = (menubar) => {
      event_registry('app:menubar:file', menubar);
      event_registry('app:menubar:data', menubar);
      event_registry('app:menubar:plots', menubar);
      event_registry('app:menubar:export', menubar);
      event_registry('app:menubar:help', menubar);
      event_registry('app:load:csv', menubar);
      event_registry('app:load:xls', menubar);
      event_registry('app:settings:modeldates', menubar);
      event_registry('app:settings:changepoint', menubar);
      event_registry('app:view:summary', menubar);
      event_registry('app:view:units', menubar);
      event_registry('app:view:rawtimeseries', menubar);
      event_registry('app:view:estimatedtimeseries', menubar);
      event_registry('app:view:modelbeforechangepoint', menubar);

      event_registry('app:report:view:full', menubar);
      event_registry('app:export:report:docx', menubar);
      event_registry('app:export:report:pdf', menubar);
      event_registry('app:export:tables:regression', menubar);
      event_registry('app:export:tables:analysis', menubar);
      event_registry('app:export:tables:inference', menubar);
      event_registry('app:help:model', menubar);
      event_registry('app:help:program', menubar);

      self.root.querySelector(".load-test-data").onclick = () => process_data(read_test_data_source);
      self.root.querySelector(".open-home-page").onclick = () => require_("electron").remote.shell.openExternal(
        "https://biostats.kaust.edu.sa");

    };

    const register_events_in_leftpanel_panel = () => {
      event_registry('app:request:current:plot:allunits', self.refs.window.refs.panels.refs.leftpanel.refs
        .plots);
      event_registry('app:request:current:unit:allplots', self.refs.window.refs.panels.refs.leftpanel.refs
        .units);
    };

    const register_events_in_date_settings_panel = () => {

    };

    const link_visual_interface_behaviour = () => {
      connect_event_with_component_action('app:view:rawtimeseries', self.refs.window.refs.panels.refs
        .leftpanel, 'view_list_unit_names');
      connect_event_with_component_action('app:view:estimatedtimeseries', self.refs.window.refs.panels.refs
        .leftpanel, 'view_list_unit_names');
      connect_event_with_component_action('app:view:modelbeforechangepoint', self.refs.window.refs.panels.refs
        .leftpanel, 'view_list_unit_names');

      connect_event_with_component_action('app:view:units', self.refs.window.refs.panels.refs.leftpanel,
        'view_list_plot_types');
    };

    const register_events_in_main_app = () => {
      event_registry('app:request:update:dataset', self.refs.window.refs.panels.refs.date_configuration);
      let current_dataset_info = null;
      self.on('app:request:update:dataset', (dataset_info) => {
        model_data(dataset_info);
        current_dataset_info = dataset_info;
      });
      self.on('app:settings:modeldates', () => {
        self.refs.window.refs.panels.show_data_configuration(null, 1, current_dataset_info);
      })
      self.on('app:settings:changepoint', () => {
        self.refs.window.refs.panels.show_data_configuration(null, 2, current_dataset_info);
      })
    }

    const create_interactive_thumbnails = (plot_thumbnails, models, selector_cb) => {
      for (let index = 0; index < models.length; index++) {
        let element = self.root.querySelector(selector_cb(index));
        plot_thumbnails.getChart(element, "snapshoot", models[index])
      }
    };
    const create_thumbnails = (data_source) => {
      data_source.units.forEach((unit_name, k) => {
        const data_values = data_source.valuesOfUnit(unit_name);
        const max_data_values = np.max(data_values);
        const min_data_values = np.min(data_values);
        const L = 100.0;
        const N = data_values.length;
        const m = 50;
        let x_y = [
          [0, 100]
        ];
        for (let _t = 0; _t < m; _t++) {
          let d = L - L * (data_values[Number.parseInt(_t * N / m)] - min_data_values) / (
            max_data_values - min_data_values);
          x_y.push([_t * L / m, d - 6]);
        }
        x_y.push([100, 100]);

        var elements = SVG.select('.unit-miniplot.unit-' + (k + 1)).fill('rgb(58, 139, 199)')
        elements.members.forEach((element) => {
          element.clear();
          let polyline = element.polyline().style({
            "stroke": "#1979ca",
            "stroke-width": "5px",
            "fill": "rgba(25, 121, 202, 0.6)",
          });
          polyline.plot(x_y);
        });

      });
    };

    const check_data_source = () => {
      return true;
    };

    const check_model = (model_parameters) => {
      return true;
    };

    const load_model = (global_state, data_source) => {
      let {
        change_point_date,
        change_point_start_date,
        change_point_end_date
      } = global_state();
      let model_parameters = {};
      if (!!!change_point_date) {
        model_parameters.theoretical_change_point = Number.parseInt(data_source.__data[Object.keys(
          data_source.__data)[0]].length / 2);
        model_parameters.candidates_before = 5;
        model_parameters.candidates_after = 5;
      } else {
        model_parameters.theoretical_change_point = data_source.closeIndexTo(change_point_date);
        model_parameters.candidates_before = model_parameters.theoretical_executive_time_point - data_source
          .closeIndexTo(change_point_start_date) + 1;
        model_parameters.candidates_after = data_source.closeIndexTo(change_point_end_date) -
          model_parameters.theoretical_executive_time_point + 1;
      }
      return model_parameters;
    };

    const read_file_data_source = (data) => {
      return () => {
        let data_source = new RTSModel.ModelDataSourceJson();
        data_source.from(data);
        return data_source;
      }
    };

    const read_test_data_source = () => {

      let data_source = new RTSModel.ModelDataSourceLargeTest();
      return data_source;
    };

    const _process_data = (app_state, global_state, data_reader) => {
      Object.keys(global_state()).forEach((k) => delete global_state()[k])
      let data_source = data_reader();
      if (!check_data_source(data_source)) {
        throw new Error("Inconsistencies in the data source. Please check the input file.");
      }
      global_state().data_source = data_source;
      view_data_source(app_state, global_state, data_source);
    };

    const process_data = (data_reader) => {
      self.refs.window.refs.panels.show_data_configuration(data_reader());
    };

    const model_data = (modeling_info) => {

      const model_parameters = {
        theoretical_change_point: modeling_info.change_point.index.theoretical,
        candidates_before: modeling_info.change_point.index.candidates_before,
        candidates_after: modeling_info.change_point.index.candidates_after,
      }
      const data_source = modeling_info.data_source;
      const number_units = data_source.units;
      let infomodels = [];
      data_source.units.forEach((unit_name, idx) => {
        let infomodel = RTSModel.fit_model(data_source.datesOfUnit(unit_name),
          data_source.valuesOfUnit(unit_name),
          model_parameters.theoretical_change_point,
          model_parameters.candidates_before,
          model_parameters.candidates_after);
        infomodel.unit_name = unit_name;
        infomodels.push(infomodel)
        self.trigger("app:processing:models", {current: idx, total: number_units});
      });
      self.refs.window.refs.panels.refs.plot_collection.opts.models = infomodels;
      self.refs.window.refs.panels.refs.plot_collection.update();

      update_unit_names_in_outline(data_source.units);
      update_models_in_summary_reports(infomodels);
      create_thumbnails(data_source);
      self.refs.window.refs.menubar.enable_model_related_buttons();
      self.trigger('app:view:summary');
    };

    const view_data_source = (app_state, global_state, data_source) => {
      var t0
      t0 = new Date();
      self.refs.window.refs.panels.show_main_container();
      console.log("========= A", new Date() - t0);
      t0 = new Date();

      Object.keys(app_state()).forEach((k) => delete app_state()[k])
      let model_parameters = load_model(global_state, data_source);
      if (!check_model(model_parameters)) {
        throw new Error("Inconsistencies in the model. Please check the parameters.");
      }
      console.log("========= B", new Date() - t0);
      t0 = new Date();
      let infomodels = [];
      data_source.units.forEach((unit_name, idx) => {
        let infomodel = RTSModel.fit_model(data_source.datesOfUnit(unit_name),
          data_source.valuesOfUnit(unit_name),
          model_parameters.theoretical_change_point,
          model_parameters.candidates_before,
          model_parameters.candidates_after);
        infomodel.unit_name = unit_name;
        infomodels.push(infomodel)
      });
      console.log("========= C", new Date() - t0);
      t0 = new Date();
      self.refs.window.refs.panels.refs.plot_collection.opts.models = infomodels;
      self.refs.window.refs.panels.refs.plot_collection.update();

      update_unit_names_in_outline(data_source.units);
      update_models_in_summary_reports(infomodels);
      t0 = new Date();
      create_thumbnails(data_source);

      console.log("========= D", new Date() - t0);
      self.refs.window.refs.menubar.enable_model_related_buttons();
      self.trigger('app:view:summary');
      console.log("========= H", new Date() - t0);
      t0 = new Date();

      app_state().model_parameters = model_parameters;
      app_state().data_source = data_source;
      app_state().models = infomodels;
    }

    const visualizing_dataset_events = (RTSdata, selector_other_panels = "") => {
      const change_title_unit_options = (title) => {
        self.refs.window.refs.panels.change_title_unit_choose(title)
      };
      let current_plot_type = null;
      const base_categories = ["time-series", "likelihood", "time-series-estimation",
        "pre-event-residuals-histogram", "pre-event-residuals-acf", "post-event-residuals-histogram",
        "post-event-residuals-acf"
      ];

      self.on('app:view:rawtimeseries', () => {
        current_plot_type = ["plain"];
        change_title_unit_options("Original series: Choose a unit")

      });
      self.on('app:view:estimatedtimeseries', () => {
        current_plot_type = ["estimation", "loglikelihood"];
        change_title_unit_options("Estimated series: Choose a unit")

      });
      self.on('app:view:modelbeforechangepoint', () => {
        current_plot_type = ["box-plot-residuals", "combined-change-point-residuals",
          "combined-change-point-autocorrelation"
        ];
        change_title_unit_options("Model before change-point: Choose a unit")

      });

      self.on('app:request:current:unit:allplots', (params) => {
        const {
          unitindex,
          unitname
        } = params;
        self.refs.window.refs.panels.show_filtered_plots(unitindex, unitname, current_plot_type);

      });
    };

    const register_help_reports = () => {
      self.on('app:help:model', () => {
        self.refs.window.refs.panels.show_paper_model_1();
      });
    };

    const visualizing_executive_summary = () => {
      self.on('app:view:summary', () => {
        Array.from(self.root.querySelectorAll(".main-panel-container")).forEach((el) => el.classList.add(
          "hidden"));
        self.root.querySelector(".content-container").classList.remove("hidden");
        self.root.querySelector(".executive-summary").classList.remove("hidden");
        self.refs.window.refs.panels.close_left_panel();
        self.refs.window.refs.panels.refs.data_summary.update();

      });
    };

    const open_file = () => {
      const file_filter = [{
          name: 'Data source (.CSV)',
          extensions: ['csv']
        },
        {
          name: 'All Files',
          extensions: ['*']
        }
      ];
      const {
        dialog
      } = require_("electron").remote;
      const fs = require_('fs');
      const paths = dialog.showOpenDialog({
        filters: file_filter,
        properties: ['openFile'],
      });
      if (!paths) {
        return;
      }
      fs.readFile(paths[0], 'utf8', (err, data) => {
        if (err) {
          return console.error(err);
        }
        papaparse.parse(data, {
          keepEmptyRows: false,
          skipEmptyLines: true,
          complete: function (results) {
            if (results.errors.length > 0)
              throw new Error(
                "CSV file cannot be processed. Please contact the administrator."
              );
            process_data(read_file_data_source(results.data));

          }
        });
      });

    };

    const modify_model_events = (app_state, global_state) => {

      let filter_end_date = new Date("01 Jan 2100");
      let filter_start_date = new Date("01 Jan 1900");

      let change_point_date = new Date("01 Jan 2000");
      let change_point_start_date = new Date("01 Jan 1900");
      let change_point_end_date = new Date("01 Jan 2100");

      const date_errors = (message) => {
        Metro.infobox.create(`<h3>Date error</h3>
                    <p>${message}</p>`, "warning");
      };

      const update_model = () => {
        global_state().change_point_date = change_point_date
        global_state().change_point_start_date = change_point_start_date
        global_state().change_point_end_date = change_point_end_date
        view_data_source(app_state, global_state, app_state().data_source.clone());
      };

      const apply_changes = () => {
        let data_source = global_state().data_source.clone();
        data_source.filter(filter_start_date, filter_end_date);
        try {
          if (data_source.__data[Object.keys(data_source.__data)[0]].length == 0) {
            throw new Error("Invalid date!");
          }
          view_data_source(app_state, global_state, data_source);
        } catch (e) {
          Metro.infobox.create(`<h3>Error processing the data</h3>
                    <p>Please check the parameters of the model. </p>`, "warning");
        }
      };

      const check_value_in_dataset = (date) => {
        let ds = global_state().data_source;
        let units = ds.units[0];
        let min_date = ds.__data[units][0][0];
        let max_date = ds.__data[units][ds.__data[units].length - 1][0];
        return date >= min_date && date <= max_date;
      };

      self.on("app:request:change:model:date:start", (event_data) => {
        console.log("#1", 123456, filter_start_date, filter_end_date);
        const {
          date
        } = event_data;
        if (date - filter_end_date >= 0) {
          date_errors("Start date must be before than the finishing date.");
        } else {
          filter_start_date = date;
          apply_changes();
        }
      });
      self.on("app:request:change:model:date:end", (event_data) => {
        console.log("#2", 123456, filter_start_date, filter_end_date);
        const {
          date
        } = event_data;
        if (filter_start_date - date >= 0) {
          date_errors("Finish date must be after than the starting date.");
        } else {
          filter_end_date = date;
          apply_changes();
        }
      });
      self.on("app:request:change:model:changepoint", (event_data) => {
        const {
          reference,
          start,
          end
        } = event_data;
        const conditions = (
          (reference >= filter_start_date && reference <= filter_end_date) &&
          (start >= filter_start_date && start <= filter_end_date) &&
          (end >= filter_start_date && end <= filter_end_date) &&
          check_value_in_dataset(reference) &&
          check_value_in_dataset(start) &&
          check_value_in_dataset(end)
        )
        if (!conditions) {
          date_errors("Change-point candidates must be in the dataset range.");
        } else {
          change_point_date = reference;
          change_point_start_date = start;
          change_point_end_date = end;
          update_model();
        }
      });
    };

    const register_events_reports = () => {
      self.on("app:report:view:full", () => {
        self.refs.window.refs.panels.close_left_panel();
        self.refs.window.refs.panels.refs.data_summary.update();
        self.refs.window.refs.panels.show_full_report();
      });
      self.on("app:export:report:docx", () => {
        self.refs.window.refs.panels.close_left_panel();
        self.refs.window.refs.panels.refs.data_summary.update();
        self.refs.window.refs.panels.show_full_report();
        self.refs.window.refs.panels.save_full_report_as_docx();
      });
      self.on("app:export:report:pdf", () => {
        self.refs.window.refs.panels.close_left_panel();
        self.refs.window.refs.panels.refs.data_summary.update();
        self.refs.window.refs.panels.show_full_report();
        self.refs.window.refs.panels.save_full_report_as_pdf();
      });
    }

    const self = this;
    const config = opts;

    self.on("mount", () => {

      jQuery(document).bind('keydown', 'ctrl+shift+r', () => {
        window.location.reload(true);
      });
      jQuery(document).bind('keydown', 'ctrl+shift+e', () => {
        require_("electron").remote.getCurrentWindow().toggleDevTools();
      });

      register_events_in_toolbar(self.refs.window.refs.menubar);
      register_events_in_leftpanel_panel();
      register_events_in_date_settings_panel();
      register_events_in_main_app();
      register_events_reports();

      link_visual_interface_behaviour();

      modify_model_events(app_state, global_state);

      _global.__gRTSdata = null;
      _global.__gRTSdata = {
        data_source: null,
      };
      _global.__RTSdata = null;
      _global.__RTSdata = {
        model_parameters: null,
        data_source: null,
        models: null,
      }
      visualizing_dataset_events(_global.__RTSdata, ".executive-summary, .paper-pdf");
      visualizing_executive_summary(".executive-summary", ".plot_container > *");
      register_help_reports();

      self.on('app:load:csv', open_file);
    });
});