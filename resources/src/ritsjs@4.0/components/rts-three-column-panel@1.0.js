riot.tag2('rts-three-column-panel', '<div class="split-parent h-100 w-100"> <div class="leftpanel-container panel-container h-100"> <rts-outline-panel class="h-100 w-100" ref="leftpanel"></rts-outline-panel> </div> <div class="plotpanel-container panel-container h-100"> <div class="welcome-splash w-100 main-panel-container"> <img class="icon" src="../resources/images/app/app.png"> <div class="title"> Robust Interrupted Time Series Toolbox </div> <div class="research-group"> KAUST Biostatistics group </div> <div> To start choose a dataset (in CSV format) and select the model parameters. </div> <div style="padding-top: 40px; font-size:90%; "> <img src="../resources/images/app/KAUST-logo.png" style="height:100px;"> <div>Development and research supported by KAUST research grant to the KAUST Biostatistics research group (biostats.kaust.edu.sa)</div> </div> </div> <div class="date-configuration main-panel-container hidden"> <rts-date-settings ref="date_configuration"></rts-date-settings> </div> <div class="plot-container main-panel-container paper-style hidden"> <rts-plot-collection ref="plot_collection"></rts-plot-collection> </div> <div class="content-container main-panel-container hidden"> <rts-data-summary ref="data_summary"></rts-data-summary> </div> <div class="paper-model-1 paper-pdf w-100 h-100 main-panel-container hidden"> <iframe class="w-100 h-100" src="./lib/pdf.js@2.0/adapted/web/viewer.html?file=../../../../../stats-model/papers/1707.01861 - A Robust Interrupted Time Series Model for Analyzing Complex Healthcare Intervention Data.pdf#zoom=page-width&"></iframe> </div> <div class="report-container main-panel-container paper-style hidden"> <rts-report ref="full_data_report"></rts-report> </div> </div> <div class="rightpanel-container panel-container h-100"> </div> </div>', 'rts-three-column-panel .welcome-splash.hidden,[data-is="rts-three-column-panel"] .welcome-splash.hidden{ display: none; } rts-three-column-panel .welcome-splash,[data-is="rts-three-column-panel"] .welcome-splash{ text-align: center; margin-top: 50px; font-size: 12pt; font-weight: 200; } rts-three-column-panel .welcome-splash .icon,[data-is="rts-three-column-panel"] .welcome-splash .icon{ max-width: 250px; opacity: 0.5; } rts-three-column-panel .welcome-splash .title,[data-is="rts-three-column-panel"] .welcome-splash .title{ font-size: 18pt; font-weight: 200; } rts-three-column-panel .welcome-splash .research-group,[data-is="rts-three-column-panel"] .welcome-splash .research-group{ font-size: 12pt; font-weight: normal; margin-bottom: 30px; } rts-three-column-panel .plot_container [class^=\'apexcharts\'],[data-is="rts-three-column-panel"] .plot_container [class^=\'apexcharts\']{ cursor: default !important; } rts-three-column-panel .plot_container .apexcharts-menu.open,[data-is="rts-three-column-panel"] .plot_container .apexcharts-menu.open{ background-color: #dadbdc; background-color: #f5f6f7; border: 1px solid rgba(0, 0, 0, 0.2); } rts-three-column-panel .plot_container .apexcharts-menu.open .apexcharts-menu-item,[data-is="rts-three-column-panel"] .plot_container .apexcharts-menu.open .apexcharts-menu-item{ margin: 1px; border: 1px solid transparent; } rts-three-column-panel .plot_container .apexcharts-menu.open .apexcharts-menu-item:hover,[data-is="rts-three-column-panel"] .plot_container .apexcharts-menu.open .apexcharts-menu-item:hover{ background-color: rgba(164, 206, 249, 0.2); border: 1px solid #a4cef9; } rts-three-column-panel .plot_container [class^=\'apexcharts\'].selected,[data-is="rts-three-column-panel"] .plot_container [class^=\'apexcharts\'].selected{ border: none; box-shadow: none !important; } rts-three-column-panel .plot_container [class^=\'apexcharts\'].selected::before,[data-is="rts-three-column-panel"] .plot_container [class^=\'apexcharts\'].selected::before,rts-three-column-panel .plot_container [class^=\'apexcharts\'].selected::after,[data-is="rts-three-column-panel"] .plot_container [class^=\'apexcharts\'].selected::after{ display: none; } rts-three-column-panel .leftpanel-container,[data-is="rts-three-column-panel"] .leftpanel-container{ overflow: hidden; } rts-three-column-panel .panel-container,[data-is="rts-three-column-panel"] .panel-container{ transitio1n: 10ms linear width; min-width: 0px !important; overflow-x: hidden; } rts-three-column-panel .panel-container.hidden,[data-is="rts-three-column-panel"] .panel-container.hidden{ } rts-three-column-panel .panel-container,[data-is="rts-three-column-panel"] .panel-container{ } rts-three-column-panel .plotpanel-container,[data-is="rts-three-column-panel"] .plotpanel-container{ overflow: auto; height: 100% !important; } rts-three-column-panel .hidden,[data-is="rts-three-column-panel"] .hidden{ display: none; } rts-three-column-panel .content-container,[data-is="rts-three-column-panel"] .content-container{} rts-three-column-panel .paper-style,[data-is="rts-three-column-panel"] .paper-style{ background-color: white; border: 1px solid #bbbec2; box-shadow: 0 10px 16px 0 rgba(0, 0, 0, 0.2), 0 6px 20px 0 rgba(0, 0, 0, 0.19) !important; min-width: 700px !important; width: 80%; position: relative; margin: 15px auto; margin-bottom: 40px; padding: 30px; } rts-three-column-panel .main-viewport,[data-is="rts-three-column-panel"] .main-viewport{ overflow: auto !important; display: block !important; padding: 10px; } rts-three-column-panel .outline-panel.hidden-panel,[data-is="rts-three-column-panel"] .outline-panel.hidden-panel{ width: 1px !important; max-width: 1px !important; min-width: 1px !important; transition: 600ms linear width; } rts-three-column-panel .outline-panel,[data-is="rts-three-column-panel"] .outline-panel{ min-width: 100px !important; transition: 600ms linear width; } rts-three-column-panel .main-viewport>*,[data-is="rts-three-column-panel"] .main-viewport>*{ min-width: 500px; } rts-three-column-panel .main-viewport .apexcharts-canvas .apexcharts-toolbar .hidden,[data-is="rts-three-column-panel"] .main-viewport .apexcharts-canvas .apexcharts-toolbar .hidden{ display: none; } rts-three-column-panel .main-viewport .apexcharts-canvas .apexcharts-toolbar,[data-is="rts-three-column-panel"] .main-viewport .apexcharts-canvas .apexcharts-toolbar{ font-family: unset; font-size: 40pt; } rts-three-column-panel .split-parent,[data-is="rts-three-column-panel"] .split-parent{ display: flex; } rts-three-column-panel .gutter,[data-is="rts-three-column-panel"] .gutter{ background-color: #dadbdc; cursor: ew-resize; background: white; border: 1px solid #dadbdc; } rts-three-column-panel .gutter.gutter-horizontal .gutter-close,[data-is="rts-three-column-panel"] .gutter.gutter-horizontal .gutter-close{ background-color: rgba(58, 139, 199, 1) !important; } rts-three-column-panel .gutter.gutter-horizontal.gutter-2,[data-is="rts-three-column-panel"] .gutter.gutter-horizontal.gutter-2{ display: none; } rts-three-column-panel .plotpanel-container.panel-container,[data-is="rts-three-column-panel"] .plotpanel-container.panel-container{ margin-left: 5px; } rts-three-column-panel .leftpanel-container [ref="leftpanel"],[data-is="rts-three-column-panel"] .leftpanel-container [ref="leftpanel"],rts-three-column-panel .leftpanel-container [ref="rightpanel"],[data-is="rts-three-column-panel"] .leftpanel-container [ref="rightpanel"]{ position: absolute; top: 0; left: 0; }', '', function(opts) {


    const self = this;
    const config = opts;

    let splitter;

    const create_splitter = (ratio_width) => {
      const PO = self.root.querySelector(".leftpanel-container");
      const PP = self.root.querySelector(".plotpanel-container");
      const PC = self.root.querySelector(".rightpanel-container");
      const containers = [PO, PP, PC];

      ratio_width[1] += ratio_width[2];
      ratio_width[2] = 0;

      splitter = Split(containers, {
        sizes: ratio_width,
        minSize: [250, 250, 255].map((c, i) => Math.abs(ratio_width[i] - 0) < 1e-3 ? 0 : c),
        expanedToMin: true,
        gutter: function (index, direction) {
          let gutter = document.createElement('div')
          gutter.className = 'gutter gutter-' + direction + ' gutter-' + index
          gutter.style.height = '100%'
          let gutter_close = document.createElement('div')
          gutter_close.className = 'gutter-close'
          gutter_close.style.cssText = `
                        height: 10%;
                        width: 200%;
                        background-color: rgba(0, 0, 0, 0.5);
                        position: absolute;
                        top: 45%;
                        right: 0px;
                        cursor: cell;
                    `;
          gutter.appendChild(gutter_close);
          gutter_close.onclick = (index == 1) ? self.toggle_left_panel : self
            .toggle_right_panel;
          return gutter
        },
        gutterSize: 4.5,
      });
    };

    const update_panel_widths = (index, action) => {
      action = !!action ? action : "show-hide";
      if (index == 1) return;
      const PO = self.root.querySelector(".leftpanel-container");
      const PP = self.root.querySelector(".plotpanel-container");
      const PC = self.root.querySelector(".rightpanel-container");
      const containers = [PO, PP, PC];
      let widths = containers.map((e) => e.clientWidth);
      const total_width = widths.reduce((a, b) => a + b)

      if (widths[index] > 0 && action.contains("hide")) {
        widths[1] += widths[index];
        widths[index] = 0;
      } else if (action.contains("show")) {
        widths[1] += widths[index] - 150 ;
        widths[index] = 150;
      }
      const ratio_width = widths.map((e) => e * 100.0 / total_width);
      splitter.destroy();
      create_splitter(ratio_width);
      window.dispatchEvent(new Event('resize'));
    };

    const panel_visibility = (selector, action) => {
      self.root.querySelector(selector).classList[action]("hidden");
    };

    self.close_left_panel = (set_blank=false) => {
      if(set_blank){
        self.refs.leftpanel.set_blank();
      }
      update_panel_widths(0, "hide");
    };

    self.show_left_panel = () => {
      update_panel_widths(0, "show");
    };

    self.close_right_panel = () => {
      update_panel_widths(2, "hide");
    };

    self.show_right_panel = () => {
      update_panel_widths(2, "show");
    };

    self.toggle_right_panel = () => {
      update_panel_widths(2);
    };

    self.toggle_left_panel = () => {
      update_panel_widths(0);
    };

    self.change_unit_names = (unit_names) => {
      self.refs.leftpanel.opts.unit_names = unit_names;
      self.refs.leftpanel.update();
    };

    self.change_summary_models = (models) => {
      self.refs.data_summary.opts.models = models;
      self.refs.data_summary.update();
    };

    self.change_report_models = (models) => {
      self.refs.full_data_report.opts.models = models;
      self.refs.full_data_report.update();
    };

    self.change_title_unit_choose = (title = "Units") => {
      self.refs.leftpanel.opts.unit_title = title;
      self.refs.leftpanel.update();
    };

    self.change_title_plot_choose = (title = "Plots") => {
      self.refs.leftpanel.opts.plot_title = title;
      self.refs.leftpanel.update();
    };

    self.show_main_panel = (selector) => {
      Array.from(self.root.querySelectorAll(".main-panel-container")).forEach(el => el.classList.add("hidden"));
      self.root.querySelector(selector).classList.remove("hidden");
    };

    self.show_filtered_plots = (unitindex, unitname, plot_types) => {
      self.show_main_panel(".plot-container");
      self.refs.plot_collection.update()
      self.refs.plot_collection.opts.plots = plot_types.map((plot_type) => ({
        index: unitindex - 1,
        type: plot_type,
        title: null,
        subtitle: `Unit: '${unitname}'`,
      }));
      self.refs.plot_collection.update()
    };

    self.show_main_container = () => {
      self.show_main_panel(".content-container");
    };

    self.show_paper_model_1 = () => {
      self.show_main_panel(".paper-model-1");
      iframe = self.root.querySelector(".paper-model-1 iframe").contentWindow.PDFViewerApplication.pdfViewer
        .scrollPageIntoView({
          pageNumber: 1,
          destArray: [null, {
            name: "XYZ"
          }, null, null, 'page-width'],
          allowNegativeOffset: true
        });
    };

    self.show_full_report = () => {
      self.show_main_panel(".report-container");
      self.refs.full_data_report.update();
    };

    self.save_full_report_as_docx = () => {
      self.show_main_panel(".report-container");
      self.refs.full_data_report.update();
      self.refs.full_data_report.save_as_docx();
    };

    self.save_full_report_as_pdf = () => {
      self.show_main_panel(".report-container");
      self.refs.full_data_report.update();
      console.error("Not implemented!");
    };

    self.on("mount", () => {
      create_splitter([15, 70, 15]);
      self.close_left_panel();
    });

    self.show_data_configuration = (data_source=null, scroll_into=null, current_dataset_info=null) => {
      self.show_main_panel(".date-configuration");
      if(!!data_source){
        self.refs.date_configuration.opts.data_source = data_source;
      }
      if(!!current_dataset_info){
        self.refs.date_configuration.update_configuration(
          current_dataset_info.change_point,
          current_dataset_info.data_range);
      }else{
        self.refs.date_configuration.update();
      }
      if(scroll_into == 1){
        self.refs.date_configuration.scroll_into_data_range_settings();
      }else if(scroll_into == 2){
        self.refs.date_configuration.scroll_into_theoretical_settings();
      }
    };
});
