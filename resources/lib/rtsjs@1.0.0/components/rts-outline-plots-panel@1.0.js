riot.tag2('rts-outline-plots-panel', '<div class="panel-title"> {opts.title} </div> <ul class="sidenav-counter sidenav-counter-expand-fs h-100 w-100"> <li><a href="#" onclick="{_trigger_event(\'RawTimeSeries\');}"> <span class="icon"> <div class="sf-icon--icon96 sf-icon--changepoint-raw-chart"></div> </span> <span class="title">Time series</span> <span class="counter">Original</span> </a></li> <li><a href="#" onclick="{_trigger_event(\'EstimatedTimeSeries\');}"> <span class="icon"> <div class="sf-icon--icon96 sf-icon--changepoint-chart"></div> </span> <span class="title">Time series</span> <span class="counter">Estimated change point</span> </a></li> <li><a href="#" onclick="{_trigger_event(\'Loglikelihood\');}"> <span class="icon"> <div class="sf-icon--icon96 sf-icon--changepoint-chart"></div> </span> <span class="title">Model</span> <span class="counter">Loglikelihood plot</span> </a></li> <li><a href="#" onclick="{_trigger_event(\'ResidualsBeforeChangePoint\');}"> <span class="icon"> <div class="sf-icon--icon96 sf-icon--changepoint-pre-acf"></div> </span> <span class="title">Model</span> <span class="counter">ACF/residuals before change point</span> </a></li> <li><a href="#" onclick="{_trigger_event(\'ResidualsAfterChangePoint\');}"> <span class="icon"> <div class="sf-icon--icon96 sf-icon--changepoint-post-acf"></div> </span> <span class="title">Model</span> <span class="counter">ACF/residuals after change point</span> </a></li> </ul>', 'rts-outline-plots-panel .icon,[data-is="rts-outline-plots-panel"] .icon{ zoom: 0.4; margin-top: -8px; } rts-outline-plots-panel *,[data-is="rts-outline-plots-panel"] *{ -webkit-user-drag: none !important; user-select: none; cursor: default; }', '', function(opts) {


        const self = this;
        const config = opts;
        config.title = !!config.title ? config.title : "";

        self._trigger_event = (name) => {

            return () => self.trigger('app:request:current:plot:allunits', {
                plottype: name.toLowerCase()
            });
        };

        self.on("mount", () => {

        });
});