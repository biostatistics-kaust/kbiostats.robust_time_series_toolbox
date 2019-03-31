
riot.tag2('sortable-file-list-input', '<div class="trigger ui labeled icon button"> <i class="right arrow icon"></i> {opts.label} </div> <div ref="modal" class="ui special modal"> <div class="header"> {opts.modal_title} </div> <div class="content scrolling content"> <div class="ui pointing below label"> {opts.modal_description} </div> <div> <csv-file-input ref="file_list_input" class="file_list_input" placeholder="Clust_Mat_8Clust_HCC01.txt" allow_multi_selections="{true}" load_content="{false}" file_filter="{opts.file_filter}" dialog_properties="{[\'openFile\',                     \'multiSelections\' ]}"></csv-file-input> </div> <virtual if="{opts.autosaved_description_valid}"> <div class="ui segment"> {opts.autosaved_description} </div> </virtual> <div class="file-input-list list-group col ui list"> <div class="list-group-item item" each="{filedata, index in filelist}"> <div class="file-buttons left floated content"> <div class="ui button" style="padding-right: 10px"><i class="file icon"></i></div> <div class="ui button remove-item" style="padding-right: 10px" onclick="{removeIndex}"><i class="delete icon"></i></div> </div> <div class="content"> <div class="header">[{index}]{filedata.replace(/^([^\\/\\\\]*[\\/\\\\])*/g, \'\')}</div> <div class="description">{filedata}</div> </div> </div> </div> </div> </div>', '', '', function(opts) {

        function default_file_filter() {
            return [{
                    name: 'Data source (.csv, .txt)',
                    extensions: ['csv', 'txt']
                },
                {
                    name: 'All Files',
                    extensions: ['*']
                }
            ];
        }

        const self = this;
        const self$ = jQuery(this.root);
        const config = opts;
        config.file_filter = config.file_filter !== undefined ? config.file_filter : default_file_filter();
        config.label = config.label !== undefined ? config.label : "Choose files";
        config.modal_title = config.modal_title !== undefined ? config.modal_title : "Files";
        config.modal_description = config.modal_description !== undefined ? config.modal_description :
            "Choose files";
        config.autosaved_description = config.autosaved_description !== undefined ? config.autosaved_description :
            "Changes are automatically saved";
        config.autosaved_description_valid = !(config.autosaved_description == null || config.autosaved_description.trim() ==
            "");

        self.filelist = []
        self.removeIndex = (event) => {
            console.log(event)
            self.filelist.splice(event.item.index, 1);
            self.update();
        }

        self.on("mount", () => {
            jQuery(self.refs.modal).modal({context: self$});
            self$.find(".trigger").on("click", () => {
                jQuery(self.refs.modal).modal("show");
            });
            self.refs.file_list_input.on("obtained-file-names", (data) => {
                self.filelist = self.filelist.concat(data);
                self.update();
            });
            new Sortable(self$.find(".file-input-list")[0], {
                animation: 150,
                onEnd: function (evt) {
                    self.filelist.splice(evt.newIndex, 0, self.filelist.splice(evt.oldIndex, 1)[0]);
                    node_changed = this.el.removeChild(this.el.children[evt.newIndex]);
                    this.el.insertBefore(node_changed, this.el.children[evt.oldIndex])
                    self.update();
                },
            });
        });
});