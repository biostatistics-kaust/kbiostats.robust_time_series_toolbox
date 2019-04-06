riot.tag2('rts-modal-window', '<window-decorator ref="window"> <yield to="title"> <span class="btn-custom" role="button"><img src="./images/app/app.png" style="width: 23px;height: 23px;position: relative;top: 7px;display: block;"></span> <span class="label">Robust Time Series Toolbox - KAUST Biostatistics Group</span> </yield> <yield to="left-title"> </yield> <yield to="right-title"> </yield> <yield to="menu-bar"> </yield> <yield to="content"> <div class="container h-100 w-100"> <div class="w-25 h-100"> <span class="big-icon-container"> <span class="mif-security mif-5x"></span> </span> </div> <div class="w-75 h-100"> <h1 class="modal-title"> <virtual if="{opts.title}"> {opts.title} </virtual> <virtual if="{!opts.title}"> Disclaimer </virtual> </h1> <div class="modal-content"> <virtual if="{opts.content}"> {opts.content} </virtual> <virtual if="{!opts.content}"> <h2>Privacy policy</h2> <p>This software do not collect any data. Once installed this toolbox do not need or use any Internet connection. Therefore,your files, settings, and all other personal data remain on-device unless you explicitly share them with a third-party app or service. </p> <h2>Data policy</h2> <p>This toolbox is made freely, and open-source, available. The authors and their organizations assume no liability for the use of and results obtained from this toolbox and, in addition, do not guarantee the safety of your data. Your use of the software is at your own risk. </p> </virtual> </div> </div> </div> </yield> </window-decorator> </div>', 'rts-modal-window h1,[data-is="rts-modal-window"] h1{ font-size: 16pt; font-weight: 400; } rts-modal-window h2,[data-is="rts-modal-window"] h2{ font-size: 14pt; margin-top: 15px; } rts-modal-window .modal-content>*:nth-child(1),[data-is="rts-modal-window"] .modal-content>*:nth-child(1){ margin-top: 0 !important; padding: 0 !important; } rts-modal-window .modal-content,[data-is="rts-modal-window"] .modal-content{ background-color: rgba(255, 255, 255, 0.8); padding: 20px; height: 90%; overflow: auto; } rts-modal-window .modal-content,[data-is="rts-modal-window"] .modal-content,rts-modal-window p,[data-is="rts-modal-window"] p{ font-size: 10pt; } rts-modal-window .big-icon-container,[data-is="rts-modal-window"] .big-icon-container{ zoom: 4; margin: auto; } rts-modal-window .big-icon-container [class*=mif-],[data-is="rts-modal-window"] .big-icon-container [class*=mif-]{ background: -webkit-gradient(linear, left top, left bottom, from(transparent), to(#158bc5)); -webkit-background-clip: text; -webkit-text-fill-color: transparent; } rts-modal-window .container,[data-is="rts-modal-window"] .container{ display: flex; max-width: 100% !important; padding: 20px; }', '', function(opts) {


    const require_ = require;

    const self = this;
    const config = opts;

    self.on("mount", () => {
      jQuery(document).bind('keydown', 'ctrl+shift+r', () => {
        window.location.reload(true);
      });
      jQuery(document).bind('keydown', 'ctrl+shift+e', () => {
        require_("electron").remote.getCurrentWindow().toggleDevTools();
      });
    });
});