'use strict';
const {
  app,
  BrowserWindow
} = require('electron');
const electron = require('electron');
///////////////////////////////////////////////////////////////////////////////
// Helpers
///////////////////////////////////////////////////////////////////////////////
const OS = () => {
  if (process.platform.indexOf("darwin") >= 0 || process.platform.indexOf("mac") >= 0) {
    return "OSX";
  } else if (process.platform.indexOf("win") >= 0) {
    return "WIN";
  } else {
    return "";
  }
}

const getIconPath = () => {
  //let mainWindowIcon = __dirname.toString() + '/resources/images/app/icon/icon_256x256';
  var path = require('path')
  let mainWindowIcon = path.join(__dirname.toString(), '/resources/images/app/icon/icon_256x256');
  if (OS() == "OSX") {
    mainWindowIcon = mainWindowIcon + ".icns";
  } else if (OS() == "WIN") {
    mainWindowIcon = mainWindowIcon + ".ico";
  } else {
    mainWindowIcon = mainWindowIcon + ".png";
  }
  mainWindowIcon = path.join(__dirname.toString(), '/resources/images/app/app.png');
  return mainWindowIcon;
}

///////////////////////////////////////////////////////////////////////////////
// Settings
///////////////////////////////////////////////////////////////////////////////
let mainWindowSettings = {
  title: 'Robust Time Series Toolbox - KAUST Biostatistics Group',
  width: 1024,
  height: 760,
  minWidth: 400,
  minHeight: 300,
  //backgroundColor: 'rgb(224, 224, 224)',
  backgroundColor: 'rgb(15,52,86)',
  icon: getIconPath(),
  acceptFirstMouse: true,
  //titleBarStyle: 'customButtonsOnHover',
  titleBarStyle: 'hiddenInset',
  //titleBarStyle: 'hidden',
  //autoHideMenuBar: !true,
  fullscreenWindowTitle: true,
  frame: false,
  allowRendererProcessReuse: true,
  webPreferences: {
/////    webSecurity: true,
    enableRemoteModule: true,
    nodeIntegration: true,
    //nodeIntegrationInWorker: true,
    //contextIsolation: false,
  },
  nodeIntegration: true,
}
let modalWindowSettings = {
  title: 'Robust Time Series Toolbox - KAUST Biostatistics Group',
  width: 600,
  height: 450,
  minWidth: 400,
  minHeight: 300,
  //backgroundColor: 'rgb(224, 224, 224)',
  backgroundColor: 'rgb(15,52,86)',
  icon: getIconPath(),
  acceptFirstMouse: true,
  //titleBarStyle: 'customButtonsOnHover',
  titleBarStyle: 'hiddenInset',
  //titleBarStyle: 'hidden',
  //autoHideMenuBar: !true,
  fullscreenWindowTitle: true,
  frame: false,
  modal: true,
  //allowRendererProcessReuse: true,
  webPreferences: {
/////    webSecurity: true,
    //nodeIntegration: true,
    //nodeIntegrationInWorker: true,
    //contextIsolation: false,
  },
}


///////////////////////////////////////////////////////////////////////////////
// Main events
///////////////////////////////////////////////////////////////////////////////
let mainWindow, modalWindow;

app.on('window-all-closed', function () {
  app.quit();
  //if(process.platform != 'darwin') {
  //  app.quit();
  //}
});

app.on('ready', function () {
  if (OS() == "OSX") {
    mainWindowSettings.frame = true;
  }
  mainWindow = new BrowserWindow(mainWindowSettings);
  modalWindowSettings.parent = mainWindow
  modalWindow = new BrowserWindow(modalWindowSettings);
  //mainWindow.setMenu(electron.Menu.buildFromTemplate(menuTemplate));
  mainWindow.setMenu(null);
  modalWindow.setMenu(null);
  mainWindow.loadURL('file://' + __dirname + '/resources/main.html');
  modalWindow.loadURL('file://' + __dirname + '/resources/modal.html');
  mainWindow.on('closed', function () {
    mainWindow = null;
  });
  modalWindow.on('closed', function () {
    
  });
  //mainWindow.maximize();
  //
  //mainWindow.toggleDevTools();
  //modalWindow.toggleDevTools();

  
});

app.on('error', function () {
  console.log("Error starting the application")
});