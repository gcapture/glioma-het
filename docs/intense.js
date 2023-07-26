// Based on the original intense.js library: https://github.com/tholman/intense-images/blob/master/intense.js
// Code for the image zoom based on: https://www.w3schools.com/howto/howto_js_image_zoom.asp
// Heavy thanks to Bhuvic Patel for the help setting this up: https://scholar.google.com/citations?user=W3Ci8_8AAAAJ&hl=en

window.requestAnimFrame = (function() {
  return (
    window.requestAnimationFrame ||
      window.webkitRequestAnimationFrame ||
      window.mozRequestAnimationFrame ||
      function(callback) {
        window.setTimeout(callback, 1000 / 60);
      }
  );
})();

window.cancelRequestAnimFrame = (function() {
  return (
    window.cancelAnimationFrame ||
      window.webkitCancelRequestAnimationFrame ||
      window.mozCancelRequestAnimationFrame ||
      window.oCancelRequestAnimationFrame ||
      window.msCancelRequestAnimationFrame ||
      clearTimeout
  );
})();

// Wrapper for the Insense figure.
var Intense = (function() {
  "use strict";

  // Define global variables.
  var KEYCODE_ESC = 27;  // Escape key keycode.

  var invertInteractionDirection = false;

  // Holds the animation frame id.
  var looper;

  // Single image
  var image;

  var closebutton;      // Will store the close button.
  var downloadbutton;   // Will store the download button.
  var magnifyingbutton; // Will store the magnifying button.
  var buttondiv;        // Will store the buttons in a div.

  var close_icon = document.createElement("i");
  var download_icon = document.createElement("i");
  var magnify_icon = document.createElement("i");

  var sourceDimensions, target;
  var targetDimensions = { w: 0, h: 0 };

  var container; // Will store the duplicated intense image.
  var containerDimensions = { w: 0, h: 0 };

  // Overflow variable before screen is locked.
  var overflowValue;

  var active = false;

  var magnifyingGlass;  // Will store the magnifying glass.
  var lens;             // Will store the lens.

  var zoomFactor = 0.5; // Stores the zoom factor.

  var cx, cy;

  /* -------------------------
    /*          UTILS
  /* -------------------------*/

    // Soft object augmentation
  function extend(target, source) {
    for (var key in source) if (!(key in target)) target[key] = source[key];

    return target;
  }

  // Applys a dict of css properties to an element
  function applyProperties(target, properties) {
    for (var key in properties) {
      target.style[key] = properties[key];
    }
  }

  // Returns whether target a vertical or horizontal fit in the page.
  // As well as the right fitting width/height of the image.
  function getFit(source) {
    if (source.w/source.h > 0.98*window.innerWidth/(0.98*window.innerHeight)) {
      var widthRatio = 0.98*window.innerWidth / source.w;
      return {
        w: source.w * widthRatio,
        h: source.h * widthRatio
      };
    } else {
      var heightRatio = 0.98*window.innerHeight / source.h;
      return { w: source.w * heightRatio, h: source.h * heightRatio};
    }
  }

  // Gets the fit of the magnifying glass when it is activated.
  function getFitMag(source) {
    if (source.w/source.h > 0.98*window.innerWidth/(0.65*window.innerHeight)) {
      var widthRatio = 0.98*window.innerWidth / source.w;
      return {
        w: source.w * widthRatio,
        h: source.h * widthRatio
      };
    } else {
      var heightRatio = 0.65*window.innerHeight / source.h;
      return { w: source.w * heightRatio, h: source.h * heightRatio};
    }
  }

  // Define function to move the lens.
  function moveLens(e) {
    var pos, x, y, a;
    // Prevent any other actions that may occur when moving over the image
    e.preventDefault();
    // Get the cursor's x and y positions:
    pos = getCursorPos(e);
    // Calculate the position of the lens:
    x = pos.x - (lens.offsetWidth / 2);
    y = pos.y - (lens.offsetHeight / 2);
    // Get the image position
    a = {left: (1*window.innerWidth - target.width)/2,
         top: (0.67*window.innerHeight - target.height)/2
    };
    // Prevent the lens from being positioned outside the image:
    if (x > target.width - lens.offsetWidth) {
      x = target.width - lens.offsetWidth;
    }
    if (x < 0) {
      x = 0;
    }
    if (y > target.height - lens.offsetHeight) {
      y = target.height - lens.offsetHeight;
    }
    if (y < 0) {
      y = 0;
    }

    // Set the position of the lens:
    lens.style.left = x + a.left + "px";
    lens.style.top = y + a.top + "px";
    // Display what the lens "sees":
    magnifyingGlass.style.backgroundPosition = "-" + (x * cx) + "px -" + (y * cy) + "px";
  }

  // Gets the cursor position.
  function getCursorPos(e) {
    var a, x = 0, y = 0;
    e = e || window.event;
    // Get the x and y positions of the image:
    a = {left: (1*window.innerWidth - target.width)/2,
         top: (0.67*window.innerHeight - target.height)/2
    };
    // Calculate the cursor's x and y coordinates, relative to the image:
      x = e.pageX - a.left;
    y = e.pageY - a.top;
    // Consider any page scrolling:
      x = x - window.pageXOffset;
    y = y - window.pageYOffset;
    return {x : x, y : y};
  }


  // Performs the image zoom for the magnifying glass.
  function imageZoom() {

    cx = magnifyingGlass.offsetWidth / lens.offsetWidth;
    cy = magnifyingGlass.offsetHeight / lens.offsetHeight;

    magnifyingGlass.style.backgroundImage = "url('" + target.src + "')";
    magnifyingGlass.style.backgroundSize = (target.width * cx) + "px " + (target.height * cy) + "px";
    magnifyingGlass.style.backgroundRepeat = "no-repeat";

    // Execute a function when someone moves the cursor over the image, or the lens:
      lens.addEventListener("mousemove", moveLens);
    target.addEventListener("mousemove", moveLens);

    // And also for touch screens:
      lens.addEventListener("touchmove", moveLens);
    target.addEventListener("touchmove", moveLens);
  }




  /* -------------------------
    /*          APP
  /* -------------------------*/

    // Start tracking.
  function startTracking(passedElements) {
    var i;
    // If passed an array of elements, assign tracking to all.
    if (passedElements.length) {
      // Loop and assign
      for (i = 0; i < passedElements.length; i++) {
        track(passedElements[i]);
      }
    } else {
      track(passedElements);
    }
  }

  // Track function.
  function track(element) {
    // Element needs a src at minumun.
    if (element.getAttribute("data-image") || element.src || element.href) {
      element.addEventListener(
        "click",
        function(e) {
          if (element.tagName === "A") {
            e.preventDefault();
          }
          if (!active) {
            init(this);
          }
        },
        false
      );
    }
  }

  // Start function.
  function start() {
    loop();
  }

  // Stop function.
  function stop() {
    cancelRequestAnimFrame(looper);
  }

  // Loop function.
  function loop() {
    looper = requestAnimFrame(loop);
  }

  // Lock scroll on the document body.
  function lockBody() {
    overflowValue = document.body.style.overflow;
    document.body.style.overflow = "hidden";
  }

  // Unlock scroll on the document body.
  function unlockBody() {
    document.body.style.overflow = overflowValue;
  }

  // Set state function.
  function setState(element, newClassName) {
    if (element) {
      element.className = element.className.replace("intense--loading", "");
      element.className = element.className.replace("intense--viewing", "");
      element.className += " " + newClassName;
    } else {
      // Remove element with class .view
      var elems = document.querySelectorAll(".intense--viewing");
      [].forEach.call(elems, function(el) {
        el.className = el.className.replace("intense--viewing", "").trim();
      });
    }
  }

  // Function to download an image from source.
  // From: https://stackoverflow.com/a/51076777

  async function imageDownload(source){
    // Transforms the src into an URL object.
    function toDataURL(url) {
      return fetch(url).then((response) => {return response.blob();}).then(blob => {return URL.createObjectURL(blob);});
    }

    // Makes the URL into a download link and clicks it.
    async function download(source) {
      const a = document.createElement("a");
      a.href = await toDataURL(source);
      a.download = "image.png";
      document.body.appendChild(a);
      a.click();
      document.body.removeChild(a);
    }

    // Call to the function.
    download(source)
  }

  // Creates a viewer.
  function createViewer() {
    // Container properties.
    var containerProperties = {
      backgroundColor: "rgba(0,0,0, 0.90)",
      width: "100%",
      height: "100%",
      position: "fixed",
      top: "0px",
      left: "0px",
      overflow: "auto",
      zIndex: "999999",
      margin: "0px",
      webkitTransition: "opacity 150ms cubic-bezier( 0, 0, .26, 1 )",
      MozTransition: "opacity 150ms cubic-bezier( 0, 0, .26, 1 )",
      transition: "opacity 150ms cubic-bezier( 0, 0, .26, 1 )",
      webkitBackfaceVisibility: "hidden",
      opacity: "0"
    };
    container = document.createElement("figure");
    target.id = "intense-img"
    container.appendChild(target);
    container.className = "intense-zoom";
    applyProperties(container, containerProperties);

    // Change cursor over image.
    var imageProperties = {
      cursor: 'zoom-out'
    };
    applyProperties(target, imageProperties);


    // Make magnifying button
    var magnifyingbutton;
    magnifyingbutton = document.createElement("button");
    magnifyingbutton.setAttribute('id', 'intense-magnifying');
    magnifyingbutton.setAttribute('class', 'btn btn-primary');
    magnifyingbutton.style["margin-right"] = "0px";
    magnifyingbutton.style["background-color"] = "#ffffff00";
    magnifyingbutton.style["border"] = "0px solid white";
    magnifyingbutton.onclick = function(){createViewerMag()};
    magnify_icon.setAttribute('class', "bi bi-zoom-in");
    magnify_icon.style["color"] = "black";
    magnifyingbutton.appendChild(magnify_icon);

    // Make download button
    downloadbutton = document.createElement("button");
    downloadbutton.setAttribute('id', 'intense-download');
    downloadbutton.setAttribute('class', 'btn btn-primary');
    downloadbutton.style["margin-right"] = "0px";
    downloadbutton.style["background-color"] = "#ffffff00";
    downloadbutton.style["border"] = "0px solid white";
    downloadbutton.onclick = function(){imageDownload(target.src)};
    download_icon.setAttribute('class', "bi bi-download");
    download_icon.style["color"] = "black";
    downloadbutton.appendChild(download_icon);

    // Make close button
    closebutton = document.createElement("button");
    closebutton.setAttribute('id', 'intense-close');
    closebutton.setAttribute('class', 'btn btn-primary');
    closebutton.style["margin-right"] = "0px";
    closebutton.style["background-color"] = "#ffffff00";
    closebutton.style["border"] = "0px solid white";
    closebutton.onclick = function(){removeViewer()};
    close_icon.setAttribute('class', "bi bi-x-lg");
    close_icon.style["color"] = "black";
    closebutton.appendChild(close_icon);

    // Make a div to hold the buttons
    buttondiv = document.createElement("div");
    buttondiv.style['position'] = 'absolute';
    buttondiv.style['width'] = 'fit-content';

    // Add the buttons
    buttondiv.appendChild(magnifyingbutton);
    buttondiv.appendChild(downloadbutton);
    buttondiv.appendChild(closebutton);
    container.appendChild(buttondiv);

    // Set the dimensions of the viewer.
    setDimensions();

    // Add the viewer to the document.
    document.body.appendChild(container);
    setTimeout(function() {
      container.style["opacity"] = "1";
    }, 10);
  }

  // Creates the viewport when magnifying glass is activated.
  function createViewerMag(){
    // Remove all listeners.
    unbindEvents();

    // Add the new ones suited for magnifying glass.
    bindEventsMag();

    // Assign the close function to the close button again.
    closebutton.onclick = function(){removeViewerMag()};

    close_icon.style["color"] = "white";
    download_icon.style["color"] = "white";


    // Remove the magnifying glass button once this is activated to avoid redundancies.
    document.getElementById("intense-magnifying").remove();

    //Make magnifying glass
    var magnifyingGlassProperties = {
      border: "1px solid #d4d4d4",
    };
    magnifyingGlass = document.createElement("div");
    magnifyingGlass.id ="intense-magnifying-glass";
    magnifyingGlass.className = "intense-magnifying-glass";
    applyProperties(magnifyingGlass, magnifyingGlassProperties);
    container.appendChild(magnifyingGlass);

    // Create lens:
      lens = document.createElement("div");
    lens.style['position'] = "absolute";
    lens.style['width'] = "fit-content";
    //lens.style['border'] = "3px green solid";
    container.appendChild(lens);

    // When magnifying glass is set, change the cursor to a cross.
    var cursorProperties = {
      cursor: 'crosshair'
    };
    applyProperties(lens, cursorProperties)

    // Set the dimensions of the image when magnifying glass is activated.
    setDimensionsMag();

    // Add the container to the document.
    document.body.appendChild(container);
    setTimeout(function() {
      container.style["opacity"] = "1";
    }, 10);
  }


  // Remove the viewers.
  function removeViewer() {
    unlockBody();
    unbindEvents();
    stop();
    document.body.removeChild(container);
    active = false;
    setState(false);
  }

  // Remove the viewers when magnified glass is activated.
  function removeViewerMag() {
    zoomFactor = 0.5;
    unlockBody();
    unbindEventsMag();
    stop();
    document.body.removeChild(container);
    active = false;
    setState(false);
  }


  // Set the dimensions.
  function setDimensions() {
    var imageDimensions = getFit(sourceDimensions);
    target.width = imageDimensions.w;
    target.height = imageDimensions.h;

    targetDimensions = { w: target.width, h: target.height };
    containerDimensions = { w: window.innerWidth, h: window.innerHeight };
    target.style['position'] = "absolute";

    // Center the image location vertically and horizontally
    target.style['left'] = (window.innerWidth - target.width)/2 + "px";
    target.style['top'] = (window.innerHeight - target.height)/2 + "px";

    // Use the image location to set position of buttons
    buttondiv.style['right'] = (window.innerWidth - target.width)/2 + "px";
    buttondiv.style['top'] = (window.innerHeight - target.height)/2 + "px";
  }

  // Sets the dimensions of the image when the magnifying glass is activated.
  function setDimensionsMag(){
    var imageDimensions = getFitMag(sourceDimensions);

    target.width = imageDimensions.w;
    target.height = imageDimensions.h;

    targetDimensions = { w: target.width, h: target.height };
    containerDimensions = { w: window.innerWidth, h: window.innerHeight };

    magnifyingGlass.style['min-width'] = String(window.innerHeight*0.32).concat("px");
    magnifyingGlass.style['min-height'] = String(window.innerHeight*0.32).concat("px");
    magnifyingGlass.style['position'] = 'absolute';
    magnifyingGlass.style['background-repeat'] = 'no-repeat';
    magnifyingGlass.style['top'] = "67%";
    magnifyingGlass.style['left'] = String(window.innerWidth/2 - window.innerHeight*0.32/2).concat("px");

    lens.style['width'] = String(window.innerHeight*0.32*zoomFactor).concat("px");
    lens.style['height'] = String(window.innerHeight*0.32*zoomFactor).concat("px");

    target.style['position'] = "absolute";
    target.style['left'] = (window.innerWidth - target.width)/2 + "px";
    target.style['top'] = (0.67*window.innerHeight - target.height)/2 + "px";

    // Get the image location to set position of buttons
    buttondiv.style['right'] = (1.0*window.innerWidth - target.width)/2 + "px";
    buttondiv.style['top'] = "67%";

    imageZoom();
  }


  function init(element) {
    setState(element, "intense--loading");
    var imageSource =
      element.getAttribute("data-image") || element.src || element.href;

    // Clear old onload message
    if (image) {
      image.onload = null;
    }

    image = new Image();
    image.onload = function() {
      sourceDimensions = { w: image.width, h: image.height }; // Save original dimensions for later.
      target = this;
      createViewer();
      lockBody();
      bindEvents();
      loop();

      setState(element, "intense--viewing");
    };

    image.src = imageSource;
  }


  function bindEvents() {
    window.addEventListener("resize", setDimensions, false);
    window.addEventListener("keyup", onKeyUp, false);
    target.addEventListener("click", removeViewer, false);
  }

  function unbindEvents() {
    window.removeEventListener("resize", setDimensions, false);
    window.removeEventListener("keyup", onKeyUp, false);
    target.removeEventListener("click", removeViewer, false);
  }


  function bindEventsMag() {
    window.addEventListener("resize", setDimensionsMag, false);
    window.addEventListener("keyup", onKeyUp, false);
    target.addEventListener("click", removeViewerMag, false);
    window.addEventListener("wheel", onScroll, {passive:false});
    //window.addEventListener("wheel", (event) => {onScroll(event)}, false);
  }

  function unbindEventsMag() {
    window.removeEventListener("resize", setDimensionsMag, false);
    window.removeEventListener("keyup", onKeyUp, false);
    target.removeEventListener("click", removeViewerMag, false);
    window.removeEventListener("wheel", onScroll, {passive:false});
  }


  // Exit on excape key pressed;
  function onKeyUp(event) {
    event.preventDefault();
    if (event.keyCode === KEYCODE_ESC) {
      removeViewer();
    }
  }

  function onScroll(event) {
    event.preventDefault();
    if (event.deltaY < 0){
      zoomFactor -= 0.05
    } else {
      zoomFactor += 0.05
    }

    if (window.innerHeight*0.32*zoomFactor >= target.height){
      zoomFactor = target.height/(window.innerHeight*0.32);
    }

    if (window.innerHeight*0.32*zoomFactor >= target.width){
      zoomFactor = target.width/(window.innerHeight*0.32);
    }

    if (zoomFactor < 0.05){
      zoomFactor = 0.05
    }

    setDimensionsMag();
    moveLens(event);
  }

  function config(options) {
    if ("invertInteractionDirection" in options)
      invertInteractionDirection = options.invertInteractionDirection;
  }

  function main(element, configOptions) {
    // Parse arguments
    if (!element) {
      throw "You need to pass an element!";
    }

    // If they have a config, use it!
      if (configOptions) {
        config(configOptions);
      }

    startTracking(element);
  }

  return extend(main, {
    resize: setDimensions,
    start: start,
    stop: stop,
    config: config
  });
})();

if (typeof module !== "undefined" && module.exports) {
  module.exports = Intense;
}
