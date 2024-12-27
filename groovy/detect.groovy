// Set image type
setImageType('BRIGHTFIELD_H_DAB');

// set image resolution
setPixelSizeMicrons(0.82, 0.82) 

// Set color deconvolution stains
setColorDeconvolutionStains('{"Name" : "H-DAB default", "Stain 1" : "Hematoxylin", "Values 1" : "0.65111 0.70119 0.29049", "Stain 2" : "DAB", "Values 2" : "0.26917 0.56824 0.77759", "Background" : " 255 255 255"}');

// Create full image annotation
createFullImageAnnotation(true);

// Run positive cell detection
runPlugin('qupath.imagej.detect.cells.PositiveCellDetection', '{"detectionImageBrightfield":"Optical density sum","requestedPixelSizeMicrons":0.0,"backgroundRadiusMicrons":8.0,"backgroundByReconstruction":true,"medianRadiusMicrons":3.0,"sigmaMicrons":3.0,"minAreaMicrons":25.0,"maxAreaMicrons":400.0,"threshold":0.1,"maxBackground":2.0,"watershedPostProcess":true,"excludeDAB":false,"cellExpansionMicrons":5.0,"includeNuclei":true,"smoothBoundaries":true,"makeMeasurements":true,"thresholdCompartment":"Cell: DAB OD mean","thresholdPositive1":0.2,"thresholdPositive2":0.4,"thresholdPositive3":0.6000000000000001,"singleThreshold":true}')


