input_path = getDirectory("Choose image folder"); 
fileList = getFileList(input_path); 
setBatchMode(true);
  for (f=0; f<fileList.length; f++) {
	image = input_path+fileList[f];
	if (startsWith(fileList[f], "Segmented")) {
		if (startsWith(fileList[f], "Segmented_G")) {
		open(image);
		//run("Invert");
		getDimensions(width, height, channels, slices, frames);
		makeRectangle(0, 0, width, height-0.078*height);
		run("Crop");
		run("Analyze Particles...", "size=50-Infinity summarize");
		run("Close All");
		} else {
		open(image);
		//run("Invert");
		getDimensions(width, height, channels, slices, frames);
		makeRectangle(0, 0, width, height-0.078*height);
		run("Crop");
		run("Analyze Particles...", "size=6-Infinity summarize");
		run("Close All");
		}
	}
  }
setBatchMode(false);