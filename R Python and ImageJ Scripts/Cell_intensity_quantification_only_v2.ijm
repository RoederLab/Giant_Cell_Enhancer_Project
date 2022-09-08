//batch processing
input_path = getDirectory("Choose image folder"); 
fileList = getFileList(input_path); 
  for (f=0; f<fileList.length; f++) {
	original_image_path = substring(input_path,0,lastIndexOf(input_path, "Results")) + substring(fileList[f], 21,lastIndexOf(fileList[f], "f")+1);
	original_image = substring(fileList[f], 21,lastIndexOf(fileList[f], "f")+1);
	image = input_path+fileList[f];
		if (startsWith(fileList[f], "Segmented_G")) {
		setBatchMode(true);
		open(original_image_path);
		run("Split Channels");
		selectWindow(original_image + " (blue)");
		close();
		selectWindow(original_image + " (red)");
		close();
		selectWindow(original_image + " (green)");
		rename(original_image);
		run("Set Measurements...", "area mean perimeter integrated median display redirect=" + original_image + " decimal=3");
		open(image);
		run("Analyze Particles...", "size=50-Infinity circularity=0.50-1.00 display summarize");
		run("Close All");
	}
	if (startsWith(fileList[f], "Segmented_S")) {
		setBatchMode(true);
		
		open(original_image_path);
		run("Split Channels");
		selectWindow(original_image + " (blue)");
		close();
		selectWindow(original_image + " (red)");
		close();
		selectWindow(original_image + " (green)");
		rename(original_image);
		run("Set Measurements...", "area mean perimeter integrated median display redirect=" + original_image + " decimal=3");
		open(image);
		run("Analyze Particles...", "size=9-Infinity circularity=0.50-1.00 display summarize");
        run("Close All");
	}
  }

