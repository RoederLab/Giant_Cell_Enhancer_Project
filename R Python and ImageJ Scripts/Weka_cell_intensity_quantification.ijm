//batch processing
input_path = getDirectory("Choose image folder"); 
fileList = getFileList(input_path); 
  for (f=0; f<fileList.length; f++) {
	image = input_path+fileList[f];
	if (endsWith(image, ".tif")) {
		//open files in Weka and get title
		setBatchMode(true);
		open(image);
        title = getTitle();
        close(image);
        setBatchMode(false);
        run("Trainable Weka Segmentation", "open="+image);
        
        //Wait for Weka to load. Code breaks if you don't wait.
		wait(3000);
		
		//Load training data and run Weka
		call("trainableSegmentation.Weka_Segmentation.loadClassifier", "/mnt/roedernas/Lilan/Confocal/Confocal-2018/4-23-2018,pAR261/pAR261-T1/WEKA_Files/classifier_GCE_nuclei_T6.model");
		call("trainableSegmentation.Weka_Segmentation.getResult");

		//Start batch mode. FYI Weka doesn't work with batch mode on.
		setBatchMode(true);

		//Reopen original image and isolate green pixel values
		open(image);
		run("Split Channels");
		selectWindow(title + " (blue)");
		close();
		selectWindow(title + " (red)");
		close();
		selectWindow(title + " (green)");
		rename(title)
		//Crop out scale bar
		getDimensions(width, height, channels, slices, frames);
		makeRectangle(0, 0, width, height-0.05*height);
		run("Crop");

		//Split the classified image into three RGB images
		selectImage("Classified image");
		saveAs("tiff", input_path + "/Results/" +  "Classified_" + title);
		run("RGB Color");
		run("Split Channels");
		
		//close the green and blue channels. We just want the red channel
		selectWindow("Classified_" + title + " (green)");
		close();
		selectWindow("Classified_" + title + " (blue)");
		close();
		selectWindow("Classified_" + title + " (red)");
		
		//make a copy of the red channel
		run("Duplicate...", " ");
		
		//Threshold out the giant cell nuclei, crop out scale bar
		selectWindow("Classified_" + title + " (red)");
		getDimensions(width, height, channels, slices, frames);
		makeRectangle(0, 0, width, height-0.05*height);
		run("Crop");
		rename("Segmented_GiantCells_" + title);
		run("Threshold...");
		setThreshold(70, 89);
		setOption("BlackBackground", false);
		run("Convert to Mask");
		saveAs("tiff", input_path + "/Results/" +  "Segmented_GiantCells_" + title);
		selectWindow("Segmented_GiantCells_" + title);
		run("Close");
		open(input_path + "/Results/" +  "Segmented_GiantCells_" + title);
		run("Set Measurements...", "area mean integrated display redirect=" + title + " decimal=3");
		run("Analyze Particles...", "size=50-Infinity summarize");
		
		//Threshold out the small cell nuclei, crop out scale bar
		selectWindow("Classified_" + title + " (red)-1");
		getDimensions(width, height, channels, slices, frames);
		makeRectangle(0, 0, width, height-0.05*height);
		run("Crop");
		rename("Segmented_SmallCells_" + title);
		run("Threshold...");
		setThreshold(194, 202);
		setOption("BlackBackground", false);
		run("Convert to Mask");
		saveAs("tiff", input_path + "/Results/" +  "Segmented_SmallCells_" + title);
		selectWindow("Segmented_SmallCells_" + title);
		run("Close");
		open(input_path + "/Results/" +  "Segmented_SmallCells_" + title);
		run("Set Measurements...", "area mean integrated display redirect=" + title + " decimal=3");
		run("Analyze Particles...", "size=6-Infinity summarize");
        //run("Set Scale...", "distance=1 known=1.18609 unit=um");
        run("Close All");
        setBatchMode(false);
	}
  }

