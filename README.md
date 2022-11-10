# Giant_Cell_Enhancer_Project
Files related to the 2022 Giant Cell Enhancer Paper (Plant Cell) Hong et al. 

Welcome to the repository of data and scripts related to the Hong et al. paper on the Giant Cell Enhancer. Use this file to navigate through and understand how to use the files in this GitHub.
#### FIMO Analysis ####
This folder contains several .tsv files summarizing the results from the Find Individual Motif Occurrences (FIMO) search used to find potential TCP, HD-ZIP, and Dof binding motifs on the giant cell enhancer. There are two files for the TCPs, since with the default threshold for calling a motif there were very few results, we were able to identify more putative binding motifs by lowering the threshold to 0.001 (a less stringent threshold).
The python file can be used for visualizing the results of the FIMO analysis. It basically counts how many motifs FIMO identifies at each nucleotide position on the Giant Cell Enhancer. You can then graph these counts in excel with a line graph to see how many motifs FIMO identified for each transcription factor family at each position on the enhancer. 
#### GCE Image Processing Data ####
The csv files in this folder contained the parsed image processing data that can be directly loaded into R to generate the graphs and statistics used in this study. 
GCE_Image_processing_data_RATIOs_GCthresh50_T6.csv has small:giant cell ratios (and other measurements not used for this study) for all the genotypes for which these ratios were used. 
GCE_Image_processing_data_pBRs_jawcrosses_T6_parsed.csv contains data for lines used in figures 3 and 6, of primary interest being raw counts of the number of giant cells and small cells. For these samples we chose not to use small:giant cell ratios to quantify the expression patterns because many of the lines from pBR63,65,57, and 69 did not have any giant cells expressing the giant cell enhancer, which would lead to zeros in the denominator if we tried to calculate small:giant cell ratio. 
GCE_Image_processing_data_trihelix_mutations_parsed.csv contains data needed for supplementary Figure S5, where we found that mutating the trihelix motif on the enhancer did not alter the expression pattern of the enhancer.
Intensity_data_ATML1_crosses_Giantv3_Smallv4_parsed.csv contains the intensity data for graphs in figure 6. Each intensity value is the mean intensity (in arbitrary units ranging from 0-255) across all small cells or all giant cells expressing the enhancer across an individual sepal. 
pAR_lazer_intensity.csv contains data on the laser intensity values used to image the lines in figure 4C. 
#### R Python and ImageJ Scripts ####
This folder just contains ImageJ scripts used to process the imaged and generate data that is not-parsed (to load it into R it helps to reformat data manually or using a script, such as the python script described later). 
Weka_Macro_script_v2.ijm is the primary script used to process imaging data. When you run the script it will ask you to choose a directory, then it will loop through all of the .tif files in that directory and process them with WEKA trainable segmentation to segment giant and small cells analyze particles to count them. Make sure you use maximum intensity projections of sepals saved as tifs and make sure there is no PI stain, rather only the channel with the giant cell enhancer signal. You will also need a WEKA classifier to perform the segmentation. If you are running the dataset from this paper you will only need to download the classifier included in a different, zipped folder in this repository and update the file path to whatever path you unzipped the classifier to. If you are attempting to use this script for a new dataset you generated, I highly recommend you use WEKA trainable segmentation to generate a new classifier trained to your dataset. Make sure before you run the program to create a new folder called “Results” in the folder containing the images you are analyzing – this is where the output from WEKA segmentation will get saved. This script includes a section where it crops off the bottom 7.8% of the segmented image because in my dataset many images had a scale bar at the bottom that would get segmented; if your images do not have scale bars you can delete or comment out that part. You may want to adjust the thresholds or counting parameters for your particular dataset. The WEKA segmentation takes quite awhile, so be sure to run this script on a fast computer (and even then a dataset of 20 images may take 20-30 minutes to process). Below is a quick summary of the way this script works:
Open directory => classify images with WEKA trainable segmentation (slow) => use thresholds isolate giant cells and small cells from classified image => count cells with Analyze Particles
Copy and paste the results from analyze particles into a text or csv file. 
Cell_counting_script.ijm Sometimes I would forget to copy and save the results from the above script or would want to try different parameters for counting Analyze Particles. This script contains just the cell counting portions of the previous script, so you don’t need to go through the very slow WEKA segmentation if you’ve already done the segmentation once before and still have the classified images saved. This script is super fast in comparison. Delete or comment out the cropping if you don’t have scale bars. Sometimes for reasons unknown ImageJ will decide to make this script not count properly (you’ll get a count of only 1 really big cell or something screwy like that), in which case de-comment the run(“invert”); line, which for some reason fixes this issue. 
Weka_cell_intensity_quantification.ijm Use this script to analyze the intensity of signal from giant cells and small cells on a sepal (used in Figure 6). It basically works the same as the first script I described, except instead of simply counting the cells it relates the segmented cell nuclei back to the original tif files to calculate signal intensity (pixel brightness) on an arbitrary scale from 0-255. It also uses a different classifier (included in a different zip folder) since I took images with different settings so as to not saturate the giant cells so we could detect differences in expression levels (whereas for simple counting we tried to saturate the signal from giant cells by using a higher laser power). 
Cell_intensity_quantification_only_v2.ijm Same as the cell counting script, except for intensity quantification where you have already saved classified images from WEKA segmentation.
#### R and Python Scripts ####
These python scripts can be used to parse the raw data output saved from the imageJ scripts in the previously described folder so that they can be easily loaded into R. However, they only work if you save your tif files with a very specific naming scheme: Genotype-ecotype-Generation#plantnumber-objective-lasersetting_Maximumintensityprojection.tif. See sample images for an example of this naming scheme. Small:giant cell ratios still need to be calculated by hand in excel (or code added to the script to calculate this). You can also manually parse your data if you want. 
The R scripts can be used to recreate the graphs found in this paper. The name of the data file you need to manually load into each script can be found commented out at the top of each script. The data files can all be found in the GCE Image Processing Data folder. The scripts can then pretty much just be run as is.
#### Cell Counting Weka Classifiers.zip ####
Classifiers and training data files needed for WEKA segmentation of datasets where we saturated signal intensity in giant cells (for counting giant and small cells). 
#### Intensity Analysis Weka Classifiers.zip ####
Classifiers and training data files needed for WEKA segmentation of datasets where we did not saturate signal intensity in giant cells (so we can compare the signal intensity of giant cells from different lines/genotypes).
#### Sample Images ####
In this folder are several sample images from selected genotypes. The genotype is the first part of each file name, including everything before the first “-“. You can use these to test the scripts and classifiers used in this study. The pAR111xCol0 image is meant to be processed with the scripts and classifier for intensity quantification, while the rest of the images are meant to be processed with the scripts and classifier for cell counting. 
