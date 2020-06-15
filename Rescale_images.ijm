dirsa = getDirectory("choose folder please");	// The first argument is the directory

dir=dirsa+"zproject\\nuclei_DNA"+ File.separator;
filename = getFileList(dir); // get the list of files
newDirm = dirsa + "Montage" + File.separator; 
// make sub directory for the first channel
newDir = dirsa + "nuclei_DNA_resizzed" + File.separator; 
File.makeDirectory(newDir); 

for (i=0; i<filename.length; i++) { 
    path=dir+filename[i]; // set path 
	open(path); //open the image
  	run("32-bit");
	getRawStatistics(nPixels, mean, min, max, std, histogram);
	setAutoThreshold("Default dark");
	getThreshold(lower, upper);
	run( "Subtract...", "value=[lower]" );
	div=(max-lower);
	run("Divide...", "value=[div]");
	saveAs("tiff", newDir + getTitle); // save the image
    run("Close All"); //close the images
} 


run("Image Sequence...", "open=newDir sort"); // open the projected images as a stack
run("Scale Bar...", "width=5 font=14 color=White background=None location=[Lower Right] bold hide label"); // add scale bar of 50 microns
run("Make Montage...", "scale=1 border=3"); // Make a montage
saveAs("Tiff", newDirm + "Montage_nuclei_DNA_rescaled"); // save the image of the montage
