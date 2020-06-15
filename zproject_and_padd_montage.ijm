/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////  SCRIPT TO PROJECT ALL THE IMAGES AND MAKE A MONTANGE FOR EACH IMAGE
///////  WRITTEN BY: SARADHA VENKATACHALAPATHY                                                                                                                   
///////  ASSUMPTIONS: The input image is a confocal zstack and the channel 1 contains the nucleus. 
///////  DESCRIPTION: The script accepts the directory to a folder containing rawimages folder and the number of channels in the raw image that need to be analysed and the taget of the stains. 
///////				  It creates 2 subfolders: one for storing the z-projected images. It opens each image, does a max intensity projection, stores the image. Then open all projected images 
///////				  and create a montage with a scale bar of 50 microns
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

args = split(getArgument()," "); 	// read in the array of argurments that are separated by a space
dirsa = args[0];	// The first argument is the directory
nchannels=args[1];	// The second argument is the number of channels
ch2_name=args[2];	// The third argument is the name of the protein in channel 2
ch3_name=args[3];	// The fourth argument is the name of the protein in channel 3
ch4_name=args[4];	// The fifth argument is the name of the protein in channel 4
dirp=args[5];	// The sixth argument is the directory of the program

setBatchMode(true);	// run in batch mode to save time and memory

dirp1= dirp + "padd_to_resize.ijm"+ File.separator; // path to rezize function
// Make the new directories for projected images and the montage
newDir1 = dirsa + "zproject" + File.separator; 
newDirm = dirsa + "Montage" + File.separator; 
File.makeDirectory(newDir1); 
File.makeDirectory(newDirm); 
 
/////// CHANNEL 1
dir=dirsa + "indivisual_nuclei_DNA"+ File.separator;
filename = getFileList(dir); // get the list of files

// make sub directory for the first channel
newDir = newDir1 + "nuclei_DNA" + File.separator; 
File.makeDirectory(newDir); 

for (i=0; i<filename.length; i++) { 
      path=dir+filename[i]; // set path 
	  open(path); //open the image
  
      run("Z Project...", "projection=[Max Intensity]"); // project images
      saveAs("tiff", newDir + getTitle); // save the image
      run("Close All"); //close the images
} 
runMacro(dirp1,newDir)

run("Image Sequence...", "open=newDir sort"); // open the projected images as a stack
run("Scale Bar...", "width=5 font=14 color=White background=None location=[Lower Right] bold hide label"); // add scale bar of 50 microns
run("Make Montage...", "scale=1 border=3"); // Make a montage
saveAs("Tiff", newDirm + "Montage_nuclei_DNA"); // save the image of the montage

/////// CHANNEL 2
if(nchannels>1){ 
	// make sub directory for the first channel
	newDir = newDir1 + "nuclei_"+ch2_name + File.separator; 
	File.makeDirectory(newDir); 

	for (i=0; i<filename.length; i++) { 
     	 path=dir+filename[i]; // set path 
	  open(path); //open the image
  
      run("Z Project...", "projection=[Max Intensity]"); // project images
      saveAs("tiff", newDir + getTitle); // save the image
      run("Close All"); //close the images
	} 
	
runMacro(dirp1,newDir)

run("Image Sequence...", "open=newDir sort"); // open the projected images as a stack
run("Scale Bar...", "width=5 font=14 color=White background=None location=[Lower Right] bold hide label"); // add scale bar of 50 microns
run("Make Montage...", "scale=1 border=3"); // Make a montage
saveAs("Tiff", newDirm + "Montage_nuclei_"+ch2_name); // save the image of the montage

}

/////// CHANNEL 3
if(nchannels>2){ 
	// make sub directory for the first channel
	newDir = newDir1 + "nuclei_"+ch3_name + File.separator; 
	File.makeDirectory(newDir); 

	for (i=0; i<filename.length; i++) { 
     	 path=dir+filename[i]; // set path 
	  open(path); //open the image
  
      run("Z Project...", "projection=[Max Intensity]"); // project images
      saveAs("tiff", newDir + getTitle); // save the image
      run("Close All"); //close the images
	} 
	
runMacro(dirp1,newDir)

run("Image Sequence...", "open=newDir sort"); // open the projected images as a stack
run("Scale Bar...", "width=5 font=14 color=White background=None location=[Lower Right] bold hide label"); // add scale bar of 50 microns
run("Make Montage...", "scale=1 border=3"); // Make a montage
saveAs("Tiff", newDirm + "Montage_nuclei_"+ch3_name); // save the image of the montage

}
/////// CHANNEL 4
if(nchannels>3){ 
	// make sub directory for the first channel
	newDir = newDir1 + "nuclei_"+ch4_name + File.separator; 
	File.makeDirectory(newDir); 

	for (i=0; i<filename.length; i++) { 
     	 path=dir+filename[i]; // set path 
	  open(path); //open the image
  
      run("Z Project...", "projection=[Max Intensity]"); // project images
      saveAs("tiff", newDir + getTitle); // save the image
      run("Close All"); //close the images
	} 
	
runMacro(dirp1,newDir)

run("Image Sequence...", "open=newDir sort"); // open the projected images as a stack
run("Scale Bar...", "width=5 font=14 color=White background=None location=[Lower Right] bold hide label"); // add scale bar of 50 microns
run("Make Montage...", "scale=1 border=3"); // Make a montage
saveAs("Tiff", newDirm + "Montage_nuclei_"+ch4_name); // save the image of the montage

}