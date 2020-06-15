/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////  SCRIPT TO IDENTIFY SEGMENT NUCLEI IN 3D AND MEASURE THE INTENSITY IN CHANNEL 1 AND 2         											                    /////////////
///////  WRITTEN BY: SARADHA VENKATACHALAPATHY                                                                                                                   /////////////
///////  ASSUMPTIONS: The input image is a confocal zstack and the channel 1 contains the nucleus 
///////  DESCRIPTION: Two dialog box opens where the user inputs the source (where the folder containing images are strored). The program creates 2 subfolders where
///////               cropped 3D nuclei are stored. The program opens nucleus channel, smoothens the stack and the thresholds the image and identifies objects in 3D. 
///////               The program then opens the raw image and crops the image in the all channels individually. 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

args = split(getArgument()," "); 	// read in the array of argurments that are separated by a space
dirsa = args[0];	// The first argument is the directory
nchannels=args[1];	// The second argument is the number of channels
ch2_name=args[2];	// The third argument is the name of the protein in channel 2
ch3_name=args[3];	// The fourth argument is the name of the protein in channel 3
ch4_name=args[4];	// The fifth argument is the name of the protein in channel 4

setBatchMode(true);	// run in batch mode to save time and memory

///creating new subdirectories for storing the processed data//
dirb= dirsa + "3d obj nuclei"+ File.separator;	// define the path to a new folder that will contain the image afer 3d segmentation to indentify spheroids
dir= dirsa + "indivisual_nuclei_DNA"+ File.separator;		// define the path to a new folder that will contain the first channel of the raw image cropped to represent each spheroid
dirc2= dirsa + "indivisual_nuclei_"+ch2_name+ File.separator; 	// define the path to a new folder that will contain the second channel of the raw image cropped to represent each spheroid
dirc3= dirsa + "indivisual_nuclei_"+ch3_name+ File.separator; 	// define the path to a new folder that will contain the third channel of the raw image cropped to represent each spheroid
dirc4= dirsa + "indivisual_nuclei_"+ch4_name+ File.separator; 	// define the path to a new folder that will contain the fourth channel of the raw image cropped to represent each spheroid
dir1= dirsa + "data"+ File.separator;
// create the paths defined above 
File.makeDirectory(dirb);
File.makeDirectory(dir); 

if(nchannels>1){ 
	File.makeDirectory(dirc2); 
}
if(nchannels>2){ 
	File.makeDirectory(dirc3); 
}
if(nchannels>3){ 
	File.makeDirectory(dirc4); 
}


//dirsa=getDirectory("Please choose the source directory");
dirw= dirsa + "after watershed"+ File.separator;
dirraw= dirsa + "rawimages"+ File.separator;


run("3D OC Options", "volume surface nb_of_obj._voxels nb_of_surf._voxels integrated_density mean_gray_value std_dev_gray_value median_gray_value minimum_gray_value maximum_gray_value centroid mean_distance_to_surface std_dev_distance_to_surface median_distance_to_surface centre_of_mass bounding_box dots_size=5 font_size=20 redirect_to=none");
run("3D Manager Options", "volume surface compactness fit_ellipse 3d_moments integrated_density mean_grey_value std_dev_grey_value mode_grey_value feret minimum_grey_value maximum_grey_value centroid_(pix) centroid_(unit) distance_to_surface centre_of_mass_(pix) centre_of_mass_(unit) bounding_box radial_distance surface_contact closest use distance_between_centers=10 distance_max_contact=1.80");

filenames=getFileList(dirw);
baseName=newArray(filenames.length);

for(f=0;f<filenames.length;f++){

	
	open(dirw+filenames[f]);
	imgName=getTitle(); 
	baseNameEnd=indexOf(imgName, ".tiff"); 
	baseName[f]=substring(imgName, 0, baseNameEnd); 
	
	getVoxelSize(width, height, depth, unit);
	a=200/(width*height*depth);// maximum volume is 1500 cu.microns
	a_1=50/(width*height*depth);// minimum volume is 200 cu.microns
	
	run("3D Objects Counter", "threshold=128 slice=12 min.=a_1 max.=a objects"); //Size filter
	saveAs("Tiff", dirb+baseName[f]+".tiff"); 
	
	run("3D Manager");
	Ext.Manager3D_AddImage();
	Ext.Manager3D_Count(nb_obj); 
	
	for(p=0;p<nb_obj;p++){
			t=p+1;
			Ext.Manager3D_Bounding3D(p,x0,x1,y0,y1,z0,z1);
			zlows=z0+1;
			zhighs=z1+1;
			wids=x1-x0;
			heigs=y1-y0;

			path=dirraw+baseName[f]+".nd2";
			
			//Channel1
			run("Bio-Formats", "open=path color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT z_begin=zlows z_end=zhighs c_begin=1 c_end=1 c_step=1");
			
			makeRectangle(x0, y0, wids, heigs);
			run("Duplicate...", "duplicate");
			names=baseName[f]+"_nucleus_"+x0+"_"+y0+"_"+z0+"_"+z1;
			saveAs("tiff",dir+names);
			
			run("Close All");
			
			if(nchannels>1){
				//Channel2
				run("Bio-Formats", "open=path color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT z_begin=zlows z_end=zhighs c_begin=2 c_end=2 c_step=1");
				makeRectangle(x0, y0, wids, heigs);
				run("Duplicate...", "duplicate");
				names=baseName[f]+"_nucleus_"+x0+"_"+y0+"_"+z0+"_"+z1;
				saveAs("tiff",dirc2+names);
				
				run("Close All");
			}

			if(nchannels>2){
				//Channel3
				run("Bio-Formats", "open=path color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT z_begin=zlows z_end=zhighs c_begin=3 c_end=3 c_step=1");
				makeRectangle(x0, y0, wids, heigs);
				run("Duplicate...", "duplicate");
				names=baseName[f]+"_nucleus_"+x0+"_"+y0+"_"+z0+"_"+z1;
				saveAs("tiff",dirc3+names);
				
				run("Close All");
			}
		
			if(nchannels>3){
				//Channel4
				run("Bio-Formats", "open=path color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT z_begin=zlows z_end=zhighs c_begin=4 c_end=4 c_step=1");
				makeRectangle(x0, y0, wids, heigs);
				run("Duplicate...", "duplicate");
				names=baseName[f]+"_nucleus_"+x0+"_"+y0+"_"+z0+"_"+z1;
				saveAs("tiff",dirc4+names);
				run("Close All");
			}
		
		
		//run("Select None");
     	}
	
	
	run("Close All");
    Ext.Manager3D_Close();
   	run("Clear Results");
}



