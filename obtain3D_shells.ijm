/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////  SCRIPT TO MAKE UNIDISTANT SHELLS OF SIZE EQUIVALENT TO THE RESOLUTION ALONG Z         											                    /////////////
///////  WRITTEN BY: SARADHA VENKATACHALAPATHY                                                                                                                   /////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
dir= getArgument();


//dirsa=getDirectory("Please choose the source directory");

dirsa = dir+"indivisual_nuclei_DNA"+File.separator;

filenames1=getFileList(dirsa);

dir2=dir+"Shells"+File.separator;
newDir113 = dir2 + "geometric_measures"+ File.separator;
newDir114= dir2 + "intensity_shells"+ File.separator;
File.makeDirectory(dir2); 
File.makeDirectory(newDir113); 
File.makeDirectory(newDir114); 
setBatchMode(true);
	
run("3D OC Options", "volume surface nb_of_obj._voxels nb_of_surf._voxels centroid std_dev_distance_to_surface bounding_box show_masked_image_(redirection_requiered) dots_size=5 font_size=10 redirect_to=none");
run("3D Manager Options", "volume surface compactness fit_ellipse 3d_moments integrated_density mean_grey_value std_dev_grey_value mode_grey_value feret minimum_grey_value maximum_grey_value centroid_(pix) centroid_(unit) distance_to_surface centre_of_mass_(pix) centre_of_mass_(unit) bounding_box radial_distance surface_contact closest use distance_between_centers=10 distance_max_contact=1.80");
run("Set Measurements...", "area mean standard modal min centroid center perimeter bounding fit shape feret's integrated median skewness kurtosis area_fraction limit display redirect=None decimal=3");

for(i=0; i<filenames1.length; i++){
		//open the nucleus can convert it to 8 bit image
		path3=dirsa+filenames1[i];
		open(path3);
		imgName1=getTitle(); 
		run("8-bit");
		
		//bin the image so that the pixel size is the similar across xy and z
		getVoxelSize(x, y, h, micron);
		f=round(h/x);
		run("Bin...", "x=f y=f z=1 bin=Average");
		getVoxelSize(x, y, h, micron);
		print(imgName1,",",x,",",y,",",h);
		
		//get filename, generate 3d object,  get the volume, and add the objects to the 3dmanager name the object
		label=File.getName(path3);
		selectWindow(imgName1);

		run("Gaussian Blur...", "sigma=3 stack");
		setAutoThreshold("Huang dark");
		run("Convert to Mask", "method=Default background=Default");
		run("Invert LUT");
		run("Fill Holes", "stack");
		run("Erode", "stack");
		run("Dilate", "stack");
		run("3D Objects Counter", "threshold=1 min.=0 max.=1000000000000 statistics objects");
		selectWindow(imgName1);
		Vol=0;
		for(l=0;l<nResults;l++){
			Vol=Vol+getResult("Volume (micron^3)",0); // get the volume of the spheroid
		}
		if (Vol >50){
			run("3D Manager");
			Ext.Manager3D_AddImage();
			Ext.Manager3D_Count(nb_obj); 
			if(nb_obj>1){
				Ext.Manager3D_SelectAll();
				Ext.Manager3D_Merge();
			}
			Ext.Manager3D_Select(0);
			Ext.Manager3D_Rename("nucleus_0");
		
			//erode the object till the volume is is less than or equal to 10 cu.microns. 
			p=0;
			print(getTitle(),",",Vol,",",p);
			while(Vol>10){
				selectWindow(imgName1);
				run("Erode (3D)", "iso=255");
				run("Z Project...", "projection=[Max Intensity]");
				//run("Threshold...");
				setThreshold(1, 4294967296);
				run("Measure");
				a=getResult("Area",0);
			
				if(a>0 ){
					run("Close");
					selectWindow(imgName1);
					run("3D Objects Counter", "threshold=1 min.=0 max.=1000000000000 statistics objects");
					Vol=0;
					for(l=0;l<nResults;l++){
						Vol=Vol+getResult("Volume (micron^3)",0);
					}
					Ext.Manager3D_AddImage();
					p=p+1;
					print(getTitle(),",",Vol,",",p);
					Ext.Manager3D_Count(nb_obj); 
					if(nb_obj>(p+1)){
						Ext.Manager3D_MultiSelect();
						for(k=p;k<nb_obj;k++){
							Ext.Manager3D_Select(k);
						}
					Ext.Manager3D_Merge();
					Ext.Manager3D_DeselectAll();		
					}
					Ext.Manager3D_Count(nb_obj); 
					Ext.Manager3D_Select((nb_obj-1));
					Ext.Manager3D_Rename("nucleus_"+p);
					Ext.Manager3D_DeselectAll();
				}
				else{
					Vol=0;
				}
			}
			Ext.Manager3D_Count(p); 
			run("Close All");
		
			path3=dirsa+filenames1[i];
			open(path3);
			run("8-bit");

			//bin the image so that the pixel size is the similar across xy and z
			getVoxelSize(x, y, h, micron);
			f=round(h/x);
			run("Bin...", "x=f y=f z=1 bin=Average");
			
			Ext.Manager3D_DeselectAll();
			Ext.Manager3D_SelectAll();
    		Ext.Manager3D_Measure();
			Ext.Manager3D_SaveResult("M",newDir113+filenames1[i]+".tsv");
			Ext.Manager3D_CloseResult("M");

			Ext.Manager3D_DeselectAll();
			Ext.Manager3D_SelectAll();
			Ext.Manager3D_Quantif();
			Ext.Manager3D_SaveQuantif(newDir114+filenames1[i]+".tsv");
			Ext.Manager3D_CloseResult("Q");
	
			run("Close All");
    		Ext.Manager3D_Close();
		}

		
    	run("Clear Results");
		
		run("Close All");
    	run("Collect Garbage");
		run("Collect Garbage");
		run("Collect Garbage");
		


}






	