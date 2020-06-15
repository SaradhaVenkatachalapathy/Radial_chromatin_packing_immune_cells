dirsa = getArgument();
//dirsa=getDirectory("choose folder");
dirb= dirsa + "3d obj nuclei"+ File.separator;
dirraw= dirsa + "rawimages"+ File.separator;

dirbc= dirsa + "3d_objects_cell_2microns"+ File.separator;
dirc= dirsa + "cell_2microns_measure"+ File.separator;
dirc1= dirc + "ch1"+ File.separator;
dirc2= dirc + "ch2"+ File.separator;
dirc3= dirc + "ch3"+ File.separator;

dirn= dirsa + "nuc_microns_measure"+ File.separator;
dirn1= dirn + "ch1"+ File.separator;
dirn2= dirn + "ch2"+ File.separator;
dirn3= dirn + "ch3"+ File.separator;


File.makeDirectory(dirbc); 
File.makeDirectory(dirc);
 File.makeDirectory(dirc1);
 File.makeDirectory(dirc2);
 File.makeDirectory(dirc3);
 File.makeDirectory(dirn); 
 File.makeDirectory(dirn1); 
 File.makeDirectory(dirn2); 
 File.makeDirectory(dirn3); 
 
run("3D OC Options", "volume surface nb_of_obj._voxels nb_of_surf._voxels integrated_density mean_gray_value std_dev_gray_value median_gray_value minimum_gray_value maximum_gray_value centroid mean_distance_to_surface std_dev_distance_to_surface median_distance_to_surface centre_of_mass bounding_box dots_size=5 font_size=20 redirect_to=none");
run("3D Manager Options", "volume surface compactness fit_ellipse 3d_moments integrated_density mean_grey_value std_dev_grey_value mode_grey_value feret minimum_grey_value maximum_grey_value centroid_(pix) centroid_(unit) distance_to_surface centre_of_mass_(pix) centre_of_mass_(unit) bounding_box radial_distance surface_contact closest use distance_between_centers=10 distance_max_contact=1.80");


setBatchMode(true);
filenames=getFileList(dirb);
filenamesraw=getFileList(dirraw);
baseName=newArray(filenames.length);

for(f=0;f<filenames.length;f++){

	//Open segmented image
	open(dirb+filenames[f]);
	imgName=getTitle(); 
	baseNameEnd=indexOf(imgName, ".tiff"); 
	baseName[f]=substring(imgName, 0, baseNameEnd); 

	// calculate the number of dilations to get 2 microns
	getPixelSize(unit, pixelWidth, pixelHeight);
	num_dia=round(2/pixelWidth);
	
	//Add to the 3D Manager
	run("3D Manager");
	Ext.Manager3D_AddImage();
	Ext.Manager3D_Count(nb_obj); 

	names=newArray(nb_obj);
	print("\\Clear"); 
	
	//figureout the name of the nucleus and dilate the nucleus upto 2 microns from edge
	for(p=0;p<nb_obj;p++){
		t=p+1;
		//figureout the name of the nucleus 
		Ext.Manager3D_Bounding3D(p,x0,x1,y0,y1,z0,z1);
		zlows=z0+1;
		zhighs=z1+1;
		wids=x1-x0;
		heigs=y1-y0;
		Ext.Manager3D_Select(p);
		names[p]=baseName[f]+"_nucleus_"+x0+"_"+y0+"_"+z0+"_"+z1;
		Ext.Manager3D_Rename(names[p]);
		
	}
	Ext.Manager3D_SelectAll();
	//Nuclear levels of the protiens
	path=dirraw+filenamesraw[f];
	run("Bio-Formats", "open=path color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT c_begin=1 c_end=1 c_step=1");
	Ext.Manager3D_Quantif();
	Ext.Manager3D_SaveQuantif(dirn1 +baseName[f]+".tsv");
	Ext.Manager3D_CloseResult("Q");
	run("Close All");
	Ext.Manager3D_SelectAll();
	run("Bio-Formats", "open=path color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT c_begin=2 c_end=2 c_step=1");
	Ext.Manager3D_Quantif();
	Ext.Manager3D_SaveQuantif(dirn2 +baseName[f]+".tsv");
	Ext.Manager3D_CloseResult("Q");
	run("Close All");
	Ext.Manager3D_SelectAll();
	run("Bio-Formats", "open=path color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT c_begin=3 c_end=3 c_step=1");
	Ext.Manager3D_Quantif();
	Ext.Manager3D_SaveQuantif(dirn3 +baseName[f]+".tsv");
	Ext.Manager3D_CloseResult("Q");
	run("Close All");
	Ext.Manager3D_Close();

	open(dirb+filenames[f]);
	
	// dilate the nucleus upto 2 microns from edge
	for(p=0;p<nb_obj;p++){
		t=p+1;
		for(k=1; k<=num_dia;k++){
			run("Dilate (3D)", "iso=t");
		}
	}
	saveAs("tiff",dirbc+baseName[f]);
	tit=getTitle();
	
	Ext.Manager3D_SelectAll();
	path=dirraw+filenamesraw[f];
	Ext.Manager3D_SelectAll();
	run("Bio-Formats", "open=path color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT c_begin=1 c_end=1 c_step=1");
	Ext.Manager3D_Quantif();
	Ext.Manager3D_SaveQuantif(dirc1 +baseName[f]+".tsv");
	Ext.Manager3D_CloseResult("Q");
	run("Close All");
	Ext.Manager3D_SelectAll();
	run("Bio-Formats", "open=path color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT c_begin=2 c_end=2 c_step=1");
	Ext.Manager3D_Quantif();
	Ext.Manager3D_SaveQuantif(dirc2 +baseName[f]+".tsv");
	Ext.Manager3D_CloseResult("Q");
	run("Close All");
	Ext.Manager3D_SelectAll();
	run("Bio-Formats", "open=path color_mode=Grayscale specify_range view=Hyperstack stack_order=XYCZT c_begin=3 c_end=3 c_step=1");
	Ext.Manager3D_Quantif();
	Ext.Manager3D_SaveQuantif(dirc3 +baseName[f]+".tsv");
	Ext.Manager3D_CloseResult("Q");
	run("Close All");
	
		
     run("Close All");
		
	Ext.Manager3D_Close();
    run("Close All");
   	run("Clear Results");

   	run("Collect Garbage");

	
}


