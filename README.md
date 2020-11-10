# Radial_chromatin_packing_immune_cells
Images were processed and further analyzed using custom programs in Fiji and R. 

### Image Analyis
The nuclear boundaries were segmented in 3D using the DAPI channels to identify individual nuclei. This is accomplished in *_Binary.ijm_* and *_3dcrop.ijm_*. 

These nuclei were eroded by 0.5 microns in x,y, and z iteratively till the volume of the eroded nucleus was less than 10 cubic microns. Then the mean intensity of each 3D ring (width 0.5 microns) in the nucleus was computed for all cells. This is accomplished in *_obtain3D_shells.ijm_*. 
The intensity fraction was calculated by normalizing the mean ring intensity for each nucleus(maximum=1). Linear interpolation was then used to compute the intensity fraction of rings that occupy  0-10% to 90-100% volume fraction of the nucleus. This is done in *_R_script_to_obtain_data.R_*

In order to calculate the cellular levels of proteins, the 3D nuclear object was dilated by 2 microns in x,y and z to obtain the corresponding cell. This was efficient as the cells were all spherical with high karyoplasmic index. The total intensity in the 3D cellular object was computed for each protein channel and their ratio was obtained for each cell. This is done in *_nuc_ring_2_micron.ijm_*

To visualize the images, zprojected montages of segemented nulei is visualized using *_zproject_and_padd_montage.ijm_* and *_Rescale_images.ijm_*

### Cluster Identification
Hierarchical clustering was performed on the dissimilarity matrix obtained from (1-spearman’s correlation) from the radial DNA intensity profile data. R libraries used: gplots, RColorBrewer and dendextend. 

[![DOI](https://zenodo.org/badge/272499677.svg)](https://zenodo.org/badge/latestdoi/272499677)
