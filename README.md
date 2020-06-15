# Radial_chromatin_packing_immune_cells
Images were processed and further analyzed using custom programs in Fiji and R. 

The nuclear boundaries were segmented in 3D using the DAPI channels to identify individual nuclei. These nuclei were eroded by 0.5 microns in x,y, and z iteratively till the volume of the eroded nucleus was less than 10 cubic microns. Then the mean intensity of each 3D ring (width 0.5 microns) in the nucleus was computed for all cells. The intensity fraction was calculated by normalizing the mean ring intensity for each nucleus(maximum=1). This is visualized in FigSA. Linear interpolation was then used to compute the intensity fraction of rings that occupy  0-10% to 90-100% volume fraction of the nucleus. The heatmap in Fig XA was visualized using functions from gplots, RColorBrewer and dendextend. 

In order to calculate the cellular levels of proteins, the 3D nuclear object was dilated by 2 microns in x,y and z to obtain the corresponding cell. This was efficient as the cells were all spherical with high karyoplasmic index. The total intensity in the 3D cellular object was computed for each protein channel and their ratio was obtained for each cell. 

Hierarchical clustering was performed on the dissimilarity matrix obtained from (1-spearman’s correlation). 
