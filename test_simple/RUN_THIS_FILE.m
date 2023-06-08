%% Run the file to see the results directly
% Just change the file name and run it with the current unit of the file
% note that there are text in the exemplified file, so the function importdata take only the data part if it's a pure text file 
% the function will not work if there are no text in the file
% By default the function will output a plot with matched spikes and a summary distribution of the
% the output figures may look difference it's dur to the centroid numbers decided by the algorithum has bias with human eyes
% But the matched results a re similar 
[SpikeFeatures, SpikeLocation] = NIE_analysis_Simp("1a_AuNpsNIE.txt","1b_AuNpsBlank.txt","A","Oxi");