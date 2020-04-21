# Wide-Field-Deconvolution
Inferring Neural Population Spiking Rate from Wide-Field Calcium Imaging.

The files attached include deconvolution algorithm implementations for recovering spiking rate in a pixel (or any ROI) from wide-field recordings using the DF/F fluorescence trace of the pixel (or ROI). 
Each algorithm requires a single parameter tuning that fits all pixels/ROIs from the same set of experiments.    

Details can be found in the accompanying manuscript https://www.biorxiv.org/content/10.1101/2020.02.01.930040v1

If you hold wide-field fluorescence recordings, after DF/F, and wish to use this code to find the spiking rate behind each pixel/region, follow the deconv_Dff_data_all_steps_using_continuously_varying.mat. 

Father support on the following topics, and more, can be found at Helpful_Tips file: 
*How to speed it up
*How sensitive are the results to the parameter values (spoiler, they are not)
*How to estimate the calcium decay rate
*The calcium baseline 
