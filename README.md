# B3AM Userguide
version 1.0
by Katrin Löer
k.loer@tudelft.nl

[![View B3AM on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://uk.mathworks.com/matlabcentral/fileexchange/128489-b3am)

Python version: https://github.com/cl-finger/B3Ampy
----------
B3AM is a toolbox for easy and fast beamforming analysis of three-component array data providing estimates of surface wave dispersion curves, frequency-dependent wavefield composition, and the direction of arrival for different wave types and frequencies from ambient seismic noise. 

B3AM performs beamforming analysis on short time-frequency windows of the provided data and identifies maxima in the beam response of each window. Each detected maximum is characterised by its wavenumber, direction of arrival (azimuth), and polarisation, i.e., wave type. A summary of the results of all windows is provided in histograms that show, for example, wavenumber as a function of frequency for different wave types.

Changes since v0.1:
- Example added 
- b3am_param.m: frequency range now to be provided manually
- plot_b3am.m: bug fixed in direction of arrival plots; export_fig.m not required anymore
- f_extrema24.m: error fixed in computation of standard deviation
- b3am_convert_iris.m/b3am_convert_mseed.m: added option to zero-pad incomplete data, improved functionality
- date2doy.m function provided for convenience
- batlowS.mat categorical colour map (Crameri, 2021) provided for convenience
- plots provided in the Folder 'Figures' have been created with this version

Parkfield Example
=================

The folder _Example_Parkfield_ contains example output data and figures as returned by B3AM for one day of ambient noise recorded at the Parkfield array, California, US (Thurber and Roecker, 2000). The data are publicly available from the Seismological Facility for the Advancement of Geoscience (SAGE, former IRIS), and can be downloaded directly into MATLAB. Go to the SAGE homepage to download the MATLAB script irisFetch.m (http://ds.iris.edu/ds/nodes/dmc/software/downloads/irisfetch.m/) and the required Java library (http://ds.iris.edu/ds/nodes/dmc/software/downloads/IRIS-WS/2-20-1/#Download). You can then use the script iris_getrawdata_example.m provided with the B3AM package to download data from the Parkfield (or another) array. In the script, specify the path to your irisFetch.m script and the Java library as both will be used in **iris_getrawdata_example.m**. Further parameters you need to define are the start and end date, the network code, names of stations in the network, channels, and storage location. **The default values in the script correspond to those used for the example**. Expect the download to take up to a few minutes per station for a single day of data depending on your network speed (here, it took around 25 minutes to download data from 34 stations).

Follow the steps below to reprodoce the figures in Example_Parkfield/Figures. Again, default values provided in the code should produce the example output.

0) What you need
-----------------

- Matlab R2020b or newer
- your 3-component array data (in mseed or .mat format)
- a txt-file containing information about the station locations in three columns:
	stationname	 longitude in degree	latitude in degree

Additional functions required NOT provided with this package - please download from FileExchange:
- rdmseed.m (https://uk.mathworks.com/matlabcentral/fileexchange/28803-rdmseed-and-mkmseed-read-and-write-miniseed-files)
- Crameri Scientific Colormaps (https://uk.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps)

Functions by other authors INCLUDED in this toolbox for convenience:
- deg2utm.m by Rafael Palacios (https://www.mathworks.com/matlabcentral/fileexchange/10915-deg2utm)
- date2doy.m by Anthony Kendall (https://www.mathworks.com/matlabcentral/fileexchange/18316-date-to-decimal-day-of-year)
- extrema.m and extrema2.m by Carlos Adrin Vargas Aguilera (https://www.mathworks.com/matlabcentral/fileexchange/12275-extrema-m-extrema2-m)
- FK3C package by Nima Riahi (https://github.com/nimariahi/fk3c)

Let's get started...
====================

1) Bring your data into the correct shape
------------------------------------------

B3AM will handle one file per day that contains seismic traces from all stations and all components.
Your data need to be sorted by component in the order E, N, Z.
If your data is already a .mat file (for example after downloading from IRIS directly into Matlab) you can use the script
**b3am_convert_iris.m** to bring the traces into the correct order. If your data comes in miniseed format, please use **b3am_convert_mseed.m** and consider the guidelines provided in the script. 

The output from either conversion script will be a (or multiple) file(s) called DAT_NN_yyyyddd.mat
where NN is the network code, yyyy is the year, and ddd the day of the year (between 001 and 365).
If you are working with multiple days of ambient noise data, all resulting DAT files will be stored in the same folder.

Eventually, B3AM requires information about the station location in form of a txt file that contains three columns:
stationname	 x-ccordinate in m	y-coordinate in m
When running the provided conversion files for IRIS or mseed data, the station file is created automatically in the same folder as the rearranged data.

Alternatively, this file can be created from any DAT file using the script **mk_stationfile.m**
DAT contains a variable called DAT.h.coords with lat and lon already converted to meter. 
The file returned from mk_stationfile.m will be called stations_utm_NN_yyyyddd.txt

The information in the stationfile will be used to compute theoretical minimum and maximum wavenumber values based on the largest and smallest station spacing within the array.
Note, however, that these limits will be the same for all days processed, so that the txt file should contain location data for the complete array, and not just the stations active on a particular day.

Before you proceed it is recommended to do a quality check on the rearranged data to see if everything is as expected and to delete stations that seem unfeasible, for example. The function 'wiggle' (available from FileExchange) can be used to display traces from multiple stations at a time. Make sure to not only delete one trace but all three channels of a dodgy station and also the corresponding information in DAT.h.channel, DAT.h.stations and DAT.h.coords. (You might find it easiest to redo the conversion and exclude the respective station from the input folder containing the mseed data.)

2) Define processing parameters to perform beamforming
-------------------------------------------------------

Open the script **b3am_param.m** and fill in the required information line by line. 
Comments and examples are provided to help you with the correct format etc.
Essential parameters concern:
- in- and output data locations
- pre-processing preferences
- frequency range
- wavenumber range
- beamforming method
- parallel computing

Run the script **b3am_check.m** which will print in the command line
- min/max station spacing
- min/max wavelength
- min/max wavenumber
- min/max frequency
and display the array response function indicating the minimum (black solid line) and maximum wavenumber (maximum of x-axis). If you are happy with these values you can proceed to the next step, otherwise make appropriate changes in b3am_param.m.

Wathelet et al. (2008) suggest to limit the wavenumber range to values for which the ARF is < 0.5 of the global maximum to obtain more robust results. So, if you observe significant sidelobes in the ARF with normalised amplitudes > 0.5 you should decrease your maximum wavenumber (kmax) accordingly. Equally, if the central peak shows amplitudes > 0.5 for values > kmin (plotted as a thick vertical line) you should increase kmin.

3) Run the beamformer
----------------------

Once all processing parameters are defined in b3am_param.m you can start the beamformer by running **b3am.m**.
Note that this script can partly run in parallel so you might want to make use of parallel computing facilities on a remote cluster.

B3AM performs the four major steps successively:
i) pre-processing (can include filtering, spectral whitening, one-bit normalisation)
ii) Fourier transformation (stored in temporary folder tmpFT/)
iii) frequency-wavenumber analysis (beamforming)
iv) identification of extrema (maxima) in the beam responses

The script provides output in the command line documenting its progress.
(When running on a remote cluster you should specify a file where this information is saved to.)

4) Retrieve the results
------------------------

In your output folder the file procpars.mat appears that contains all processing parameters used in the beamforming process, such as the resolution of the wavenumber grid (kgrid), the frequency range (freqs), and so on. These values will be taken from this file when plotting.

The beamforming results are stored in the output directory that you specified in b3am_param.m
one filed is saved for each frequency ffff, called kmax_NN_yyyyddd_ffff.mat
The information stored in each such file refer to the maxima in the beam responses and are
a_all		% beam power at all extrema
kr_all		% radial wavenumber at all extrema
kth_all		% azimuth of all extrema
pola_all	% polarization parameters for all extrema
pola_ind	% polarization indeces for all extrema
wave_ind	% wave type indeces for all extrema (0: P, 1: SV, 2: SH/Love, 3: retro. Rayleigh, 4: pro. Rayleigh)
The size of all these variables is [nwin x nwin], where nwin denotes the number of time windows and nmax the maximum number of maxima detected in any time window.
pola_all will be [nwin x nmax x 4] because 4 polarization parameters are stored (azimuth, dip, ellipticity, tilt).
When plotting, you can choose to only consider the first maximum in each window (option 'MAX1'). 

Note that beam response maps are not plotted (as you would create one for each time window processed)!

5) Plot the results
--------------------

To get a first overview of the beamforming results, you can use the script **plot_b3am.m**.
Provide the location of the beamforming results, i.e., the max files ('dir_in'), and a directory to save the figures in ('dir_out').
If all plot options are set to true you will obtain the following 8 Figures:
1) Wavefield composition: absolute contribution per frequency
2) Wavefield composition: relative contribution per frequency
3) Wavefield composition: amplitude variation with frequency
4)-6) f-k histogram plots of pro- and retrograde Rayleigh waves and Love waves, respectively
7) Surface wave dispersion curves
8) Direction of arrival as a function of frequency for all surface and body waves
Set the save options separately to decide if figures are to be saved.
The following parameters also need to be defined:
- SNR: A value for an acceptable SNR is required for the picking of dispersion curves from histograms. Experiment with this value as it will depend on the data quality, length of recording, number of stations, etc.
- maxflag: choose if you want to consider only the first/largest maximum ('MAX1') or all maxima ('NOMAX') detected in each time window 
- countflag: to plot wavefield composition decide if you want to consider the number of waves counted ('count') or the number weighted by beam power amplitude ('amp')

Note that, prior to plotting, this script performs essential analysis steps. In the section "Wave type analysis" results are sorted with respect to their detected wave type before they can be plotted accordingly.

B3AM in the Literature
======================

Finger, C., & Löer, K. (2024). Depth of sudden velocity changes derived from multi‐mode Rayleigh waves. Journal of Geophysical Research: Solid Earth, 129(3), e2023JB028322.

Löer, K., Finger, C., Obiri, E., and Kennedy, H.: A comprehensive beamforming toolbox to characterise surface and body waves in three-component ambient noise wavefields, EGU General Assembly 2023, Vienna, Austria, 24–28 Apr 2023, EGU23-5670, https://doi.org/10.5194/egusphere-egu23-5670, 2023.

Löer, K., Toledo, T., Norini, G., Zhang, X., Curtis, A., Saenger, E.H.: Imaging the Deep Structures of Los Humeros Geothermal Field, Mexico, Using Three‐Component Seismic Noise Beamforming, Seismological Research Letters (2020) 91 (6): 3269–3277.

Löer, K., Riahi, N., & Saenger, E. H. , Three-component ambient noise beamforming in the Parkfield area, Geophysical Journal International, Volume 213, Issue 3, June 2018, Pages 1478–1491, https://doi.org/10.1093/gji/ggy058

Obiri, E., Löer, K., & Finger, C. (2023, June). Wavefield composition analysis from three-component beamforming improves thickness estimates of sedimentary layers. In 84th EAGE Annual Conference & Exhibition (Vol. 2023, No. 1, pp. 1-5). European Association of Geoscientists & Engineers.

Riahi, N., and Saenger, E. H. (2014), Rayleigh and Love wave anisotropy in Southern California using seismic noise, Geophys. Res. Lett., 41, 363– 369, doi:10.1002/2013GL058518.

Riahi, N., Bokelmann, G., Sala, P., and Saenger, E. H. (2013), Time-lapse analysis of ambient surface wave anisotropy: A three-component array study above an underground gas storage, J. Geophys. Res. Solid Earth, 118, 5339– 5351, doi:10.1002/jgrb.50375.

REFERENCES
===========

Crameri, Fabio. (2021). Scientific colour maps (7.0.1). Zenodo. https://doi.org/10.5281/zenodo.5501399

Thurber, C. and Roecker, S. Parkfield Passive Seismic Array [Data set]. International Federation of Digital Seismograph Networks, 2000. http://doi.org/10.7914/SN/XN_2000

Wathelet, M., Jongmans, D., Ohrnberger, M., & Bonnefoy-Claudet, S. (2008). Array performances for ambient vibrations on a shallow structure and consequences over V s inversion. Journal of Seismology, 12, 1-19.

