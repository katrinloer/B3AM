# B3AM Userguide
(version 0.1; under development - please test and provide feedback!)

B3AM is a toolbox for easy and fast beamforming analysis of three-component array data providing estimates of surface wave dispersion curves, frequency-dependent wavefield composition, and the direction of arrival for different wave types and frequencies from ambient seismic noise. 

B3AM performs beamforming analysis on short individual time-frequency windows of the provided data and identifies maxima in the beam response of each window. Each detected maximum is characterised by its wavenumber, direction of arrival (azimuth), and polarisation, i.e., wave type. A summary of the results of all windows is provided in histograms that show, for example, wavenumber as a function of frequency for different wave types.

-----------------
0) What you need
- Matlab R2020b or newer
- your array data (in mseed or .mat format)
- a txt-file containing information about the station locations in three columns:
	stationname	 longitude in degree	latitude in degree

Download all files and folders in this repository in a directory of your choice.

Alternatively, if you are familiar with Matlab add-ons, you can download and install just the toolbox file
- B3AM.mltbx

Additional functions required NOT provided with this package - please download separately:
- rdmseed.m (https://uk.mathworks.com/matlabcentral/fileexchange/28803-rdmseed-and-mkmseed-read-and-write-miniseed-files)
- export_fig (https://uk.mathworks.com/matlabcentral/fileexchange/23629-export_fig)
- Crameri Scientific Colormaps (https://uk.mathworks.com/matlabcentral/fileexchange/68546-crameri-perceptually-uniform-scientific-colormaps)
- Categorical Palettes by Crameri (https://www.fabiocrameri.ch/categorical-colour-maps/)

Functions by other authors INCLUDED in this toolbox for convenience:
- deg2utm.m by Rafael Palacios
- extrema.m and extrema2.m by Carlos Adrin Vargas Aguilera
- FK3C package by Nima Riahi (https://github.com/nimariahi/fk3c)

Let's get started...

------------------------------------------
1) Bring your data into the correct shape

B3AM will handle one file per day that contains seismic traces from all stations and all components.
Your data need to be sorted by component in the order E, N, Z.
If your data is already a .mat file (for example after downloading from IRIS directly into Matlab) you can use the script
- b3am_convert_iris.m

to bring the traces into the correct order.
If your data comes in miniseed format, please use
- b3am_convert_mseed.m

and consider the guidelines provided in the script. 

The output from either conversion script will be a (or multiple) file(s) called
- DAT_NN_yyyyddd.mat

where NN is the network code, yyyy is the year, and ddd the day of the year.
If you are working with multiple days of ambient noise data, all resulting DAT files will be stored in the same folder.

Eventually, B3AM requires information about the station location in form of a txt file that contains three columns:
stationname	 x-ccordinate in m	y-coordinate in m
When running b3am_convert_mseed.m the station file is created automatically in the same folder as the rearranged data.

Alternatively, this file can be created from any DAT file using the script
- mk_stationfile.m

DAT contains a variable called DAT.h.coords with lat and lon already converted to meter. 
The file returned from mk_stationfile.m will be called
- stations_utm_yyyyddd.txt

where yyyy is the year (e.g., 2023) and ddd the day of year (between 001 and 365).

The information in the stationfile will be used to compute theoretical minimum and maximum wavenumber values
based on the largest and smallest station spacing within the array.
Note, however, that these limits will be the same for all days processed, so that the txt file should contain
location data for the complete array, and not just the stations active on a particular day.

Before you proceed it is recommended to do a quality check on the rearranged data to see if everything is as expected and to delete stations that seem unfeasible, for example. The function 'wiggle' (available from FileExchange) can be used to display traces from multiple stations at a time. Make sure to not only delete one trace but all three channels of a dodgy station and also the corresponding information in DAT.h.channel, DAT.h.stations and DAT.h.coords. (You might find it easiest to redo the conversion and exclude the respective station from the input folder containing the mseed data.)

-------------------------------------------------------
2) Define processing parameters to perform beamforming

Open the script
- b3am_param.m

and fill in the required information line by line. 
Comments and examples are provided to help you with the correct format etc.

Per default, the resolvable frequency range is based on conservative estimates of linearised average theoretical dispersion curves from literature values for various frequency ranges (e.g., Löer et al. 2018). 
However, the user can choose to adapt these values if, for example, information from theoretical dispersion curves is available. Note that frequency limits will be different for different modes of surface waves due to their different velocities (typically higher frequencies for higher modes), default values are adapted for the fundamental mode.

Run the script
- b3am_check.m

which will print in the command line
- min/max station spacing
- min/max wavelength
- min/max wavenumber
- min/max frequency

and display the array response function indicating the minimum (black solid line) and maximum wavenumber (maximum of x-axis). If you are happy with these values you can proceed to the next step, otherwise make appropriate changes in b3am_param.m.

Wathelet et al. (2008) suggest to limit the wavenumber range to values for which the ARF is < 0.5 of the global maximum to obtain more robust results. So, if you observe significant sidelobes in the ARF with normalised amplitudes > 0.5 you should decrease your maximum wavenumber (kmax) accordingly. Equally, if the central peak shows amplitudes > 0.5 for values > kmin (plotted as a thick vertical line) you should increase kmin.

----------------------
3) Run the beamformer

Once all processing parameters are defined in b3am_param.m you can start the beamformer by running
- b3am.m

Note that this script can partly run in parallel so you might want to make use of parallel computing facilities on a remote cluster. 

B3AM performs the four major steps successively:
- pre-processing (can include filtering, spectral whitening, one-bit normalisation, running-average-mean normalistion)
- Fourier transformation (stored in temporary folder tmpFT/)
- frequency-wavenumber analysis (beamforming)
- identification of maxima in the beam responses

The script provides output in the command line documenting its progress.
(When running on a remote cluster you should specify a file where this information is saved to.)

------------------------
4) Retrieve the results

In your output folder the file procpars.m appears that contains all processing parameters used in the beamforming process, such as the resolution of the wavenumber grid (kgrid), the frequency range (frees), and so on. These values will be taken from this file when plotting.

The beamforming results are stored in the output directory that you specified in b3am_param.m
one filed is saved for each frequency ffff, called
- kmax_NN_yyyyddd_ffff.mat

The information stored in each such file refer to the maxima in the beam responses and are
- a_all: beam power at all extrema
- kr_all: radial wavenumber at all extrema
- kth_all: azimuth of all extrema
- pola_all: polarization parameters for all extrema
- pola_ind: polarization indeces for all extrema
- wave_ind: wave type indeces for all extrema (0: P, 1: SV, 2: SH/Love, 3: retro. Rayleigh, 4: pro. Rayleigh)

The size of all these variables is [nwin x nwin], where nwin denotes the number of time windows and nmax the maximum number of maxima detected in any time window.
pola_all will be [nwin x nmax x 4] because 4 polarization parameters are stored (azimuth, dip, ellipticity, tilt).
When plotting, you can choose to only consider the first maximum in each window (option 'MAX1'). 

Note that beam response maps are not plotted (as you would create one for each time window processed)!

--------------------
5) Plot the results

To get a first overview of the beamforming results, you can use the script
- plot_b3am.m

Provide the location of the beamforming results, i.e., the max files ('dir_in'), and a directory to save the figures in ('dir_out').
If all plot options are set to true you will obtain the following 8 Figures:
- Wavefield composition: absolute contribution per frequency
- Wavefield composition: relative contribution per frequency
- Wavefield composition: amplitude variation with frequency
- f-k histogram plots of pro- and retrograde Rayleigh waves and Love waves, respectively
- Surface wave dispersion curves
- Direction of arrival as a function of frequency for all surface and body waves

Set the save options separately to decide if figures are to be saved.
The following parameters also need to be defined:
- SNR: A value for an acceptable SNR is required for the picking of dispersion curves from histograms. Experiment with this value as it will depend on the data quality, length of recording, number of stations, etc.
- maxflag: choose if you want to consider only the first/largest maximum ('MAX1') or all maxima ('NOMAX') detected in each time window 
- countflag: to plot wavefield composition decide if you want to consider the number of waves counted ('count') or the number weighted by beam power amplitude ('amp')

Note that, prior to plotting, this script performs essential analysis steps. In the section "Wave type analysis" results are sorted with respect to their detected wave type before they can be plotted accordingly.

-------------------------
REFERENCES & PUBLICATIONS

Finger, C. and Löer, K.: Depth of sudden velocity increases from multi-mode Rayleigh waves derived with three-component ambient noise beamforming, EGU General Assembly 2023, Vienna, Austria, 24–28 Apr 2023, EGU23-12396, https://doi.org/10.5194/egusphere-egu23-12396, 2023.

Löer, K., Finger, C., Obiri, E., and Kennedy, H.: A comprehensive beamforming toolbox to characterise surface and body waves in three-component ambient noise wavefields, EGU General Assembly 2023, Vienna, Austria, 24–28 Apr 2023, EGU23-5670, https://doi.org/10.5194/egusphere-egu23-5670, 2023.

Löer, K., Toledo, T., Norini, G., Zhang, X., Curtis, A., Saenger, E.H.: Imaging the Deep Structures of Los Humeros Geothermal Field, Mexico, Using Three‐Component Seismic Noise Beamforming, Seismological Research Letters (2020) 91 (6): 3269–3277.

Löer, K., Riahi, N., & Saenger, E. H. , Three-component ambient noise beamforming in the Parkfield area, Geophysical Journal International, Volume 213, Issue 3, June 2018, Pages 1478–1491, https://doi.org/10.1093/gji/ggy058

Riahi, N., and Saenger, E. H. (2014), Rayleigh and Love wave anisotropy in Southern California using seismic noise, Geophys. Res. Lett., 41, 363– 369, doi:10.1002/2013GL058518.

Riahi, N., Bokelmann, G., Sala, P., and Saenger, E. H. (2013), Time-lapse analysis of ambient surface wave anisotropy: A three-component array study above an underground gas storage, J. Geophys. Res. Solid Earth, 118, 5339– 5351, doi:10.1002/jgrb.50375.

Wathelet, M., Jongmans, D., Ohrnberger, M., & Bonnefoy-Claudet, S. (2008). Array performances for ambient vibrations on a shallow structure and consequences over V s inversion. Journal of Seismology, 12, 1-19.

