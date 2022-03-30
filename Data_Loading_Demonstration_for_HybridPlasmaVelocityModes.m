%% Preamble:
%This MATLAB program reads in patterns of dominant spatial and temporal
% horizontal ionospheric plasma velocity variation covering the northern
% polar region, for user-speficied calendar months within the range 1997.0
% to 2009.0. the program demonstrates how to combine the two models
% provided within each month, in order to give a continuous reanalysis
% record of ionospheric plasma velocity variation.
%
%The data are described in the associated README, as well as in the paper
% 'A reanalysis of SuperDARN plasma velocity variability for solar cycle
% 23', by Shore et al., intended for publication in JGR Space Physics.
%
%Author: Robert Shore, March 2022.
% Email: robore@bas.ac.uk, ORCID: orcid.org/0000-0002-8386-1425.
%
%Point of contact for dataset: Mervyn Freeman, mpf@bas.ac.uk.

%% Options, to be set by the user: 

%State which directory the plasma velocity mode netcdf files are stored in:
file_directory = '/users/reader/specify_your_own_data_storage_location/';

%Specify year and month of the monthly data file to load in.  These could
% be iterated over to cover all 144 available months.
file_year = 2001;%1997 to 2008 inclusive.
file_month = 2;%months 1 to 12.

%% Load spatial and temporal positional data.
disp('Loading spatial and temporal positional data...')

%Load bin location data: the centroids are the locations of the centre of
% each bin, and the limits of the bins give the region over which the
% EOF prediction at the centroid is assumed to apply. There are 559 bins
% in the northern polar region which is the area of focus for this
% analysis.  They are ordered approximately by latitude, then longitude,
% but this does not always apply near the 0/360 degree longitude boundary.
bin_centroids_colatitude = ncread([file_directory 'HybridPlasmaVelocityModes-ds01.nc'],'bin_centroids_colatitude');%size [559 (i.e. number of bins) by 1].
bin_centroids_longitude =  ncread([file_directory 'HybridPlasmaVelocityModes-ds01.nc'],'bin_centroids_longitude');%size [559 (i.e. number of bins) by 1].
bin_limits_colatitude =    ncread([file_directory 'HybridPlasmaVelocityModes-ds01.nc'],'bin_limits_colatitude');%size [559 (i.e. number of bins) by 2].
bin_limits_longitude =     ncread([file_directory 'HybridPlasmaVelocityModes-ds01.nc'],'bin_limits_longitude');%size [559 (i.e. number of bins) by 2].

%Load bin times for the specified month. These are the means of the
% 5-minute spans used to capture data for each bin.
% Column format: year, month, day, hour, minute, second.
bin_times = ncread([file_directory 'HybridPlasmaVelocityModes-ds02.nc'], ['bin_times_' num2str(file_year) num2str(file_month,'%2.2d')]);%size [(number of 5-minute epochs in month) by 6].

disp('... data loaded.')

%% Load the 'background mean' field: the temporal mean for each spatial bin location.
disp('Loading background mean field data...')

%Load monthly mean field at bin locations:
bin_means_north = ncread([file_directory 'HybridPlasmaVelocityModes-ds03.nc'], ['bin_means_' num2str(file_year) num2str(file_month,'%2.2d') '_north']);%size [559 (i.e. number of bins) by 1].
bin_means_east = ncread([file_directory 'HybridPlasmaVelocityModes-ds03.nc'], ['bin_means_' num2str(file_year) num2str(file_month,'%2.2d') '_east']);%size [559 (i.e. number of bins) by 1].

disp('... data loaded.')

%% Load the spatial and temporal eigenvectors for 'model 1', and reconstruct that model for each location and time in the specified month:
disp('Loading model 1 data...')

%Load spatial and temporal eigenvectors for the specified month:
%Preallocate storage:
model_1_spatial_eigenvectors_north_component =  NaN(size(bin_centroids_colatitude,1),10);%size [559 (i.e. number of bins) by 10 (number of modes)].
model_1_spatial_eigenvectors_east_component =  NaN(size(bin_centroids_colatitude,1),10);%size [559 (i.e. number of bins) by 10 (number of modes)].
model_1_temporal_eigenvectors = NaN(size(bin_times,1),10);%size [(number of 5-minute epochs in month) by 10 (number of modes)].
for mode_index = 1:10
    %Load spatial eigenvector components:
    model_1_spatial_eigenvectors_north_component(:,mode_index) = ncread([file_directory 'HybridPlasmaVelocityModes-ds04.nc'], ['model_1_eig_s_' num2str(file_year) num2str(file_month,'%2.2d') '_mode' num2str(mode_index,'%2.2d') '_north']);%size of extraction: [559 (i.e. number of bins) by 1].
    model_1_spatial_eigenvectors_east_component(:,mode_index) = ncread([file_directory 'HybridPlasmaVelocityModes-ds04.nc'], ['model_1_eig_s_' num2str(file_year) num2str(file_month,'%2.2d') '_mode' num2str(mode_index,'%2.2d') '_east']);%size of extraction: [559 (i.e. number of bins) by 1].
    
    %Load temporal eigenvector:
    model_1_temporal_eigenvectors(:,mode_index) = ncread([file_directory 'HybridPlasmaVelocityModes-ds06.nc'], ['model_1_eig_t_' num2str(file_year) num2str(file_month,'%2.2d') '_mode' num2str(mode_index,'%2.2d')]);%size of extraction: [(number of 5-minute epochs in month) by 1].
end%Loop over each mode, store corresponding eigenvectors.

%Reconstruct the plasma velocity matrix for all time and locations,
% and each component.  This sums all modes:
model_1_plasma_velocity_north_component = model_1_temporal_eigenvectors * model_1_spatial_eigenvectors_north_component';%size of result: [(number of 5-minute epochs in month) by (number of bins)].
model_1_plasma_velocity_east_component = model_1_temporal_eigenvectors * model_1_spatial_eigenvectors_east_component';%size of result: [(number of 5-minute epochs in month) by (number of bins)].

disp('... data loaded.')

%% Load the spatial and temporal eigenvectors for 'model 2', and reconstruct that model for each location and time in the specified month:
disp('Loading model 2 data...')

%Load spatial and temporal eigenvectors for the specified month:
%Preallocate storage:
model_2_spatial_eigenvectors_north_component =  NaN(size(bin_centroids_colatitude,1),10);%size [559 (i.e. number of bins) by 10 (number of modes)].
model_2_spatial_eigenvectors_east_component =  NaN(size(bin_centroids_colatitude,1),10);%size [559 (i.e. number of bins) by 10 (number of modes)].
model_2_temporal_eigenvectors = NaN(size(bin_times,1),10);%size [(number of 5-minute epochs in month) by 10 (number of modes)].
for mode_index = 1:10
    %Load spatial eigenvector components:
    model_2_spatial_eigenvectors_north_component(:,mode_index) = ncread([file_directory 'HybridPlasmaVelocityModes-ds05.nc'], ['model_2_eig_s_' num2str(file_year) num2str(file_month,'%2.2d') '_mode' num2str(mode_index,'%2.2d') '_north']);%size of extraction: [559 (i.e. number of bins) by 1].
    model_2_spatial_eigenvectors_east_component(:,mode_index) = ncread([file_directory 'HybridPlasmaVelocityModes-ds05.nc'], ['model_2_eig_s_' num2str(file_year) num2str(file_month,'%2.2d') '_mode' num2str(mode_index,'%2.2d') '_east']);%size of extraction: [559 (i.e. number of bins) by 1].
    
    %Load temporal eigenvector:
    model_2_temporal_eigenvectors(:,mode_index) = ncread([file_directory 'HybridPlasmaVelocityModes-ds07.nc'], ['model_2_eig_t_' num2str(file_year) num2str(file_month,'%2.2d') '_mode' num2str(mode_index,'%2.2d')]);%size of extraction: [(number of 5-minute epochs in month) by 1].
end%Loop over each mode, store corresponding eigenvectors.

%Reconstruct the plasma velocity matrix for all time and locations,
% and each component. This sums all modes:
model_2_plasma_velocity_north_component = model_2_temporal_eigenvectors * model_2_spatial_eigenvectors_north_component';%size of result: [(number of 5-minute epochs in month) by (number of bins)].
model_2_plasma_velocity_east_component = model_2_temporal_eigenvectors * model_2_spatial_eigenvectors_east_component';%size of result: [(number of 5-minute epochs in month) by (number of bins)].

disp('... data loaded.')

%% Load the hybrid model choice values for this month, and use them to construct the hybrid model infill values:
disp('Constructing hybrid model...')

hybrid_model_choice_values = ncread([file_directory 'HybridPlasmaVelocityModes-ds08.nc'], ['hybrid_model_choice_' num2str(file_year) num2str(file_month,'%2.2d')]);%size of extraction: [559 (i.e. number of bins) by 1].

%Preallocate storage:
hybrid_model_plasma_velocity_north_component = NaN(size(model_1_plasma_velocity_north_component));%size [(number of 5-minute epochs in month) by (number of bins)].
hybrid_model_plasma_velocity_east_component = NaN(size(model_1_plasma_velocity_north_component));%size [(number of 5-minute epochs in month) by (number of bins)].
%Loop over bins, and assign the model values based on the model choice
% value:
for i_bin = 1:559
    if(hybrid_model_choice_values(i_bin) == 1)
        hybrid_model_plasma_velocity_north_component(:,i_bin) = model_1_plasma_velocity_north_component(:,i_bin);
        hybrid_model_plasma_velocity_east_component(:,i_bin) = model_1_plasma_velocity_east_component(:,i_bin);
    elseif(hybrid_model_choice_values(i_bin) == 2)
        hybrid_model_plasma_velocity_north_component(:,i_bin) = model_2_plasma_velocity_north_component(:,i_bin);
        hybrid_model_plasma_velocity_east_component(:,i_bin) = model_2_plasma_velocity_east_component(:,i_bin);
    end%Conditional: choose model values based on earlier computations of prediction efficiency.
end%Loop over each spatial bin.


%Add the background mean which was removed prior to performing the EOF
% analysis:
hybrid_model_plasma_velocity_north_component = hybrid_model_plasma_velocity_north_component + repmat(bin_means_north',[size(bin_times,1) 1]);%size [(number of 5-minute epochs in month) by (number of bins)].
hybrid_model_plasma_velocity_east_component = hybrid_model_plasma_velocity_east_component + repmat(bin_means_east',[size(bin_times,1) 1]);%size [(number of 5-minute epochs in month) by (number of bins)].

%The rows of hybrid_model_plasma_velocity_*_component now give a plasma
% velocity prediction at the locations given by the vectors
% bin_centroids_colatitude and bin_centroids_longitude, for the times given
% by the vector bin_times.  The units are metres per second.  Should the
% user wish to display the limits of the cell over which the original data
% were captured, refer to the first and second columns (again, of the
% corresponding rows) of the matrices bin_limits_colatitude and
% bin_limits_longitude.  Alternatively, the centroid values can be linearly
% interpolated between in order to give a prediction at any point inside
% the region of coverage. 

disp('... model constructed. Program complete.')
