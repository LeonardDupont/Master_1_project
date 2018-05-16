%
% This function is useful to align calcium data with electrophysiology
% data. It uses other functions to extract only fluorescence parts of the
% calcium signal and align it on the Vm one of the patched cell using the
% fire input. The electrophysiology data should have been opened with
% 'opening_abf.m' beforehand. 
%
%   INPUT
%
%   electrophysiology_data    The electrophsyiology file containing the
%                             trace to align data on
%
%   channel_number            The channel number to be used for the plot
%                             and alignment
%   
%   calcium_file              The raw calcium file number, '*...*' 
%
%   ephys_sweeps              The number of sweeps in the ephys file, h.lActualEpisodes
%
%
%   This programme notably calls in zprofile for the ROIs. 



%importing the calcium file

calcium_data = dir(calcium_file);
%filtered_data = filter_hilolike_from_file_par(calcium_data.name);
filtered_data = calcium_data.name;

[stats,datasum]= zprofile(filtered_data,'fps',18.38); %here we draw ROIs that we'll be using
%frame_time = linspace(1,length(stats{1}.time,1)); %generating a time-vector based on frames
[loc_ONSET,loc_OFFSET] = fluorescence_variations(calcium_file);

individual_fluorescence_sweeps = zeros(loc_OFFSET(1)-loc_ONSET(1)-1,length(loc_ONSET),length(stats));
%here we create a big 3D matrix with D1 = points (value in frame) minus 2 because we start at loc_ONSET +1 and stop at loc_OFFSET -1, D2 =
%sweeps, D3 = ROIs


for region=1:length(stats) %let us now fill in this matrix
    for sweep=1:length(loc_ONSET)
        c = 1;
        for value=(loc_ONSET(sweep)+1):(loc_OFFSET(sweep)-1) %we crop the region a bit to be sure
            individual_fluorescence_sweeps(c,sweep,region) = stats{region}.intensraw(value); %we are storing the value of the raw intensity
            c = c+1;
        end
    end
end

[fires, ~] = fire_detection(1,electrophysiology_data,h);

start_sweep = fires(loc_ONSET(1)+1) ;
stop_sweep = fires(loc_OFFSET(1)-1) ; %we want to plot traces (electrophy and calcium) which span the same time duration, 
%so we need to find - based on the frames - from which point to start and
%at which point to stop in the ephys trace. We still crop it a bit to be
%sure not to have weird fluorescence traces

%We now need to interpolate the calcium data to be able to plot it against
%the electrophysiology one

interpolated_calcium_data = zeros(stop_sweep-start_sweep+1,length(loc_ONSET),length(stats)); %This matrix is going to replace the individual fluorescence one. Our goal is to generate missing points between
%fire points in the electrophysiology dataset which are the only locations to which a calcium value is currently attributed. 
current_fire_points = fires(loc_ONSET(1)+1:loc_OFFSET(1)-1); %These are the current electrophysiology points for which there already exists a fluorescence value
ultimate_resolution = linspace(0,stop_sweep-start_sweep+1,stop_sweep-start_sweep+1); %This is the list of points for which we eventually need a value

for region=1:length(stats)
    for sweep=1:length(loc_ONSET)
        interpolated_calcium_data(:,sweep,region) = interp1(current_fire_points,individual_fluorescence_sweeps(:,sweep,region),ultimate_resolution); %Hence we interpolate
    end
end


spikes = zeros(ephys_sweeps,1); %we create a list that will contain the location of the first spike for each
%sweep, it is just for plotting purposes
for k=1:ephys_sweeps
    [~,LOC]=findpeaks(electrophysiology_data(:,channel,k),'MinPeakHeight',20); %finding the first peak of the k-ieth sweep
    spikes(k)= LOC(1,1); 
end



for region=1:length(stats) %first we choose a ROI
    
    y1 = min(electrophysiology_data(:,channel,1));
    y2 = max(interpolated_calcium_data(:,1,region));
    diff = y2 - y1;
    diffa = abs(diff);
    if diff == diffa
     shift = 30 + diff;
    else
     shift = 30 - diff; %with this small calculation we are able to align the y axes of both traces, then to shift the calcium one by 30 units... (*)
    end
    
    figure('Name',num2str(region)); hold on %creating a figure specifically for this region and naming it correctly
    plot(squeeze(electrophysiology_data(start_sweep:stop_sweep,channel,1))) %We start by the ephys trace, sweep1: it is cropped to the fluorescence_ON period
    for n=1:length(loc_ONSET) %and then plot all fluorescence traces recorded in this ROI while light was on
        plot(interpolated_calcium_data(:,n,region)-(shift+(n-1)*30)) %...(*) and then each time we plot a new sweep, it is shifted by 30 units down again
        star_position = max(interpolated_calcium_data(:,n,region))- min(interpolated_calcium_data(:,n,region));
        plot(spikes(n),star_position - (shift+(n-1)*30), '*') %we add a start at the location of the spiking onset, it has to be shifted along the sweeps
    end
    axis tight
end


 %in addition we find the location of the first spike of the first sweep

amplitude_deltaF = zeros(length(stats),1); %we want to calculate the variation amplitude of fluorescence for each ROI, we still perform this calculation on the non-interpolated values 

for region=1:length(stats)
    amplitude_deltaF(region) = max(mean(individual_fluorescence_sweeps(:,:,region))) - min(mean(individual_fluorescence_sweeps(:,:,region)));
%... so we do so using the mean, store everything in a list
end

red = [1, 0, 0]; %now here we go a bit fancy, we create a color gradient based on amplitude, defining pink (min) and red (max
pink = [255, 192, 203]/255;
colors_p = [linspace(red(1),pink(1),length(stats))', linspace(red(2),pink(2),length(stats))', linspace(red(3),pink(3),length(stats))'];

region_ranking = zeros(length(stats),1);
for region=1:length(stats)
    mini = min(amplitude_deltaF);
    mini_r = find(amplitude_deltaF==mini);
    region_ranking(region)=mini_r;
    amplitude_deltaF(mini_r)=10e6;
end

first_spike = spikes(1);

axh=[]; figure('Name','Average traces'), hold on  %This time we create one only figure for all ROIs with average traces 
axh(1) = subplot(5*length(stats),1,1);
plot(squeeze(electrophysiology_data(start_sweep:stop_sweep,channel,1)))
for region=1:length(stats)
    x= (region-1)*5;
    if x==0
        x=1;
    end
    axh(2) = subplot(x:region*5,1,2);
    for sweep=1:length(loc_ONSET)
        plot(individual_fluorescence_sweeps(:,sweep,region),'Color',[.5 .5 .5])
    end
    plot(squeeze(mean(individual_fluorescence_sweeps(:,:,region),2)),'r','Linewidth',2)
    star_position = max(mean(individual_fluorescence_sweeps(:,:,region),2))- min(mean(individual_fluorescence_sweeps(:,:,region),2));
    plot(first_spike,star_position), '*')
end
linkaxes(axh,x)



figure
imagesc(datasum)
axis equal
hold on

for region=1:length(stats)
    plot(stats{region}.boundary{1}(:,2),stats{region}.boundary{1}(:,1),'Color',colors_p(find(region_ranking==region),:))
end



    
    






