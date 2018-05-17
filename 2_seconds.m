[shortened_calcium_data,start_point,start_frame,stop_point,stop_frame] = shorten_2s_protocol(calcium_file,electrophysiology_data);

[stats,datasum] = zprofile();


[fires, fire_number] = fire_detection(1,electrophysiology_data,h)

loc_ONSET = zeros(h.lActualEpisodes,1);
loc_OFFSET = zeros(h.lActualEpisodes,1);
for k=1:length(loc_ONSET)
    if k==1:
        loc_ONSET(k)=1;
        loc_OFFSET(k) = k*(stop_frame-start_frame+1);
    else
        loc_ONSET(k) = (k-1)*(stop_frame-start_frame+1); %the next calcium sweep is actually just a shift of the window by fire_number
        loc_OFFSET(k) = k*(stop_frame-start_frame+1);
    end
end

individual_fluorescence_sweeps = zeros(start_frame - stop_frame,length(loc_ONSET),length(stats));
%here we create a big 3D matrix with D1 = points (value in frame) minus 2 because we start at loc_ONSET +1 and stop at loc_OFFSET -1, D2 =
%sweeps, D3 = ROIs
for region=1:length(stats) %let us now fill in this matrix
    for sweep=1:length(loc_ONSET)
        c=1;
        for value=loc_ONSET(sweep):loc_OFFSET(sweep) %we crop the region a bit to be sure
            individual_fluorescence_sweeps(c,sweep,region) = stats{region}.intensraw(value); %we are storing the value of the raw intensity
            c=c+1;
        end
    end
end

start_sweep = start_point;
stop_sweep = stop_point; %just for notation simplicity purposes

%now we can interpolate the calcium data with the electrophy one
interpolated_calcium_data = zeros(stop_sweep-start_sweep+1,length(loc_ONSET),length(stats)); %This matrix is going to replace the individual fluorescence one. Our goal is to generate missing points between
%fire points in the electrophysiology dataset which are the only locations to which a calcium value is currently attributed. 
nb_fire_points = stop_frame-start_frame+1; %These are the current electrophysiology points for which there already exists a fluorescence value
current_fire_points = linspace(0,stop_sweep-start_sweep+1,nb_fire_points);
ultimate_resolution = linspace(0,stop_sweep-start_sweep+1,stop_sweep-start_sweep+1); %This is the list of points for which we eventually need a value

for region=1:length(stats)
    for sweep=1:length(loc_ONSET)
        interpolated_calcium_data(:,sweep,region) = interp1(current_fire_points,individual_fluorescence_sweeps(:,sweep,region),ultimate_resolution); %Hence we interpolate
    end
end
