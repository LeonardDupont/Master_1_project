function [shortened_calcium_data,start_point,start_frame,stop_point,stop_frame] = shorten_2s_protocol(calcium_file,electrophysiology_data,h,si,channel)
%The 2s_protocol consists in 2s empty + 2s current injection + 6s empty
%which can generate huge calcium_imaging files. This is a problem for the
%filtering time is exploding. This function takes the original, far too big
%calcium movie (unfiltered) and the corresponding electrophysiology trace
%and crops it from 0.5s before current injection onset to 2s after current
%injection offset. 
%
%   INPUT
%
%   calcium_file      the calcium file number '*...*' in the cd
%
%   electrophy_data   coming from the opening_abf function
%   si
%   h
%
%   channel           the channel of the voltage-clamp with APs from 
%                     the 2s_protocol
%
%
%   OUTPUT
%
%   shortened_ca_data    the array of shortened_calcium_data
%   
%   start_point          the point (ephys) of start, 0.5s before current
%                        onset
%   stop_point           same for stop, 2s after current offset
%
%   stop/start_frame     same principle, but with frames from the imaging





calcium = dir(calcium_file);
calcium_data = imread_universal(calcium.name);

[LOC,~] = findpeaks(electrophysiology_data(:,channel,1), 'MinPeakHeight',20); 
first_spike = LOC(1,1); %corresponds to the beginning of the current injection

%now, we want to keep points up to 0.5s before the beginning and 1s after
%the end of the current injection

start_sweep = first_spike - (0.5/si)*1e6;
stop_sweep = first_spike + (3/si)*1e6; %3 and not 2 because the current injection itself spans 2s

[fires, ~] = fire_detection(1,electrophysiology_data,h);
%now we need to find the closest fire corresponding to these two points

start_point = 1e100;
stop_point = 1e100;
for k=1:length(fires) %we go through the list and find the minimal difference between the two points and a fire point
    diff_start = abs(fires(k)-start_sweep);
    diff_stop = abs(fires(k)-stop_sweep);
    if diff_start < start_point 
        start_point = fires(k);
        start_frame = k; %we also store the frame number
    end
    if diff_stop < stop_point
        stop_point = fires(k);
        stop_frame = k;
    end
end




shortened_calcium_data = zeros(512,512,(stop_frame-start_frame)*h.lActualEpisodes);
for k=1:(stop_frame-start_frame)*h.lActualEpisodes
    shortened_calcium_data(:,:,k) = calcium_data(:,:,k);
end

prompt = 'Would you like to write a shortened .tiff file? (0 or 1)'
write = input(prompt);

if write %then we write a shortened tiff file
    cropped_name = strsplit(calcium.name,'.sif');
    full_name = strcat(cropped_name, '_shortened_2s_file','.tiff');
    util.imwritetiff(shortened_calcium_data,full_name,'close',true)
end


end

