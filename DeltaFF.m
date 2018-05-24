function [DeltaFF_data] = DeltaFF(region,sweep,interpolated_calcium_data,electrophysiology_data,channel,si,start_sweep,stop_sweep)
% This function takes a raw intensity and gives back the deltaF over F
% values. It uses the points before the spiking time and is hence to
% be applied on calcium-imaging data only. It should also be applied on
%interpolated calcium data to ease things a bit.
%
%   INPUT 
%   
%   region      the ROI in question0
%
%   sweep       the fluorscence sweep in question
%
%   interpol... the interpolated calcium data, 3D matrix created before
%
%   electroph.. the electrophysiology data matrix imported earlier w/
%               opening abf
%   
%   channel     also defined earlier in the script; the channel w/ spikes
%
%   si          the sampling interval in us
%
%   OUTPUT
%   The function returns a vector with deltaFF values that can be stored in
%   a new matrix later on


[~,LOC]=findpeaks(electrophysiology_data(start_sweep:stop_sweep,channel,sweep),'MinPeakProminence',12); %finding the first peak of the k-ieth sweep
spike_1 = LOC(1,1); %this is the location of the first spike, we can take points before

F0 = mean(interpolated_calcium_data(1:(spike_1-(10^4)/si),sweep,region)); % we take all points from the beginning up to 10ms before the first spike

DeltaFF_data = zeros(length(interpolated_calcium_data),1);
for k=1:length(interpolated_calcium_data)
    DeltaFF_data(k) = (interpolated_calcium_data(k,sweep,region) - F0)/F0; %we now calculate the actual \Delta F/F values
end



    