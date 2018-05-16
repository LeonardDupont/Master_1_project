function [fires, fire_number] = fire_detection(sweep,data,h)
% This function takes a sweep corresponding to the fire input of the Andor
% camera and detects the number and mainly position of fires within the
% sweep. This should be useful to align calcium imaging and
% electrophysiology data. One should open the electrophysiology file with
% opening_abf.m for this programme to run properly.
%
%
% INPUT
%
%   sweep       The sweep number to be analysed (1 to h.lActualEpisodes)
%   data, h     Both coming from [data,si,h]=opening_abf('*...*')
%
% OUTPUT
%
%   fires       Is a vector containing the index of fire starts
%               (discontinuous region, ramp up)
%
%   fire_number The number of 'fire' detected
%
%At the end, a plot can be realised with stars above each fire to verify
%that they all were detected
%

indexed_fire_starts = zeros(h.sweepLengthInPts,1); %creating a too big vector
for k=1:(h.sweepLengthInPts-1)
    detection = data(k+1,5,sweep) - data(k,5,sweep);
    if detection>1.7
        indexed_fire_starts(k) = 1;
    end
end
fires = find(indexed_fire_starts==1);
fire_number = length(fires);
figure; hold on
plot(data(:,5,3))
for j=1:(h.sweepLengthInPts-1)
    if indexed_fire_starts(j)==1
        plot(j+1,data(j+1,5,sweep)+0.5,'*')
    end
end


    