function [aligned_data] = aligning_traces(datafile,ch_pks,ch_trc)
%This function takes a .abf file corresponding to a recorded pair and
%aligns all the sweeps on a common point. Indeed, testing the connection
%between two cells requires the injection of spikes in one of them, but
%jitter exists and needs to be removed to properly average the response
%INPUT
%
%   data        the name of the abf file, with hyphens
%   ch_pks      the channel number on which the alignment has to be based
%               (voltage trace with the spikes)
%   ch_trc      the channel number of the sweeps that need to be aligned
%
%OUTPUT
%
%   jitter      returns the maximum jitter accumulated during the
%               experiment
%   plots       

%Characterisation of the recording, new matrix preparation
[data, ~ , h] = abfload(datafile);
points = h.sweepLengthInPts;
sweeps = h.lActualEpisodes;

T_sweep = h.si * points;
time = linspace(0,T_sweep,h.si); %a list with time values 




spikes = zeros(sweeps,1);
for k=1:sweeps
    [~,LOC]=findpeaks(data(:,ch_pks,k),'MinPeakHeight',20); %finding the first peak of the k-ieth sweep
    spikes(k)= LOC(1,1); 
end

max_jitter = max(spikes) - min(spikes); %how big is the maximum jitter? 
aligned_data = zeros(points+max_jitter,2,sweeps); % creating a bigger matrix to host all centered sweeps

%Alignement of all sweeps on the spike which time of occurrence is closest
%to 0 (found with min(spikes))
for j=1:sweeps
    start = max(spikes) - spikes(j); %for every sweep, we determine where we start filling the new matrix to align
    for i=1:points
        aligned_data(i+start,1,j) = data(i,ch_trc,j);
    end
end
local_max = max(aligned_data((min(spikes)-200):(min(spikes)+200),1,:)); %finding a local maximum... 
MAX = max(local_max);
figure; hold on
plot(min(spikes),MAX+2,'*') %... such that we can plot a star to indicate where the spike occured
aligned_data = aligned_data((max_jitter+1):(points-1),1,:); 

%a problem is that we have some remaining 0s when data was aligned, since
%every sweep only occupies 'points' entries of the aligned_data matrix.
%Therefore, we curtail the matrix to x:y, where x is the first compulsorily
%filled value for all sweeps (max_jitter+1) and y is the last compulsorily
%filled value for all sweeps (nb of points -1). 



for k=1:sweeps
    plot(smooth(squeeze(aligned_data(:,1,k))',10),'-','color',[.5 .5 .5]);
end
    plot(smooth(mean(squeeze(aligned_data(:,1,:))'),10),'k-','linewidth',2)
    xlabel('t (ms)')
    ylabel('V_{m} (mV)')
    title('Aligned responses of the cell on the first spike of its pair')

disp('The max jitter is') ; disp(max_jitter) ; disp('ms')



end 
