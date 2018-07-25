% This script was used to analyse the set of sweeps from 2018-06-27
% (recording 25) showing a cell with little recurrent inhibition but
% showing strong, well-timed EPSPs after a 2s current injection

% First we load and prepare the data
electrophy = load('20180627_rec25.mat');
for k=1:40
    electrophysiology(:,k) = electrophy.(['Trace_1_25_',num2str(k),'_1'])(:,2);
end

%%


% Then we go through the sweep and detect peaks (id est EPSPs outside of the 2s current injection)
peak_times = zeros(40,5); %maximum 5 EPSPs detected in response here
peak_heights = zeros(40,5);
for sweep=1:40
    [pks,locs,width,prominence] = findpeaks(electrophysiology(:,sweep),'MinPeakProminence',4); %then we take the peaks that are at least 4mV local amplitude
    c=1;
    for peak=1:length(pks)
        if locs(peak)<20000 || locs(peak)>40011
            peak_times(sweep,c) = locs(peak);
            peak_heights(sweep,c) = prominence(peak);
            c=c+1;
        end
    end
end

%%

%To remove the zeros
for k=1:40
    for j=1:5
        if peak_times(k,j)==0
        peak_times(k,j)=NaN;
        end
        if peak_heights(k,j)==0
        peak_heights(k,j)=NaN;
        end
    end
end
%%

%Now we can make means and histograms without too many problems. 
for k=1:40
    figure
    plot(electrophysiology(:,k))
    title(num2str(k))  %also plotting all sweeps in different figures quickly to check if some detected EPSPs are not fake
end 

%%
figure; histogram(peak_heights,50), title('Histogram of EPSPs amplitude')
figure; histogram(peak_times,50), title('Histogram of EPSPs time')

recurrent_sweeps= [1,6,7,8,9,10,11,14,16,17,23,24,25,29,31,33,34];
back_to_zero = zeros(130000,40);
for sweep=1:40
    back_to_zero(:,sweep) = electrophysiology(:,sweep) - min(electrophysiology(:,sweep));
end
figure
plot(back_to_zero(:,recurrent_sweeps)), title('Corresponding (selected) ephys sweeps')