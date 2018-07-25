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
    [pks,locs,width,prominence] = findpeaks(smooth(back_to_zero(:,sweep)),'MinPeakProminence',3.5,'MaxPeakWidth',5000); %then we take the peaks that are at least 4mV local amplitude
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
[x1,x2] = size(peak_times);
for k=1:40
    for j=1:x2
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

recurrent_sweeps= [1,6,7,8,9,10,11,14,16,17,23,24,25,28,29]; %this should be modified depending on the recording = the list of sweeps showing recurrent EPSPs
%back_to_zero = zeros(130000,40);
%for sweep=1:40
  %  back_to_zero(:,sweep) = (electrophysiology(:,sweep) - min(electrophysiology(:,sweep)))*1000; %to put the unit back to mV
%end
figure; hold on
plot(back_to_zero(:,recurrent_sweeps)), title('Corresponding (selected) ephys sweeps'); 
axis tight
xlabel('t (s)')
ylabel('V_{m} (mV)')

%%
c=1;
for i=1:40
for j=1:5
if isnan(peak_heights(i,j))==0
amplitude(c) = peak_heights(i,j); %amplitude is a list of all non-NaN EPSP amplitudes
time(c) = peak_times(i,j); %and time is a list of corresponding times of these very EPSPs
c=c+1;
end
end
end

time_sorted = sort(time); %so then we sort time to get the time distribution of EPSPs
for k=1:length(time)
amplitude_sorted(k) = amplitude(find(time==time_sorted(k))); %but we want amplitude to be sorted accordingly toom so we use the find function
end

%%
clear amplitude_sp
clear time_sp
clear c
c=1;
[x1,x2] = size(peak_heights_sp);
for i=1:x1
    for j=1:x2
        if isnan(peak_heights_sp(i,j))==0
            amplitude_sp(c) = peak_heights_sp(i,j);
            time_sp(c) = peak_times_sp(i,j);
            c=c+1;
        end
   end
end

time_sorted_sp = sort(time_sp);

for k=1:length(time_sp)
    amplitude_sorted_sp(k) = amplitude_sp(find(time_sp==time_sorted_sp(k)));
end

%%
clear categorised_epsps
clear stop
bin_nb = 3000;
stop=bin_nb;

%this part of the code is supposed to categorise epsps based on their time
%of occurence (it packs epsps together if they're in the same time interval)
c=1;
cat=1;
index=1;
for k=time_sorted %we go through the sorted times
    if k<stop %while the times of occurence are smaller than the upper boundary of the first category
        categorised_epsps(c,cat)=amplitude_sorted(index); %the we add the epsp to category 1
        c=c+1;
    else %otherwise it means we reached another category
        while k>stop
            stop = stop + bin_nb; %so until we aren't in the right one, we shift
            cat=cat+1;
            c=1;
        end
        categorised_epsps(c,cat)=amplitude_sorted(index); %once it is reached, we can add the epsp
    end
    index=index+1; %and then we try the next one
end

[x1,x2] = size(categorised_epsps);
time_bins = [1:bin_nb:x2*bin_nb];
means_epsps = mean(categorised_epsps);
figure; hold on
bar(time_bins,means_epsps), xlabel('time category'), ylabel('Mean EPSP amplitude (mV)'), box off

%%
clear categorised_epsps_sp
clear stop
bin_nb = 5000;
time_bins = [1:bin_nb:130000];
stop=bin_nb;


c=1;
cat=1;
index=1;
for k=time_sorted_sp
    if k<stop
        categorised_epsps_sp(c,cat)=amplitude_sorted_sp(index);
        c=c+1;
    else 
        while k>stop
            stop = stop + bin_nb;
            cat=cat+1;
            c=1;
        end
        categorised_epsps_sp(c,cat)=amplitude_sorted(index);
    end
    index=index+1;
end

[x1,x2] = size(categorised_epsps_sp)
time_bins = [1:bin_nb:x2*bin_nb];
means_epsps_sp = mean(categorised_epsps_sp);
bar(time_bins,means_epsps_sp)  
