electrophy = load('20180705.mat');
%%
%in this file, we are interested in the following series:
%   - 19; spontaneous before any stimulation (actually just IV and adaptation)
%   - 21; 2s protocole
%   - 22; spontaneous after the 2s stimulation

sp_before = zeros(130000,15);
stimulation = zeros(130000,40);
sp_after = zeros(130000,40);
for k=1:15
    sp_before(:,k) = electrophy.data.series(19).trace(1).data(:,k);
end
for k=1:40
    stimulation(:,k) = electrophy.data.series(21).trace(1).data(:,k);
    sp_after(:,k) = electrophy.data.series(22).trace(1).data(:,k);
end
%%
%% this can be useful to plot all sweeps from a series
%sp_before
for k=1:length(sp_before)
    figure, hold on
    plot(smooth(sp_before(:,k))), title(num2str(k));
    hold off
end

%%
%stimulation
for k=1:length(stimulation)
    figure, hold on
    plot(smooth(stimulation(:,k))), title(num2str(k));
    hold off
end

%%
%sp_after
for k=1:length(sp_after)
    figure, hold on
    plot(smooth(sp_after(:,k))), title(num2str(k));
    hold off
end

%%
%% Now we want to start to detect EPSPs, take their time and amplitude.

peak_times_before = zeros(15,5); %maximum 5 EPSPs detected in response here
peak_heights_before = zeros(15,5);
peak_times_stim = zeros(40,5);
peak_heights_stim = zeros(40,5);
peak_times_after = zeros(40,5);
peak_heights_after = zeros(40,5);

%%
back_to_zero_stim = zeros(130000,40);
for k=1:40
    back_to_zero_stim(:,k)= stimulation(:,k) - min(stimulation(:,k));
end
%% Here we are only going to take EPSPs > 3.5mV into account; quite restrictive
%before
for sweep=1:15
    [pks,locs,width,prominence] = findpeaks(smooth(sp_before(:,sweep)),'MinPeakProminence',3.5,'MaxPeakWidth',5000); %then we take the peaks that are at least 4mV local amplitude
    c=1;
    for peak=1:length(pks)
        peak_times_before(sweep,c) = locs(peak);
        peak_heights_before(sweep,c) = prominence(peak);
        c=c+1;
    end
end
%%
%stim
for sweep=1:40
    [pks,locs,width,prominence] = findpeaks(smooth(back_to_zero_stim(:,sweep)),'MinPeakProminence',3.5,'MaxPeakWidth',5000); %then we take the peaks that are at least 4mV local amplitude
    c=1;
    for peak=1:length(pks) %this part allows us to remove the EPSPs detected between 2 and 4s: these are actually elicited spikes from the current injection
        if locs(peak)<20000 || locs(peak)>40011
            peak_times_stim(sweep,c) = locs(peak);
            peak_heights_stim(sweep,c) = prominence(peak);
            c=c+1;
        end
    end
end
%%
%after
for sweep=1:40
    [pks,locs,width,prominence] = findpeaks(smooth(sp_after(:,sweep)),'MinPeakProminence',3.5,'MaxPeakWidth',5000); %then we take the peaks that are at least 4mV local amplitude
    c=1;
    for peak=1:length(pks)
        peak_times_after(sweep,c) = locs(peak);
        peak_heights_after(sweep,c) = prominence(peak);
        c=c+1;
    end
end
%%
%Now we want to remove the 0s from our tables
[x1,x2] = size(peak_times_stim);
[x3,x4] = size(peak_times_before);
[x5,x6] = size(peak_times_after);
for k=1:x3
    for j=1:x4
        if peak_times_before(k,j)==0
        peak_times_before(k,j)=NaN;
        end
        if peak_heights_before(k,j)==0
        peak_heights_before(k,j)=NaN;
        end
    end
end
for k=1:x1
    for j=1:x2
        if peak_times_stim(k,j)==0
        peak_times_stim(k,j)=NaN;
        end
        if peak_heights_stim(k,j)==0
        peak_heights_stim(k,j)=NaN;
        end
    end
end
for k=1:x5
    for j=1:x6
        if peak_times_after(k,j)==0
        peak_times_after(k,j)=NaN;
        end
        if peak_heights_after(k,j)==0
        peak_heights_after(k,j)=NaN;
        end
    end
end

%% Now we can start doing some real analysis
% We'll start of with the most simple feat of data: the number of EPSP per
% trial (close to the notion of probability)
clear mag_before
clear mag_stim
clear mag_after
clear c

nb_before = zeros(x3,1); %this list will store the number of epsps per trial

nb_stim = zeros(x1,1);

nb_after = zeros(x5,1);


%before
c=1;
for k=1:x3
    nb_before(k) = length(peak_times_before(k,:)) - sum(isnan(peak_times_before(k,:)));
    if nb_before(k)~=0
        mag_before(c) = nanmean(peak_heights_before(k,:));
        c=c+1;
    end
end
c=1;
for k=1:x1
    nb_stim(k) = length(peak_times_stim(k,:)) - sum(isnan(peak_times_stim(k,:)));
    if nb_stim(k)~=0
        mag_stim(c) = nanmean(peak_heights_stim(k,:));
        c=c+1;
    end
end
c=1;
for k=1:x5
    nb_after(k) = length(peak_times_after(k,:)) - sum(isnan(peak_times_after(k,:)));
    if nb_after(k)~=0
        mag_after(c) = nanmean(peak_heights_after(k,:));
        c=c+1;
    end
end
%stim

mean_nbs = [mean(nb_before), mean(nb_stim), mean(nb_after)];
std_nbs = [std(nb_before),std(nb_stim),std(nb_after)];
mean_mags = [mean(mag_before),mean(mag_stim),mean(mag_after)];
std_mags = [std(mag_before),std(mag_stim),std(mag_after)];

%% amplitude only 
b = barwitherr(std_mags,mean_mags,'FaceColor','flat'); 

b.CData(1,:) = [0.8 0.898 1];
b.CData(2,:) = [.4 .4 1];
b.CData(3,:) = [0 .4 .8];
xticks([1,2,3])
set(gca,'XTickLabel',{'spontaneous <','stim','spontaneous >'})
ylabel('EPSP amplitude (mV)')
box off

p = ranksum(mag_after,mag_stim) %Mann Whitney U 
[h,p] = ttest2(mag_after,mag_stim) %ttest 2 samples 

%% number of EPSPs

b = barwitherr(std_nbs,mean_nbs,'FaceColor','flat');

b.CData(1,:) = [0.8 0.898 1];
b.CData(2,:) = [.4 .4 1];
b.CData(3,:) = [0 .4 .8];
xticks([1,2,3])
set(gca,'XTickLabel',{'spontaneous <','stim','spontaneous >'})
ylabel('Number of EPSPs (>3.5mV) per sweep')
box off

%%
%% Now we might want to categorise our EPSPs in time categories to see if there is a time-distribution, or even a time distribution of amplitude...
clear c
clear amplitude
clear time
c=1;
for i=1:40
    for j=1:5
      if isnan(peak_heights_stim(i,j))==0
        amplitude(c) = peak_heights_stim(i,j); %amplitude is a list of all non-NaN EPSP amplitudes
         time(c) = peak_times_stim(i,j); %and time is a list of corresponding times of these very EPSPs
          c=c+1;
      end
    end
end

time_sorted = sort(time); %so then we sort time to get the time distribution of EPSPs
for k=1:length(time)
    amplitude_sorted(k) = amplitude(find(time==time_sorted(k))); %but we want amplitude to be sorted accordingly toom so we use the find function
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

%%
[x1,x2] = size(categorised_epsps);
time_bins_2 = [1:bin_nb:x2*bin_nb];
means_epsps = mean(categorised_epsps);
figure; hold on
bar(time_bins_2,means_epsps), xlabel('time category'), ylabel('Mean EPSP amplitude (mV)'), box off

%%
%% plotting all the sweeps showing recurrent EPSPs during the stimulation protocole
selected_sweeps = [2,3,4,9,11,12,13,14,15,22,23,24,27,28,31,33,36,38];
figure; hold on 
box off
for k = selected_sweeps
    plot(smooth(back_to_zero_stim(:,k)))
end
xlabel('t (s)')
ylabel('V_{m} (mV)')
axis tight

%% Same but with the after to see if there is a peculiar time distribution

selected_sweeps_after = [1,2,3,4,5,8,11,17,18,19,20,21,22,24,25,30,31];
figure; hold on 
box off
for k = selected_sweeps_after
    plot(smooth(sp_after(:,k)))
end
xlabel('t (s)')
ylabel('V_{m} (mV)')
axis tight


%%
clear stop
clear counter
clear cat
clear time_cat_stim

stop = 5000;
time_bins = [1:stop:130000]; 
time_cat_stim = zeros(130000/stop,1);
counter=0;
cat=1;

for k=time_sorted
    if k<stop
        counter=counter+1;
    else
        while k>stop
            time_cat_stim(cat)=counter;
            stop=stop+5000;
            cat=cat+1;
            counter=0;
        end
        counter=counter+1;
    end
end

 bar(time_bins,time_cat_stim,'FaceColor',[.8 .8 1])      
            
%% same but with a subplot and the time distribution under it (stimulation )

axh=[]; figure, hold on
axh(1)=subplot(26,1,1:15); hold on
selected_sweeps = [2,3,4,9,11,12,13,14,15,22,23,24,27,28,31,33,36,38];
for k = selected_sweeps
    plot(smooth(back_to_zero_stim(:,k)))
    axis tight
end
set(gca,'TickLength',[0,0])
xlabel('t (s)')
ylabel('V_{m} (mV)')


axh(2)=subplot(26,1,16:20); hold on
bar(time_bins,time_cat_stim,'FaceColor',[.8 .8 1]) 
ylabel('Number of EPSPs')
axis tight


axh(3)=subplot(26,1,21:26); hold on
bar(time_bins_2,means_epsps,'FaceColor',[.8 0 .4])
ylabel('EPSP amplitude')
axis tight

linkaxes(axh,'x')



