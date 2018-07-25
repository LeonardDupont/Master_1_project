fluorescence_2 = csvread('BrightnessOverTime-06202018-1340-066_Cycle00001-botData.csv',1,0);
[x1,x2] = size(fluorescence_2);

Maxes = zeros(x2,1);
Mines = zeros(x2,1);
figure; hold on
title('Fluorescence traces of a few defined ROIs')
ylabel('Intensity')
xlabel('t (s)')
h = 1:x2-1;
for region=1:x2-1
    Back_to_zero = fluorescence_2(:,region+1) - min(fluorescence_2(:,region+1));
    h(region)=plot(fluorescence_2(:,1),Back_to_zero+sum(Maxes),'Linewidth',1.5);
    Maxes(region) = max(Back_to_zero);
end
axis tight
ROIs = {};
for k=1:x2-1
    ROIs{k} = strcat('ROI',num2str(k));
end
legend(h,ROIs)
%% concatenating the ephys sweeps

concatenated_ephys = zeros(1300000,1);
waiting=zeros(1300000,1);
start = 1;
stop = 100000;
for sweep=1:10
concatenated_ephys(start:stop) = ephys_data(:,sweep,2);
waiting(start:stop)=NaN;
last = ephys_data(end,sweep,2);
for k=1:30000
waiting(stop+k) = last;
concatenated_ephys(stop+k)=NaN;
end
start = start + 130000;
stop = stop + 130000;
end

%%
time_ephys = linspace(0,13,1300000);
time_calcium = linspace(0,13,253);

axh=[];figure; hold on
ax(1) = subplot(20,1,1:9); hold on
plot(time_ephys,1000*concatenated_ephys,'Color','k');
set(gca,'Fontsize',12,'Ticklength',[0.003 0.003])
ylabel('V_{m} (mV)')
axis tight
ax(2) = subplot(20,1,10:20);
plot(time_calcium,smooth(corrected_fluorescence(1:253)),'color',[0 0 0.6],'Linewidth',1.5);
ylabel('Raw intensity')
set(gca,'Fontsize',12,'Ticklength',[0.003 0.003])
xlabel('t (s)')
box off
axis tight
linkaxes(ax,'x')