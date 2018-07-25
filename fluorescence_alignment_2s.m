function [extracted_fluorescence_sweeps] = fluorescence_alignment_2s(filename,nb_sweeps)
%This function takes the whole fluorescence trace from the 2P imaging, cuts
%it in small pieces according to the 2s protocol waveform and gives back a
%plot with aligned fluorescence sweeps
%   INPUT
%   filename      the .csv filename in the cd, or just '*...*'
%
%   nb_sweeps     the number of sweeps ran in the electrophysiology
%   protocol

%% Creating a matrix with individual fluorescence sweeps

sweep_length_s = 10; %specific to the 2s protocol

fluorescence = dir(filename);
fluorescence_data = csvread(fluorescence.name,1,0); %now we have our matrix

[~,x2] = size(fluorescence_data);
sampling_rate = fluorescence_data(2,1) - fluorescence_data(1,1); %this is the sampling rate of the 2P imaging
nb_points = floor(sweep_length_s/sampling_rate);

extracted_fluorescence_sweeps = zeros(nb_points,nb_sweeps,x2-1); % we create a matrix where we'll store the information

for region=1:x2-1
    start=1;
    stop=nb_points;
    for sweep=1:nb_sweeps
        extracted_fluorescence_sweeps(:,sweep,region) = fluorescence_data(region+1,start:stop);
    end
end

time_vector = fluorescence_data(1,1:nb_points); % just creating a time vector in s
%% Plotting individual fluorescence sweeps

y_names = 1:nb_sweeps;

for region=1:x2-1
    figure; hold on
    ylabel('Trial')
    xlabel('t (s)')
    Maxes = zeros(x2,1);
    y_positions = zeros(nb_sweeps,1); %we'll use this list to have the right x_ticks
    for sweep=1:nb_sweeps
        Back_to_zero = extracted_fluorescence_sweeps(:,sweep,region) - min(extracted_fluorescence_sweeps(:,sweep,region)); %we bring each sweep back to 0
        plot(time_vector,Back_to_zero+sum(Maxes),'Linewidth',1.5) %so we plot the sweep with the adapted positive offset
        y_positions(sweep) = sum(Maxes);
        Maxes(region-1) = max(Back_to_zero); %and store the max value here for the next round of the loop
    end
    
    maximulus = max(Maxes);
    current_trace = zeros(nb_points,1);
    c=1;
    for k=time_vector
        if k<2
            current_trace(c) = 0;
            c=c+1;
        elseif k>4
            current_trace(c) = 0;
            c=c+1;
        elseif (4>k) && (k>2)
            current_trace(c)=maximulus/2;
            c=c+1;
        end
    end % this section here is completely specific to the 2s pulse and is used to plot a fake I trace above. 
    
    yticks(flipud(y_positions));
    yticklabels(fliplr(y_names))
    axis tight
    plot(time_vector,current_trace + sum(Maxes),'Linewidth',1.5,'Color','k') %plotting this all above
    title(['ROI',num2str(region)])    
end

