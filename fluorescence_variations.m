function [loc_ONSET,loc_OFFSET] = fluorescence_variations(calcium_file)
%This function is used on .sif files coming from calcium imaging. To limit
%bleaching, fluorescence is not 'ON' during the whole sweep/whole movie.
%Hence there are some dark periods that are a pain for the fluorescence
%signal AND the visual information. This function is part of a skeleton
%allowing to concentrate only on fluorescence_ON periods. 
%
%   INPUT
%
%   calcium_data    the previously imported calcium data from the sif file
%
%   OUTPUT
%
%   loc_ONSET       list with all frames indices at which fluorescence was 
%                   turned on
%   loc_OFFSET      same, but turned off   
%
%

movie =  dir(calcium_file);
calcium_data = imread_universal(movie.name);


summed_frames = sumframe(calcium_data); %summing over all frames the value of each pixel
fluorescence_ON = summed_frames > 2.5e6;%finding the frames where the fluorescence was ON, threshold is bad
fluorescence_OFF = summed_frames < 2.5e6;
nb_frames = length(summed_frames);

ONSET = zeros(nb_frames,1); 
OFFSET = zeros(nb_frames,1); %creating two too-long lists for the following steps

for k=2:nb_frames
    switching = fluorescence_ON(k) - fluorescence_ON(k-1); %we are looking for switches of the fluorescence
    
    if switching==1 %then it was switched on
        ONSET(k)=1; %adding a 1 at the right position
        fluorescence_OFF(k) = 1; %this way we make sure we loose discontinuities
        
    elseif switching==-1 %then it was switched off
        OFFSET(k)=1;     %ditto
        fluorescence_OFF(k-1)=1;
    end
end
loc_ONSET = find(ONSET==1); %finally, we create two shorter lists only keeping the indices
loc_OFFSET = find(OFFSET==1);


cropped_data = calcium_data; %just for visual input, could be used later to visualise data without gaps
cropped_data(:,:,fluorescence_OFF(:)) = [];
sf = sumframe(cropped_data);
figure;
plot(sf) 
    

