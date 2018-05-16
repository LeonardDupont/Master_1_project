function [] = subtract_mean_from_file(calcium_file)
%This function should be used to subtract the mean of a given pixel to
%this pixel in each frame of the calcium movie. It should be applied on
%hilolike_filtered data. 

whichfile = dir(calcium_file);
calcium_data = imread_universal(whichfile.name);

averaged_data = mean(calcium_data,3); %averaging it over the third dimension
averaged_rep = repmat(averaged_data,1,1,length(calcium_data));

subtracted_mean = double(calcium_data) - averaged_rep;
subtracted_uint = uint8(subtracted_mean); %otherwise it's a double, not easy to write a tif from this

split_name = strsplit(whichfile.name,'.tif'); %splitting the string at a given delimiter position
new_name = strcat(split_name{1},'_subtracted_mean','.tif'); %used to concatenate strings

util.imwritetiff(subtracted_uint,new_name,'close',true) %important to close the file 

end

