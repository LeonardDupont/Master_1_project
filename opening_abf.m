function [data,si,h] = opening_abf(file_number)
  
    %The purpose of this function is to open a .abf file and prepare it for
    %further analysis. 
    %
    %   INPUT
    %
    %   file_number        corresponds to the recording number of the file
    %                      with stars around it. Example = 0012 -->
    %                      '*0012*'
    %
    %   OUTPUT
    %   The given output is in the typicalformat of the abffile
    %   
    %   data               the value matrix
    %   si                 sampling interval in us
    %   h                  the header
    % 
    %   A plot with the different I-V channels as subplots and channel
    %   number as a subplot title is also given, allowing the simplified
    %   use of other programmes. 

    data_file = dir(file_number); %this gives a struct, first part is
    %a char
    [data, si, h] = abfload(data_file.name);
    channels = h.nADCNumChannels;
    split = channels / 2;
    figure;
    for k=1:channels
        subplot(split,2,k)
        plot(data(:,k,1))
        title(num2str(k))
    end  

end
