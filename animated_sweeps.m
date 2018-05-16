function [] = animated_sweeps(electrophysiology_data,spikes,channel,start_sweep,stop_sweep,interpolated_calcium_data,region,loc_ONSET)
%This function is strictly used in the context of electrophysiology data
%against calcium data alignment. It is quite fitting in that one sweep at
%a time is being displayed, but one can navigate through sweeps using left
%and right arrow keys. All the input variables can be properly generated
%using the calcium analysis script.

 fig=figure('Name',['ROI',num2str(region)]);
 set(fig, 'KeyPressFcn',@keypress)
 
 sweep_max = length(loc_ONSET);
 sweep = 1;
    
    function redraw()
        clf
        axh=[];  %This time we create one only figure for all ROIs with average traces 
        axh(1) = subplot(15,1,1:5);hold on
        plot(squeeze(electrophysiology_data(start_sweep:stop_sweep,channel,sweep)))
        title(['sweep ', num2str(sweep)])
        axh(2) = subplot(15,1,6:15); hold on
        plot(interpolated_calcium_data(:,sweep,region),'b','Linewidth',2)
        star_position = max(interpolated_calcium_data(:,sweep,region))-20;
        plot(spikes(sweep),star_position, '*')
        linkaxes(axh,'x')
    end

redraw();

    function keypress(~,evnt)
        switch lower(evnt.Key)
            case 'rightarrow'
                sweep = min(sweep+1,sweep_max);
            case 'leftarrow'
                sweep = max(1 , sweep-1);
            otherwise
                return
        end
        redraw()
    end
end


