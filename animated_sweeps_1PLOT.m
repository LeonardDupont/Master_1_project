function [] = animated_sweeps_1PLOT(normalised_sweeps,start_sweep,stop_sweep,interpolated_calcium_data_DeltaFF,region,ephys_sweeps)
%This function is strictly used in the context of electrophysiology data
%against calcium data alignment. It is quite fitting in that one sweep at
%a time is being displayed, but one can navigate through sweeps using left
%and right arrow keys. All the input variables can be properly generated
%using the calcium analysis script.

 fig=figure('Name',['ROI',num2str(region)]); 
 set(fig, 'KeyPressFcn',@keypress)
 
 sweep_max = ephys_sweeps;
 sweep = 1;
    
    function redraw()
        clf
        hold on
        plot(squeeze(normalised_sweeps(start_sweep:stop_sweep,sweep)))
        title(['sweep ', num2str(sweep)])
        xlabel('t (10^{-4} s)')
        ylabel('V_{m}')
        min_elec = min(normalised_sweeps(start_sweep:stop_sweep,sweep));
        max_calcium = max(interpolated_calcium_data_DeltaFF(:,sweep,region));
        x = -.5 - min_elec + max_calcium;
        plot(interpolated_calcium_data_DeltaFF(:,sweep,region)-x,'b','Linewidth',2)
        [~,LOC] = findpeaks(interpolated_calcium_data_DeltaFF(:,sweep,region),'MinPeakProminence',.1);
        vline(LOC,'-.')
        axis tight
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


