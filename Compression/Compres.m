function [CCG_shiftloc,CCG_shiftglobal,Peak_Distance] = Compres(CellN1_SpikeTime,CellN2_SpikeTime,CellN1_1Drate,CellN2_1Drate,smoothing_range_local,smoothing_range_global, ...
    ccg_binsize,track_length,CellN1_PositionPhase,CellN2_PositionPhase)
% Input:
% CellN1_SpikeTime : SpikeTime of cell number 1
% CellN2_SpikeTime : SpikeTime of cell number 2
% CellN1_rows : 1D firing rate of cell number 1 in bins
% CellN2_rows : 1D firing rate of cell number 2
% smoothing_range_local : smoothing range for the local small peaks in CCG
% smoothing_range_global : smoothing range for the whole CCG
% ccg_binsize : time bin size of the CCG
% CellN1_PositionPhase : (Optional for plot) Position and phase of the spikes for the cell number 1
% CellN2_PositionPhase : (Optional for plot) Position and phase of the spikes for the cell number 2


% Written by Farnaz Sharif
% Buzsaki lab, NYU Neuroscience Institute, New York University, Langone Medical Center
% Dec 2019


%% Distance between the peaks
[Peakrate1,PeakPosition1]=max(CellN1_1Drate);
[Peakrate2,PeakPosition2]=max(CellN2_1Drate);
Peak_Distance=(PeakPosition2-PeakPosition1).*(track_length/length(CellN1_1Drate));


%% Fast theta-cycle difference between place fields

[ccg,ccg_time] = CCG([CellN1_SpikeTime;CellN2_SpikeTime], [ones(size(CellN1_SpikeTime)) ; 2*ones(size(CellN2_SpikeTime))], 'duration',2,'binsize',ccg_binsize);

ccg_time = ccg_time * 1000;
ccgNorm = ccg(:,1,2)'./max(ccg(:,1,2));

ccg_smooth_local = nanconv(ccgNorm,gausswin(smoothing_range_local)','edge');

    [pks,locs] = findpeaks(ccg_smooth_local);

    if isempty(locs)~=1
        [~,J]=sort(abs(ccg_time(locs)),'ascend');
        CCG_shiftloc(:,1)=ccg_time(locs(J));
        CCG_shiftloc(:,2)=pks(J);
    else
        CCG_shiftloc=[nan nan];
    end
    
    ccg_smooth_global = nanconv(ccgNorm,gausswin(smoothing_range_global)','edge');
    [CCG_shiftglobal(1,2),nd_max]=max(ccg_smooth_global);
    CCG_shiftglobal(1,1)=ccg_time(nd_max);

%% plot compression

if nargin >8
    %%
    figure('position',[200 200 1600 300])
    subplot(1,4,1)
    plot(ccg_time,ccgNorm);
    hold on
    plot(ccg_time,ccg_smooth_local,'b','linewidth',2)
    hold on
    plot([ccg_time==0 ccg_time==0],[0 1],'g','linewidth',1)
    hold on
    plot(CCG_shiftloc(1,1),CCG_shiftloc(1,2),'ok','markersize',5,'linewidth',2)
    hold on
    plot([0 CCG_shiftloc(1,1)],[CCG_shiftloc(1,2),CCG_shiftloc(1,2)],'k:','linewidth',2)
    hold on
    plot(ccg_time,ccg_smooth_global,'r','linewidth',2)
    hold on
    plot(CCG_shiftglobal(1,1),CCG_shiftglobal(1,2),'ok','markersize',5,'linewidth',2)
    hold on
    plot([0 CCG_shiftglobal(1,1)],[CCG_shiftglobal(1,2) CCG_shiftglobal(1,2)],'k:','linewidth',2)
    ylabel('Normalized CCG')
    xlabel('Time (ms)')
    title(['Local \Deltat= ' num2str(CCG_shiftloc(1,1)) 'ms' '   Global \Deltat= ' num2str(CCG_shiftglobal(1,1)) 'ms'])
    
    subplot(1,4,2)
    plot(CellN1_1Drate,'linewidth',2)
    hold on
    plot(CellN2_1Drate,'linewidth',2)
    hold on
    plot([PeakPosition1 PeakPosition1],[0 max([Peakrate1;Peakrate2])],'k','linewidth',2)
    hold on
    plot([PeakPosition2 PeakPosition2],[0 max([Peakrate1;Peakrate2])],'k','linewidth',2)
    hold on
    plot([PeakPosition1 PeakPosition2],[max([Peakrate1;Peakrate2]) max([Peakrate1;Peakrate2])],'k:','linewidth',2)
    
    title(['Peak distance=' num2str(Peak_Distance) ' cm'])
    ylabel('Firingrate (Hz)')
    xlabel('Position (bin)')
    
    subplot(1,4,3)
    plot(CellN1_PositionPhase(:,1),CellN1_SpikeTime,'.')
    hold on
    plot(CellN2_PositionPhase(:,1),CellN2_SpikeTime,'.')
    % xlim([0 track_length])
    ylabel('Time (ms)')
    xlabel('Position (cm)')
    title('Spike raster')
    subplot(1,4,4)
    plot([CellN1_PositionPhase(:,1); CellN1_PositionPhase(:,1)],radtodeg([CellN1_PositionPhase(:,2); CellN1_PositionPhase(:,2)+2*pi]),'.')
    hold on
    plot([CellN2_PositionPhase(:,1); CellN2_PositionPhase(:,1)],radtodeg([CellN2_PositionPhase(:,2); CellN2_PositionPhase(:,2)+2*pi]),'.')
    ylim([0 720])
    % xlim([0 track_length])
    ylabel('Phase (degree)')
    xlabel('Position (cm)')
    title('Phase precession')
end

end

