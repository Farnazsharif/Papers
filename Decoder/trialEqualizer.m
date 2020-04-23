% Written by Farnaz Sharif
% Buzsaki lab, NYU Neuroscience Institute, New York University, Langone Medical Center
% Sept 2019

function [Bin_tr_Rate,AllMaps]=trialEqualizer(rate,tr_imresize,smoothingRange)
% equalize and normalize the trial numbers and generatting psudo-trials#### 
% normalize each map trial by trial
smoothingRange=10;

xbinNumber=max(rate(:,1));



bin_tr=[];
for tr=1:tr_imresize
    binTr(:,1)=1:xbinNumber;
    binTr(:,2)=ones(xbinNumber,1)*tr;
    bin_tr=[bin_tr;binTr];
end
Bin_tr_Rate(:,1:2)=bin_tr;


CellN_all=1:size(rate,2)-4;
    for CellN=CellN_all;
        
         smoothTime=smooth(rate(:,2),smoothingRange);
         smoothSpikes=smooth(rate(:,4+CellN),smoothingRange);
         smoothRate=smoothSpikes./smoothTime;
         Ratemap=reshape(smoothRate,xbinNumber,length(smoothRate)/xbinNumber)';
         RatemapNormal=matnorm(imresize(Ratemap,[tr_imresize xbinNumber],'lanczos3'),2);
         Bin_tr_Rate(:,2+CellN)=reshape(RatemapNormal',[xbinNumber*tr_imresize 1 ]);
         AllMaps{CellN}=RatemapNormal;
         
    end