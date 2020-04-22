%   Compression tutorial is relatted to Figure 4 of the fowllowing paper: 
%   "Subcircuits of deep and superficial CA1 place cells support efficient spatial
%    coding across heterogeneous environments", doi: https://doi.org/10.1101/2020.04.17.047399

%    Dataset gathered in Sebastien Royer lab (2017), animal is running head-fixed on a non motorized treadmill


% Written by Farnaz Sharif
% Buzsaki lab, NYU Neuroscience Institute, New York University, Langone Medical Center
% Dec 2019

%% Part 1 load data
clear
clc

cd Data
load('CellN1_SpikeTime','CellN1_SpikeTime')
load('CellN2_SpikeTime','CellN2_SpikeTime')
load('CellN1_1Drate','CellN1_1Drate')
load('CellN2_1Drate','CellN2_1Drate')
load('CellN1_PositionPhase','CellN1_PositionPhase')
load('CellN2_PositionPhase','CellN2_PositionPhase')
cd ..

%% Part 2 calculate and plot compression

smoothing_range_local=10;  % for fitting a line to all the local maxima of the ccg
smoothing_range_global=40; % for fitting a line to the whole ccg
track_length=230;          % linear track length
ccg_binsize=0.01;

[CCG_shiftloc,CCG_shiftglobal,Peak_Distance] = Compres(CellN1_SpikeTime,CellN2_SpikeTime,CellN1_1Drate,CellN2_1Drate,smoothing_range_local,smoothing_range_global, ...
    ccg_binsize,track_length,CellN1_PositionPhase,CellN2_PositionPhase);



