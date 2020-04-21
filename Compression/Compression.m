clc
clear

filename='SFA4_S3_TRD2';
load([filename '.mat'])
load([filename 'PlaceField.mat'])

%% find spikes and phases
Maingroup=[6 10 12 13];
CellVectoreTut=Cell_group.C(:,3);
% Pairs_txvtlph_1=[];
Pairs_txvtlph_1{1,1}=Cellinfo.Cell_txvtlph{CellVectoreTut(6)};
Pairs_txvtlph_1{1,2}=Cellinfo.Cell_txvtlph{CellVectoreTut(7)};
Pairs_txvtlph_1{2,2}=Cellinfo.Cell_txvtlph{CellVectoreTut(9)};
Pairs_txvtlph_1{3,2}=Cellinfo.Cell_txvtlph{CellVectoreTut(10)};
Pairs_txvtlph_1{4,2}=Cellinfo.Cell_txvtlph{CellVectoreTut(13)};

Pairs_txvtlph_2{1,1}=Cellinfo.Cell_txvtlph{CellVectoreTut(10)};
Pairs_txvtlph_2{1,2}=Cellinfo.Cell_txvtlph{CellVectoreTut(6)};
Pairs_txvtlph_2{2,2}=Cellinfo.Cell_txvtlph{CellVectoreTut(7)};
Pairs_txvtlph_2{3,2}=Cellinfo.Cell_txvtlph{CellVectoreTut(8)};
Pairs_txvtlph_2{4,2}=Cellinfo.Cell_txvtlph{CellVectoreTut(9)};
Pairs_txvtlph_2{5,2}=Cellinfo.Cell_txvtlph{CellVectoreTut(12)};
Pairs_txvtlph_2{6,2}=Cellinfo.Cell_txvtlph{CellVectoreTut(13)};

Pairs_txvtlph_3{1,1}=Cellinfo.Cell_txvtlph{CellVectoreTut(12)};
Pairs_txvtlph_3{1,2}=Cellinfo.Cell_txvtlph{CellVectoreTut(7)};
Pairs_txvtlph_3{2,2}=Cellinfo.Cell_txvtlph{CellVectoreTut(8)};
Pairs_txvtlph_3{3,2}=Cellinfo.Cell_txvtlph{CellVectoreTut(9)};
Pairs_txvtlph_3{4,2}=Cellinfo.Cell_txvtlph{CellVectoreTut(10)};

Pairs_txvtlph_4{1,1}=Cellinfo.Cell_txvtlph{CellVectoreTut(13)};
Pairs_txvtlph_4{1,2}=Cellinfo.Cell_txvtlph{CellVectoreTut(6)};
Pairs_txvtlph_4{2,2}=Cellinfo.Cell_txvtlph{CellVectoreTut(7)};
Pairs_txvtlph_4{3,2}=Cellinfo.Cell_txvtlph{CellVectoreTut(9)};
Pairs_txvtlph_4{4,2}=Cellinfo.Cell_txvtlph{CellVectoreTut(10)};
Pairs_txvtlph_4{5,2}=Cellinfo.Cell_txvtlph{CellVectoreTut(12)};

%% find rates

Matrix=xttsc;
smooth=10;delta=0.10;TRX=[[1 20 40]; [1 20 40]+45];
tr_1=[1,fix(length(Matrix)/200),length(Matrix)/100];
[Field_edge,rows,PF] = Feild_width(filename,CellVectoreTut,smooth,delta,Matrix,tr_1,TRX);

Pairs_rate_1{1,1}(:,:)=PF(6,:,:);
Pairs_rate_1{1,2}(:,:)=PF(7,:,:);
Pairs_rate_1{2,2}(:,:)=PF(9,:,:);
Pairs_rate_1{3,2}(:,:)=PF(10,:,:);
Pairs_rate_1{4,2}(:,:)=PF(13,:,:);

Pairs_rate_2{1,1}(:,:)=PF(10,:,:);
Pairs_rate_2{1,2}(:,:)=PF(6,:,:);
Pairs_rate_2{2,2}(:,:)=PF(7,:,:);
Pairs_rate_2{3,2}(:,:)=PF(8,:,:);
Pairs_rate_2{4,2}(:,:)=PF(9,:,:);
Pairs_rate_2{5,2}(:,:)=PF(12,:,:);
Pairs_rate_2{6,2}(:,:)=PF(13,:,:);

Pairs_rate_3{1,1}(:,:)=PF(12,:,:);
Pairs_rate_3{1,2}(:,:)=PF(7,:,:);
Pairs_rate_3{2,2}(:,:)=PF(8,:,:);
Pairs_rate_3{3,2}(:,:)=PF(9,:,:);
Pairs_rate_3{4,2}(:,:)=PF(10,:,:);

Pairs_rate_4{1,1}(:,:)=PF(13,:,:);
Pairs_rate_4{1,2}(:,:)=PF(6,:,:);
Pairs_rate_4{2,2}(:,:)=PF(7,:,:);
Pairs_rate_4{3,2}(:,:)=PF(9,:,:);
Pairs_rate_4{4,2}(:,:)=PF(10,:,:);
Pairs_rate_4{5,2}(:,:)=PF(12,:,:);
%%

Field_edge_1{1,1}=Field_edge(6,:);
Field_edge_1{1,2}=Field_edge(7,:);
Field_edge_1{2,2}=Field_edge(9,:);
Field_edge_1{3,2}=Field_edge(10,:);
Field_edge_1{4,2}=Field_edge(13,:);

Field_edge_2{1,1}=Field_edge(10,:);
Field_edge_2{1,2}=Field_edge(6,:);
Field_edge_2{2,2}=Field_edge(7,:);
Field_edge_2{3,2}=Field_edge(8,:);
Field_edge_2{4,2}=Field_edge(9,:);
Field_edge_2{5,2}=Field_edge(12,:);
Field_edge_2{6,2}=Field_edge(13,:);

Field_edge_3{1,1}=Field_edge(12,:);
Field_edge_3{1,2}=Field_edge(7,:);
Field_edge_3{2,2}=Field_edge(8,:);
Field_edge_3{3,2}=Field_edge(9,:);
Field_edge_3{4,2}=Field_edge(10,:);

Field_edge_4{1,1}=Field_edge(13,:);
Field_edge_4{1,2}=Field_edge(6,:);
Field_edge_4{2,2}=Field_edge(7,:);
Field_edge_4{3,2}=Field_edge(9,:);
Field_edge_4{4,2}=Field_edge(10,:);
Field_edge_4{5,2}=Field_edge(12,:);

%%
save('Pairs_txvtlph_1','Pairs_txvtlph_1')
save('Pairs_txvtlph_2','Pairs_txvtlph_2')
save('Pairs_txvtlph_3','Pairs_txvtlph_3')
save('Pairs_txvtlph_4','Pairs_txvtlph_4')

%% Part 1 load data
clear
clc
cd SingleCellData
load('Pairs_txvtlph_4.mat');
load('Pairs_txvtlph_3.mat');
load('Pairs_txvtlph_2.mat');
load('Pairs_txvtlph_1.mat');
cd ..

%% inputs
CN=2;
speed_thr=5;

% spike information
txvtlph_1=Pairs_txvtlph_1{1, 1};
txvtlph_2=Pairs_txvtlph_1{CN, 2};

% select spikes
% All running spikes

nd=find(txvtlph_1(:,3)>speed_thr);
txvtlph_1=txvtlph_1(nd,:);
nd=find(txvtlph_2(:,3)>speed_thr);
txvtlph_2=txvtlph_2(nd,:);
%***********************************
CellN1_SpikeTime=txvtlph_1(:,1); 
CellN2_SpikeTime=txvtlph_2(:,1); 
length(CellN1_SpikeTime)
length(CellN2_SpikeTime)

%% Select Infiled running spikes
clc
track_length=230;
binN=100;

Field_edge=Field_edge_1{1, 1};
TX=txvtlph_1(:,1:2);
Boundries=[Field_edge(1,1)*track_length/binN Field_edge(1,3)*track_length/binN];
rows=mean(Pairs_rate_1{1,1});
Spike_nd=[];
[Spike_nd] = InfieldspikesNdex_circular(TX,track_length,Boundries,rows);
txvtlph_1=txvtlph_1(Spike_nd,:);

Field_edge=Field_edge_1{CN, 2};
TX=txvtlph_2(:,1:2);
Boundries=[Field_edge(1,1)*track_length/binN Field_edge(1,3)*track_length/binN];
rows=mean(Pairs_rate_1{CN,2});
Spike_nd=[];
[Spike_nd] = InfieldspikesNdex_circular(TX,track_length,Boundries,rows);
txvtlph_2=txvtlph_2(Spike_nd,:);
%***********************************

CellN1_SpikeTime=txvtlph_1(:,1); 
CellN2_SpikeTime=txvtlph_2(:,1); 
length(CellN1_SpikeTime)
length(CellN2_SpikeTime)

%% INPUTS 
CellN1_SpikeTime=txvtlph_1(:,1); 
CellN2_SpikeTime=txvtlph_2(:,1); 
% mean rate
CellN1_1Drate=mean(Pairs_rate_1{1,1});
CellN2_1Drate=mean(Pairs_rate_1{CN,2});

CellN1_PositionPhase=txvtlph_1(:,[2,6]);
CellN2_PositionPhase=txvtlph_2(:,[2,6]);

smoothing_range_local=10;
smoothing_range_global=40;
ccg_binsize=0.01;

[CCG_shiftloc,CCG_shiftglobal,Peak_Distance] = compresion(CellN1_SpikeTime,CellN2_SpikeTime,CellN1_1Drate,CellN2_1Drate,smoothing_range_local,...
                                                       smoothing_range_global,ccg_binsize,CellN1_PositionPhase,CellN2_PositionPhase,track_length);























