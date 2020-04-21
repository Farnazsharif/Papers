clc
clear

%%%% CA3
FN_TRD1={'SFA4_S5','SFA4_S6','SFA4_S7','SFA4_S9','SFA5_S5','SFA5_S6','SFA5_S7','SFA5_S9','SFA3_S5','SFA3_S6'};
FN_TRD2={'SFA4_S5','SFA4_S6','SFA4_S7','SFA4_S9','SFA5_S5','SFA5_S6','SFA5_S7','SFA5_S9','SFA3_S5','SFA3_S6'};
dir='E:\Farnaz\chronic\SFA5\SFA5_S5';
cd (dir)
TRD= 1;
FN=eval(['FN_TRD' num2str(TRD)]);

h=5;
SpeedThreshold=5;
smoothingRange=10;
xbinNumber=100;

filename=([FN{h} '_TRD' num2str(TRD)]);
load([filename])

%%
% single unit

Training_Percentile=.70;
SpeedThreshold=5;
smoothingRange=1;
xbinNumber=100;
CellVector=Cell_group.P;
CellVector((isnan(CellVector(:,1))==1),:)=[];

% gsub=[CellVector(:,4);G(46);G(1)];
gsub=[CellVector(:,4)];
[phase , rate]=makexttsc_PhS(gsub,xbinNumber,SpeedThreshold,filename);
Matrix=rate;
phase(isnan(phase)==1)=0;
position_train=[];spk_trains=[];phase_trains=[];
trial_numbers=unique(rate(:,3));

%%
cd C:\Users\aza\Dropbox\PS\Decoding
save(['Matrix_P'],'Matrix')
cd (dir)

%% generate shuffeld trials

% first we ranomize the trials 
random_tr_Number = randperm(length(trial_numbers));
Shuffled_trial(:,1)=trial_numbers(random_tr_Number);
% next we circshift the trials and repeat for all posible combinations
for iter = 1:length(Shuffled_trial)-1
    Shuffled_trial(:,iter+1)=circshift(Shuffled_trial(:,iter),1);
end

training_trails=round(Training_Percentile*size(Shuffled_trial,1));
nd_training_trails=1:training_trails*xbinNumber;
nd_test_trails=nd_training_trails(end)+1:size(rate,1);

%%
smoothingRange= 10;
positionDecodingMaxCorr = table;
CellN=10;
Bin_tr_Rate=[];
Bin_tr_Rate=rate(:,[1,3,4+CellN]);
Bin_tr_Rate(:,3)=smooth(Bin_tr_Rate(:,3),smoothingRange);

y_test=[];
X_test=[];
correct_class_decision_values=[]; normalized_rank_results=[]; rank_confusion_matrix=[];

for iter = 1:size(Shuffled_trial,1)
    
    % find Shuffle trials in the rate matrix
    Shuf_rate=[];
    for i=1:size(Shuffled_trial,1)
        nd=find(Bin_tr_Rate(:,2)==Shuffled_trial(i,iter));
        Shuf_rate=[Shuf_rate;Bin_tr_Rate(nd,:)];
    end
   
    % train and test
    % select 80% of the trials for training the model and 20% for test
    rate_train=Shuf_rate(nd_training_trails,3)';
    position_train=Shuf_rate(nd_training_trails,1)';
    
    rate_test=Shuf_rate(nd_test_trails,3)';
    position_test=Shuf_rate(nd_test_trails,1)';
    
    %% rate coding model
    ytest=[];cl=[];
    cl = max_correlation_coefficient_CL;
    cl = train(cl,rate_train,position_train);
%     yfit_rate=test(cl,rate_test);
    [ytest decision_values]=test(cl,rate_test);
    mse_rate = mean((ytest-position_test).^2);
    create_rank_confusion_matrix=1;
    [correct_class_decision_values(:,iter), normalized_rank_results(:,iter), rank_confusion_matrix(:,:,iter)] = ...
    get_rank_and_decision_value_results_Farnaz(ytest, cl.labels, decision_values, create_rank_confusion_matrix);

    Cltraining(iter,:)=cl.templates;
    y_test(iter,:)=ytest;
    X_test(iter,:)=position_test;
    
    %David test method
    ytest=[];
    for ts = 1:length(position_test)
        ytest(ts) = test(cl,rate_test(ts));
    end
    mse_rate_Dav = mean((ytest-position_test).^2);
    
    
    % estimate chance
    ytest=[];
    rr = randperm(length(rate_train));
    rrr = randperm(length(rate_test));
    cl = max_correlation_coefficient_CL;
    cl = train(cl,rate_train(rr),position_train);
    ytest=test(cl,rate_test(rrr));
    mse_chance_rate = mean((ytest-position_test).^2);
    
    % estimate real field decoding eror
    ytest=[];
    cl = max_correlation_coefficient_CL;
    cl = train(cl,Bin_tr_Rate(nd_training_trails,3)',Bin_tr_Rate(nd_training_trails,1)');
    ytest=test(cl,Bin_tr_Rate(nd_test_trails,3)');
    mse_rate_real= mean((ytest-Bin_tr_Rate(nd_test_trails,1)').^2);
    
    % store results
    struct.mse_rate=mse_rate;
    struct.mse_rate_Dav=mse_rate_Dav;
    struct.mse_rate_real=mse_rate_real;
    struct.mse_chance_rate=mse_chance_rate;
    % struct.condition = cond;
    struct.smoothingRange = smoothingRange;
    struct.iter = iter;
    struct.Training_Percentile = Training_Percentile;
    struct.trialOrder = Shuffled_trial(:,iter)';
    positionDecodingMaxCorr = [positionDecodingMaxCorr;struct2table(struct)];
    
    clear var yfit_rate struct
end



% close all
figure('position',[100 400 1000 200])
subplot(1,3,1)
plot(correct_class_decision_values,'.')
title([num2str(mean(mean(correct_class_decision_values)))])
xlabel('Position-data for test')
ylabel('correct-class-decision-values')
subplot(1,3,2)
plot(normalized_rank_results,'.')
title([num2str(mean(mean(normalized_rank_results)))])
xlabel('normalized-rank-results')
ylabel('correct-class-decision-values')
subplot(1,3,3)
imagesc(mean(rank_confusion_matrix,3))
xlabel('real position')
ylabel('Predictted Position')

figure('position',[100 400 1000 200])

% PLOT Place fileds
subplot(1,5,1)
nPFbin=100;
Matrix=rate;
smoothT=smooth1D(Matrix(:,2),smoothingRange,1);
smoothC=smooth1D(Matrix(:,CellN+4),smoothingRange,1);
ncell=length(smoothC(1,:));
rate2=smoothC./repmat(smoothT,1,ncell);
xtsr=[Matrix(:,[1 3 4]) rate2];
rr=xtsr(:,4);
matC=reshape(rr,nPFbin,length(rr)/nPFbin)';
imagesc(matnorm(matC,2))
xlabel('Bin')
ylabel('Trial number')
title( 'RateMap')

subplot(1,5,2)
plot(mean(matC),'k','linewidth',2)
xlabel('bin ')
ylabel('Rate')
title( 'Rate curve')

%plot eror VS iteratin
subplot(1,5,3)
Colmn=1;
U1=table2array(positionDecodingMaxCorr(:,Colmn));
Colmn=2;
U2=table2array(positionDecodingMaxCorr(:,Colmn));
Colmn=3;
U3=table2array(positionDecodingMaxCorr(:,Colmn));
plot(U1,'linewidth',2)
hold on
plot(U2,'linewidth',2)
hold on
plot(U3,'linewidth',2)
hold on
plot([0 iter],[mean(U1) mean(U1)],'r')
hold on
plot([0 iter],[median(U1) median(U1)],'k')
xlabel('iteration number')
ylabel('MSE')
% legend('mse','mseD','mseR')
title( 'eror with iteratin')

subplot(1,5,4)
plot(cl.labels, cl.templates,'k','linewidth',2)
xlabel('cl.labels')
ylabel('cl.templates')
title( 'Model templates')

subplot(1,5,5)
h=histogram2(X_test,y_test,'DisplayStyle','tile','ShowEmptyBins','on');
hold on
plot((mean(matC)./max(mean(matC)))*100,'w','linewidth',1)
hold on
plot(1:xbinNumber,'w','linewidth',1)
colormap jet
xlabel('position')
ylabel('Predicted position')
title('Model performace')


figure('position',[200 400 1000 200])
for i=1:size(y_test,1)
    plot(y_test(i,:),'.b','linewidth',.5)
    hold on
    plot(X_test(i,:),'k','linewidth',2)
    hold on
end
hold on
M=mean(matC)./max(mean(matC))*100;
M1=repmat(M,1,6);
plot(M1,'r','linewidth',2)
legend('Prediction','Position','Rate')
xlabel('Position-data for test= 20% (6 trial) ')
ylabel('normalized-rank-results')
title('data for test= 20% (6 trial)')


figure('position',[200 400 1000 200])
for i=1:size(correct_class_decision_values,2)
%     plot(correct_class_decision_values(:,i),'.b','linewidth',.5)
    plot(normalized_rank_results(:,i),'.b','linewidth',.5)

    
%     hold on
%     plot(X_test(i,:),'k','linewidth',2)
    hold on
end
hold on
M=(mean(matC)./max(mean(matC)));
M1=repmat(M,1,6);
plot(M1,'r','linewidth',2)
legend('Prediction','Position','Rate')
xlabel('Position-data for test= 20% (6 trial) ')
ylabel('Prediction for all iteration (22 trail)')
title('data for test= 20% (6 trial)')



%%
M=mean(matC)./max(mean(matC))*100;
M1=repmat(M,1,6)
figure
plot(M1)

%%
figure
plot(Bin_tr_Rate(nd_test_trails,1)')
figure
plot(Bin_tr_Rate(nd_training_trails,1)')

%%
clc
cl=[];
cl = max_correlation_coefficient_CL;
ytest=[];
cl = train(cl,Bin_tr_Rate(nd_training_trails,3)',Bin_tr_Rate(nd_training_trails,1)');
ytest=test(cl,Bin_tr_Rate(nd_test_trails,3)');
mse_rate_real = mean((ytest-Bin_tr_Rate(nd_test_trails,1)').^2)

%%
% close all
figure
plot(cl.templates)
figure
plot(ytest,'.-')
hold on
plot(Shuf_rate(nd_test_trails,1),'.')
figure
plot((ytest-Shuf_rate(nd_test_trails,1)').^2)
hold on
plot(Shuf_rate(nd_test_trails,3)'.^2*1000)

%%
for ts = 1:length(position_test)
    ytest(ts) = test(cl,[round(rates_trains_smooth_test(:,ts));theta_test(ts)]);
end
struct.mse_rate = mean((ytest-position_test).^2);


%%
trial_numbers=unique(Matrix(:,3))%#change 5 or 4
for t = 1:length(trial_numbers)
    ndx_trial=[];
    ndx_trial=find(Matrix(:,3)==trial_numbers(t));%#change
    phase_trains_C{h}{t} = phase(ndx_trial,5:end)';%#change
    spk_trains_C{h}{t} = rate(ndx_trial,5:end)';%#change
    position_C{h}{t} = Matrix(ndx_trial,1)';%#change
end


%%
positionDecodingMaxCorr_C=Dav_positionDecodingMaxCorr(phase_trains_C,spk_trains_C,position_C,smoothingRange);


%%
figure
% plot(circ_smoothTS(phase_trains{cond}{r(t)}(cell,:),wind,'method','mean','exclude',0))
plot(ytest,'.')
figure
% plot(circ_smoothTS(phase_trains{cond}{r(t)}(cell,:),wind,'method','mean','exclude',0))
plot(position_test,'.')

















