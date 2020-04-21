


%% single unit
close all
clear
clc
CellNumber=1:9;
load('Matrix_C.mat')
rate=Matrix(:,[1:4 CellNumber+4]);
Training_Percentile=.70;
SpeedThreshold=5;
smoothingRange=10;
xbinNumber=100;
tr_imresize=28;

[Bin_tr_Rate,AllMaps]=trialEqualizer(rate,tr_imresize,smoothingRange);


%% generate shuffeld trials
trial_numbers=unique(Bin_tr_Rate(:,2));

% first we ranomize the trials 
random_tr_Number = randperm(length(trial_numbers));
Shuffled_trial(:,1)=trial_numbers(random_tr_Number);
% next we circshift the trials and repeat for all posible combinations
for iter = 1:length(Shuffled_trial)-1
    Shuffled_trial(:,iter+1)=circshift(Shuffled_trial(:,iter),1);
end

training_trails=round(Training_Percentile*size(Shuffled_trial,1));
nd_training_trails=1:training_trails*xbinNumber;
nd_test_trails=nd_training_trails(end)+1:size(Bin_tr_Rate,1);

%%
positionDecodingMaxCorr = table;

y_test=[];
X_test=[]; 
Cl_temp=[];Cl_X=[];
CellN=[3:size(Bin_tr_Rate,2)];

for iter = 1:size(Shuffled_trial,1)
    
    % find Shuffle trials in the rate matrix
    Shuf_rate=[];
    for i=1:size(Shuffled_trial,1)
        nd=find(Bin_tr_Rate(:,2)==Shuffled_trial(i,iter));
        Shuf_rate=[Shuf_rate;Bin_tr_Rate(nd,:)];
    end
   
    % train and test
    rate_train=Shuf_rate(nd_training_trails,CellN)';
    position_train=Shuf_rate(nd_training_trails,1)';
    
    rate_test=Shuf_rate(nd_test_trails,CellN)';
    position_test=Shuf_rate(nd_test_trails,1)';
    
    %% rate coding model
    ytest_rate=[];cl=[];
    cl = max_correlation_coefficient_CL;
    cl = train(cl,rate_train,[position_train]);
    [ytest_rate decision_values]=test(cl,rate_test);
    mse_rate = mean((ytest_rate-position_test).^2);

    create_rank_confusion_matrix=1;
    [correct_class_decision_values(:,iter), normalized_rank_results(:,iter), rank_confusion_matrix(:,:,iter)] = ...
    get_rank_and_decision_value_results_Farnaz(ytest_rate, cl.labels, decision_values, create_rank_confusion_matrix);
    Cl_temp=[Cl_temp;cl.templates];
    Cl_X=[Cl_X;cl.labels];
    y_test(iter,:)=ytest_rate;
    X_test(iter,:)=position_test;
        
    % estimate chance
    ytest_rate=[];
    rr = randperm(length(rate_train));
    rrr = randperm(length(rate_test));
    cl_chance = max_correlation_coefficient_CL;
    cl_chance  = train(cl_chance ,rate_train(rr),position_train);
    ytest_rate=test(cl_chance ,rate_test(rrr));
    mse_chance_rate = mean((ytest_rate-position_test).^2);
   
    % store results
    struct.mse_rate=mse_rate;
    struct.mse_chance_rate=mse_chance_rate;
    struct.smoothingRange = smoothingRange;
    struct.iter = iter;
    struct.Training_Percentile = Training_Percentile;
    struct.trialorder = Shuffled_trial(:,iter)';
    positionDecodingMaxCorr = [positionDecodingMaxCorr;struct2table(struct)];
    
    clear var yfit_rate struct 
end

Cl_X=[Cl_X;Cl_X];
%%
Maps=zeros(tr_imresize,xbinNumber);
rows=zeros(1,xbinNumber);
for i=1:size(AllMaps,2)
    Maps=Maps+AllMaps{:,i};
    rows(i,:)=mean(AllMaps{:,i});
end
Maps=Maps./size(AllMaps,2);

% close all
figure('position',[100 400 1000 200])

% Plot all Place maps
subplot(1,5,1)
imagesc(Maps)
xlabel('Bin')
ylabel('Trial number')
title( 'RateMap')

% Plot all firing rate curves
subplot(1,5,2)
plot(1:xbinNumber,rows,'k','linewidth',2)
hold on
plot(1:xbinNumber,mean(rows,1),'r','linewidth',2)
xlabel('bin ')
ylabel('Rate')
title( 'Rate curve')

%plot eror VS iteratin
subplot(1,5,3)
Colmn=1;
U1=table2array(positionDecodingMaxCorr(:,Colmn));
Colmn=2;
U2=table2array(positionDecodingMaxCorr(:,Colmn));

plot(U1,'b','linewidth',2)
hold on
plot(U2,'r','linewidth',2)
hold on
plot([0 iter],[mean(U1) mean(U1)],'b')
hold on
plot([0 iter],[median(U1) median(U1)],':b')
hold on
plot([0 iter],[mean(U2) mean(U2)],'r')
xlabel('iteration number')
ylabel('MSE')
title( 'eror with iteratin')


subplot(1,5,4)
plot(1:size(Cl_temp,2),(Cl_temp),'k','linewidth',2)
hold on
plot(mean(Cl_X),mean(Cl_temp),'r','linewidth',2)
xlabel('cl.labels')
ylabel('cl.templates')
title( 'Model templates')

subplot(1,5,5)
h=histogram2(X_test,y_test,'DisplayStyle','tile','ShowEmptyBins','on');
hold on
plot(1:xbinNumber,mean(rows,1)./max(mean(rows,1))*xbinNumber,'w','linewidth',1)
hold on
plot(1:xbinNumber,'w','linewidth',1)
colormap jet
xlabel('position')
ylabel('Predicted position')
title('Model performace')

figure('position',[200 400 1000 200])
plot(1:size(y_test,2),y_test,'.b','linewidth',.5)
hold on
plot(1:size(X_test,2),X_test,'k','linewidth',2)
hold on
M=[];
M=mean(rows,1)./max(mean(rows,1))*xbinNumber;
M=repmat(M,1,size(y_test,2)./xbinNumber);
plot(M,'r','linewidth',2)
legend('Prediction','Position','Rate')
xlabel('Position-data for test= 30%  ')
ylabel('Prediction for all iteration trail')
title('data for test= 30% ( trial)')
%%
figure('position',[100 400 1000 200])
subplot(1,4,1)
imagesc(matnorm(mean(rank_confusion_matrix,3),1))

axis xy
xlabel('real position')
ylabel('Predictted Position')

subplot(1,4,2)
plot(normalized_rank_results,'.')
title([num2str(mean(mean(normalized_rank_results)))])
xlabel('normalized-rank-results')
ylabel('correct-class-decision-values')

subplot(1,4,3)
M=[];
M=reshape(correct_class_decision_values,xbinNumber,tr_imresize,size(correct_class_decision_values,1)./xbinNumber);
M=correct_class_decision_values;
plot(correct_class_decision_values,'.')
title([num2str(mean(mean(correct_class_decision_values)))])
xlabel('Position-data for test')
ylabel('correct-class-decision-values')

subplot(1,4,4)
plot(1:800,decision_values,'.')


%%
decisionVal=-.2;MSE=1000;
decisionVal=-.25;MSE=1100;


figure
h=histogram2(X_test,normalized_rank_results','DisplayStyle','tile','ShowEmptyBins','on');
%%
clc
%A=normalized_rank_results;
A=correct_class_decision_values;
M=zeros(xbinNumber,tr_imresize);
M=[];
for i=1:size(A,1)./xbinNumber
    N(i,1)=[xbinNumber*(i-1)+1];
    N(i,2)=[xbinNumber*i];
    M(:,:,i)=A(N(i,1):N(i,2),:);
end

M=nanmean(M,3);


% close all
figure
plot(mean(M,2))
figure
imagesc(M')
axis xy

%%
figure
M=[];
M=AllMaps{1,4};
imagesc(M)



