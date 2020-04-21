clear
% % cd E:\Farnaz\chronic\N_LFP_CA1
dir= 'E:\Farnaz\chronic\N_LFP_CA1'
cd (dir)
load('rate_C_TRD1.mat')
% load('rate_C_TRD2.mat')
load('rate_P_TRD1.mat')
% load('rate_P_TRD2.mat')
cd ..
%%
tr_imresize=40;
Training_Percentile=.70;
SpeedThreshold=5;
smoothingRange=10;
xbinNumber=100;
Big_M=rate_P;
rank_confusion_matrix_all=[];
CellN_all=0;
bin_tr=[];
for tr=1:tr_imresize
    binTr(:,1)=1:xbinNumber;
    binTr(:,2)=ones(xbinNumber,1)*tr;
    bin_tr=[bin_tr;binTr];
 end
 tic        
for Session=1:size(Big_M,2)
    Session
    rate=[];
    rate=Big_M{1, Session} ;
    trial_numbers=unique(bin_tr(:,2));
%     trial_numbers=unique(rate(:,3)); % Original 1
   
    
    
%  generate shuffeld trials ###############################################
    
    % first we ranomize the trials
    random_tr_Number = randperm(length(trial_numbers));
    Shuffled_trial(:,1)=trial_numbers(random_tr_Number);
    % next we circshift the trials and repeat for all posible combinations
    for iter = 1:length(Shuffled_trial)-1
        Shuffled_trial(:,iter+1)=circshift(Shuffled_trial(:,iter),1);
    end
    
    training_trails=round(Training_Percentile*length(trial_numbers));
    nd_training_trails=1:training_trails*xbinNumber;
    nd_test_trails=nd_training_trails(end)+1:length(trial_numbers)*xbinNumber;
%##########################################################################
    
    
    for CellN=1:size(rate,2)-4;
        positionDecodingMaxCorr=[];positionDecodingMaxCorr = table;       
        CellN_all=CellN_all+1;
% equalize and normalize the trial numbers and generatting psudo-trials#### 
         Bin_tr_Rate=[];
         smoothT=smooth1D(rate(:,2),smoothingRange,1);
         smoothC=smooth(rate(:,4+CellN),smoothingRange);
         rr=smoothC./smoothT;
         matC=reshape(rr,xbinNumber,length(rr)/xbinNumber)';
         C2=imresize(matC,[tr_imresize xbinNumber],'lanczos3');
         C2=matnorm(imresize(matC,[tr_imresize xbinNumber],'lanczos3'),2);
         Bin_tr_Rate(:,1:2)=bin_tr;
         Bin_tr_Rate(:,3)=reshape(C2',[xbinNumber*tr_imresize 1 ]);
%        Bin_tr_Rate=rate(:,[1,3,4+CellN]); % Original 2
%##########################################################################


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
            [correct_class_decision_values, normalized_rank_results, rank_confusion_matrix] = ...
            get_rank_and_decision_value_results_Farnaz(ytest, cl.labels, decision_values, create_rank_confusion_matrix);
            
  
            % estimate chance
            ytest=[];
            rr = randperm(length(rate_train));
            rrr = randperm(length(rate_test));
            cl = max_correlation_coefficient_CL;
            cl = train(cl,rate_train(rr),position_train);
            ytest=test(cl,rate_test(rrr));
            mse_chance_rate = mean((ytest-position_test).^2);
            
            % store results
            struct.mse_rate=mse_rate;
            struct.mse_chance_rate=mse_chance_rate;
            struct.CCDV_Mean=mean(correct_class_decision_values);
            struct.CCDV_Med=median(correct_class_decision_values);
            struct.NRR_Mean= mean(normalized_rank_results);
            struct.NRR_Med= median(normalized_rank_results);
            struct.smoothingRange = smoothingRange;
            struct.iter = iter;
            struct.Training_Percentile = Training_Percentile;
            struct.trialOrder = Shuffled_trial(:,iter)';
            
            positionDecodingMaxCorr = [positionDecodingMaxCorr;struct2table(struct)];
            rank_confusion_matrix_all=cat(3,rank_confusion_matrix_all,rank_confusion_matrix);
            
            
            clear var yfit_rate struct
            
        end
        DecodingMaxCorr_all{CellN_all} = positionDecodingMaxCorr;
    end
 
end
toc
%%
cd (dir)
% save(['DecodingSingleCell_C_TRD1.mat'],'DecodingMaxCorr_all')
% save(['rank_confusion_matrix_C_TRD1.mat'],'rank_confusion_matrix_all')
save(['DecodingSingleCell_P_TRD1.mat'],'DecodingMaxCorr_all')
save(['rank_confusion_matrix_P_TRD1.mat'],'rank_confusion_matrix_all')
cd ..

%%
M=mean(rank_confusion_matrix_all,3);
size(M)
figure
imagesc(matnorm(M,1))
colormap jet
%%
clc
clear
cd E:\Farnaz\chronic\N_LFP_CA1
load('DecodingSingleCell_P_TRD1.mat')
Table_P1=DecodingMaxCorr_all;
load('DecodingSingleCell_C_TRD1.mat')
Table_C1=DecodingMaxCorr_all;
cd ..
cd E:\Farnaz\chronic\N_LFP_CA3
load('DecodingSingleCell_P_TRD1.mat')
Table_P3=DecodingMaxCorr_all;
load('DecodingSingleCell_C_TRD1.mat')
Table_C3=DecodingMaxCorr_all;
cd ..

%%
    U_C1=[];
    colmn1=1; colmn2=2;colmn3=3;colmn4=4;colmn5=5;colmn6=6;
    Mt=Table_C1;
    for k=1:size(Mt,2)  
%             U_C1=[U_C1;mean(table2array(Mt{1,k}(:,colmn1))) mean(table2array(Mt{1,k}(:,colmn2)))...
%                     mean(table2array(Mt{1,k}(:,colmn3))) mean(table2array(Mt{1,k}(:,colmn4)))];
              U_C1=[U_C1;table2array(Mt{1,k}(:,colmn1)) table2array(Mt{1,k}(:,colmn2))...
                    table2array(Mt{1,k}(:,colmn3)) table2array(Mt{1,k}(:,colmn4)) table2array(Mt{1,k}(:,colmn5)) table2array(Mt{1,k}(:,colmn6))];
                
    end
     
    U_P1=[];
    Mt=Table_P1;
    for k=1:size(Mt,2)  
              U_P1=[U_P1;table2array(Mt{1,k}(:,colmn1)) table2array(Mt{1,k}(:,colmn2))...
                    table2array(Mt{1,k}(:,colmn3)) table2array(Mt{1,k}(:,colmn4)) table2array(Mt{1,k}(:,colmn5)) table2array(Mt{1,k}(:,colmn6))];
    end
    
    U_C3=[];
        Mt=Table_C3;
    for k=1:size(Mt,2)  
              U_C3=[U_C3;table2array(Mt{1,k}(:,colmn1)) table2array(Mt{1,k}(:,colmn2))...
                    table2array(Mt{1,k}(:,colmn3)) table2array(Mt{1,k}(:,colmn4)) table2array(Mt{1,k}(:,colmn5)) table2array(Mt{1,k}(:,colmn6))];
    end
    
    U_P3=[];
    Mt=Table_P3;
    for k=1:size(Mt,2)  
              U_P3=[U_P3;table2array(Mt{1,k}(:,colmn1)) table2array(Mt{1,k}(:,colmn2))...
                    table2array(Mt{1,k}(:,colmn3)) table2array(Mt{1,k}(:,colmn4)) table2array(Mt{1,k}(:,colmn5)) table2array(Mt{1,k}(:,colmn6))];
    end
    
%%

figure
plot(U_C1(:,4),'r')
hold on
plot(U_P1(:,4),'b')

%%
 Q={'mse-rate','mse-rate-chance','CCDV_Mean','CCDV_Med','NRR_Mean','NRR_Med','mse-phase-all-chance'};

Col_C=1;Col_P=1;

    CA1.C_1=U_C1(:,Col_C);%'mse-rate'
%     CA1.C_1(CA1.C_1==0)=[];
    
    CA1.P_1=U_P1(:,Col_P);%'mse-rate'
%     CA1.P_1(CA1.P_1==0)=[];
    
    CA1.C_2=U_C1(:,2)+500;%'mse-rate-chance'
%     CA1.C_2(CA1.C_2==0)=[];
    
    CA1.P_2=U_P1(:,2)+500;%'mse-rate-chance'
%     CA1.P_2(CA1.P_2==0)=[];

    CA3.C_1=U_C3(:,Col_C)-100;%'mse-phase'
%     CA3.C_1(CA3.C_1==0)=[];
    
    CA3.P_1=U_P3(:,Col_P)-100;%'mse-phase'
%     CA3.P_1(CA3.P_1==0)=[];
    
    CA3.C_2=U_C3(:,2)+500;%'mse-phase-all-chance'
%     CA3.C_2(CA3.C_2==0)=[];
    
    CA3.P_2=U_P3(:,2)+500;%'mse-phase-all-chance'
%     CA3.P_2(CA3.P_2==0)=[];

%     barplot_f(CA3,CA1,Q{Col_C})
Boxplot_f(CA3,CA1,Q{Col_C})
[Pa(1,1),Pk(1,1)] = Stest( CA1.C_1',CA1.P_1')%
[Pa(1,1),Pk(1,1)] = Stest( CA3.C_1',CA3.P_1')%
[Pa(1,1),Pk(1,1)] = Stest( CA1.C_1',CA3.C_1')%
[Pa(1,1),Pk(1,1)] = Stest( CA1.P_1',CA3.P_1')%
ylim([0 4500])










