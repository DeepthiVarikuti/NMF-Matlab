clear all;close all;

indir = '/Users/dvarikuti/Multistate/Projects/VBM/Generalisation_15Sep/FZJ1000_55to76AgeRanged_AllOther_mainanalysis';



components_change = [25,50,75,150,175,200,225,250,275,300,325,350,375,400,425,450,500,575,600,625,690];

for i_numcomp = 1:length(components_change)
    K = components_change(i_numcomp);
    
    load(fullfile(indir,num2str(K),'Preddat.mat'));
    AE.oneR(i_numcomp,:) = abs(mean(predmat.one.raw)-predmat.one.age');
    AE.twoR(i_numcomp,:) = abs(mean(predmat.two.raw)-predmat.two.age');
    AE.threeR(i_numcomp,:) = abs(predmat.three.raw-predmat.three.age);
    AE.fourR(i_numcomp,:) = abs(mean(predmat.four.raw)-predmat.four.age');
    AE.fiveR(i_numcomp,:) = abs(mean(predmat.five.raw)-predmat.five.age');
    AE.sixR(i_numcomp,:) = abs(predmat.six.raw-predmat.six.age);
    
    
    
    
    MAE.one(3,i_numcomp) = corr(mean(predmat.one.raw)',predmat.one.age);
    MAE.two(3,i_numcomp) = corr(mean(predmat.two.raw)',predmat.two.age);
    MAE.three(3,i_numcomp) = corr(predmat.three.raw,predmat.three.age);
    MAE.four(3,i_numcomp) = corr(mean(predmat.four.raw)',predmat.four.age);
    MAE.five(3,i_numcomp) = corr(mean(predmat.five.raw)',predmat.five.age);
    MAE.six(3,i_numcomp) = corr(predmat.six.raw,predmat.six.age);
    
    outAE(i_numcomp).oneR = find(AE.oneR(i_numcomp,:)>10);
    outAE(i_numcomp).twoR = find(AE.twoR(i_numcomp,:)>10);
    outAE(i_numcomp).threeR = find(AE.threeR(i_numcomp,:)>10);
    outAE(i_numcomp).fourR = find(AE.fourR(i_numcomp,:)>10);
    outAE(i_numcomp).fiveR = find(AE.fiveR(i_numcomp,:)>10);
    outAE(i_numcomp).sixR = find(AE.sixR(i_numcomp,:)>10);
    
    
    eachrep(i_numcomp).oneR = abs(bsxfun(@minus,predmat.one.raw',predmat.one.age));
    eachrep(i_numcomp).twoR = abs(bsxfun(@minus,predmat.two.raw',predmat.two.age));
    eachrep(i_numcomp).fourR = abs(bsxfun(@minus,predmat.four.raw',predmat.four.age));
    eachrep(i_numcomp).fiveR = abs(bsxfun(@minus,predmat.five.raw',predmat.five.age));
    
    
    
    minoutAE(i_numcomp).oneR = find(sum( eachrep(i_numcomp).oneR>4) == min(sum( eachrep(i_numcomp).oneR>4)));
     minoutAE(i_numcomp).twoR = find(sum( eachrep(i_numcomp).twoR>4) == min(sum( eachrep(i_numcomp).twoR>4)));
     minoutAE(i_numcomp).fourR = find(sum(eachrep(i_numcomp).fourR>4) == min(sum(eachrep(i_numcomp).fourR>4)));
     minoutAE(i_numcomp).fiveR = find(sum( eachrep(i_numcomp).fiveR>4) == min(sum( eachrep(i_numcomp).fiveR>4)));
     
    
    
    
    
    outAEAge(i_numcomp).oneR = predmat.one.age(outAE(i_numcomp).oneR);
    outAEAge(i_numcomp).twoR = predmat.two.age(outAE(i_numcomp).twoR);
    outAEAge(i_numcomp).threeR = predmat.three.age(outAE(i_numcomp).threeR);
    outAEAge(i_numcomp).fourR = predmat.four.age(outAE(i_numcomp).fourR);
    outAEAge(i_numcomp).fiveR = predmat.five.age(outAE(i_numcomp).fiveR);
    outAEAge(i_numcomp).sixR = predmat.six.age(outAE(i_numcomp).sixR);
    
    
        
%     outAEagerange(1).valR(1,i_numcomp)= min(outAEAge(i_numcomp).oneR);
%      outAEagerange(1).valR(2,i_numcomp)= max(outAEAge(i_numcomp).oneR);
%      
%       outAEagerange(2).valR(1,i_numcomp)= min(outAEAge(i_numcomp).twoR);
%      outAEagerange(2).valR(2,i_numcomp)= max(outAEAge(i_numcomp).twoR);
%      
%        outAEagerange(3).valR(1,i_numcomp)= min(outAEAge(i_numcomp).threeR);
%      outAEagerange(3).valR(2,i_numcomp)= max(outAEAge(i_numcomp).threeR);
%   
%        outAEagerange(4).valR(1,i_numcomp)= min(outAEAge(i_numcomp).fourR);
%      outAEagerange(4).valR(2,i_numcomp)= max(outAEAge(i_numcomp).fourR);
%      
%        outAEagerange(5).valR(1,i_numcomp)= min(outAEAge(i_numcomp).fiveR);
%      outAEagerange(5).valR(2,i_numcomp)= max(outAEAge(i_numcomp).fiveR);
%      
%        outAEagerange(6).valR(1,i_numcomp)= min(outAEAge(i_numcomp).sixR);
%      outAEagerange(6).valR(2,i_numcomp)= max(outAEAge(i_numcomp).sixR);
  
   
    
     outAElenR(1,i_numcomp) = length(outAE(i_numcomp).oneR)/length(predmat.one.age) * 100;
    outAElenR(2,i_numcomp) = length(outAE(i_numcomp).twoR)/length(predmat.two.age) * 100;
    outAElenR(3,i_numcomp) = length(outAE(i_numcomp).threeR)/length(predmat.three.age) * 100;
    outAElenR(4,i_numcomp) = length(outAE(i_numcomp).fourR)/length(predmat.four.age) * 100;
    outAElenR(5,i_numcomp) = length(outAE(i_numcomp).fiveR)/length(predmat.five.age) * 100;
    outAElenR(6,i_numcomp) = length(outAE(i_numcomp).sixR)/length(predmat.six.age) * 100;
    
    
    AE.oneA(i_numcomp,:) = abs(mean(predmat.one.adj)-predmat.one.age');
    AE.twoA(i_numcomp,:)= abs(mean(predmat.two.adj)-predmat.two.age');
    AE.threeA(i_numcomp,:) = abs(predmat.three.adj-predmat.three.age);
    AE.fourA(i_numcomp,:) = abs(mean(predmat.four.adj)-predmat.four.age');
    AE.fiveA(i_numcomp,:) = abs(mean(predmat.five.adj)-predmat.five.age');
    AE.sixA(i_numcomp,:) = abs(predmat.six.adj-predmat.six.age);
    
    
    MAE.one(4,i_numcomp) = corr(mean(predmat.one.adj)',predmat.one.age);
    MAE.two(4,i_numcomp) = corr(mean(predmat.two.adj)',predmat.two.age);
    MAE.three(4,i_numcomp) = corr(predmat.three.adj,predmat.three.age);
    MAE.four(4,i_numcomp) = corr(mean(predmat.four.adj)',predmat.four.age);
    MAE.five(4,i_numcomp) = corr(mean(predmat.five.adj)',predmat.five.age);
    MAE.six(4,i_numcomp) = corr(predmat.six.adj,predmat.six.age);
    
    outAE(i_numcomp).oneA = find(AE.oneA(i_numcomp,:)>10);
    outAE(i_numcomp).twoA = find(AE.twoA(i_numcomp,:)>10);
    outAE(i_numcomp).threeA = find(AE.threeA(i_numcomp,:)>10);
    outAE(i_numcomp).fourA = find(AE.fourA(i_numcomp,:)>10);
    outAE(i_numcomp).fiveA = find(AE.fiveA(i_numcomp,:)>10);
    outAE(i_numcomp).sixA = find(AE.sixA(i_numcomp,:)>10);
    
     outAEAge(i_numcomp).oneA = predmat.one.age(outAE(i_numcomp).oneA);
    outAEAge(i_numcomp).twoA = predmat.two.age(outAE(i_numcomp).twoA);
    outAEAge(i_numcomp).threeA = predmat.three.age(outAE(i_numcomp).threeA);
    outAEAge(i_numcomp).fourA = predmat.four.age(outAE(i_numcomp).fourA);
    outAEAge(i_numcomp).fiveA = predmat.five.age(outAE(i_numcomp).fiveA);
    outAEAge(i_numcomp).sixA = predmat.six.age(outAE(i_numcomp).sixA);
   
    
%     outAEagerange(1).val(1,i_numcomp)= min(outAEAge(i_numcomp).oneA);
%      outAEagerange(1).val(2,i_numcomp)= max(outAEAge(i_numcomp).oneA);
%      
%       outAEagerange(2).val(1,i_numcomp)= min(outAEAge(i_numcomp).twoA);
%      outAEagerange(2).val(2,i_numcomp)= max(outAEAge(i_numcomp).twoA);
%      
%        outAEagerange(3).val(1,i_numcomp)= min(outAEAge(i_numcomp).threeA);
%      outAEagerange(3).val(2,i_numcomp)= max(outAEAge(i_numcomp).threeA);
%   
%        outAEagerange(4).val(1,i_numcomp)= min(outAEAge(i_numcomp).fourA);
%      outAEagerange(4).val(2,i_numcomp)= max(outAEAge(i_numcomp).fourA);
%      
%        outAEagerange(5).val(1,i_numcomp)= min(outAEAge(i_numcomp).fiveA);
%      outAEagerange(5).val(2,i_numcomp)= max(outAEAge(i_numcomp).fiveA);
%      
%        outAEagerange(6).val(1,i_numcomp)= min(outAEAge(i_numcomp).sixA);
%      outAEagerange(6).val(2,i_numcomp)= max(outAEAge(i_numcomp).sixA);
  
  
  
  
    
      outAElenA(1,i_numcomp) = length(outAE(i_numcomp).oneA)/length(predmat.one.age) * 100;
    outAElenA(2,i_numcomp) = length(outAE(i_numcomp).twoA)/length(predmat.two.age) * 100;
    outAElenA(3,i_numcomp) = length(outAE(i_numcomp).threeA)/length(predmat.three.age) * 100;
    outAElenA(4,i_numcomp) = length(outAE(i_numcomp).fourA)/length(predmat.four.age) * 100;
    outAElenA(5,i_numcomp) = length(outAE(i_numcomp).fiveA)/length(predmat.five.age) * 100;
    outAElenA(6,i_numcomp) = length(outAE(i_numcomp).sixA)/length(predmat.six.age) * 100;
    
    MAE.one(1,i_numcomp) = median(AE.oneR(i_numcomp,:));
    MAE.two(1,i_numcomp) = median(AE.twoR(i_numcomp,:));
    MAE.three(1,i_numcomp) = median(AE.threeR(i_numcomp,:));
    MAE.four(1,i_numcomp) = median(AE.fourR(i_numcomp,:));
    MAE.five(1,i_numcomp) = median(AE.fiveR(i_numcomp,:));
    MAE.six(1,i_numcomp) = median(AE.sixR(i_numcomp,:));
    
    MAE.one(2,i_numcomp) = median(AE.oneA(i_numcomp,:));
    MAE.two(2,i_numcomp) = median(AE.twoA(i_numcomp,:));
    MAE.three(2,i_numcomp) = median(AE.threeA(i_numcomp,:));
    MAE.four(2,i_numcomp) = median(AE.fourA(i_numcomp,:));
    MAE.five(2,i_numcomp) = median(AE.fiveA(i_numcomp,:));
    MAE.six(2,i_numcomp) = median(AE.sixA(i_numcomp,:));
    
    
     %%%%%%%%% root mean square error
    
    rmsAE.one(1,i_numcomp) = sqrt(mean(AE.oneR(i_numcomp,:).^2));
    rmsAE.two(1,i_numcomp) = sqrt(mean(AE.twoR(i_numcomp,:).^2));
    rmsAE.three(1,i_numcomp) = sqrt(mean(AE.threeR(i_numcomp,:).^2));
    rmsAE.four(1,i_numcomp) = sqrt(mean(AE.fourR(i_numcomp,:).^2));
    rmsAE.five(1,i_numcomp) = sqrt(mean(AE.fiveR(i_numcomp,:).^2));
    rmsAE.six(1,i_numcomp) = sqrt(mean(AE.sixR(i_numcomp,:).^2));
    
    rmsAE.one(2,i_numcomp) = sqrt(mean(AE.oneA(i_numcomp,:).^2));
    rmsAE.two(2,i_numcomp) = sqrt(mean(AE.twoA(i_numcomp,:).^2));
    rmsAE.three(2,i_numcomp) = sqrt(mean(AE.threeA(i_numcomp,:).^2));
    rmsAE.four(2,i_numcomp) = sqrt(mean(AE.fourA(i_numcomp,:).^2));
    rmsAE.five(2,i_numcomp) = sqrt(mean(AE.fiveA(i_numcomp,:).^2));
    rmsAE.six(2,i_numcomp) = sqrt(mean(AE.sixA(i_numcomp,:).^2));
    
    %%%%%%%%%%%%%%%%%%% mean abs error
    
     MeanAE.one(1,i_numcomp) = mean(AE.oneR(i_numcomp,:));
    MeanAE.two(1,i_numcomp) = mean(AE.twoR(i_numcomp,:));
    MeanAE.three(1,i_numcomp) = mean(AE.threeR(i_numcomp,:));
    MeanAE.four(1,i_numcomp) = mean(AE.fourR(i_numcomp,:));
    MeanAE.five(1,i_numcomp) = mean(AE.fiveR(i_numcomp,:));
    MeanAE.six(1,i_numcomp) = mean(AE.sixR(i_numcomp,:));
    
    MeanAE.one(2,i_numcomp) = mean(AE.oneA(i_numcomp,:));
    MeanAE.two(2,i_numcomp) = mean(AE.twoA(i_numcomp,:));
    MeanAE.three(2,i_numcomp) = mean(AE.threeA(i_numcomp,:));
    MeanAE.four(2,i_numcomp) = mean(AE.fourA(i_numcomp,:));
    MeanAE.five(2,i_numcomp) = mean(AE.fiveA(i_numcomp,:));
    MeanAE.six(2,i_numcomp) = mean(AE.sixA(i_numcomp,:));
   
    
    disp(K);
    
    %          for i_rep = 1:length(predmat.one.repimat.GlmBeta(1,1,:))
    %              for i_fold = 1:length(predmat.one.repimat.GlmBeta(1,:,1))
    %                  logifunc = predmat.one.repimat.GlmBeta(:,i_fold,i_rep)~=0;
    %                  Repeatedcomp(i_rep).fold(i_fold).ind = find(logifunc==1);
    %                  %Repeatedcomp(i_rep,i_fold).ind = find(logifunc==1);
    %
    %                  %predmat.one.repimat.GlmBeta(:,)
    %              end
    %          end
    
    clear logifunc;
    for i_rep = 1:length(predmat.one.repimat.GlmBeta(1,1,:))
        logifunc.one(:,i_rep) = sum(predmat.one.repimat.GlmBeta(:,:,i_rep)~=0,2);
        logifunc.two(:,i_rep) = sum(predmat.two.repimat.GlmBeta(:,:,i_rep)~=0,2);
        logifunc.four(:,i_rep) = sum(predmat.four.repimat.GlmBeta(:,:,i_rep)~=0,2);
        logifunc.five(:,i_rep) = sum(predmat.five.repimat.GlmBeta(:,:,i_rep)~=0,2);
    end
    
    logifunc.three = find((predmat.three.repimat.GlmBeta ~=0) == 1);
    logifunc.six = find((predmat.six.repimat.GlmBeta~=0) ==1);
    
    thrshold = (length(predmat.one.repimat.GlmBeta(1,1,:))-1);
    
    Cmn.one(i_numcomp).ind = find(sum(logifunc.one>8,2)>thrshold);
    Cmn.two(i_numcomp).ind = find(sum(logifunc.two>8,2)>thrshold);
    Cmn.three(i_numcomp).ind = logifunc.three;
    Cmn.four(i_numcomp).ind = find(sum(logifunc.four>8,2)>thrshold);
    Cmn.five(i_numcomp).ind = find(sum(logifunc.five>8,2)>thrshold);
    Cmn.six(i_numcomp).ind= logifunc.six;
    
    
    
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% save componenent in nifti image FZJ100
%     max_iter = 50000; initMeth = 1;
%     localdir = '/Users/dvarikuti/Multistate';
%     MainprjtWord = 'OPNMF_693_outliersremoved_NoScaling';%'OPNMF_allOther'; %'OPNMF_693_outliersremoved_NoScaling';
%     subkeyword = 'FZJ1000_55to76AgeRanged_Dec2015_RSandVBMandDWI';%'AllOther_VBM';%'FZJ1000_55to76AgeRanged_Dec2015_RSandVBMandDWI'
%     
%     cd(fullfile(localdir,'Projects','VBM',MainprjtWord,subkeyword));
%     
%     load NNMF2VBM
%     project_word = ['result_',num2str(K),'_',num2str(max_iter),'_',num2str(initMeth)];
%     projectDir = fullfile('/Users/dvarikuti/Multistate','Projects','VBM',MainprjtWord,subkeyword,project_word);
%     filename = 'Components_used_atleat8fold_25repi.nii';
%     cd(projectDir);
%     load NNMF_W_H
%     
%    
%     NewW = zeros(size(W));
%     NewW(:,Cmn.one(:,i_numcomp).ind) = W(:,Cmn.one(:,i_numcomp).ind);
%     new_avg_RsErr = zeros(121,145,121);
%     [~, MPM] = yael_kmax (single(NewW)',1);
%     
%     %rmdir(fullfile(indir,num2str(K),predmat.one.name),'s');
%     mkdir(fullfile(indir,num2str(K),predmat.one.name));
%     
% %     MPM(MPM~=1)=0;
%     MPM(MPM==1)=0;
%     idx = sub2ind(size(new_avg_RsErr),XYZ(:,1),XYZ(:,2),XYZ(:,3));
%     new_XX = zeros(121,145,121);
%     new_XX(idx) = MPM;
%     Vos = rmfield(Vgm,'pinfo');
%     Vos.fname = fullfile(indir,num2str(K),predmat.one.name,filename);
%     Vos = spm_write_vol(Vos,new_XX);
%     
%     
%     NewW = zeros(size(W));
%     NewW(:,Cmn.three(:,i_numcomp).ind) = W(:,Cmn.three(:,i_numcomp).ind);
%     new_avg_RsErr = zeros(121,145,121);
%     [~, MPM] = yael_kmax (single(NewW)',1);
%     
%     %rmdir(fullfile(indir,num2str(K),predmat.three.name),'s');
%     mkdir(fullfile(indir,num2str(K),predmat.three.name));
%     
% %     MPM(MPM~=1)=0;
%     MPM(MPM==1)=0;
%     idx = sub2ind(size(new_avg_RsErr),XYZ(:,1),XYZ(:,2),XYZ(:,3));
%     new_XX = zeros(121,145,121);
%     new_XX(idx) = MPM;
%     Vos = rmfield(Vgm,'pinfo');
%     Vos.fname = fullfile(indir,num2str(K),predmat.three.name,filename);
%     Vos = spm_write_vol(Vos,new_XX);
%     
%     NewW = zeros(size(W));
%     NewW(:,Cmn.five(:,i_numcomp).ind) = W(:,Cmn.five(:,i_numcomp).ind);
%     new_avg_RsErr = zeros(121,145,121);
%     [~, MPM] = yael_kmax (single(NewW)',1);
%     
%     %rmdir(fullfile(indir,num2str(K),predmat.five.name),'s');
%     mkdir(fullfile(indir,num2str(K),predmat.five.name));
%     
%     MPM(MPM==1)=0;
% %     MPM(MPM~=1)=0;
%     idx = sub2ind(size(new_avg_RsErr),XYZ(:,1),XYZ(:,2),XYZ(:,3));
%     new_XX = zeros(121,145,121);
%     new_XX(idx) = MPM;
%     Vos = rmfield(Vgm,'pinfo');
%     Vos.fname = fullfile(indir,num2str(K),predmat.five.name,filename);
%     Vos = spm_write_vol(Vos,new_XX);
%     
%     
%      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% save componenent in nifti image FZJ100
%     max_iter = 50000; initMeth = 1;
%     localdir = '/Users/dvarikuti/Multistate';
%     MainprjtWord = 'OPNMF_allOther';%'OPNMF_allOther'; %'OPNMF_693_outliersremoved_NoScaling';
%     subkeyword = 'AllOther_VBM';%'AllOther_VBM';%'FZJ1000_55to76AgeRanged_Dec2015_RSandVBMandDWI'
%     
%     cd(fullfile(localdir,'Projects','VBM',MainprjtWord,subkeyword));
%     
%     load NNMF2VBM
%     project_word = ['result_',num2str(K),'_',num2str(max_iter),'_',num2str(initMeth)];
%     projectDir = fullfile(localdir,'Projects','VBM',MainprjtWord,subkeyword,project_word);
%     cd(projectDir);
%     load NNMF_W_H
%     
%    
%     NewW = zeros(size(W));
%     NewW(:,Cmn.two(:,i_numcomp).ind) = W(:,Cmn.two(:,i_numcomp).ind);
%     new_avg_RsErr = zeros(121,145,121);
%     [~, MPM] = yael_kmax (single(NewW)',1);
%     
% %     rmdir(fullfile(indir,num2str(K),predmat.two.name),'s');
%     mkdir(fullfile(indir,num2str(K),predmat.two.name));
%     
%     MPM(MPM==1)=0;
% %     MPM(MPM~=1)=0;
%     idx = sub2ind(size(new_avg_RsErr),XYZ(:,1),XYZ(:,2),XYZ(:,3));
%     new_XX = zeros(121,145,121);
%     new_XX(idx) = MPM;
%     Vos = rmfield(Vgm,'pinfo');
%     Vos.fname = fullfile(indir,num2str(K),predmat.two.name,filename);
%     Vos = spm_write_vol(Vos,new_XX);
%     
%     
%     NewW = zeros(size(W));
%     NewW(:,Cmn.four(:,i_numcomp).ind) = W(:,Cmn.four(:,i_numcomp).ind);
%     new_avg_RsErr = zeros(121,145,121);
%     [~, MPM] = yael_kmax (single(NewW)',1);
%     
% %     rmdir(fullfile(indir,num2str(K),predmat.four.name),'s');
%     mkdir(fullfile(indir,num2str(K),predmat.four.name));
%     
%     MPM(MPM==1)=0;
% %     MPM(MPM~=1)=0;
%     idx = sub2ind(size(new_avg_RsErr),XYZ(:,1),XYZ(:,2),XYZ(:,3));
%     new_XX = zeros(121,145,121);
%     new_XX(idx) = MPM;
%     Vos = rmfield(Vgm,'pinfo');
%     Vos.fname = fullfile(indir,num2str(K),predmat.four.name,filename);
%     Vos = spm_write_vol(Vos,new_XX);
%     
%     NewW = zeros(size(W));
%     NewW(:,Cmn.six(:,i_numcomp).ind) = W(:,Cmn.six(:,i_numcomp).ind);
%     new_avg_RsErr = zeros(121,145,121);
%     [~, MPM] = yael_kmax (single(NewW)',1);
%     
% %    rmdir(fullfile(indir,num2str(K),predmat.six.name),'s');
%     mkdir(fullfile(indir,num2str(K),predmat.six.name));
%     
%     MPM(MPM==1)=0;
% %     MPM(MPM~=1)=0;
%     idx = sub2ind(size(new_avg_RsErr),XYZ(:,1),XYZ(:,2),XYZ(:,3));
%     new_XX = zeros(121,145,121);
%     new_XX(idx) = MPM;
%     Vos = rmfield(Vgm,'pinfo');
%     Vos.fname = fullfile(indir,num2str(K),predmat.six.name,filename);
%     Vos = spm_write_vol(Vos,new_XX);
%     
    
    
end

figure;
subplot(1,2,1);plot(components_change,outAElenR); xlabel('Components'); ylabel('percentage of subjects with prediction error >10 years'); title('Without slope adjustment');
ylim([0 70]);
legend('FZJ1000 components test n train FZJ100','FZJ1000 components test and train ALl other','FZJ1000 components test FZJ1000 and train Allother','ALLother components test and train all other','ALLother components test and train FZJ1000','ALLother components test all other and train FZJ1000');
subplot(1,2,2);plot(components_change,outAElenA); xlabel('Components'); ylabel('percentage of subjects with prediction error > 10 years'); title('With slope adjustment');
ylim([0 70]);
%legend('FZJ1000 components test n train FZJ100','FZJ1000 components test and train ALl other','FZJ1000 components test FZJ1000 and train Allother','ALLother components test and train all other','ALLother components test and train FZJ1000','ALLother components test all other and train FZJ1000');

maxval = max([max(max(MAE.one(1:2,:))),max(max(MAE.two(1:2,:))),max(max(MAE.three(1:2,:))),max(max(MAE.four(1:2,:))),max(max(MAE.five(1:2,:))),max(max(MAE.six(1:2,:)))]);
minval = min([min(min(MAE.one(1:2,:))),min(min(MAE.two(1:2,:))),min(min(MAE.three(1:2,:))),min(min(MAE.four(1:2,:))),min(min(MAE.five(1:2,:))),min(min(MAE.six(1:2,:)))]);


% maxval = max([max(MAE.oneR),max(MAE.twoR),max(MAE.threeR),max(MAE.fourR),max(MAE.fiveR),max(MAE.sixR),max(MAE.oneA),max(MAE.twoA),max(MAE.threeA),max(MAE.fourA),max(MAE.fiveA),max(MAE.sixA)]);
% minval = min([min(MAE.oneR),min(MAE.twoR),min(MAE.threeR),min(MAE.fourR),min(MAE.fiveR),min(MAE.sixR),min(MAE.oneA),min(MAE.twoA),min(MAE.threeA),min(MAE.fourA),min(MAE.fiveA),min(MAE.sixA)]);

H = figure; clf;
subplot(2,3,1); hold on;plot(components_change,MAE.one(1,:),'b.-'); plot(components_change,MAE.one(2,:),'r.-'); title('compdiv: FZJ100, FZ1000 (train & test same data)'); xlabel('NMF Components'); ylabel('Median absolute error'); ylim([minval-1 maxval+1]);
subplot(2,3,2); hold on; plot(components_change,MAE.two(1,:),'b.-');plot(components_change,MAE.two(2,:),'r.-'); title('compdiv: FZJ100, AllOther (train & test same data)'); xlabel('NMF Components'); ylabel('Median absolute error');ylim([minval-1 maxval+1]);
subplot(2,3,3); hold on; plot(components_change,MAE.three(1,:),'b.-');plot(components_change,MAE.three(2,:),'r.-');title('compdiv: FZJ100, train on FZJ100 & test on AllOther'); xlabel('NMF Components'); ylabel('Median absolute error'); ylim([minval-1 maxval+1]);
subplot(2,3,5); hold on; plot(components_change,MAE.four(1,:),'b.-');plot(components_change,MAE.four(2,:),'r.-');title('compdiv: AllOther , AllOther  (train & test same data)'); xlabel('NMF Components'); ylabel('Median absolute error'); ylim([minval-1 maxval+1]);
subplot(2,3,4); hold on; plot(components_change,MAE.five(1,:),'b.-');plot(components_change,MAE.five(2,:),'r.-');title('compdiv: AllOther , FZ1000 (train & test same data)'); xlabel('NMF Components'); ylabel('Median absolute error');ylim([minval-1 maxval+1]);
subplot(2,3,6); hold on; plot(components_change,MAE.six(1,:),'b.-');plot(components_change,MAE.six(2,:),'r.-');title('compdiv: AllOther, train on AllOther & test on FZJ1000'); xlabel('NMF Components'); ylabel('Median absolute error'); ylim([minval-1 maxval+1]);legend('without slope adj', 'with slope adj');


maxval = max([max(max(MAE.one(3:4,:))),max(max(MAE.two(3:4,:))),max(max(MAE.three(3:4,:))),max(max(MAE.four(3:4,:))),max(max(MAE.five(3:4,:))),max(max(MAE.six(3:4,:)))]);
minval = min([min(min(MAE.one(3:4,:))),min(min(MAE.two(3:4,:))),min(min(MAE.three(3:4,:))),min(min(MAE.four(3:4,:))),min(min(MAE.five(3:4,:))),min(min(MAE.six(3:4,:)))]);


H = figure; clf;
subplot(2,3,1); hold on;plot(components_change,MAE.one(3,:),'b.-'); plot(components_change,MAE.one(4,:),'r.-'); title('compdiv: FZJ100, FZ1000 (train & test same data)'); xlabel('NMF Components'); ylabel('Correlations'); ylim([minval-0.1 maxval+0.1]);
subplot(2,3,2); hold on; plot(components_change,MAE.two(3,:),'b.-');plot(components_change,MAE.two(4,:),'r.-'); title('compdiv: FZJ100, AllOther (train & test same data)'); xlabel('NMF Components'); ylabel('Correlations');ylim([minval-0.1 maxval+0.1]);
subplot(2,3,3); hold on; plot(components_change,MAE.three(3,:),'b.-');plot(components_change,MAE.three(4,:),'r.-');title('compdiv: FZJ100, train on FZJ100 & test on AllOther'); xlabel('NMF Components'); ylabel('Correlations'); ylim([minval-0.1 maxval+0.1]);
subplot(2,3,5); hold on; plot(components_change,MAE.four(3,:),'b.-');plot(components_change,MAE.four(4,:),'r.-');title('compdiv: AllOther , AllOther  (train & test same data)'); xlabel('NMF Components'); ylabel('Correlations'); ylim([minval-0.1 maxval+0.1]);
subplot(2,3,4); hold on; plot(components_change,MAE.five(3,:),'b.-');plot(components_change,MAE.five(4,:),'r.-');title('compdiv: AllOther , FZ1000 (train & test same data)'); xlabel('NMF Components'); ylabel('Correlations');ylim([minval-0.1 maxval+0.1]);
subplot(2,3,6); hold on; plot(components_change,MAE.six(3,:),'b.-');plot(components_change,MAE.six(4,:),'r.-');title('compdiv: AllOther, train on AllOther & test on FZJ1000'); xlabel('NMF Components'); ylabel('Correlations'); ylim([minval-0.1 maxval+0.1]);legend('without slope adj', 'with slope adj');


% H = figure; clf;
% subplot(2,3,1); hold on; bar(components_change,MAE.one'); title('compdiv: FZJ100, FZ1000 (train & test same data)'); xlabel('NMF Components'); ylabel('Median absolute error'); ylim([minval maxval+1]);
% subplot(2,3,2); hold on; bar(components_change,MAE.two'); title('compdiv: FZJ100, AllOther (train & test same data)'); xlabel('NMF Components'); ylabel('Median absolute error');ylim([minval maxval+1]);
% subplot(2,3,3); hold on; bar(components_change,MAE.three');title('compdiv: FZJ100, train on FZJ100 & test on AllOther'); xlabel('NMF Components'); ylabel('Median absolute error'); ylim([minval maxval+1]);
% subplot(2,3,5); hold on; bar(components_change,MAE.four'); title('compdiv: AllOther , AllOther  (train & test same data)'); xlabel('NMF Components'); ylabel('Median absolute error'); ylim([minval maxval+1]);
% subplot(2,3,4); hold on; bar(components_change,MAE.five');title('compdiv: AllOther , FZ1000 (train & test same data)'); xlabel('NMF Components'); ylabel('Median absolute error');ylim([minval maxval+1]);
% subplot(2,3,6); hold on; bar(components_change,MAE.six'); title('compdiv: AllOther, train on AllOther & test on FZJ1000'); xlabel('NMF Components'); ylabel('Median absolute error'); ylim([minval maxval+1]);legend('without slope adj', 'with slope adj');



maxval = max([max(max(rmsAE.one(1:2,:))),max(max(rmsAE.two(1:2,:))),max(max(rmsAE.three(1:2,:))),max(max(rmsAE.four(1:2,:))),max(max(rmsAE.five(1:2,:))),max(max(rmsAE.six(1:2,:)))]);
minval = min([min(min(rmsAE.one(1:2,:))),min(min(rmsAE.two(1:2,:))),min(min(rmsAE.three(1:2,:))),min(min(rmsAE.four(1:2,:))),min(min(rmsAE.five(1:2,:))),min(min(rmsAE.six(1:2,:)))]);

H = figure; clf;
subplot(2,3,1); hold on;plot(components_change,rmsAE.one(1,:),'b.-'); plot(components_change,rmsAE.one(2,:),'r.-'); title('compdiv: FZJ100, FZ1000 (train & test same data)'); xlabel('NMF Components'); ylabel('root mean square error'); ylim([minval-1 maxval+1]);
subplot(2,3,2); hold on; plot(components_change,rmsAE.two(1,:),'b.-');plot(components_change,rmsAE.two(2,:),'r.-'); title('compdiv: FZJ100, AllOther (train & test same data)'); xlabel('NMF Components'); ylabel('root mean square error');ylim([minval-1 maxval+1]);
subplot(2,3,3); hold on; plot(components_change,rmsAE.three(1,:),'b.-');plot(components_change,rmsAE.three(2,:),'r.-');title('compdiv: FZJ100, train on FZJ100 & test on AllOther'); xlabel('NMF Components'); ylabel('root mean square error'); ylim([minval-1 maxval+1]);
subplot(2,3,5); hold on; plot(components_change,rmsAE.four(1,:),'b.-');plot(components_change,rmsAE.four(2,:),'r.-');title('compdiv: AllOther , AllOther  (train & test same data)'); xlabel('NMF Components'); ylabel('root mean square error'); ylim([minval-1 maxval+1]);
subplot(2,3,4); hold on; plot(components_change,rmsAE.five(1,:),'b.-');plot(components_change,rmsAE.five(2,:),'r.-');title('compdiv: AllOther , FZ1000 (train & test same data)'); xlabel('NMF Components'); ylabel('root mean square error');ylim([minval-1 maxval+1]);
subplot(2,3,6); hold on; plot(components_change,rmsAE.six(1,:),'b.-');plot(components_change,rmsAE.six(2,:),'r.-');title('compdiv: AllOther, train on AllOther & test on FZJ1000'); xlabel('NMF Components'); ylabel('root mean square error'); ylim([minval-1 maxval+1]);legend('without slope adj', 'with slope adj');

maxval = max([max(max(MeanAE.one(1:2,:))),max(max(MeanAE.two(1:2,:))),max(max(MeanAE.three(1:2,:))),max(max(MeanAE.four(1:2,:))),max(max(MeanAE.five(1:2,:))),max(max(MeanAE.six(1:2,:)))]);
minval = min([min(min(MeanAE.one(1:2,:))),min(min(MeanAE.two(1:2,:))),min(min(MeanAE.three(1:2,:))),min(min(MeanAE.four(1:2,:))),min(min(MeanAE.five(1:2,:))),min(min(MeanAE.six(1:2,:)))]);

H = figure; clf;
subplot(2,3,1); hold on;plot(components_change,MeanAE.one(1,:),'b.-'); plot(components_change,MeanAE.one(2,:),'r.-'); title('compdiv: FZJ100, FZ1000 (train & test same data)'); xlabel('NMF Components'); ylabel('mean abs error'); ylim([minval-1 maxval+1]);
subplot(2,3,2); hold on; plot(components_change,MeanAE.two(1,:),'b.-');plot(components_change,MeanAE.two(2,:),'r.-'); title('compdiv: FZJ100, AllOther (train & test same data)'); xlabel('NMF Components'); ylabel('mean abs error');ylim([minval-1 maxval+1]);
subplot(2,3,3); hold on; plot(components_change,MeanAE.three(1,:),'b.-');plot(components_change,MeanAE.three(2,:),'r.-');title('compdiv: FZJ100, train on FZJ100 & test on AllOther'); xlabel('NMF Components'); ylabel('mean abs error'); ylim([minval-1 maxval+1]);
subplot(2,3,5); hold on; plot(components_change,MeanAE.four(1,:),'b.-');plot(components_change,MeanAE.four(2,:),'r.-');title('compdiv: AllOther , AllOther  (train & test same data)'); xlabel('NMF Components'); ylabel('mean abs error'); ylim([minval-1 maxval+1]);
subplot(2,3,4); hold on; plot(components_change,MeanAE.five(1,:),'b.-');plot(components_change,MeanAE.five(2,:),'r.-');title('compdiv: AllOther , FZ1000 (train & test same data)'); xlabel('NMF Components'); ylabel('mean abs error');ylim([minval-1 maxval+1]);
subplot(2,3,6); hold on; plot(components_change,MeanAE.six(1,:),'b.-');plot(components_change,MeanAE.six(2,:),'r.-');title('compdiv: AllOther, train on AllOther & test on FZJ1000'); xlabel('NMF Components'); ylabel('mean abs error'); ylim([minval-1 maxval+1]);legend('without slope adj', 'with slope adj');

