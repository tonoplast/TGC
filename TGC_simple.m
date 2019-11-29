close all; clear;

eeglab;
filterfolder ='E:\MATLAB\R2018b\toolbox\signal\signal\'; %% To use Matlab filter and not Fieldtrip one
% caploc='F:\Data\EEG\standard-10-5-cap385.elp'; %path containing electrode positions

% Path (loading / saving)
inPath = 'D:\3_Frequency\thetagamma\H301\'; %% Directory of your cleaned file
outPath = 'D:\3_Frequency\thetagamma\H301\'; %% Saving path

   for h=1:5 %To run Theta-gamma coupling 5 times --> Depending on your need
       
    %% load files (cleaned data) here %%
    EEG = pop_loadset('filename','004_postmarainter_epochcorrect.set','filepath',inPath);
    
    % setting SAMPLING RATE here %
    set_srate = 250;
    label_multiplier = 1000 / set_srate; % to be used later for edge effect %
    
    EEG = pop_resample( EEG, set_srate); %% down-sampling should reduce time it takes.
    [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG,  1);


    %%
    
    frontchan = {'F3','Fz','F4'}; %% channel of interest 1
    backchan = {'P1','Pz','P3'}; %% channel of interest 2
    
%     frontchan = {'Fz'}; %% channel of interest 1
%     backchan = {'Pz'}; %% channel of interest 2
    
    %%
    %%%
    [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    THISCHAN=find(ismember({EEG.chanlocs.labels},backchan)); %% gamma
    THATCHAN=find(ismember({EEG.chanlocs.labels},frontchan)); %% theta
    [ch,pnts,eps]=size(EEG.data);
    

    %% selecting epochs and shuffling %%%
    I=1:eps;
    I= shuffle(I);
    keep=[I(1:round(eps*0.5))]; %%% keeping 10 epochs (consider changing as needed - better to have longer data)
%      keep=[I(1:30)]; 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%% PAC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    try
        %these=unique([EEG.event(:).epoch]);
        these=keep;
    catch
        disp('Not enough epochs');
        these=1:eps; % use as many as i got
    end
    
    ALLTRIALS=EEG;
    EEG = pop_select( EEG,'trial',these);
    
    %%% all trials or use above
    % EEG=ALLTRIALS;
    EEG = eeg_checkset( EEG );
    
    clear glmf mi plv vtk
%     [c,m,n]= size([EEG.data(THISCHAN,:,:)]);
    [ch,pnts,eps]=size(EEG.data([THISCHAN],:,:));
    
    %% DATA zero padding and reshaping
    
    datarange=[1:EEG.pnts];
    
    data1 = EEG.data([THISCHAN],datarange,:);
    data2 = EEG.data([THATCHAN],datarange,:);
    
    %%% meaning the data if more than 1 chan
    data1 = mean(data1,1);
    data2 = mean(data2,1);
    
    label = zeros(1,pnts,eps); %% zero padding labelling for edge effect problem
    label(:,[1:(500/label_multiplier),end-(500/label_multiplier):end],:)=1; %% edge effect removal step (taking sampling rate into account)
    
    dat1 = squeeze(data1)';
    dat2 = squeeze(data2)';
    label=squeeze(label)';
    
    
    [m,n]=size(dat1);
    data1 =reshape([dat1'],[1,m*n]);
    data2= reshape([dat2'],[1,m*n]);
    label= logical(reshape([label'],[1,m*n]));
    label_ori = label;
    
    % alternative label %
    label_alt = zeros(1,length(data1));
    label_alt(:,[1:(2000/label_multiplier),end-(2000/label_multiplier):end])=1;
    label = logical(label_alt);

    select=[];
    srate=EEG.srate;
    
        %% Filtering (adaptive - script attached)
    %srate=1000
    cfg =[];
    cfg.sr= srate;
    cfg.lo_bounds= [3 9];
    cfg.lo_step=0.1;
    cfg.lo_bandwidth=2;
    
    cfg.hi_bounds= [20 70];
    cfg.hi_step= 1;
    cfg. hi_bandwidth='adaptive';
    [gamma,theta]=setup_adaptivefilterbands(cfg);
    
    %cd(filterfolder);
    %% PAC estimation --> calculates four different methods
    
    bb=size(gamma,2);
    cc=size(theta,2);
    aa = size(data1,3);
   
    tic
    
    clear glmf mi plv vtk
    [glmf,plv,vtk,mi] = deal(zeros(bb,cc,aa)); % -- initialize output matrices
    
    for a = 1:aa
        raw_signal = data1(1,:,a);
        raw_signal2=data2(1,:,a);
        label=logical(label);
        
        ntimepoints = [length(raw_signal)]- [sum(label)];
        for b = 1:bb
            [x_gamma,y_gamma]=butter(2,[gamma(1,b) gamma(2,b)]/(srate/2),'bandpass');
            gamma_wave= filtfilt(x_gamma,y_gamma, double(raw_signal));
            gamma_z = hilbert(gamma_wave);
            gamma_amp= abs(gamma_z);
            
            for c = 1:cc
                [x_theta,y_theta]=butter(2,[theta(1,c) theta(2,c)]/(srate/2),'bandpass');
                theta_wave= filtfilt(x_theta,y_theta, double(raw_signal2));
                theta_z=hilbert(theta_wave);
                theta_phase=angle(theta_z);
                
                gamma_ampfilt= filtfilt(x_theta,y_theta, double(gamma_amp));
                
                gamma_amp_z=hilbert(gamma_ampfilt);
                gamma_amp_phase = angle(gamma_amp_z);
                %%%%%%%%%%
                
                thetaphase = theta_phase(~label);
                gammapow = gamma_amp(~label);
                gammapow1=gamma_ampfilt(~label);
                gammapowfilt =gamma_ampfilt(~label);
                gammapowphase=gamma_amp_phase(~label);
                nbins = 18;
                
                %%%% Tort's Modulation Index (Tort et al., 2010)
%                 thetaphase_bin = ceil( tiedrank( thetaphase ) / (ntimepoints / nbins) ); % -- bin the theta phase angles into nbins -- NOTE: tiedrank also exists in eeglab toolbox; when added to path, may cause conflict
%                 gammapow_bin = zeros(1,nbins);
%                 for k=1:nbins
%                     gammapow_bin(k) = squeeze(mean(gammapow(thetaphase_bin==k))); % -- compute mean gamma power in each bin
%                 end
%                 gammapow_bin = gammapow_bin ./ sum(gammapow_bin); % -- normalize
%                 
%                mi(b,c,a) = (log(nbins) + sum(gammapow_bin.*log(gammapow_bin)) ) ./ log(nbins); % -- compute MI
%                 
%                 %%%  -- Phase Locking Value (Cohen, 2008; Colgin et al 2009)
%                 plv(b,c,a) = abs(mean(exp(1i*( thetaphase - angle(hilbert(detrend(gammapow))) ))));
%                 
%                 %%%     Voytek method
%                 vtk(b,c,a) = cfc_est_voytek( thetaphase', gammapowfilt');
%                 
                % glm method with theta filtered (Used this method)
%                  glmf(b,c,a) = cfc_est_glm(thetaphase,gammapowfilt); %% using this at the end
                
                % For masking (GLM)
                [out,stat] = cfc_est_glmstats(thetaphase,gammapowfilt);
                glmf(b,c,a) =out;
                pvals=stat.stats.p;
                glmfP1(b,c,a)= pvals(1);
                glmfP2(b,c,a)= pvals(2);
                glmfP3(b,c,a)= pvals(3);
            end
            disp((b/bb)*100);
            disp('%Still going..');
        end
        
    end
    disp(toc);
    disp('%Finished!!!');
    
    %% plot results
    close all;
 
    % just for cropping of the figure
    theta=mean(theta,1); % makes axis for plots 
    gamma=mean(gamma(:,:,1),1);

    % Using glmf here  
    TGC_plotter(glmf,theta,gamma,[4,8],[30,60],1,'GLMF')

    
    
    %% Makes masking and saves
   
    figure;
    in=glmfP1;
    moo=in;
    pthresh= 0.05/(size(glmfP1,1)*size(glmfP1,2));%(31*41); % Set your threshold here.
    moo((in<pthresh))=0;
    moo((in>pthresh))=1;
%     imshow(moo)
    
    moo2(:,:,1)=[glmf];
    moo2(:,:,2)=[glmfP1];
    moo2(:,:,3)=[glmfP2];
    moo2(:,:,4)=[glmfP3];
    pval=logical(in>pthresh);
    for i =1:4
        mi =   moo2(:,:,i);     
        mi(pval)=NaN;
        moo3(:,:,i)=mi;
     end

    mkdir([outPath]);
    
    cd([outPath])
    fn=['GLMf_stat_repeat_' num2str(h)];
    save(fn,'moo3');
    save('freqs','theta','gamma');
    
    fn2=['GLMf_stat_repeat_no_mask' num2str(h)];
    save(fn2,'glmf');
        
   end

   close all;
%% Loading saved TGCs

cd([outPath])

hh=[];hh=h;
for i = 1:hh
    load(['GLMf_stat_repeat_' num2str(i) '.mat'])
    load('freqs.mat')    
    moo4(:,:,:,i)=moo3;
    clear moo3
    
    load(['GLMf_stat_repeat_no_mask' num2str(i) '.mat'])
    mooX(:,:,:,i) = glmf;
end

hh=[];hh=h;
m=[]; s=1;
for i=1:hh
        m(:,:,s) = moo4(:,:,1,i);
    s=s+1;
end

%%
% Plotting each TGC with masking
TGC_plotter(m,theta,gamma,[4,8],[30,60],1,'GLMf')

% Average of TGCs (or median) --> This would be the stimulation
% frequency
TGC_plotter(mean(m,3),theta,gamma,[4,8],[30,60],1,'GLMf')
% TGC_plotter(median(m,3),theta,gamma,[4,8],[30,60],1,'GLMf')
% TGC_plotter(mean(mooX,3),theta,gamma,[4,8],[30,60],1,'GLMf')


