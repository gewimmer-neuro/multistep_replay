% TrySeq_currother_example.m
% 
% TDLM sequenceness calculation
% 
% input MEG:  nstates x ntimes 'state evidence' matrix
% 
% script by Yunzhe Liu 2018
% adapted by Elliott Wimmer 2019-2021


%% initial set up
clear

procstartup_example

nSubj=is.nSubj;
fprintf(['\n' num2str(nSubj) ' subjects \n'])


analysisdir='/data1/multimeg/analysisr';
addpath(analysisdir);
selectedSubj = is.selectedSubj;
cd(analysisdir);

maxLag =   length(is.lgncy); % 60: evaluate cross-correlation up to 600ms
nstates =                12; % fixed at 12
timepoint =               1; % fixed at 1
L1l =       length(is.ENL1);
Ltt = length(is.whichTimes);
l1param =                 2; % 2
setcon =                  1; % use constant in regression (default yes)
dozscore =                1; % YES to across-lag z-score
ntrialexclude =           4; % fixed exclusion to allow for adjustment to start of task


% n.b.:  need to accomodate the variable # trials in s207-s209
% s207: 100 trials
% s208: 88 trials
% s209: 120 trials


%%% default settings
pretime = 'x2p5konly'; % 2.5 s pre-choice period only
prefile = '2p5konly';
tpinclude = 17:251; % exclude first 160 ms of stimulus processing from replay analysis
descriptiontime = (['\n' 'exclude first 160 ms of stimulus processing from replay analysis' '\n']);
fprintf(descriptiontime)

fprintf(['\n' 'transmat = current + other' '\n'])

input('\nConfirm the above\n')
fprintf('Ok...\n\n');


sf1 = cell(nSubj,1);
sb1 = cell(nSubj,1);
sf2 = cell(nSubj,1);
sb2 = cell(nSubj,1);
ss = cell(nSubj,1);
sc = cell(nSubj,1);


tic

%% loop through subjects
for iSj=1:nSubj % 1:nSubj
    
    % load in transition matrices and trial number
    load([is.rootBehav 's' is.fnSID{iSj} '_trial2state.mat'],'trial2state'); % trial (1:ntrials)
    load([is.rootBehav 's' is.fnSID{iSj} '_transmat.mat']); % transition matrix
    load([is.rootBehav 's' is.fnSID{iSj} '_transotherpathmat.mat']); % transition matrices for other world
    
    
    %% set transition matrix
    transmat{1} = transmat; % current
    transmat{2} = transothermat; % other both
    
    
    %% load in MEG data
    S = load([analysisdir '/classifiers/testcue/Mreds_' prefile '_s' is.fnSID{iSj}]);
    fprintf(descriptiontime)
    
    % Mreds = nshuf(1) * timepoint (1) * trials * L1 params
    Mreds = S.Mreds;
    clear S; % set this subject's preds
    
    % squeeze Mreds
    Mreds=squeeze(Mreds(1,timepoint,:,:)); % look at all trials
    % just transpose, no need for reordering
    Mreds = Mreds';
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nTr = size(Mreds,2);
    L1l = length(is.ENL1);
    Ltt = length(is.whichTimes);
    
    iShuf=1;
    
    sf1{iSj} = NaN(1,nTr, iShuf, maxLag+1);
    sb1{iSj} = NaN(1,nTr, iShuf, maxLag+1);
    sf2{iSj} = NaN(1,nTr, iShuf, maxLag+1);
    sb2{iSj} = NaN(1,nTr, iShuf, maxLag+1);
    ss{iSj} = NaN(1,nTr, iShuf, maxLag+1);
    sc{iSj} = NaN(1,nTr, iShuf, maxLag+1);
    
    %%% exclude first 4 trials to allow for adjustment to beginning of task
    if ntrialexclude>0
        trial2state(1:ntrialexclude) = NaN;
        fprintf(['\n' 'Excluding initial ' num2str(ntrialexclude) ' trials. ' '\n\n']);
    end
    
    %% loop through trials
    for iTr=1:nTr
        
        trialnow=trial2state(iTr);
        
        prSj=squeeze(Mreds(:,iTr));
        
        vec=l1param;
        
        Xtemp=prSj{vec,1};
        
        % select time period
        Xtemp = Xtemp(tpinclude,:);
        
        if any(isnan(Xtemp(:))) || any(isnan(trialnow)) || isnan(choice2state(iTr)) % if X is nan
            disp(['iSj=' is.fnSID{iSj} ' (trial ' num2str(trialnow) ') NaN'])
            continue
        end
        
        disp(['iSj=' is.fnSID{iSj} ' (trial ' num2str(trialnow) ')'])
        
        X = Xtemp;
        clear Xtemp
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % regression
        nbins=maxLag+1;
        warning off
        dm=[toeplitz(X(:,1),[zeros(nbins,1)])];
        dm=dm(:,2:end);
        
        for kk=2:nstates
            temp=toeplitz(X(:,kk),[zeros(nbins,1)]);
            temp=temp(:,2:end);
            dm=[dm temp];
        end
        
        warning on
        
        Y=X;
        
        betas = NaN(nstates*maxLag, nstates);
        
        % coact: control for coactivation of other states
        for ilag=1:maxLag
            zinds = (1:maxLag:nstates*maxLag) + ilag - 1;
            temp = pinv([dm(:,zinds) ones(length(dm(:,zinds)),1)])*Y;
            betas(zinds,:)=temp(1:end-1,:);
        end
        premodel = 'coact_';
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%
        betasnbins64=reshape(betas,[maxLag nstates^2]);
        
        TFcurr = transmat{1};
        TFcurr = TFcurr{iTr};
        
        TFother = transmat{2};
        TFother = TFother{iTr};
        
        % set forward trans matrices for regression, and backwards
        T1f = TFcurr;
        T1b = T1f'; % backwards is transpose of forwards
        
        T2f = TFother;
        T2b = T2f';
        
        bbb=pinv([T1f(:) T1b(:) T2f(:) T2b(:) squash(eye(nstates)) squash(ones(nstates))])*(betasnbins64'); % basic with constant
        preregtype = 'basic';
        colcons = 6;
        
        sf1{iSj}(1,iTr,iShuf,2:end)=zscore(bbb(1,:));
        sb1{iSj}(1,iTr,iShuf,2:end)=zscore(bbb(2,:));
            
        sf2{iSj}(1,iTr,iShuf,2:end)=zscore(bbb(3,:));
        sb2{iSj}(1,iTr,iShuf,2:end)=zscore(bbb(4,:));
        
        ss{iSj}(1,iTr,iShuf,2:end)=zscore(bbb(colcons-1,:)); % autocorrelation
        sc{iSj}(1,iTr,iShuf,2:end)=zscore(bbb(colcons,:)); % constant
        
    end
    
    
    clear X prSj betas betasnbins64 Mreds trans*mat transmat* bbb dm T1* T2* trial2state
end

% report
fprintf(['\n' 'transmat = current and other' '\n'])

toc


%% concatenate
if iscell(sf1) && iscell(sb1)
    sf1cell = sf1;
    sb1cell = sb1;
    sf2cell = sf2;
    sb2cell = sb2;
    
    sscell = ss;
    sccell = sc;
    
    sf1 = cell2mat(sf1);
    sb1 = cell2mat(sb1);
    sf2 = cell2mat(sf2);
    sb2 = cell2mat(sb2);
    ss = cell2mat(ss);
    sc = cell2mat(sc);
    
    sf1 = squeeze(sf1);
    sb1 = squeeze(sb1);
    sf2 = squeeze(sf2);
    sb2 = squeeze(sb2);
    ss = squeeze(ss);
    sc = squeeze(sc);
end


nSubj = size(sf1,1);
maxLag = 60;
cTime = 0:10:maxLag*10;

% get means
sf1m = [];
sb1m = [];
sf2m = [];
sb2m = [];
ssm = [];
scm = [];
for iSj = 1:size(sf1,1)
    sf1m(iSj,:,:) = squeeze(nanmean(sf1(iSj,:,:),2));
    sb1m(iSj,:,:) = squeeze(nanmean(sb1(iSj,:,:),2));
    sf2m(iSj,:,:) = squeeze(nanmean(sf2(iSj,:,:),2));
    sb2m(iSj,:,:) = squeeze(nanmean(sb2(iSj,:,:),2));
    ssm(iSj,:,:) = squeeze(nanmean(ss(iSj,:,:),2));
    scm(iSj,:,:) = squeeze(nanmean(sc(iSj,:,:),2));
end



%%% save
save(['seqcue_currothertx_' preregtype num2str(ntrialexclude) 'x_' pretime '_' premodel 'n' num2str(size(sf1,1))])

