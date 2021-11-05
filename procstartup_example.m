% procstartup.m
% 
% 

%%% specify subjects
specifysubjects_example


%%% general stuff
basedir = '/data1/multimeg/';

%%% set is values
is.fs = filesep;
is.rootBehav = [basedir 'regfiles/'];
is.rootMEG = [basedir 'optdata/'];
is.networkPath = basedir; % where to pull the unpreprocessed MEG data from
is.OPTPath = [basedir 'optdata/'];
is.writePath = '/data/multimeg/analysisl';
is.smoothFact = 6;
is.msPerSample = 1000 * is.smoothFact / 1200; % 1200
is.lgncy = 1:60;  % the latencies to consider cross-correlation at, in units of samples
is.Ltimes = is.lgncy*is.msPerSample/1000;
is.nShuf = 20; %20;%20;%29;
is.ENL1 = 0.001:0.001:0.01; % 0.001:0.001:0.01   % lasso regularization L1=lambda*alpha
is.selectedSubj = 1:is.nSubj;
is.highpass = 0.5; % decide which hiph pass spec is used for the further analysis
is.bandpass = [];
is.whichSubj = is.selectedSubj;