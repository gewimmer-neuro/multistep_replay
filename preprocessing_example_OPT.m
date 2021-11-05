% preprocessing_example_OPT.m
%
% by Yunzhe Liu 2018
% adapted by Elliott Wimmer 2019-2021
% 
% provide input "is" struct with base paths, subjects, run labels per subject, and other variables


procstartup_example


%% get reward onset events
codeset = 1; % cue = 1; reward onset = 3;

epochtimerange = [-0.5 2.5]; % pre-choice period and 500 ms baseline

%% loop through subjects
for subj = 1:nSubj
    
    tic
    
    fprintf(['\n' '      preprocessing reward period s' is.fnSID{subj} '\n\n']);
    
    spL = [is.OPTPath strrep(is.fnDate{subj}, '-', '_') '/' is.fnMEG{subj} '_'];  % local base path
    
    %% set path for raw data
    spN = [is.networkPath strrep(is.fnDate{subj},'-','') '/' is.fnMEG{subj} '_']; % network base path
    
    spS = [is.networkPath is.fnDate{subj}]; % subject base path
    subjNum = is.fnDate{subj}; subjNum = str2double(subjNum(2:4));
    
    %% select runs and loop through runs
    rewInd = find(ismember(is.MEGruns{subj},'rew'));
    nsession=length(rewInd); % usually 6 reward blocks
    
    for run = 1:nsession
        
        MEGrunlabel = is.MEGruns{subj}{run};
        MEGrunlabel = MEGrunlabel(1:3);
        
        inSt = num2str(run, '%02d');
        localPath = [spL inSt '.ds'];
        mkdir(localPath);
        ds_folder = [spN inSt '.ds'];
        
        %% Import - using OPT
        S = struct;
        S.fif_file = ds_folder;
        S.spm_file = fullfile(localPath,'spmeeg.mat');
        
        S.other_channels = {'UADC001','UADC002','UADC003','UADC004'}; % added photodiode channel to 3 eyelink channelss
        
        D = osl_convert_script(S);
        
        fprintf(['\n' '      preprocessing reward period '  is.fnDate{subj} ' run ' num2str(run) '\n'])
        
        % add up other channels
        D = D.chantype(find(strcmp(D.chanlabels,'UADC001')),'EOG1');
        D = D.chantype(find(strcmp(D.chanlabels,'UADC002')),'EOG2');
        D = D.chantype(find(strcmp(D.chanlabels,'UADC003')),'EOG3');
        D = D.chantype(find(strcmp(D.chanlabels,'UADC004')),'PHDI'); % add photodiode channel
        
        D.save()
        % Crop unused data - Note: spm_eeg_crop uses a 1e-3*max notation
        S = struct;
        S.D = fullfile(localPath,'spmeeg.mat');
        S.prefix='';
        event = ft_read_event(ds_folder);
        
        fprintf(['\n' '      preprocessing reward period '  is.fnDate{subj} ' run ' num2str(run) '\n'])
        
        % crop the data using slightly modified script
        sample = [event(find(strcmp('frontpanel trigger', {event.type}))).sample]';
        S.timewin = [0,round(sample(end)/1200*1000)+40000]; % for 1200/sec sampling rate
        D=spm_eeg_crop_YL(S);
        D.save()
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %% OPT Preprocessing Pipeline!
        
        %% Phase 1 - Filter
        opt=[];
        opt.maxfilter.do=0;
        
        spm_files{1}=fullfile(localPath,'spmeeg.mat');
        structural_files{1}=[]; % leave empty if no .nii structural file available
        
        opt.spm_files=spm_files;
        opt.datatype='ctf';
        
        % HIGHPASS
        if is.highpass==0
            opt.highpass.do=0;
            opt.dirname=fullfile(localPath,['nohighpassrew_',num2str(is.highpass)]);
        else
            opt.highpass.cutoff = is.highpass; % create different folder for different frenquency band (0.5; 1; 20)
            opt.dirname=fullfile(localPath,['highpassrew_',num2str(opt.highpass.cutoff)]);
            opt.highpass.do=1;
        end
        
        % Notch filter settings
        opt.mains.do=1;
        
        % DOWNSAMPLING
        opt.downsample.do=0;
        
        % IDENTIFYING BAD SEGMENTS
        opt.bad_segments.do=0;
        
        % Set to 0 for now
        opt.africa.do=0;
        opt.epoch.do=0;
        opt.outliers.do=0;
        opt.coreg.do=0;
        
        %%%%%%%%%%%%%%%%%%%%%
        opt = osl_run_opt(opt);
        
        
        %% Phase 2 - downsample + bad segment
        opt2=[];
        opt2.maxfilter.do=0;
        
        opt2.spm_files = opt.results.spm_files;
        opt2.datatype='ctf';
        
        % optional inputs
        opt2.dirname=opt.dirname; % directory opt settings and results will be stored
        opt2.convert.spm_files_basenames = opt.results.spm_files_basenames;
        
        % DOWNSAMPLING
        opt2.downsample.do=1;
        opt2.downsample.freq=100;
        
        % IDENTIFYING BAD SEGMENTS
        opt2.bad_segments.do=1;
        
        % Set to 0 for now
        opt2.africa.do=0;
        opt2.epoch.do=0;
        opt2.outliers.do=0;
        opt2.coreg.do=0;
        opt2.highpass.do=0;
        opt2.mains.do=0;
        %%%%%%%%%%%%%%%%%%%%%
        opt2 = osl_run_opt(opt2);
        
        
        %% Phase 3 - ICA
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        opt3 = [];
        % required inputs
        opt3.spm_files = opt2.results.spm_files;
        opt3.datatype='ctf';
        
        % optional inputs
        opt3.dirname=opt2.dirname; % directory opt settings and results will be stored
        
        % africa settings
        opt3.africa.do=1;
        opt3.africa.ident.artefact_chans = {'EOG1','EOG2','EOG3'};
        opt3.africa.precompute_topos = false;
        opt3.africa.ident.mains_kurt_thresh = 0.5;
        %opt3.africa.ident.func = @identify_artefactual_components_auto;
        %opt3.africa.ident.max_num_artefact_comps = 50;
        opt3.africa.ident.do_kurt = true;
        opt3.africa.ident.do_cardiac = true;
        opt3.africa.ident.do_plots = true;
        opt3.africa.ident.do_mains = false;
        
        % turn the remaining options off
        opt3.maxfilter.do=0;
        opt3.convert.spm_files_basenames = opt2.results.spm_files_basenames;
        opt3.downsample.do=0;
        opt3.highpass.do=0;
        opt3.bad_segments.do=0;
        opt3.epoch.do=0;
        opt3.outliers.do=0;
        opt3.coreg.do=0;
        
        opt3 = osl_run_opt(opt3);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %% Phase 4 - Epoch + coreg
        opt4 = [];
        % required inputs
        opt4.spm_files = opt3.results.spm_files;
        opt4.datatype='ctf';
        
        % optional inputs
        opt4.dirname=opt3.dirname; % directory opt settings and results will be stored
        
        % Epoching settings
        opt4.epoch.do=1;
        
        opt4.epoch.time_range = epochtimerange;
        
        allevent=[event.value];
        eventtally = [];
        
        m = 1;
        
        % find event==codeset (here = 1)
        eventnum = find(allevent==codeset);
        eventnum = eventnum(1);
        
        % only keep event type codeset
        for ii=eventnum:eventnum
            
            if ~(ismember(allevent(ii),eventtally))
                opt4.epoch.trialdef(m).eventvalue = allevent(ii);
                opt4.epoch.trialdef(m).eventtype = 'frontpanel trigger';
                
                opt4.epoch.trialdef(m).conditionlabel = ['All' num2str(codeset)];
                
                m = m+1;
                eventtally = [eventtally allevent(ii)];
            end
        end
        
        opt4.bad_segments.do=0;
        opt4.outliers.do=1;
        
        
        % coreg for subsequent source analysis - already done
        opt4.coreg.do=1;
        opt4.coreg.use_rhino=0; % unless have head point data
        opt4.coreg.useheadshape=0;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % turn the remaining options off
        opt4.maxfilter.do=0;
        opt4.convert.spm_files_basenames = opt3.results.spm_files_basenames;
        opt4.downsample.do=0;
        opt4.africa.todo.ica=0;
        opt4.africa.todo.ident=0;
        opt4.africa.todo.remove=0;
        
        cd([is.OPTPath is.fnDate{subj}])
        
        opt4=osl_run_opt(opt4);
        
        %% display results and cleanup
        opt4 = osl_load_opt(opt4.dirname);
        
        
        close all
        fclose('all');
    end
    
    toc
end



