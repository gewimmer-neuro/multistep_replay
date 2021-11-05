% specifysubjects_example.m
% 
% 
% 
% MEGruns: 'loc' = localizer; 'struct' = structure learning; 'rew' = reward learning; 'post' = post-localizer


is.nSubj = 0;


%% 


is.fnDate{is.nSubj} = 's200'; is.fnSID{is.nSubj} = '200';
is.fnBehav{is.nSubj} = {'001'};
is.fnMEG{is.nSubj} = 'MG09999_Elliott_20200101';
is.MEGruns{is.nSubj} = {'loc' 'loc' 'loc' 'loc' 'struct' 'struct' 'rew' 'rew' 'rew' 'rew' 'rew' 'rew' 'post'};
is.taskversion{is.nSubj}=1;
is.nSubj = is.nSubj + 1;



