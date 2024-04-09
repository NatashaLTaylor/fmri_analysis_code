%% Timeseries Analysis
%% Functional Connectivity
% load timeseries
sub = sprintf('%s%d','sub-s',test_grp(ID,1));
filename = sprintf('%s%s',sub,'_ses-2_timeseries.mat'); %'sub-s558_ses-2_timeseries.mat'
load(filename)
%functional connectivity
ts_corr = corr(ts); %ROIs x ROIs output matrix

nPairs = 502; %number of unique pairs from template (ie., ROIs)
nTime = size(ts,1); %number of time points

%% Design Matrix from event.tsv files
%Design Matrix, hrf of the events of interest
%example of of loading tsv from different dataset
    cd ('C:\Users\natas\Documents\PhD\fMRI_stop_task\fmri_stop_data\');
    sub = sprintf('%s%d','sub-s',test_grp(ID,1));
    cd(sub) 
%load in file, loads each column as separate variable
    events = dir('*_run-1_events.tsv');
    tdfread(events.name);
    %binarise the 'trial_type' column
    bin_trial_type = zeros(length(trial_type),1);
    for i=1:length(trial_type)
        if trial_type(i,:) == 'stop_success';
        bin_trial_type(i,1)=1;
        elseif trial_type(i,:) == 'stop_failure';
        bin_trial_type(i,1)=2;
        else  trial_type(i,:) == 'go          ';
        bin_trial_type(i,1)=0;
        end
    end
loc_sf = find(bin_trial_type==2);
trial_prior_after_sf= onset(loc_sf);
for n=1:length(loc_sf)
    if n==1
        a = loc_sf(n,1);
        b(n,1) = a-1;
        b(n,2) = a+1;
        %c = vertcat(b,loc_sf);
    else
        a = loc_sf(n,1);
        b(n,1) = a-1;
        b(n,2) = a+1;
        %d = [c(1:n,:); b; c(n+1:end,:)]
    end
end

a= [0,1,2];
a= a.';
%1st column = go trials, 2nd column = ss, 3rd column = sf
for i=1:3
    trial_type_onsets(:,i) = onset.*(bin_trial_type==(a(i,1)));
    trial_type_duration(:,i) = block_duration.*(bin_trial_type==(a(i,1)));
end


%% dsmtx by trial type
dsmtx_go = dsmtx_hrf_epoch(trial_type_onsets(:,1),trial_type_duration(:,1),0.68);
if length(dsmtx_go)~=size(ts,1) %if length of dsmtx not same as ts
    dsmtx_go = vertcat(dsmtx_go,zeros(size(ts,1)-length(dsmtx_go),1));
end
dsmtx_sf = dsmtx_hrf_epoch(trial_type_onsets(:,3),trial_type_duration(:,3),0.68);
if length(dsmtx_sf)~=size(ts,1)
    dsmtx_sf = vertcat(dsmtx_sf,zeros(size(ts,1)-length(dsmtx_sf),1));
end
dsmtx_ss = dsmtx_hrf_epoch(trial_type_onsets(:,2),trial_type_duration(:,2),0.68);
if length(dsmtx_ss)~=size(ts,1)
    dsmtx_ss = vertcat(dsmtx_ss,zeros(size(ts,1)-length(dsmtx_ss),1));
end


dsmtx_all = zeros(size(ts,1),3);
for i=1:3
    if i==1
        dsmtx_all(:,i) = dsmtx_go(1:size(ts,1),:);
    elseif i==2
        dsmtx_all(:,i) = dsmtx_ss(1:size(ts,1),:);
    else 
        dsmtx_all(:,i) = dsmtx_sf(1:size(ts,1),:);
    end
end

%% Gitting GLM to Time-series
%fit glm to time-series
for qq = 1:nPairs
    beta_ts(:,qq) = glmfit(All_dsmtx(ID).dsmtx,ts(:,qq)); 
    % first row of beta is y-intercept
    %second row of beta corresponds to dsmtx column 1 (eg., go) etc..
%sprintf('%d',qq)
end

%% Dynamic Functional Connectivity
% coupling comes from - http://github.com/macshine/coupling/
%previous loading timeseries from above
window = 10; %variably change window from 5 - 20
mtd = coupling(ts,window,0,0); %time x nodes
% flatten out approach for MTD
for nn=1:nROI
template = find(tril(ones(nROI))-eye(nROI)); %try taking lower triangle instead tril
end
nTime = size(mtd,3);
for tt = 1:nTime
temp = mtd(:,:,tt);
mtd_flat(:,tt) = temp(template); %flattens mtd across timepoints
end

% MTD fitted with glm to events
nPairs = size(template,1);
for qq = 1:nPairs
beta_mtd(:,qq) = glmfit(dsmtx,mtd_flat(qq,:)'); 
end


% running through each ROI of the MTD
for jj = 1:nROI
    for kk = 1:nROI
        beta_mtd_full(jj,kk,:) = glmfit(dsmtx,mtd(jj,kk,:));
    end
end



%% Example loop through subjects
for ID=1:49
sub = sprintf('%s%d','sub-s',test_grp(ID,1));
% standard functional connectivity
ts_corr = corr(All_ts_test(ID).timeseries'); %calculate the Pearson's correlation between regions
nPairs = 502; %number of unique pairs from template
nTime = size(All_ts_test(ID).timeseries,2); %number of time points

if size(All_ts_test(ID).timeseries,2)==size(All_dsmtx(ID).dsmtx,1);
    sprintf('ts and dsmtx match')
else 
    sprintf('ts and dsmtx no match')
end

%fit glm to time-series
for qq = 1:nPairs
beta_ts(:,qq) = glmfit(All_dsmtx(ID).dsmtx,All_ts_test(ID).timeseries(qq,:)); 
%sprintf('%d',qq)
end

All_glm_fc(ID) = struct('sub',{sub},'fc',{ts_corr},'glm_timeseries',{beta_ts});
sprintf('%d',ID)
clear ts_corr beta_ts nTime sub
end



%% Graph Theory Analysis Example
%requires installation of brain connectivity toolbox https://sites.google.com/site/bctnet/
%Output:     
%              ci       time-resolved community assignment
%              q        time-resolved modularity
%              p        time-resolved participation coefficient
%              z        time-resolved module-degree z-score
%              hc       cartographic profile
%              f        flexibility

filename=dir('*mtd.mat');
for i=1:length(filename)
    %get 'sub-...' from the name
    subject_file= filename(i).name; %filename of time-series data
    split = strsplit(subject_file,'_');
    subnum = split(1); %sub-.. section
    subnum =cell2mat(subnum);
    ses = split(2); %session ses-..
    ses = cell2mat(ses);
    
    load([subject_file]); %load in mtd

%data size is 502 x 502 x 509
    gamma = 1;
    beta = 0.75;

    [ci,q,part,modz,hc,f] = integration_plus5(mtd,gamma,beta);

    %saving files
    save([subnum '_' ses '_part.mat'],"part");
    save([subnum '_' ses '_flex.mat'],"f");
    save([subnum '_' ses '_modz.mat'],"modz");
    save([subnum '_' ses '_modularity.mat'],"q");
    save([subnum '_' ses '_community.mat'],"ci");
    save([subnum '_' ses '_cartographic.mat'],"hc");

end


