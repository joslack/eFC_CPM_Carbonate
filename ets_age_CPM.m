%% CPM predicting age using code from:
%Shen, X., Finn, E. S., Scheinost, D., Rosenberg, M. D., Chun, M. M.,
%          Papademetris, X., & Constable, R. T. (2017). Using connectome-based
%          predictive modeling to predict individual behavior from brain connectivity.
%          Nature Protocols, 12(3), 506.

addpath(genpath('CPM'));
clear all
close all
clc




n = 214;
dropframes = 50;

% load in activity data
load('nki_peaks');
subs = unique(subjid);


% load in FC
ts = zeros(896,n,length(subs));
[ntime,nnodes] = size(squeeze(ts(:,:,1)));
nedges = nnodes*(nnodes-1)/2;
[u,v] = find(triu(ones(nnodes),1));
idx = (v - 1)*nnodes + u;
filenames = dir(fullfile('../rest/*645*.mat'));
for i = 1:length(subs)        % loop over all subjects
    for j = 1:length(filenames)
        if contains(filenames(j).name,subs(i))
            filename = filenames(j).name
            break;
        end
    end
    if ~isempty(filename)       % if file exists
        % take first file in list
        load(fullfile('../rest',filename));                % load file
       
        
        z = ...                                         % z-score time series
            zscore(parcel_time{1}(dropframes + 1:end - dropframes,:));
        %             ets = fcn_edgets(z);                            % calculate ets
        %             r = sum(ets.^2,2).^0.5;                         % calculate rss
        %             R(:,i) = r;                                     % store rss
        ets = z(:,u).*z(:,v);
        efc = triu(corr(ets));
        At = efc.';
        m  = (1:size(At,1)).' >= (1:size(At,2));
        vefc(i,:)  = At(m);
    end % calculate fc
    
    
    %         end
end


% loop over clusters

%% generate cpm network

% keep = ~isnan(sum(centroids,1));

% this all comes from figs 3 and 4 in shen et al.

% i renamed variables so that they are identical to what's in the shen
% paper.

all_behav = a;
all_mats = vefc;

% all_behav = all_behav(keep);
% all_mats = all_mats(:,:,keep);

no_sub = size(all_mats,3);
no_node = size(all_mats,1);

behav_pred_pos = zeros(no_sub,1);
behav_pred_neg = zeros(no_sub,1);

thresh = 0.01; % 0.001, 0.05, or adjust for multiple comparisons

for leftout = 1:no_sub
    
    train_mats = all_mats;
    train_mats(:,:,leftout) = [];
    train_vcts = reshape(train_mats,[],size(train_mats,3));
    
    train_behav = all_behav;
    train_behav(leftout) = [];
    
    [r_mat,p_mat] = corr(train_vcts',train_behav); % corr(x,y,'type','spearman')
    r_mat = reshape(r_mat,no_node,no_node);
    p_mat = reshape(p_mat,no_node,no_node);
    
    pos_mask = zeros(no_node,no_node);
    neg_mask = pos_mask;
    
    pos_edges = find(r_mat > 0 & p_mat < thresh);
    neg_edges = find(r_mat < 0 & p_mat < thresh);
    
    pos_mask(pos_edges) = 1;
    neg_mask(neg_edges) = 1;
    
    train_sumpos = zeros(no_sub-1,1);
    train_sumneg = train_sumpos;
    
    for ss = 1:size(train_sumpos)
        train_sumpos(ss) = sum(sum(train_mats(:,:,ss).*pos_mask))/2;
        train_sumneg(ss) = sum(sum(train_mats(:,:,ss).*neg_mask))/2;
    end
    
    fit_pos = polyfit(train_sumpos,train_behav,1); % robustfit
    fit_neg = polyfit(train_sumneg,train_behav,1);
    
    test_mat = all_mats(:,:,leftout);
    test_sumpos = sum(sum(test_mat.*pos_mask))/2;
    test_sumneg = sum(sum(test_mat.*neg_mask))/2;
    
    behav_pred_pos(leftout) = fit_pos(1)*test_sumpos + fit_pos(2);
    behav_pred_neg(leftout) = fit_neg(1)*test_sumneg + fit_neg(2);
    
    disp(leftout)
    
    plot(all_behav(1:leftout),behav_pred_pos(1:leftout),'r.'); lsline;
    drawnow;
    
    
end

[R_pos, P_pos] = corr(behav_pred_pos,all_behav);
[R_neg, P_neg] = corr(behav_pred_neg,all_behav);

% to test: modify script slightly to do CPM with traditional functional
% connectivity

figure(1); plot(behav_pred_pos,all_behav,'r.'); lsline;
figure(2); plot(behav_pred_neg,all_behav,'b.'); lsline;



