%% PERMUTATION TESTING - True Left hem vs Pseudo Left hem
% Pseudo left hem is the right hem flipped onto left hem
% Script to test lateralization of HCP variants across handedness groups
% input is variant maps 
% output is a matrix with r values for n permutations, with the option of
% creating a plot, too.
% Written by Diana P -- January 2021

clear all

%--------------------------------------------------------------------------
%% PATHS
%--------------------------------------------------------------------------

%root_dir = '/projects/p31161/lateralization_code/'; %location of scripts and needed files
root_dir = '/Users/dianaperez/Desktop/lateralization_code/'
addpath(genpath(root_dir))
%addpath(genpath('/projects/b1081/Scripts')); %add cifti-matlab-master folder to path
varmap_loc = [root_dir '/Variant_Maps/new_split_vars/reassigned/']; %location of variant maps
varmap_str = '_uniqueIDs_afterReassign.dtseries.nii'; 
template_loc = [varmap_loc '100206_border1ectopic2.dtseries.nii'];
output_dir = [root_dir 'testing_output/spatial_distribution/'];
if ~exist(output_dir, 'file')
    mkdir(output_dir)
end

%--------------------------------------------------------------------------
%% VARIABLES
%--------------------------------------------------------------------------
numperms = 1000;% number of permutations
plot_results = 1;

%initialize variables
spCorrs=zeros(numperms,1);% initialize mat of values i'm permuting
perm_diffs=zeros(numperms,1);% initialize mat of values i'm permuting

in_str = 'HCP752';
out_str = strcat(in_str, '_LHandvsRHand');
group_str = {'HCP752_LH', 'HCP752_RH'};

LH = []; MH = []; RH = []; % initialize mat of groups of subs; LH = left handers, RH = right handers, MH = 'middle' handers
load([root_dir '/needed_files/goodSubs752.mat'])
all_subs = goodSubs752;
for s = 1:length(all_subs)
    if all_subs(s,2) < -28
        LH = [LH; all_subs(s,:)];
    elseif all_subs(s,2) < 48
        MH = [MH; all_subs(s,:)];
    elseif all_subs(s,2) <= 100
        RH = [RH; all_subs(s,:)];
    end
end

subs = {LH; RH};
numSubs = length(LH) + length(RH);
numLH = length(LH);      

%--------------------------------------------------------------------------
%% BEGIN ANALYSIS
%--------------------------------------------------------------------------
left_cifti_loc = [root_dir strcat(in_str, '_true_left_handers.mat')];
right_cifti_loc = [root_dir strcat(in_str, '_true_right_handers.mat')];
if exist(left_cifti_loc, 'file')
    load(left_cifti_loc)   
else        
    for n = 1:length(LH)
%         true_left_handers = [];
        cifti = ft_read_cifti_mod([varmap_loc num2str(LH(n,1)) varmap_str]);
        true_left_handers(:,n) = cifti.data;
    end
    save(left_cifti_loc, 'true_left_handers')
end
if exist(right_cifti_loc, 'file')
    load(right_cifti_loc)
else
    for n = 1:length(RH)
        cifti = ft_read_cifti_mod([varmap_loc num2str(RH(n,1)) varmap_str]);
        true_right_handers(:,n) = cifti.data;
    end
    save(right_cifti_loc, 'true_right_handers')
end

% matrix with 1's and 0's, will determine if sub's maps will be flipped
% -- will be randomized later
flip_switch = zeros(numSubs,1);

%--------------------------------------------------------------------------
%% GET CORRELATION BETWEEN TRUE MAPS
%-------------------------------------------------------------------------- 
tic;    

[true_overlap_left] = makemaps(true_left_handers);
[true_overlap_right] = makemaps(true_right_handers);

% get correlation
true_spCorr = corr(true_overlap_left, true_overlap_right);
trueDiff = true_overlap_left - true_overlap_right;
save([output_dir out_str '_true_spatialCorr.mat'],'true_spCorr')               
save([output_dir out_str '_true_diffMap.mat'],'trueDiff')

all_data = [true_left_handers true_right_handers];
%--------------------------------------------------------------------------
%% RUN PERMUTATIONS - MAKE RANDOMIZED MAPS
%--------------------------------------------------------------------------
%flip_switch(1:length(numLH)) = 1; % add 1's
        
for p = 1:numperms % will go through this loop for each permutation    
    disp(['Permutation #' num2str(p)])
    rng('shuffle'); % need seed for randomizer

    % randomize 1's and 0's
    ind = randperm(numSubs)';
    %rand_data = all_data(:,ind);
    rand_LH = all_data(:,ind(1:numLH));
    rand_RH = all_data(:,ind(numLH+1:end));
    
    [overlap_pseudo_left] = makemaps(rand_LH);
    [overlap_pseudo_right] = makemaps(rand_RH);

    % get correlation
    spCorrs(p) = corr(overlap_pseudo_left, overlap_pseudo_right);
    perm_diff = overlap_pseudo_left - overlap_pseudo_right;
    diffmats(:,p)= perm_diff;
end 
        
save([output_dir out_str '_spCorr_' num2str(numperms) '_permutations.mat'],'spCorrs')
save([output_dir out_str '_diffMap_' num2str(numperms) '_permutations.mat'],'diffmats')

%--------------------------------------------------------------------------
%% PLOT RESULTS
%--------------------------------------------------------------------------

if plot_results == 1
    figure;
    ind = zeros(length(perm_diffs),1);
    ind(1,1) = 1;
    handles = plotSpread(spCorrs, 'categoryMarkers', {'o', '.'}, 'categoryLabels', {'Permuted Correlation', 'True Correlation'}, 'categoryColors', {'k', 'r'}, 'categoryIdx', [ind])
    scatter(1,spCorrs(1,:), 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red', 'SizeData', 50)
    ax = gca;
    ax.FontSize = 12;
    axis([0.5, 1.5, min(spCorrs)-std(spCorrs), max(spCorrs)+std(spCorrs)])
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.3, 0.3, 0.3, 0.7]); %first and second control position on screen, third controls width, and fourth controls height
    print(gcf,[output_dir out_str '_permuted_spatialCorrs_' num2str(numperms) '_permutations_figure.jpg'],'-dpng','-r300');
    close gcf
end

    

function [overlap] = makemaps(data)
 overlap = sum(data,2)/length(data);
end


