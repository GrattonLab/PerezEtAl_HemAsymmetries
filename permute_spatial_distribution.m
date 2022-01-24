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


addpath(genpath('/Users/dianaperez/Box/Dependencies/cifti-matlab-master'));
root_dir = '/Users/dianaperez/Desktop/lateralization_code/';
addpath(genpath(root_dir))
varmap_loc = '/Volumes/RESEARCH_HD/HCP_Variants/new_split_vars/reassigned/'; %location of variant maps 
varmap_str = '_uniqueIDs_afterReassign.dtseries.nii'; 
template_loc = '/Volumes/RESEARCH_HD/HCP_Variants/new_split_vars/reassigned/100206_border1ectopic2.dtseries.nii';
output_dir = [root_dir '/testing_output/spatial_distribution/'];
if ~exist(output_dir, 'file')
    mkdir(output_dir)
end

%--------------------------------------------------------------------------
%% VARIABLES
%--------------------------------------------------------------------------
numperms = 1000;% number of permutations
% sample that is being analyzed
MSC = 0;
HCP = 1;
handedness = 1;
plot = 1;

%initialize variables
spCorrs=zeros(numperms,1);% initialize mat of values i'm permuting
perm_diffs=zeros(numperms,1);% initialize mat of values i'm permuting

if MSC
    in_str = {'MSC'};
else
    if handedness
        in_str = {'HCP752_LH', 'HCP752_RH'};
    else
        in_str = {'HCP384'};
    end
end

if HCP
    if handedness
        LH = []; MH = []; RH = []; % initialize mat of groups of subs; LH = left handers, RH = right handers, MH = 'middle' handers
        load([root_dir '/needed_files/goodSubs752.mat'])
        all_subs = goodSubs752;
        for n = 1:length(all_subs)
            if all_subs(n,2) < -28
                LH = [LH; all_subs(n,:)];
            elseif all_subs(n,2) < 48
                MH = [MH; all_subs(n,:)];
            elseif all_subs(n,2) <= 100
                RH = [RH; all_subs(n,:)];
            end
        end
        subs = {LH; RH};
        group = {'LH', 'RH'};
    else
        load([root_dir '/needed_files/goodSubs384.mat'])
        all_subs = goodSubs384;
        subs = {all_subs};
    end
elseif MSC
    all_subs = {'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07', 'MSC09', 'MSC10'};
    subs = {all_subs};
end

%--------------------------------------------------------------------------
%% BEGIN ANALYSIS
%--------------------------------------------------------------------------
template = ft_read_cifti_mod(template_loc);
brainstructure = single(template.brainstructure);

for g = 1:numel(subs)
    out_str = append(in_str{g}, '_LhemvsRhem');
    numSubs = length(subs{g});
    left_cifti_loc = [root_dir append(in_str{g}, '_true_left_hems.mat')];
    right_cifti_loc = [root_dir append(in_str{g}, '_true_right_hems.mat')];
    if exist(left_cifti_loc, 'file')
        load(left_cifti_loc)
        load(right_cifti_loc)
    else        
        true_left_hems = [];
        true_right_hems = [];
        for n = 1:length(subs{g})
            if HCP
                cifti = ft_read_cifti_mod([varmap_loc num2str(subs{g,1}(n)) varmap_str]);
            elseif MSC
                cifti = ft_read_cifti_mod([varmap_loc subs{g,1}{n} varmap_str]);
            end

            new_cifti = [];
            count = 1;

            for i = 1:length(cifti.brainstructure)
                if cifti.brainstructure(i) > 0 
                    if cifti.data(count) >0
                        new_cifti(i,1) = 1;
                    else new_cifti(i,1) = 0;
                    end
                    count = count + 1;
                else
                    new_cifti(i,1) = 0;
                end
            end    

            true_left_hems(:,n) = single(new_cifti(1:(length(new_cifti)/2)));
            true_right_hems(:,n) = single(new_cifti((length(new_cifti)/2)+1:end));
        end

        save(left_cifti_loc, 'true_left_hems')
        save(right_cifti_loc, 'true_right_hems')
    end


    % matrix with 1's and 0's, will determine if sub's maps will be flipped
    % -- will be randomized later
    flip_switch = zeros(numSubs,1);

    %--------------------------------------------------------------------------
    %% GET CORRELATION BETWEEN TRUE MAPS
    %-------------------------------------------------------------------------- 
    [true_overlap_left, true_overlap_flip_left] = makemaps(true_left_hems, true_right_hems, flip_switch, brainstructure);

    % get correlation
    true_spCorr = corr(true_overlap_left, true_overlap_flip_left);
    trueDiff = true_overlap_left - true_overlap_flip_left;
    save([output_dir out_str '_true_spatialCorr.mat'],'true_spCorr')               
    save([output_dir out_str '_true_diffMap.mat'],'trueDiff')
    % save overlap map for this group
    %% PROBLEM - it's two left hems so missing some vertices -- might be worth just including the separate function for making overlap maps
%     overlap = [true_overlap_left; true_overlap_flip_left];
%     template.data = overlap;
%     ft_write_cifti_mod([output_dir out_str '_overlap_map.dtseries.nii'], template);
    

    %--------------------------------------------------------------------------
    %% RUN PERMUTATIONS - MAKE RANDOMIZED MAPS
    %--------------------------------------------------------------------------
    flip_switch(1:(round(length(flip_switch)/2))) = 1; % add 1's

    for p = 1:numperms % will go through this loop for each permutation
    tic;    
        disp(['Permutation #' num2str(p)])
        rng('shuffle'); % need seed for randomizer

        % randomize 1's and 0's
        ind = randperm(length(flip_switch))';
        for n = 1:length(flip_switch)
            rand_flip_switch(n,1) = flip_switch(ind(n));
        end

        [overlap_left, overlap_flip_left] = makemaps(true_left_hems, true_right_hems, rand_flip_switch, brainstructure);

        % get correlation
        spCorrs(p) = corr(overlap_left, overlap_flip_left);
        permDiff = overlap_left - overlap_flip_left;
        diffmats(:,p)= permDiff;

    end 
        
    save([output_dir out_str '_permuted_spatialCorrs_' num2str(numperms) '_permutations.mat'],'spCorrs')
    save([output_dir out_str '_permuted_diffMap_' num2str(numperms) '_permutations.mat'],'diffmats')
        
    %--------------------------------------------------------------------------
    %% PLOT RESULTS
    %--------------------------------------------------------------------------
    if plot
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
end

function [overlap_left overlap_flip_left] = makemaps(left_hem, right_hem, rand_flip_switch, brainstructure)
 % initialize variables to include left hem and flipped left hem
left = [];
flip_left = [];
brainstructure = brainstructure(1:32492);
for x = 1:length(rand_flip_switch)
    if rand_flip_switch(x) == 0
        left(:,x) = left_hem(brainstructure==1,x);
        flip_left(:,x) = right_hem(brainstructure==1,x);
    elseif rand_flip_switch(x) == 1
        left(:,x) = right_hem(brainstructure==1,x);
        flip_left(:,x) = left_hem(brainstructure==1,x);
    end
end

% divide by number of subs to obtain proportion
overlap_left = sum(left,2)/length(rand_flip_switch);        
overlap_flip_left = sum(flip_left,2)/length(rand_flip_switch);
end

function [LH_files, RH_files, middle_files] = makingdatafiles(dataLoc, fileName, LH, middle, RH)
%% function to make list of paths to files
    LH_files = []; middle_files = []; RH_files = [];
    for x = 1:length(LH)
        file = [dataLoc num2str(LH(x)) fileName];
        LH_files{x,1} = file;
    end
    for y = 1:length(middle)
        file = [dataLoc num2str(middle(y)) fileName];
        middle_files{y,1} = file;
    end
    for z = 1:length(RH)
        file = [dataLoc num2str(RH(z)) fileName];
        RH_files{z,1} = file;
    end
end
