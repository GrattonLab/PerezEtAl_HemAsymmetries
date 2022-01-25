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

root_dir = '/projects/p31161/lateralization_code/'; %location of scripts and needed files
addpath(genpath(root_dir))
addpath(genpath('/projects/b1081/Scripts')); %add cifti-matlab-master folder to path
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
plot = 1;

%initialize variables
spCorrs=zeros(numperms,1);% initialize mat of values i'm permuting
perm_diffs=zeros(numperms,1);% initialize mat of values i'm permuting

in_str = 'HCP752';
out_str = strcat(in_str, '_LHandvsRHand');
group_str = {'HCP752_LH', 'HCP752_RH'};

LH = []; MH = []; RH = []; % initialize mat of groups of subs; LH = left handers, RH = right handers, MH = 'middle' handers
load([root_dir '/PerezEtAl_HemAsymmetries/needed_files/goodSubs752.mat'])
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
template = ft_read_cifti_mod(template_loc);
brainstructure = single(template.brainstructure);
clear template
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
all_subs_cifti = [true_left_handers true_right_handers];
% matrix with 1's and 0's, will determine if sub's maps will be flipped
% -- will be randomized later
flip_switch = zeros(numSubs,1);

%--------------------------------------------------------------------------
%% GET CORRELATION BETWEEN TRUE MAPS
%-------------------------------------------------------------------------- 
tic;    

[true_overlap_left, true_overlap_flip_left] = makemaps(true_left_hems, true_right_hems, flip_switch, brainstructure);

% get correlation
true_spCorr = corr(true_overlap_left, true_overlap_flip_left);
trueDiff = true_overlap_left - true_overlap_flip_left;
save([output_dir out_str '_true_spatialCorr.mat'],'true_spCorr')               
save([output_dir out_str '_true_diffMap.mat'],'trueDiff')

%     --------------------------------------------------------------------------
    %% RUN PERMUTATIONS - MAKE RANDOMIZED MAPS
%    --------------------------------------------------------------------------
        flip_switch(1:length(numLH)) = 1; % add 1's
        
        for p = 1:numperms % will go through this loop for each permutation
        tic;    
            disp(['Permutation #' num2str(p)])
            rng('shuffle'); % need seed for randomizer

            % randomize 1's and 0's
            ind = randperm(length(flip_switch))';
            for s = 1:length(flip_switch)
                rand_flip_switch(s,1) = flip_switch(ind(s));
            end
                                 
            [overlap_left, overlap_flip_left] = makemaps(groups{g}, true_left_hems, true_right_hems, rand_flip_switch, brainstructure);
            
            % get correlation
            spCorrs(p,t) = corr(overlap_left, overlap_flip_left);
            permDiff = overlap_left - overlap_flip_left;
            diffmats(:,p)= permDiff;
            disp(['Correlation between left hem and flipped left hem for group ' groups{g} ' is ' num2str(spCorrs(p,t))]) 
        toc
        end 
        
        disp(['Mean correlation between left hem and flipped left hem for group ' groups{g} ' is ' num2str(mean(spCorrs(p,t)))]) 
        save([output_dir '/VariantsvsFlippedVariants_' groups{g} '_spCorr_' num2str(numperms) '_thresh_' num2str(threshold(t)) '_permutations.mat'],'spCorrs')
        save([output_dir '/VariantsvsFlippedVariants_' groups{g} '_diffMap_' num2str(numperms) '_thresh_' num2str(threshold(t)) '_permutations.mat'],'diffmats')
        
    %--------------------------------------------------------------------------
    %% PLOT RESULTS
    %--------------------------------------------------------------------------
    
        if plot_results == 1
            figure;
            plot(1:numperms, spCorrs, 'bo')
            ylim([.4 1]);
            hold on
            plot(round(numperms/2), true_spCorr, 'ro', 'MarkerFaceColor', 'r')
            ylabel('r value');
            xlabel('Permutations');
            title(['Left Hem vs Pseudo Left Hem Variant Maps - ' groups{g}]);
            m = findobj(gca, 'Type', 'line');
            hleg1 = legend(m(1:2), 'true map', 'randomized maps', 'Location', 'SouthWest');
            hleg1.FontSize = 14;
            ax = gca;
            ax.FontSize = 14;
            set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.5, 0.7, 0.5, 0.7]);
            print(gcf,[output_dir 'VariantsvsFlippedVariants_' groups{g} '_spCorr_'  num2str(threshold(t)) '_thresh_' num2str(numperms) '_permutations_results.jpg'],'-dpng','-r300');
            close gcf
        end
        toc
    



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
