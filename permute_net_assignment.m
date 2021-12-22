%% Permutations - Differences in Variant Properties Between Handedness Groups
clear all
%--------------------------------------------------------------------------
%% PATHS 
root_dir = '/Users/dianaperez/Desktop/lateralization_code/'; % location of code directory
data_location = [root_dir 'testing_output/'];
output_dir = [root_dir 'testing_output/'];
%% VARIABLES
numperms = 1000;% number of permutations
% sample that is being analyzed
MSC = 0;
HCP = 1;
handedness = 0;
plot = 1;

%initialize variables
p_vals = [];
perm_diffs=zeros(numperms,1);% initialize mat of values i'm permuting

if MSC
    in_str = 'MSC';
    out_str = append(in_str, '_LhemvsRhem_');
else
    if handedness
        in_str = 'HCP752';
        out_str = append(in_str, '_LHandvsRHand_');
    else
        in_str = 'HCP384';
        out_str = append(in_str, '_LhemvsRhem_');
    end
end

%% LOAD files
if handedness
    LHand_info = load([data_location in_str '_LH_networksxHem.mat']);
    RHand_info = load([data_location in_str '_RH_networksxHem.mat']);
    numsubs = height(LHand_info.networksxHem.clustersLH)+height(RHand_info.networksxHem.clustersLH);
    numLhanders = height(LHand_info.networksxHem.clustersLH);
    flip_switch = zeros(numsubs,1);
    flip_switch(floor(1:numLhanders)) = 1;
else
    allSubs_info = load([data_location in_str '_networksxHem.mat']);
    numsubs = height(allSubs_info.networksxHem.clustersLH);
    flip_switch = zeros(numsubs,1);
    flip_switch(floor(1:numsubs/2)) = 1;
end

network_names = {'DMN','Vis','FP',' ','DAN',' ','VAN','Sal','CO','SMd','SMl','Aud','Tpole','MTL','PMN','PON'};

for net = 1:numel(network_names)
    if handedness       
        [true_diff perm_diffs] = network_comparisons(LHand_net_info, RHand_net_info, flip_switch, numperms, 1, 0, net);
        p = min([length(find(perm_diffs<true_diff))/numperms, length(find(perm_diffs>true_diff))/numperms]);
        disp(['Number of Variants assigned to ' network_names{net} ' (permute handedness) - ' num2str(p)])
    else
        [true_diff perm_diffs] = network_comparisons(allSubs_info.networksxHem.clustersLH(:,net), allSubs_info.networksxHem.clustersRH(:,net), flip_switch, numperms, 0, 1, net);
        p = min([length(find(perm_diffs<true_diff))/numperms, length(find(perm_diffs>true_diff))/numperms]);
        disp(['Number of Variants assigned to ' network_names{net} ' (permute hemisphere) - ' num2str(p)])
    end 
    p_vals = [p_vals; p];
    output_mat(:,net) = [true_diff; perm_diffs'];
end

save([output_dir 'network_perms_numVars.mat'], 'output_mat');
[p_fdr p_masked] = FDR(p_vals, .025);

if plot
    network_names = {'DMN','Vis','FP','DAN','VAN','Sal','CO','SMd','SMl','Aud','Tpole','MTL','PMN','PON'};
    rgb_colors = [1 0 0; %DMN   
              0 0 .6; %Vis
              1 1 0; %FP
              0 .8 0; %DAN
              0 .6 .6;%VAN
              0 0 0; % Sal
              .3 0 .6; %CON
              .2 1 1; %SMd
              1 .5 0; % SMl
              .6 .2 1; %Aud
              .2 1 .2; %Tpole
              0 .2 .4; %MTL
              0 0 1; %PMN
              .8 .8 .6]; %PON
    good_nets = [output_mat(:,1:3) output_mat(:,5) output_mat(:,7:end)];
    rgb_for_plot = {};
    for n = 1:numer(network_names)
        rgb_for_plot{n} = rgb_colors(n,:);
    end
    ind = zeros(1001,1);
    ind(1,1) = 1;
    net_inds = [ind ind ind ind ind ind ind ind ind ind ind ind ind ind];
    handles = plotSpread(good_nets, 'categoryMarkers', {'x', '.'}, 'categoryLabels', {'Permuted Differences','True Difference'}, 'distributionColors', rgb_for_plot', 'categoryIdx', net_inds, 'xNames', network_names)
    ax = gca;
    ax.FontSize = 24;
    scatter([1:14],good_nets(1,:), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k', 'SizeData', 50)
    scatter(6,good_nets(1,6), 'MarkerEdgeColor', 'w', 'MarkerFaceColor', 'w', 'SizeData', 50)
    ylabel('Difference in Number of Variants')
    xlabel('Assigned Network')
    axis([0, 15, -0.4, .7])
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.3, 0.3, 0.9, 0.7]);
    print(gcf,[output_dir out_str '_PermutationTesting_networkAssignment.jpg'],'-dpng','-r300');
    close gcf
end
    

function [true_diff perm_diffs] = network_comparisons(left, right, flip_switch, numperms, hand, hemi, net)
if hand
    lefthand_lefthem = left.networksxHem.clustersLH(:,net);
    lefthand_righthem = left.networksxHem.clustersRH(:,net);
    righthand_lefthem = right.networksxHem.clustersLH(:,net);
    righthand_righthem = right.networksxHem.clustersRH(:,net);
    
    if hemi
        lefthand_diff = lefthand_lefthem - lefthand_righthem;
        lefthand_truediff = mean(lefthand_diff);
        righthand_diff = righthand_lefthem - righthand_righthem;
        righthand_truediff = mean(righthand_diff);
        true_diff = lefthand_truediff - righthand_truediff; 
        all_subs = [lefthand_diff;righthand_diff];
    else
        left = lefthand_lefthem + lefthand_righthem;
        right = righthand_lefthem + righthand_righthem;
        true_diff = mean(left) - mean(right);
        all_subs = [left; right];
    end
    for p = 1:numperms
        rng('shuffle');
        ind = randperm(length(flip_switch))';
        rand_flip_switch = flip_switch(ind);
        pseudo_left = all_subs(rand_flip_switch==1);
        pseudo_right = all_subs(rand_flip_switch==0);
        perm_diffs(p) = mean(pseudo_left) - mean(pseudo_right);
    end
else
    % left_hem = allSubs_info.variants_info.left_hem.group_avg{:,index};
    % right_hem = allSubs_info.variants_info.right_hem.group_avg{:,index};
    true_diff = mean(left - right);
    for p = 1:numperms
        rng('shuffle');
        ind = randperm(length(flip_switch))';
        rand_flip_switch = flip_switch(ind);
        for s = 1:length(left)
            if rand_flip_switch(s) == 1
                pseudo_left(s,1) = left(s);
                pseudo_right(s,1) = right(s);
            elseif rand_flip_switch(s) == 0
                pseudo_right(s,1) = left(s);
                pseudo_left(s,1) = right(s);
            end
        end
        perm_diffs(p) = mean(pseudo_left - pseudo_right);
    end
    
end

end



