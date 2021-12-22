%% Run permutations on comparisons of network variant properties

%% Permutations - Differences in Variant Properties Between Handedness Groups
clear all
%--------------------------------------------------------------------------
%% PATHS 
root_dir = '/Users/dianaperez/Desktop/lateralization_code/'; % location of code directory
data_location = [root_dir 'testing_output/'];
output_dir = [root_dir 'testing_output/'];
%% VARIABLES
numperms = 100;% number of permutations
% sample that is being analyzed
MSC = 1;
HCP = 0;
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
    LHand_info = load([data_location in_str '_LH_variants_info.mat']);
    RHand_info = load([data_location in_str '_RH_variants_info.mat']);
    numsubs = height(LHand_info.variants_info.left_hem.group_avg)+height(RHand_info.variants_info.left_hem.group_avg);
    numLhanders = height(LHand_info.variants_info.left_hem.group_avg);
    flip_switch = zeros(numsubs,1);
    flip_switch(floor(1:numLhanders)) = 1;
else
    allSubs_info = load([data_location in_str '_variants_info.mat']);
    numsubs = height(allSubs_info.variants_info.left_hem.group_avg);
    flip_switch = zeros(numsubs,1);
    flip_switch(floor(1:numsubs/2)) = 1;
end


if handedness
    %% Number of variants        
    [true_diff perm_diffs] = permute_values(LHand_info.variants_info.both_hems{:,1}, RHand_info.variants_info.both_hems{:,1}, flip_switch, numperms, 1);  
    p = min([length(find(perm_diffs<true_diff))/numperms, length(find(perm_diffs>true_diff))/numperms]);
    disp(['Number of Variants (permute handedness) - ' num2str(p)])
    p_vals = [p_vals; p];
    output_mat = [true_diff; perm_diffs'];
    save([output_dir out_str 'num_vars.mat'], 'output_mat');
    if plot
        title_str = 'Number of Variants Across Handedness Groups';
        y_label = 'Difference in Number of Variants (left handers - right handers)';
        outfile = [output_dir out_str '_PermutationTesting_NumberofVariants.jpg'];
        plot_permtesting(output_mat, title_str, y_label, outfile); 
    end
   
    %% Number of variant vertices
    [true_diff perm_diffs] = permute_values(LHand_info.variants_info.both_hems{:,2}, RHand_info.variants_info.both_hems{:,2}, flip_switch, numperms, 1);    
    p = min([length(find(perm_diffs<true_diff))/numperms, length(find(perm_diffs>true_diff))/numperms]);
    disp(['Number of variant vertices (permute handedness) - ' num2str(p)])
    p_vals = [p_vals; p];
    output_mat = [true_diff; perm_diffs'];
    save([output_dir out_str 'var_vertices.mat'], 'output_mat');
    if plot
        title_str = 'Number of Variant Vertices Across Handedness Groups';
        y_label = 'Difference in Variant Vertices Across Handedness Groups (left handers - right handers)';
        outfile = [output_dir out_str '_PermutationTesting_numVarVerts.jpg'];
        plot_permtesting(output_mat, title_str, y_label, outfile); 
    end
    
    %% Average variant size (in vertices)
    [true_diff perm_diffs] = permute_values(LHand_info.variants_info.both_hems{:,3}, RHand_info.variants_info.both_hems{:,3}, flip_switch, numperms, 1);    
    p = min([length(find(perm_diffs<true_diff))/numperms, length(find(perm_diffs>true_diff))/numperms]);
    disp(['Average variant size in vertices (permute handedness) - ' num2str(p)])
    p_vals = [p_vals; p];
    output_mat = [true_diff; perm_diffs'];
    save([output_dir out_str 'varSize_verts.mat'], 'output_mat');
    if plot
        title_str = 'Average variant size in vertices across handedness groups';
        y_label = 'Difference in variant size in vertices (left handers - right handers)';
        outfile = [output_dir out_str '_PermutationTesting_variantSizeVerts.jpg'];
        plot_permtesting(output_mat, title_str, y_label, outfile); 
    end

    %% Average variant size (in surface area)
    [true_diff perm_diffs] = permute_values(LHand_info.variants_info.both_hems{:,4}, RHand_info.variants_info.both_hems{:,4}, flip_switch, numperms, 1);    
    p = min([length(find(perm_diffs<true_diff))/numperms, length(find(perm_diffs>true_diff))/numperms]);
    disp(['Average Variant Surface Area (permute handedness) - ' num2str(p)])
    p_vals = [p_vals; p];
    output_mat = [true_diff; perm_diffs'];
    save([output_dir out_str 'var_surfArea.mat'], 'output_mat');
    if plot
        title_str = 'Average variant surface area across handedness groups';
        y_label = 'Difference in variant surface area (left handers - right handers)';
        outfile = [output_dir out_str '_PermutationTesting_VariantSurfArea.jpg'];
        plot_permtesting(output_mat, title_str, y_label, outfile); 
    end

else
    
    %% Number of variants
    [true_diff perm_diffs] = permute_values(allSubs_info.variants_info.left_hem.group_avg{:,1}, allSubs_info.variants_info.right_hem.group_avg{:,1}, flip_switch, numperms, 0);
    p = min([length(find(perm_diffs<true_diff))/numperms, length(find(perm_diffs>true_diff))/numperms]);
    disp(['Number of Variants (permute hemispheres) - ' num2str(p)])
    p_vals = [p_vals; p];
    output_mat = [true_diff; perm_diffs'];
    save([output_dir out_str 'num_vars.mat'], 'output_mat');
    if plot
        title_str = 'Number of Variants Across Hemispheres';
        y_label = 'Difference in Number of Variants (left hem - right hem)';
        outfile = [output_dir out_str '_PermutationTesting_NumberofVariants.jpg'];
        plot_permtesting(output_mat, title_str, y_label, outfile);
    end
    
    %% Number of variant vertices
    [true_diff perm_diffs] = permute_values(allSubs_info.variants_info.left_hem.group_avg{:,2}, allSubs_info.variants_info.right_hem.group_avg{:,2}, flip_switch, numperms, 0);
    p = min([length(find(perm_diffs<true_diff))/numperms, length(find(perm_diffs>true_diff))/numperms]);
    disp(['Number of Variants (permute hemispheres) - ' num2str(p)])
    p_vals = [p_vals; p];
    output_mat = [true_diff; perm_diffs'];
    save([output_dir out_str 'num_verts.mat'], 'output_mat'); 
    if plot
        title_str = 'Number of Variant Vertices Across Hemispheres';
        y_label = 'Difference in Number of Variant Vertices (left hem - right hem)';
        outfile = [output_dir out_str '_PermutationTesting_NumberofVariantVertices.jpg'];
        plot_permtesting(output_mat, title_str, y_label, outfile);
    end

    %% Average variant size (vertices)
    [true_diff perm_diffs] = permute_values(allSubs_info.variants_info.left_hem.group_avg{:,3}, allSubs_info.variants_info.right_hem.group_avg{:,3}, flip_switch, numperms, 0);
    p = min([length(find(perm_diffs<true_diff))/numperms, length(find(perm_diffs>true_diff))/numperms]);
    disp(['Average Variant Size in Vertices (permute hemispheres) - ' num2str(p)])
    p_vals = [p_vals; p];
    output_mat = [true_diff; perm_diffs'];
    save([output_dir out_str 'var_size_verts.mat'], 'output_mat');
    if plot
        title_str = 'Average Variant Size in Vertices Across Hemispheres';
        y_label = 'Difference in Average Variant Size in Vertices (left hem - right hem)';
        outfile = [output_dir out_str '_PermutationTesting_VariantSizeVerts.jpg'];
        plot_permtesting(output_mat, title_str, y_label, outfile);
    end
    
    %% Average variant size (surface area)
    [true_diff perm_diffs] = permute_values(allSubs_info.variants_info.left_hem.group_avg{:,4}, allSubs_info.variants_info.right_hem.group_avg{:,4}, flip_switch, numperms, 0);
    p = min([length(find(perm_diffs<true_diff))/numperms, length(find(perm_diffs>true_diff))/numperms]);
    p_vals = [p_vals; p];
    disp(['Average Variant Size in Surface Area (permute hemispheres) - ' num2str(p)])
    output_mat = [true_diff; perm_diffs'];
    save([output_dir out_str 'var_size_surf_area.mat'], 'output_mat');
    if plot
        title_str = 'Average Variant Size in Surface Area Across Hemispheres';
        y_label = 'Difference in Average Variant Size in Surface Area (left hem - right hem)';
        outfile = [output_dir out_str '_PermutationTesting_VariantSizeSurfArea.jpg'];
        plot_permtesting(output_mat, title_str, y_label, outfile);
    end

end
%final_pvals = p_vals(find(p_vals(:,3)>0),3);
[p_fdr p_masked] = FDR(p_vals, .025);


function [true_diff perm_diffs] = permute_values(left, right, flip_switch, numperms, hand, hemi, index)
    if hand
        true_diff = mean(left) - mean(right);
        all_subs = [left; right];
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



function plot_permtesting(perm_diffs, title_str, y_label, outfile)
ind = zeros(length(perm_diffs),1);
ind(1,1) = 1;
handles = plotSpread(perm_diffs, 'categoryMarkers', {'o', '.'}, 'categoryLabels', {'Permuted Difference', 'True Difference'}, 'categoryColors', {'k', 'r'}, 'categoryIdx', [ind])
scatter(1,perm_diffs(1,:), 'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red', 'SizeData', 50)
ax = gca;
ax.FontSize = 12;
axis([0.5, 1.5, min(perm_diffs)-5, max(perm_diffs)+5])
title(title_str)
ylabel(y_label)
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.3, 0.3, 0.3, 0.7]); %first and second control position on screen, third controls width, and fourth controls height
print(gcf,outfile,'-dpng','-r300');
close gcf
end

