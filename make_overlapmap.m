%% Overlap Map Making Script
% input: variant maps
% output: a map of proportion of subjects that have overlapping variants at given vertex

clear all
%--------------------------------------------------------------------------
%% PATHS
addpath(genpath('/Users/dianaperez/Box/Dependencies/cifti-matlab-master'));
root_dir = '/Users/dianaperez/Desktop/lateralization_code/';
%root_dir = '/projects/p31161/lateralization_code/';
addpath(genpath(root_dir))
%varmap_loc = [root_dir '/Variant_Maps/new_split_vars/reassigned/'];
varmap_str = '_reassigned.dtseries.nii';
varmap_loc = '/Volumes/RESEARCH_HD/HCP_Variants/new_split_vars/reassigned/'; %location of variant maps 
%template_loc = '/Volumes/RESEARCH_HD/HCP_Variants/new_split_vars/reassigned/100206_border1ectopic2.dtseries.nii';
template_loc = [varmap_loc '100206_border1ectopic2.dtseries.nii'];
out_dir = [root_dir '/testing_output/spatial_distribution/'];
if ~exist(out_dir, 'file')
    mkdir(out_dir)
end
%% VARIABLES
% sample being analyzed
MSC = 0;
HCP = 1;
handedness = 0;
wholeGroup = 0;
% set net_id to 0 to create overlap map for all network variants, otherwise
% set net_id to the id number of a specific network
net_id = 1; 
network_names = {'DMN','Vis','FP','Unassigned','DAN','Unassigned','VAN','Sal','CO','SMd','SMl','Aud','Tpole','MTL','PMN','PON'};

if HCP 
 % 1 for 752 subs (includes sets of relatives) 0 for 384 subs (excludes relatives of subjects)
    if wholeGroup == 1 % 
        load('goodSubs752.mat')
        all_subs = goodSubs752;
        in_str = {'HCP752'};
    elseif wholeGroup == 0
        load('goodSubs384.mat')
        all_subs = goodSubs384;
        in_str = {'HCP384'};
    end
    if handedness
        in_str = {'HCP752_LH', 'HCP752_MH', 'HCP752_RH'};
        LH=[]; MH=[]; RH=[];
        for s = 1:length(all_subs)
            if all_subs(s,2) < -28
                LH = [LH; all_subs(s,:)];
            elseif all_subs(s,2) < 48
                MH = [MH; all_subs(s,:)];
            elseif all_subs(s,2) <= 100
                RH = [RH; all_subs(s,:)];
            end
        end
    end
elseif MSC
    in_str = {'MSC'};
    all_subs = {'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07', 'MSC09', 'MSC10'};
end

subs = {all_subs}; %{LH, MH, RH};

for g = 1:numel(subs)
    out_mat = [root_dir strcat(in_str{g}, '_variant_maps.mat')];
    if exist(out_mat, 'file')
        load(out_mat)
    else
        disp('Creating matrix with variant maps...')
        variant_map_mat = [];
        for n=1:length(subs{g})
            if HCP
                cifti = ft_read_cifti_mod([varmap_loc num2str(subs{g,1}(n)) varmap_str]);
            elseif MSC
                cifti = ft_read_cifti_mod([varmap_loc subs{g,1}{n} varmap_str]);               
            end
            variant_map_mat(:,n) = cifti.data;
        end
        save(out_mat, 'variant_map_mat');
    end
            
    disp('Creating overlap map...')
    template = ft_read_cifti_mod(template_loc);
    [overlap_map] = makemap(variant_map_mat, net_id, template);
    if net_id == 0;
        out_cifti = [out_dir in_str{g} '_overlap_map.dtseries.nii'];
    else
        out_cifti = [out_dir in_str{g} '_' network_names{net_id} '_overlap_map.dtseries.nii'];
    end
    ft_write_cifti_mod(out_cifti, overlap_map);
end

function [overlap_map] = makemap(variant_map_mat, net_id, template)
numsubs = length(variant_map_mat);
if net_id == 0
    group_map = logical(variant_map_mat);
else
    group_map = logical(variant_map_mat==net_id);
end
group_overlap = sum(group_map,2)/numsubs;
overlap_map = template;
overlap_map.data = group_overlap;
end

