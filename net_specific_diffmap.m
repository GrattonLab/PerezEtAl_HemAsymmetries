%% MAKE DIFF MAP NET SPECIFIC

clear all

%--------------------------------------------------------------------------
%% PATHS
%--------------------------------------------------------------------------

%root_dir = '/projects/p31161/lateralization_code/'; %location of scripts and needed files
root_dir = '/Users/dianaperez/Desktop/lateralization_code/';
addpath(genpath(root_dir))
%addpath(genpath('/projects/b1081/Scripts')); %add cifti-matlab-master folder to path
varmap_loc = '/Volumes/RESEARCH_HD/HCP_Variants/new_split_vars/reassigned/';
%varmap_loc = [data_dir '/Variant_Maps/new_split_vars/reassigned/']; %location of variant maps
varmap_str = '_uniqueIDs_afterReassign.dtseries.nii';
net_assign_str = '_reassigned.dtseries.nii';
template_loc = [varmap_loc '100206_border1ectopic2.dtseries.nii'];
output_dir = root_dir;
if ~exist(output_dir, 'file')
    mkdir(output_dir)
end

%--------------------------------------------------------------------------
%% VARIABLES
%--------------------------------------------------------------------------
network_id = 9;
in_str = 'HCP752';
out_str = strcat(in_str, '_LHand-RHand');
group_str = {strcat(in_str, '_LH'), strcat(in_str, '_RH')};
network_names = {'DMN','Vis','FP',' ','DAN',' ','VAN','Sal','CO','SMd','SMl','Aud','Tpole','MTL','PMN','PON'};

LH = []; MH = []; RH = []; % initialize mat of groups of subs; LH = left handers, RH = right handers, MH = 'middle' handers
load([root_dir '/needed_files/goodSubs752.mat'])
all_subs = goodSubs752;
for s = 1:length(all_subs)
    if all_subs(s,2) < -40
        LH = [LH; all_subs(s,:)];
    elseif all_subs(s,2) < 40
        MH = [MH; all_subs(s,:)];
    elseif all_subs(s,2) <= 100
        RH = [RH; all_subs(s,:)];
    end
end

subs = {LH; RH};
   

%--------------------------------------------------------------------------
%% BEGIN ANALYSIS
%--------------------------------------------------------------------------
left_cifti_loc = [root_dir strcat(in_str, 'LH_net_assignments.mat')];
right_cifti_loc = [root_dir strcat(in_str, 'RH_net_assignments.mat')];
if exist(left_cifti_loc, 'file')
    load(left_cifti_loc)   
else        
    for n = 1:length(LH)
%         true_left_handers = [];
        cifti = ft_read_cifti_mod([varmap_loc num2str(LH(n,1)) net_assign_str]);
        LHand_net_assign(:,n) = cifti.data;
    end
    save(left_cifti_loc, 'LHand_net_assign')
end
if exist(right_cifti_loc, 'file')
    load(right_cifti_loc)
else
    for n = 1:length(RH)
        cifti = ft_read_cifti_mod([varmap_loc num2str(RH(n,1)) net_assign_str]);
        Rhand_net_assign(:,n) = cifti.data;
    end
    save(right_cifti_loc, 'Rhand_net_assign')
end

data = {LHand_net_assign, Rhand_net_assign};
template = ft_read_cifti_mod(template_loc);

numsubs = size(LHand_net_assign,2);
LH_bin_data = logical(LHand_net_assign==network_id);
LH_overlap_map = sum(LH_bin_data,2)/numsubs;
template.data = LH_overlap_map;
ft_write_cifti_mod([output_dir group_str{1} '_' network_names{network_id} '_Variants_Overlap.dtseries.nii'], template);

numsubs = size(Rhand_net_assign,2);
RH_bin_data = logical(Rhand_net_assign==network_id);
RH_overlap_map = sum(RH_bin_data,2)/numsubs;
template.data = RH_overlap_map;
ft_write_cifti_mod([output_dir group_str{2} '_' network_names{network_id} '_Variants_Overlap.dtseries.nii'], template);

diff_map = LH_overlap_map - RH_overlap_map;
template.data = diff_map;
ft_write_cifti_mod([output_dir out_str '_' network_names{network_id} '_Variants_DiffMap.dtseries.nii'], template);


