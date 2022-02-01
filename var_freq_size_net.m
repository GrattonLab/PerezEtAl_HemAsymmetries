%% Variant frequency, size, and network assignment analysis

clear all

%--------------------------------------------------------------------------
%% PATHS 
root_dir = '/Users/dianaperez/Desktop/'; % location of code directory
data_location = '/Volumes/RESEARCH_HD/HCP_Variants/new_split_vars/reassigned/';
%data_location = '/Volumes/RESEARCH_HD/HCP_Variants/MSC_split/';
output_dir = '/Users/dianaperez/Desktop/lateralization_code/testing_output/';
if ~exist(output_dir)
    mkdir(output_dir)
end

% file name for relevant files (variant map, network assignment and border/ectopic labels
varmap_str = '_uniqueIDs_afterReassign.dtseries.nii';
netAssign_str = '_reassigned.dtseries.nii';

%load file with surface area measurement for each vertex
surf_areas = ft_read_cifti_mod([root_dir 'lateralization_code/needed_files/surf_areas_verts.dtseries.nii']);

%% VARIABLES
% sample that is being analyzed
MSC = 0;
HCP = 1;
if HCP
    wholeGroup = 1;
    handedness = 1;
end

plot = 1;
% surface areas for each hem
surf_areas_LHem = surf_areas.data(1:29696);
surf_areas_RHem = surf_areas.data(29697:end);
clear surf_areas

% initialize variables
variants_info.left_hem = [];
variants_info.right_hem = [];
tmp_LHem = []; tmp_LHem_group = []; tmp_RHem = []; tmp_RHem_group = []; 
networksxHem = [];

% load subs
if HCP
    if wholeGroup == 1 % 
        out_str = {'HCP752'};
        load([root_dir '/lateralization_code/needed_files/goodSubs752.mat'])
        all_subs = goodSubs752;
    elseif wholeGroup == 0
        out_str = {'HCP384'};
        load([root_dir '/lateralization_code/needed_files/goodSubs384.mat'])
        all_subs = goodSubs384;
    end
    if handedness 
        out_str = {'HCP752_LH', 'HCP752_RH', 'HCP752_MH'}; % output files will contain these strings
        LH = []; MH = []; RH = [];
        for s = 1:length(all_subs)
                if all_subs(s,2) < -40
                    LH = [LH; all_subs(s,:)];
                elseif all_subs(s,2) < 40
                    MH = [MH; all_subs(s,:)];
                elseif all_subs(s,2) <= 100
                    RH = [RH; all_subs(s,:)];
                end
        end
        subs = {LH;RH;MH};
    else
        subs = {all_subs};
    end
elseif MSC == 1
    out_str = {'MSC'};
    all_subs = {'MSC01', 'MSC02', 'MSC03', 'MSC04', 'MSC05', 'MSC06', 'MSC07', 'MSC09', 'MSC10'}';
    subs = {all_subs};
end

for g = 1:numel(subs)
    for n = 1:length(subs{g})        
        %output files
        outfile = sprintf('%s/%s_variants_info.mat', output_dir, out_str{g});
        outfile2 = sprintf('%s/%s_networksxHem.mat', output_dir, out_str{g});
        
        tmp = [];         
        %% Load variant maps with unique ID's, network assignments, and border/ectopic info
        if HCP
            net_assign = ft_read_cifti_mod([data_location num2str(subs{g,1}(n)) netAssign_str]);
            unique_IDs = ft_read_cifti_mod([data_location num2str(subs{g,1}(n)) varmap_str]);
        else
            net_assign = ft_read_cifti_mod([data_location subs{g,1}{n} netAssign_str]);
            unique_IDs = ft_read_cifti_mod([data_location subs{g,1}{n} varmap_str]);
        end
       
        %% Separate data by hemispheres
        Lhem(:,1) = unique_IDs.data(1:29696,1);
        Rhem(:,1) = unique_IDs.data(29697:end,1);
        Lhem(:,2) = net_assign.data(1:29696,1);
        Rhem(:,2) = net_assign.data(29697:end,1);

        
        %% Get unique variant IDs
        Lhem_var_IDs = unique(Lhem(:,1));
        Rhem_var_IDs = unique(Rhem(:,1));
        % get rid of 0's (non-variant)
        Lhem_var_IDs = Lhem_var_IDs(2:end);
        Rhem_var_IDs = Rhem_var_IDs(2:end);
        
        if isempty(Lhem_var_IDs) %if sub has no variants in that hemisphere
            disp([num2str(subs(n)) ': no variant info']);
            networksxHem.clustersLH(n,:) = 0;
            networksxHem.verticesLH(n,:) = 0;
            tmp_LHem_group(n,:) = [0 0 0];   
        else
            sub_info = [];

            % get info for each variant in left hem
            for v = 1:length(Lhem_var_IDs)
                verts = find(Lhem(:,1)==Lhem_var_IDs(v)); % get vertices for each variant
                num_verts = length(verts); % variant size (number of vertices)
                surf_area = sum(surf_areas_LHem(verts)); %variant size (surface area)
                net = Lhem(verts(1),2); %network assignment
                sub_info(v,:) = [Lhem_var_IDs(v) num_verts surf_area net]; %save to matrix
            end
            
            for net = 1:16 %find number of variants assigned to each network
                vars = find(sub_info(:,4)==net);
                networksxHem.clustersLH(n,net) = length(vars); % number of clusters assigned to each network
                networksxHem.verticesLH(n,net) = sum(sub_info(vars,2)); % number of variant vertices assigned to each network
                networksxHem.surfaceLH(n,net) = sum(sub_info(vars,3)); % sum of surface area of variants assigned to each network
            end
            
            tmp_LHem = [tmp_LHem; sub_info]; %save to matrix
            tmp = [tmp; sub_info];
            num_vars = size(sub_info,1); %number of variants for each subject
            num_verts = sum(sub_info(:,2)); % number of variant vertices (sum across all variants)
            avg_num_verts = mean(sub_info(:,2)); %average variant size
            avg_surf_area = mean(sub_info(:,3));
            
            tmp_LHem_group(n,:) = [num_vars num_verts avg_num_verts avg_surf_area]; %save to matrix
        end
        
        if isempty(Rhem_var_IDs)
            disp([num2str(subs(n)) ': subjects RH skipped, no variant info']);
        else
            %do the same for the right hemisphere
            sub_info = [];
            for v = 1:length(Rhem_var_IDs)
                verts = find(Rhem(:,1)==Rhem_var_IDs(v));
                num_verts = length(verts);
                surf_area = sum(surf_areas_RHem(verts));
                net = Rhem(verts(1),2);
                sub_info(v,:) = [Rhem_var_IDs(v) num_verts surf_area net];
            end
            
            for net = 1:16
                vars = find(sub_info(:,4)==net);
                networksxHem.clustersRH(n,net) = length(vars);
                networksxHem.verticesRH(n,net) = sum(sub_info(vars,2));
                networksxHem.surfaceLH(n,net) = sum(sub_info(vars,3));
            end
            tmp_RHem = [tmp_RHem; sub_info];
            tmp = [tmp; sub_info];
            num_vars = size(sub_info,1);
            num_verts = sum(sub_info(:,2));
            avg_num_verts = mean(sub_info(:,2));
            avg_surf_area = mean(sub_info(:,3));
            tmp_RHem_group(n,:) = [num_vars num_verts avg_num_verts avg_surf_area];
        end
        both_hems(n,:) = [length(tmp) sum(tmp(:,2)) mean(tmp(:,2)) mean(tmp(:,3))]; 
end
    %add labels to columns and save to structure
    table = array2table(tmp_LHem_group, 'VariableNames', {'numVars', 'numVerts', 'avgVarSize', 'avgSurfArea'});
    variants_info.left_hem.group_avg = table; clear tmp_LH_group
    table = array2table(tmp_RHem_group, 'VariableNames', {'numVars', 'numVerts', 'avgVarSize', 'avgSurfArea'});
    variants_info.right_hem.group_avg = table; clear tmp_RH_group
    table = array2table(both_hems, 'VariableNames', {'numVars', 'numVerts', 'avgVarSize', 'avgSurfArea'});
    variants_info.both_hems = table; clear both_hems
    
    %save structure to output file
    save(outfile, 'variants_info');
    save(outfile2, 'networksxHem');
    
if plot
    % plot number of vars
    leftHem = variants_info.left_hem.group_avg{:,1};
    rightHem = variants_info.right_hem.group_avg{:,1};
    [h,mu,sigma,q,notch] = al_goodplot(leftHem,0.1,.2,[ 0    0.4470    0.7410],'left');
    [h,mu,sigma,q,notch] = al_goodplot(rightHem,0.11,.2,[0.8500    0.3250    0.0980],'right');
    xticks([0.105])
    xticklabels('Number of Variant Regions')
    axis([-0.2, 0.5, min([leftHem; rightHem])-std([leftHem; rightHem]), max([leftHem; rightHem])+std([leftHem; rightHem])])
    ax = gca;
    ax.FontSize = 24;
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.3, 0.3, 0.3, 0.9]); %first and second control position on screen, third controls width, and fourth controls height
    m = findobj(gca, 'Type', 'Scatter');
    hleg1 = legend(m(1:2), 'Right Hemisphere', 'Left Hemisphere', 'Location', 'SouthOutside');
    print(gcf, [output_dir cell2mat(out_str(g)) '_numVars.jpg'], '-dpng', '-r300')
    close gcf
    % plot number of var vertices
    leftHem = variants_info.left_hem.group_avg{:,2};
    rightHem = variants_info.right_hem.group_avg{:,2};
    [h,mu,sigma,q,notch] = al_goodplot(leftHem,0.1,.2,[ 0    0.4470    0.7410],'left');
    [h,mu,sigma,q,notch] = al_goodplot(rightHem,0.11,.2,[0.8500    0.3250    0.0980],'right');
    xticks([0.105])
    xticklabels('Number of Variant Vertices')
    axis([-0.2, 0.5, min([leftHem; rightHem])-std([leftHem; rightHem]), max([leftHem; rightHem])+std([leftHem; rightHem])])
    ax = gca;
    ax.FontSize = 24;
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.3, 0.3, 0.3, 0.9]); %first and second control position on screen, third controls width, and fourth controls height
    m = findobj(gca, 'Type', 'Scatter');
    hleg1 = legend(m(1:2), 'Right Hemisphere', 'Left Hemisphere', 'Location', 'SouthOutside');
    print(gcf, [output_dir cell2mat(out_str(g)) '_numVarVerts.jpg'], '-dpng', '-r300')
    close gcf
    % Average variant size in vertices
    leftHem = variants_info.left_hem.group_avg{:,3};
    rightHem = variants_info.right_hem.group_avg{:,3};
    [h,mu,sigma,q,notch] = al_goodplot(leftHem,0.1,.2,[ 0    0.4470    0.7410],'left');
    [h,mu,sigma,q,notch] = al_goodplot(rightHem,0.11,.2,[0.8500    0.3250    0.0980],'right');
    xticks([0.105])
    xticklabels('Average Variant Size')
    ylabel('Number of vertices')
    axis([-0.2, 0.5, min([leftHem; rightHem])-std([leftHem; rightHem]), max([leftHem; rightHem])+std([leftHem; rightHem])])
    ax = gca;
    ax.FontSize = 24;
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.3, 0.3, 0.3, 0.9]); %first and second control position on screen, third controls width, and fourth controls height
    m = findobj(gca, 'Type', 'Scatter');
    hleg1 = legend(m(1:2), 'Right Hemisphere', 'Left Hemisphere', 'Location', 'SouthOutside');
    print(gcf, [output_dir cell2mat(out_str(g)) '_avgSizeVerts.jpg'], '-dpng', '-r300')
    close gcf
    % Average variant size in surface area
    leftHem = variants_info.left_hem.group_avg{:,4};
    rightHem = variants_info.right_hem.group_avg{:,4};
    [h,mu,sigma,q,notch] = al_goodplot(leftHem,0.1,.2,[ 0    0.4470    0.7410],'left');
    [h,mu,sigma,q,notch] = al_goodplot(rightHem,0.11,.2,[0.8500    0.3250    0.0980],'right');
    xticks([0.105])
    xticklabels('Average Variant Size')
    ylabel('Surface area (mm^2)')
    axis([-0.2, 0.5, min([leftHem; rightHem])-std([leftHem; rightHem]), max([leftHem; rightHem])+std([leftHem; rightHem])])
    ax = gca;
    ax.FontSize = 24;
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.3, 0.3, 0.3, 0.9]); %first and second control position on screen, third controls width, and fourth controls height
    m = findobj(gca, 'Type', 'Scatter');
    hleg1 = legend(m(1:2), 'Right Hemisphere', 'Left Hemisphere', 'Location', 'SouthOutside');
    print(gcf, [output_dir cell2mat(out_str(g)) '_avgSizeSurfArea.jpg'], '-dpng', '-r300')
    close gcf
    % Network assignment
    network_names = {'DMN'	'Vis'	'FP'	'DAN'	'VAN'	'Sal'	'CO'	'SMd'	'SMl'	'Aud'	'Tpole'	'MTL'	'PMN'	'PON'};
    good_nets_LH = [networksxHem.clustersLH(:,1:3) networksxHem.clustersLH(:,5) networksxHem.clustersLH(:,7:end)];
    good_nets_RH = [networksxHem.clustersRH(:,1:3) networksxHem.clustersRH(:,5) networksxHem.clustersRH(:,7:end)];
    good_nets_LH = mean(good_nets_LH);
    good_nets_RH = mean(good_nets_RH);
    means = [good_nets_LH' good_nets_RH'];
    bar(1:14, means)
    legend('Left Hem', 'Right Hem')
    xticks(1:14)
    xticklabels(network_names)
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.5, 0.7, 0.9, 0.7]);
    ylabel('Average Number of Variants');
    xlabel('Network');
    title('Variants Assigned to Each Network Across Hemispheres');
    ax = gca;
    ax.FontSize = 24;
    print(gcf,[output_dir cell2mat(out_str(g)) '_netAssignment.jpg'],'-dpng','-r300');
    close gcf
end
end


