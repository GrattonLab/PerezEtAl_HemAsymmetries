%% network symmetry

% load cifti
group_avg = ft_read_cifti_mod('/Users/dianaperez/Desktop/Yeo_17nets.dtseries.nii');
left_hem = group_avg.data(1:29696);
right_hem = group_avg.data(29697:end);

%load surf area file
surf_areas = ft_read_cifti_mod('/Users/dianaperez/Desktop/Thesis/surf_areas_verts.dtseries.nii');
lh_surf = surf_areas.data(1:29696);
rh_surf = surf_areas.data(29697:end);

num_nets = unique(group_avg.data);
net_size = [];

for net = 1:length(num_nets)
    lh_verts = find(left_hem==num_nets(net));
    rh_verts = find(right_hem==num_nets(net));
    net_size(net,1) = length(lh_verts);
    net_size(net,2) = length(lh_verts)/length(left_hem);
    net_size(net,3) = sum(lh_surf(lh_verts));
    net_size(net,4) = sum(lh_surf(lh_verts))/sum(lh_surf);
    net_size(net,5) = length(rh_verts);
    net_size(net,6) = length(rh_verts)/length(right_hem);
    net_size(net,7) = sum(rh_surf(rh_verts));
    net_size(net,8) = sum(rh_surf(rh_verts))/sum(rh_surf);
end

