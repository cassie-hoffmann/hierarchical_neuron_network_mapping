function [cc_nodes, node_positions, adj_matrix] = create_network_from_stacks(node_stack, edge_stack, node_stack_dil, dendrite_stack, min_volume)

unique_labels = unique(node_stack_dil);
unique_labels(unique_labels == 0) = [];  % Remove background

cc_nodes.Connectivity = 26;
cc_nodes.ImageSize = size(node_stack_dil);
cc_nodes.NumObjects = numel(unique_labels);

[rows, cols, z] = ind2sub(size(edge_stack),find(node_stack_dil));  % indices + label values
lin_idx = sub2ind(size(edge_stack),rows(:), cols(:), z(:));
[label_idx, ~, label_vals] = find(node_stack_dil);


% Get unique nonzero labels
unique_labels = unique(node_stack_dil);
unique_labels(unique_labels == 0) = [];

% prepare new connected component structure
cc_nodes = struct();
cc_nodes.Connectivity = 26;
cc_nodes.ImageSize = size(node_stack_dil);
cc_nodes.PixelIdxList = {};
cc_nodes.NumObjects = 0;

label_offset = 0; % counts how many new labels added so far

for i = 1:numel(unique_labels)
    lbl = unique_labels(i);

    % Mask for current label
    mask = node_stack_dil == lbl;

    % Find separate connected components within this label
    cc_local = bwconncomp(mask, 26);

    if cc_local.NumObjects == 0
        continue;
    end

    % For the first connected component, keep its spot (original)
    cc_nodes.NumObjects = cc_nodes.NumObjects + 1;
    cc_nodes.PixelIdxList{cc_nodes.NumObjects} = cc_local.PixelIdxList{1};

    % For remaining components (if any), add them as new unique labels
    for j = 2:cc_local.NumObjects
        cc_nodes.NumObjects = cc_nodes.NumObjects + 1;
        cc_nodes.PixelIdxList{cc_nodes.NumObjects} = cc_local.PixelIdxList{j};
    end
end

% Convert to label matrix (each object has a unique label)
node_stack_dil = labelmatrix(cc_nodes);


%%
disp('somata labeled');

stats_nodes = regionprops3(node_stack_dil, 'Centroid', 'Volume');
disp('stats nodes completed');
% filter out debris
%min_volume =  550000; %was 400 and 400000 and 4500000
valid_nodes = stats_nodes.Volume >= min_volume;
stats_nodes = stats_nodes(valid_nodes, :);

cc_nodes.PixelIdxList = cc_nodes.PixelIdxList(valid_nodes);
cc_nodes.NumObjects = numel(cc_nodes.PixelIdxList);
cc_nodes.Size = cellfun(@numel, cc_nodes.PixelIdxList);
combined_mask = edge_stack | dendrite_stack;


overlap_counts = zeros(cc_nodes.NumObjects, 1);


se = strel('sphere', 10);  % choose radius in voxels (e.g. 1–3 voxels)
cc_nodes.NeuriteOverlap = zeros(cc_nodes.NumObjects, 1);

for iii = 1:cc_nodes.NumObjects
    % Extract current node’s voxel indices
    idx = cc_nodes.PixelIdxList{iii};

    % Create a temporary mask for this node (minimal bounding box)
    tempMask = false(size(combined_mask));
    tempMask(idx) = true;

    % Dilate this node locally
    tempMaskDil = imdilate(tempMask, se);

    overlap_voxels = combined_mask & tempMaskDil;

    cc_nodes.NeuriteOverlap(iii,1) = nnz(overlap_voxels); % number of overlapping voxels
end


num_nodes = height(stats_nodes);
node_positions = stats_nodes.Centroid;



disp('node positions calculated');
fprintf('stats_nodes rows:  %d\n', height(stats_nodes));
fprintf('pixelidxlist entries   %d\n', numel(cc_nodes.PixelIdxList));
fprintf('node_positions rows:   %d\n', size(node_positions,1));



% find connections

cc_edges_axons = bwconncomp(edge_stack, 26);
cc_edges_dendrites = bwconncomp(dendrite_stack, 26);
adj_matrix = zeros(num_nodes);

se2d = strel('disk', 10, 0);



%% For each edge component
for i = 1:cc_edges_axons.NumObjects
    % Create mask for specific neurite segment
    edge_mask = false(size(edge_stack));
    dendrite_mask = false(size(dendrite_stack));
    edge_mask(cc_edges_axons.PixelIdxList{i}) = true;

    %dilate segment at relevant slices

    nonEmptySlices = find(any(any(edge_mask, 1), 2));

    for ii = 1:numel(nonEmptySlices)
        z = nonEmptySlices(ii);
        edge_mask(:,:,z) = imdilate(edge_mask(:,:,z), se2d);
    end


    is_connected = false(1, num_nodes);

    for j = 1:num_nodes
        node_pixels = cc_nodes.PixelIdxList{j};
        if any(edge_mask(node_pixels))
            is_connected(j) = true;
        end
    end

    connected_nodes = find(is_connected);

    fprintf('established connection %d of %d\n', i, cc_edges_axons.NumObjects);

    %  local adjacency matrix for this component
    local_adj = zeros(num_nodes);
    if length(connected_nodes) >= 2
        for k = 1:length(connected_nodes)-1
            for l = k+1:length(connected_nodes)
                local_adj(connected_nodes(k), connected_nodes(l)) = 1;
                local_adj(connected_nodes(l), connected_nodes(k)) = 1;
            end
        end
    end


    adj_matrix = adj_matrix + local_adj;


end

%% add dendrites

for i = 1:cc_edges_dendrites.NumObjects
    % Create mask for specific dendrite segment
    edge_mask = false(size(edge_stack));
    edge_mask(cc_edges_dendrites.PixelIdxList{i}) = true;

    % Dilate segment at relevant slices
    nonEmptySlices = find(any(any(edge_mask, 1), 2));
    for ii = 1:numel(nonEmptySlices)
        z = nonEmptySlices(ii);
        edge_mask(:,:,z) = imdilate(edge_mask(:,:,z), se2d);
    end

    % Check connectivity to nodes
    is_connected = false(1, num_nodes);
    for j = 1:num_nodes
        node_pixels = cc_nodes.PixelIdxList{j};
        if any(edge_mask(node_pixels))
            is_connected(j) = true;
        end
    end

    connected_nodes = find(is_connected);

    fprintf('established dendrite connection %d of %d\n', i, cc_edges_dendrites.NumObjects);

    % Local adjacency for this dendrite component
    local_adj = zeros(num_nodes);
    if length(connected_nodes) >= 2
        for k = 1:length(connected_nodes)-1
            for l = k+1:length(connected_nodes)
                local_adj(connected_nodes(k), connected_nodes(l)) = 1;
                local_adj(connected_nodes(l), connected_nodes(k)) = 1;
            end
        end
    end

    % Add connections to existing adjacency matrix
    adj_matrix = adj_matrix + local_adj;
end


end
