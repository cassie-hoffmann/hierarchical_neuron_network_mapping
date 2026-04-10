function cluster_network_data = microscale_neurosphere_mapping(somata_folder, neurite_folder)

addpath 2019_03_03_BCT/;

cnt = 0;

subdirs_s = dir(somata_folder);
subdirs_s = subdirs_s([subdirs_s.isdir]);
subdirs_s = subdirs_s(~ismember({subdirs_s.name}, {'.', '..'}));

subdirs_n = dir(neurite_folder);
subdirs_n = subdirs_n([subdirs_n.isdir]);
subdirs_n = subdirs_n(~ismember({subdirs_n.name}, {'.', '..'}));


somata_files   = dir(fullfile(somata_folder, '*.tiff'));
neurites_files = dir(fullfile(neurite_folder, '*.tiff'));

% ensure both folders have the same number of image stacks
num_stacks = min(length(somata_files), length(neurites_files));

% process each image stack pair
for o = 1:num_stacks
    somata_stack_path  = fullfile(somata_folder, somata_files(o).name);
    neurite_stack_path = fullfile(neurite_folder, neurites_files(o).name);

    fprintf('Processing image %s \n', somata_files(o).name);

    % read neurite stack using Bio-Formats
    reader   = bfGetReader(neurite_stack_path);
    reader_s = bfGetReader(somata_stack_path);
    omeMeta  = reader.getMetadataStore();

    % Get dimensions of neurite image
    sizeX_neurites = reader.getSizeX();
    sizeY_neurites = reader.getSizeY();
    sizeZ = reader.getSizeZ();
    sizeC = reader.getSizeC();

    somata_info = imfinfo(somata_stack_path);

    sizeX_somata = somata_info(1).Width;
    sizeY_somata = somata_info(1).Height;


    sizeX = min(sizeX_somata, sizeX_neurites);
    sizeY = min(sizeY_somata, sizeY_neurites);


    % sanity check for channel count
    if sizeC < 2
        warning('Less than 2 channels found in %s. Skipping file.', neurite_stack_path);
        reader.close();
        continue;
    end

    axon_stack = zeros(sizeY, sizeX, sizeZ);
    dendrite_stack = zeros(sizeY, sizeX, sizeZ);
    somata_stack = zeros(sizeY, sizeX, sizeZ);


    % read each Z slice for each channel
    for z = 1:sizeZ
        axon_index = reader.getIndex(z - 1, 0, 0) + 1;      % Channel 0 = axon
        dendrite_index = reader.getIndex(z - 1, 1, 0) + 1;  % Channel 1 = dendrite

        axon_plane = bfGetPlane(reader, axon_index);
        axon_plane = axon_plane(1:sizeY,1:sizeX);
        dendrite_plane = bfGetPlane(reader, dendrite_index);
        dendrite_plane = dendrite_plane(1:sizeY,1:sizeX);

        axon_stack(:, :, z) = axon_plane;
        dendrite_stack(:, :, z) = dendrite_plane;
        somata_stack(:,:,z) = imread(somata_stack_path, z); % impose neuite stack size because post-stardist upscaling adds extra y pixels
    end


    axon_thresh = otsuThreshold(axon_stack);
    dendrite_thresh = otsuThreshold(dendrite_stack);

    axon_stack = axon_stack > axon_thresh;
    dendrite_stack = dendrite_stack > dendrite_thresh;


    voxel_size_x = omeMeta.getPixelsPhysicalSizeX(0).value().doubleValue();  % in µm
    voxel_size_y = omeMeta.getPixelsPhysicalSizeY(0).value().doubleValue();  % in µm
    voxel_size_z = omeMeta.getPixelsPhysicalSizeZ(0).value().doubleValue();  % in µm

    reader.close();


    disp('image has been read');


    se = strel3d('sphere', 3); %was 9
    somata_stack_dil = imdilate(somata_stack,se);

    neuron_stack = uint8(axon_stack)-uint8(somata_stack_dil);
    neuron_stack(neuron_stack==-1) = 0;
    neuron_stack = bwareaopen(neuron_stack, 300, 26); %eliminate fragments

    dendrite_stack = uint8(dendrite_stack)-uint8(somata_stack_dil);
    dendrite_stack(dendrite_stack==-1) = 0;
    dendrite_stack = bwareaopen(dendrite_stack, 300, 26); %eliminate fragments

    % remove somata outliers
    cluster_boundary = imclose_3d(somata_stack_dil, 200); %was 2

    cc = bwconncomp(cluster_boundary, 26); % 26-connectivity for 3D

    % find the largest object by voxel count
    numPixels = cellfun(@numel, cc.PixelIdxList);
    [~, idx] = max(numPixels);

    % create new mask with only that object
    cluster_boundary = false(size(cluster_boundary));
    cluster_boundary(cc.PixelIdxList{idx}) = true;

    %close with bigger se now that extraneous objects have been removed
    cluster_boundary = imclose_3d(cluster_boundary, 30);


    proj = any(cluster_boundary, 3);   % Y x X, logical

    %calculate cluster size with convex hull

    area_px   = bwarea(proj);

    area_um2  = area_px * voxel_size_x * voxel_size_y;


    se2 = strel('disk',30);
    %erode to make up for the large structuring element in closing
    for z = 1:size(cluster_boundary, 3)
        cluster_boundary(:,:,z) = imerode(cluster_boundary(:,:,z), se2);
    end
    cleaned_somata = cluster_boundary + somata_stack_dil;
    %cleaned_somata(find(cleaned_somata == 1)) = 0; %remove convex hull
    cleaned_neurites = cluster_boundary & neuron_stack;
    cleaned_dendrites = cluster_boundary & dendrite_stack;
    neurite_stack = cleaned_neurites;


    % determine cluster size
    numSlices = size(cleaned_somata, 3);
    hullAreas = zeros(numSlices, 1);

    for i = 1:numSlices

        [y, x] = find(cleaned_somata(:,:,i));

        if numel(x) < 3  % Convex hull requires at least 3 points
            continue;
        end


        [~, minIdx] = min(x);
        [~, maxIdx] = max(x);
        k = [minIdx, maxIdx];

        hullAreas(i) = polyarea(x(k), y(k));
    end

    % Find the image with the largest convex hull
    [size_hull, largestIdx] = max(hullAreas);



    %% Coloured neurite objects

    CC = bwconncomp(neurite_stack, 26);  % 26-connectivity for 3D
    num_components = CC.NumObjects;


    adj_matrix = [];

    [node_properties, node_positions, adj_matrix_temp] = create_network_from_stacks_nodeprops(cleaned_somata, cleaned_neurites, cleaned_somata, cleaned_dendrites, 550000);

    if isempty(adj_matrix)
        adj_matrix = zeros(size(adj_matrix_temp)); % Initialize after getting size
    end

    adj_matrix = adj_matrix + adj_matrix_temp;

    network_graph = graph(adj_matrix);
    adj_matrix(adj_matrix > 1) = 1;
    num_nodes = numnodes(network_graph);
    num_edges = numedges(network_graph);

    % rank nodes

    sizeVec = node_properties.Size(:);
    overlapVec = node_properties.NeuriteOverlap(:);


    degree = sum(adj_matrix, 2); % Degrees
    mean_degree = mean(degree);
    clust_coeff_original = mean(clustering_coef_bu(adj_matrix)); % Clustering coefficient
    distance_original = distance_bin(adj_matrix);

    char_path_original = charpath(distance_original,0,0); % Characteristic path length
    [Ci,modularity_orig] = modularity_und(adj_matrix);
    density = density_und(adj_matrix);

    spl_null_emp = log(length(adj_matrix))/log(mean_degree);
    cc_null_emp = mean_degree/length(adj_matrix);


    spl_emp_norm = char_path_original./spl_null_emp;
    cc_emp_norm = clust_coeff_original./cc_null_emp;

    %% generate random graph with same number of nodes/edges

    random_adj = randmio_und(adj_matrix,1);

    random_adj(random_adj > 1) = 1;

    degree_rand = sum(random_adj, 2); % Degree
    degree_rand_mean = mean(degree_rand);
    clust_coeff_rand = mean(clustering_coef_bu(random_adj)); % Clustering coefficient
    distance_rand = distance_bin(random_adj);
    char_path_rand = charpath(distance_rand,0,0); % Characteristic path length
    [Ci,modularity_rand] = modularity_und(random_adj);
    density_rand = density_und(random_adj);


    spl_null_rand = log(length(random_adj))/log(degree_rand_mean);
    cc_null_rand = degree_rand_mean/length(random_adj);


    spl_rand_norm = char_path_rand./spl_null_rand;
    cc_rand_norm = clust_coeff_rand./cc_null_rand;


    %create struct to store cluster data
    cnt = cnt + 1;
    disp(cnt);disp(o);
    cluster_network_data(o) = struct('cluster_id', somata_files(o).name, 'cluster_connectivity_array', adj_matrix,'size', area_um2, 'age', 'w1','degree', degree,'density', density, 'mean_degree', mean_degree, 'clustering_coefficient', clust_coeff_original, 'char_path_length', char_path_original, 'modularity', modularity_orig, 'spl_emp_norm', spl_emp_norm, 'cc_emp_norm', cc_emp_norm, 'node_positions', node_positions, 'node_properties', node_properties);
    random_network_data(o) = struct('cluster_id', somata_files(o).name, 'cluster_connectivity_array', random_adj,'size', area_um2, 'age', 'w1','degree_rand', degree_rand,'density_rand', density_rand, 'mean_degree_rand', degree_rand_mean,'clustering_coefficient', clust_coeff_rand, 'char_path_length', char_path_rand, 'modularity', modularity_rand, 'spl_rand_norm', spl_rand_norm, 'cc_rand_norm', cc_rand_norm);

    % 
    % [~, lastFolder] = fileparts(somata_folder);
    % 
    % save(['cluster_network_data_size_nodeprop_w4_lot2_' lastFolder '.mat'], 'cluster_network_data','-v7.3');
    % disp('saved');


end
end

function t = otsuThreshold(img)
% compute Otsu threshold on the non-zero distribution if zeros dominate
img1 = img(:);
img1 = img1(~isnan(img1) & ~isinf(img1));
if isempty(img1)
    t = 0;
    return;
end
imgNonZero = img1(img1>0);
if isempty(imgNonZero)
    % fallback to entire range
    imgNonZero = img1;
end
% Normalize to [0,1] for graythresh
maxv = double(max(imgNonZero));
if maxv==0
    t = 0;
    return;
end
normvals = imgNonZero / maxv;
try
    level = graythresh(normvals);
    t = level * maxv;
catch
    % in case graythresh fails, use median*0.5
    t = median(imgNonZero) * 0.5;
end
end

