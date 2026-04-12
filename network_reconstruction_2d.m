function network_reconstruction_2d

    clear all; close all;
    %% user input

    input_folder = "/path/to/input/folder";
    output_folder = "/path/to/output/folder";

    av_soma_size = x; %average size of individual soma in pixels
    proportions = [1]; %replace with proportion for inter-cluster connectivity ranging from 0.1 to 1
    filetype = '*.tiff'; %file type of input images
    som_thresh_min = x; %minimum pixel value of somata objects
    som_thresh_max = x; %maximum pixel value of somata objects
    neur_thresh_min = x; %minimum pixel value of neurite objects
    neur_thresh_max = x; %maximum pixel value of neurite objects

    som_filter = x; % lower size threshold (pixels) for valid somata objects
    neurite_filter = x; % lower size threshold (pixels) for valid neurite objects
    ang_dil_90_130 = x; % angular neurite dilation - 90-130 degrees (upper right orientation)
    ang_dil_131_180 = x; % angular neurite dilation - 131-180 degrees (right orientation)
    ang_dil_70_89 = x;% angular neurite dilation - 70-89 degrees (upper left orientation)
    ang_dil_0_69 = x;% angular neurite dilation - 0-69 degrees (left orientation)
    ang_dil_filter = x; % secondary neurite size filter for each orientation mask
    visualize = 1; % 1 = visualise network, 0 = do not visualize network

    %% end of user input    
    
    %% Store files in cell array
    
    imageFiles = dir(fullfile(input_folder,filetype));
    
    % filter out files starting with "._"
    validFiles = imageFiles(~startsWith({imageFiles.name}, '._'));
    
    imagesData = cell(length(validFiles), 1);
    
    for f = 1:length(validFiles)
        imagePath = fullfile(input_folder, validFiles(f).name);
        imagesData{f} = imread(imagePath); % Read the image
    end
    
    metricsData = cell(length(imageFiles), 1);
    
    %% Calculate metrics
    
    for time = 1:length(imagesData)
        disp(time);
        node_point_coords = [];
    
        image = double(imagesData{time});
        img_som_dil = image >= som_thresh_min & image <= som_thresh_max;
        img_neur_dil = image >= neur_thresh_min & image <= neur_thresh_max;
    
        img_som_dil = bwareaopen(img_som_dil, som_filter);
        img_neur_dil = bwareaopen(img_neur_dil, neurite_filter);
        cluster_ct = 0;
    
        %parse neurite mask into angular orientation masks
        [image_b, image_g, image_p, image_r] = isolate_angular_neurites_split(img_neur_dil, ang_dil_131_180,ang_dil_90_130,ang_dil_70_89,ang_dil_0_69, ang_dil_filter);
        image_b = bwmorph(image_b, "skel", Inf); image_g = bwmorph(image_g, "skel", Inf); image_p = bwmorph(image_p, "skel", Inf);image_r = bwmorph(image_r, "skel", Inf);
    
        N=[size(image_b,1),size(image_b,2)];
        se=strel('disk',1);
    
        cc=bwconncomp(image_b);
        KK=length(cc.PixelIdxList);
        cnt=0;
        img_sep=zeros(N);
    
        for i=1:KK
            cnt=cnt+1;
            a{cnt}=cc.PixelIdxList{i};
            img_sep(a{cnt})=cnt;
        end
    
        cc=bwconncomp(image_g);
        KK=length(cc.PixelIdxList);
    
        for i=1:KK
            cnt=cnt+1;
            a{cnt}=cc.PixelIdxList{i};
            img_sep(a{cnt})=cnt;
        end
    
        cc=bwconncomp(image_p);
        KK=length(cc.PixelIdxList);
    
        for i=1:KK
            cnt=cnt+1;
            a{cnt}=cc.PixelIdxList{i};
            img_sep(a{cnt})=cnt;
        end
    
        cc=bwconncomp(image_r);
        KK=length(cc.PixelIdxList);
    
        for i=1:KK
            cnt=cnt+1;
            a{cnt}=cc.PixelIdxList{i};
            img_sep(a{cnt})=cnt;
        end
    
        C=cnt;
        SubdivideNodes=1;
        ct=0;
    
        ind=find(img_som_dil==1);
        nodes=zeros(size(img_som_dil));
        nodes(ind)=1;          %set nodes to a value of 1
        nodes=medfilt2(nodes);  %remove speckles in image
        N=size(nodes);
        cc=bwconncomp(nodes);  %give nodes a distinct label
        nNodes=cc.NumObjects;
        adj=zeros(nNodes,nNodes); %initialize connectivity matrix
    
        for i=1:nNodes
            sz(i)=length(cc.PixelIdxList{i});
        end
    
        %% Establish intra-cluster connectivity
    
        master_adj_proportions = cell(size(proportions));
        cluster_intracluster_sizes = [];
        x_means_clusters = {[]}; y_means_clusters = {[]};
    
        for p = 1:length(proportions)
            ct=0;
            for i=1:nNodes
                if sz(i)>av_soma_size && SubdivideNodes
                    cluster_ct = cluster_ct+1;
                    [row,col]=ind2sub(size(nodes),cc.PixelIdxList{i});
                    number_nodes=ceil(sz(i)/av_soma_size);
                    idx=kmeans([row,col],number_nodes);
                    first_ct = ct+1;
                    x_means_temp = 0; y_means_temp = 0;
                    for j=1:number_nodes
                        ct=ct+1;
                        nodes(cc.PixelIdxList{i}(idx==j))=ct;
                        y_means_temp(j,:) = mean(row(idx == j));
                        x_means_temp(j,:) = mean(col(idx == j));
                    end

                    x_means_clusters{cluster_ct,:} = x_means_temp; % store x-means for this node
                    y_means_clusters{cluster_ct,:} = y_means_temp; % store y-means for this node
                    intracluster_adj_size = ct-first_ct+1;
                    intracluster_adj_size = ct - first_ct + 1;
                    n = intracluster_adj_size;

                    adj(first_ct:ct, first_ct:ct) = 0;
                    num_possible_edges = n * (n - 1) / 2;
                    num_edges = round(proportions * num_possible_edges);
                    upper_idx = find(triu(ones(n), 1));
                    selected_edges = upper_idx(randperm(length(upper_idx), num_edges));
                    temp_block = zeros(n);
                    temp_block(selected_edges) = 1;
                    temp_block = temp_block + temp_block';
                    adj(first_ct:ct, first_ct:ct) = temp_block;

                else

                    ct=ct+1;
                    nodes(cc.PixelIdxList{i})=ct;
                end
                adj(logical(eye(size(adj)))) = 0;
                master_adj_proportions{p}=adj; %just clusters

            end
    
            fprintf('%d nodes, subdivided to %d nodes\n',nNodes,ct);
        end
    
    
    
        for s = 1:length(master_adj_proportions)
            adj = master_adj_proportions{s};
            img_nodes_small=nodes(1:N(1),1:N(2));
            img_small = zeros(size(img_nodes_small));
    
            img=zeros(N(1),N(2));      %image to hold axons
            img_tmp=zeros(N(1),N(2));  %image to hold a target axon temporarily
            se = strel('disk',3);
            cnt=0;
            adj2=adj;
    
            for i=1:length(a) %loop over axons
                img_tmp(a{i}(:,:))=1;
                img_tmp=imdilate(img_tmp,se); %dilate the axon slightly
                img(~~img_tmp)=cnt;
                %calculate matrix
                img_small=img_tmp(1:N(1),1:N(2));
                ind=find(img_small);
                vals=setdiff(unique(img_nodes_small(ind)),0);
                for n1=1:length(vals)
                    for n2=1:length(vals)
                        adj2(vals(n1),vals(n2))=1;
                    end
                end
                img_tmp=zeros(N(1),N(2));
            end
            adj2(logical(eye(size(adj2)))) = 0;
            master_adj_proportions{s} = adj2;
        end
    
    
        adj_size = size(adj);
        adj2_size = size(adj2);
        adj(adj2_size(1), adj2_size(2)) = 0;
        mask = adj == 0;
        adj(mask) = adj2(mask);
    
        %% Assign neighbour input node to terminal node in cluster
    
        adj3 = zeros(size(adj2));
        array1=[];
    
        for array = 1:min(size(adj, 1), size(adj2, 1))
            if any(adj(array, :) ~= 0) && any(adj2(array, :) ~= 0)
                % Find non-zero columns in adj2 for the current row
                root_nodes = find(adj2(:, array) ~= 0);
                [rows, cols] = ind2sub(size(adj2(:, array)), root_nodes);
                input_node = array(cols(1));
                cluster_number = adj(root_nodes(1), input_node);
                array1 = find(adj == cluster_number);
            end
        end
    
        for o = array1'
            [row_idx, col_idx] = ind2sub(size(adj), o);
            adj3(row_idx, input_node) = 1;
        end
    
        node_point_coords(ct,1) =  mean(cols);
        node_point_coords(ct,2) =  mean(rows);
    
        %% Visualize network graphs
    
        if visualize == 1
            for p = 1:length(master_adj_proportions)
    
                adj = master_adj_proportions{p};
    
                x = zeros(size(adj,1),1);
                y = zeros(size(adj,1),1);
    
    
                for w = 1:size(adj,1)
    
                    [row,col] = find(nodes == w);
    
                    if ~isempty(row)
                        x(w) = mean(col);
                        y(w) = mean(row);
                    else
                        x(w) = NaN;
                        y(w) = NaN;
                    end
    
                end
    
                node_point_coords(:,1) = y;
                node_point_coords(:,2) = x;
    
                hf = figure;
                hf.Color = 'w';
    
                imagesc(nodes);
                colormap colorcube;
                hold on;
                axis off;
    
                for i = 1:size(adj,1)
                    for j = i+1:size(adj,1)
                        if adj(i,j) && ~isnan(x(i)) && ~isnan(x(j))
                            line([node_point_coords(i,2), node_point_coords(j,2)], [node_point_coords(i,1),  node_point_coords(j,1)], ...
                                'Color','red','LineWidth',0.5);
                        end
                    end
                end
    
                for i = 1:size(adj,1)
                    if ~isnan(x(i))
                        text(node_point_coords(i,2), node_point_coords(i,1), num2str(i), 'Color','w');
                    end
                end
    
            end
        end
    
        %% Graph analysis
    
        [Ci,Q_overall]=modularity_und(adj,1);
        number_modules=length(unique(Ci));
        number_of_nodes = length(adj);
        dist_matrix = zeros(size(adj));
        dist_matrix = pdist2(node_point_coords, node_point_coords);
        dist_matrix(adj==0) = 0; %keep only connected distances
        mean_edge_distance = mean(dist_matrix(triu(adj > 0)));
        degree_empirical = degrees_und(adj);
        mean_degree = mean(mean(degree_empirical));
    
        %community/module detection
        n = length(Ci);
        within_mask = false(n);
        between_mask = false(n);
    
        for m = unique(Ci(:))'
            ix = (Ci == m);
            within_mask(ix, ix) = true; %fill within-module blocks
        end
    
        adj_intra = within_mask  & (adj > 0);
    
        [network_density,node_num,edge_num] = density_und(adj);
        [module_density,node_num_intra,edge_num_intra] = density_und(adj_intra);
    
    
        empirical_metrics(time,:)=table(number_of_nodes, mean_edge_distance, network_density, module_density, mean_degree, number_modules);
    
    
        parts = strsplit(input_folder, filesep);
        tag = parts{end};
        fname1 = fullfile(output_folder, ['network_metrics_' tag '.mat']);
    
         save(fname1, 'empirical_metrics', '-v7.3');
    
    
    end

end
