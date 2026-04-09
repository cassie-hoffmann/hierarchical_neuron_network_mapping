function [angular_neurites1, angular_neurites2, angular_neurites3, angular_neurites4] = isolate_angular_neurites_split(neurite_image, dilationfactor1, dilationfactor2, dilationfactor3, dilationfactor4, cleanfactor)

    bw = neurite_image;
    
    skeleton = bwmorph(bw, 'skel', Inf);
    
    branchPoints = bwmorph(skeleton, 'branchpoints');
    
    se = strel('disk', 1);
    dilatedBranchPoints = imdilate(branchPoints, se);
    
    segments = skeleton & ~dilatedBranchPoints;
    
    segments = bwareaopen(segments, 4); 

    labeledSegments = bwlabel(segments);
    
    props = regionprops(labeledSegments, 'Orientation', 'PixelIdxList');
    
    numGroups = 6;
    windowSize = 360 / numGroups;
    
    segmentGroups = cell(1, numGroups);
    seg_orientations = cell(1, numGroups);
    
    for i = 1:numGroups
        segmentGroups{i} = {};
        seg_orientations{i} = [];
    end
    
    % sort segments into angular groups
    for k = 1:length(props)
        % normalize orientation to [0, 360)
        orientation = mod(props(k).Orientation, 360);
        
        % determine the group based on the orientation, ensuring 1 to numGroups
        groupIdx = floor(orientation / windowSize) + 1;
        groupIdx = min(groupIdx, numGroups); % Cap at numGroups for boundary cases
        
        % store the segment in the appropriate group
        segmentGroups{groupIdx}{end+1}(:,1) = props(k).PixelIdxList;
        seg_orientations{groupIdx}(end+1) = mean(orientation);
    end
    
    angular_neurites_masks = cell(1, numGroups);

          for g = 1
        angular_neurites_masks{g} = false(size(bw));
        segs = false(size(bw));

        for j = 1:length(segmentGroups{g})
            segment = false(size(bw));
            segment(segmentGroups{g}{j}) = 1;
            segs(segmentGroups{g}{j}) = 1;
            angle = seg_orientations{g}(j); % get the orientation for this segment
        
            dilatedSegment = dilateInDirection(segment, angle, dilationfactor1);
            angular_neurites_masks{g} = angular_neurites_masks{g} | dilatedSegment;
        end

        % clean up small objects in each mask
        angular_neurites_masks{g} = bwareaopen(angular_neurites_masks{g}, cleanfactor);

    end

    for g = 2
        angular_neurites_masks{g} = false(size(bw));
        segs = false(size(bw));

        for j = 1:length(segmentGroups{g})
            segment = false(size(bw));
            segment(segmentGroups{g}{j}) = 1;
            segs(segmentGroups{g}{j}) = 1;
            angle = seg_orientations{g}(j); % Get the orientation for this segment
            dilatedSegment = dilateInDirection(segment, angle, dilationfactor2);
            angular_neurites_masks{g} = angular_neurites_masks{g} | dilatedSegment;
        end

    end

         for g = 5
        angular_neurites_masks{g} = false(size(bw));
        segs = false(size(bw));

        for j = 1:length(segmentGroups{g})
            segment = false(size(bw));
            segment(segmentGroups{g}{j}) = 1;
            segs(segmentGroups{g}{j}) = 1;
            angle = seg_orientations{g}(j); % Get the orientation for this segment
            dilatedSegment = dilateInDirection(segment, angle, dilationfactor3);
            angular_neurites_masks{g} = angular_neurites_masks{g} | dilatedSegment;
        end

    end

 for g = 6
        angular_neurites_masks{g} = false(size(bw));
        segs = false(size(bw));

        for j = 1:length(segmentGroups{g})
            segment = false(size(bw));
            segment(segmentGroups{g}{j}) = 1;
            segs(segmentGroups{g}{j}) = 1;
            angle = seg_orientations{g}(j); % Get the orientation for this segment
            dilatedSegment = dilateInDirection(segment, angle, dilationfactor4);
            angular_neurites_masks{g} = angular_neurites_masks{g} | dilatedSegment;
        end

    end



    angular_neurites1 = angular_neurites_masks{1}; 
    angular_neurites2 = angular_neurites_masks{2};
    angular_neurites3 = angular_neurites_masks{5};
    angular_neurites4 = angular_neurites_masks{6};
end

% function to dilate segment in a specified direction
function dilatedSegment = dilateInDirection(segment, angle, numDilations)
    se = strel('line', numDilations, angle); % create line structuring element
    dilatedSegment = imdilate(segment, se);
end







