function [strokePlaneVectorTop, strokePlaneVectorBottom, y_axis_projected_on_strokePlane, strokePlaneAngle] = ...
    strokePlaneVectors(input_points)

if size(input_points,2) == 1
    input_points = reshape(input_points, size(input_points,1), size(input_points,3));
end

if(sum(isnan(input_points),'all') > 0)
    n_nanpoints = char(num2str(sum(isnan(input_points),'all')/3));
    for i=1:3
        input_points(:,i) = fillmissing(input_points(:,i),"linear");
    end
    warning(['Input points contained ', n_nanpoints, ' nan points. These have been filled in using linear interpolation of neighboring, nonmissing values.'])
end

% Step 1: find the mean of the points
average = mean(input_points, 1, 'omitnan');

% Step 2: subtract the mean from all points
subtracted = bsxfun(@minus, input_points, average);

% Step 3: perform SVD
[~, ~, V] = svd(subtracted);

% Step 4: find the direction vector 
% which is the right singular vector corresponding to the largest singular value
direction = V(:, 1);

% Line is 'avg' and 'direction'
p0 = average;
d = direction;

% Parametric equation: P = p0 + t*d
get_line_point = @(t) [p0(1) + t*d(1), p0(2) + t*d(2), p0(3) + t*d(3)];

% the plane is defined by any two vectors on it, so we make two whatever
% ones now, but we want the output ones to correspond to top and bottom
% wingposition, so we will refine them to that end later
strokePlaneVector1_first = get_line_point(0);
strokePlaneVector2_first = get_line_point(1);

% get the angle between y axis and vector p projected on yz plane
angle_to_y_axis = @(p) acosd(dot(p([2 3]), sign(average(2))*[1 0])/(norm(p([2 3]))*norm([1 0])));

% get the same angle but now by giving a value for t in the parametric
% equation P = p0 + t*d
angle_to_y_axis_from_t = @(t) angle_to_y_axis(get_line_point(t));

% minimize angle between line from equation and y axis in yz plane
y_axis_t = fminsearch(angle_to_y_axis_from_t,-10);

% get y_axis projected on stroke plane
y_axis_projected_on_strokePlane = get_line_point(y_axis_t);

flap_angle_of_input_points = flapAngle(...
    strokePlaneVector1_first, ...
    strokePlaneVector2_first, ...
    y_axis_projected_on_strokePlane, ...
    input_points);

[~, index_of_max_angle] = max(flap_angle_of_input_points);
[~, index_of_min_angle] = min(flap_angle_of_input_points);

normal_of_strokePlane = cross(strokePlaneVector1_first, strokePlaneVector2_first);
normal_of_strokePlane = normal_of_strokePlane/norm(normal_of_strokePlane);

nstrpl = normal_of_strokePlane;

strokePlaneVectorTop_short = ...
    input_points(index_of_max_angle,:) - ...
    dot(input_points(index_of_max_angle,:), nstrpl) * nstrpl;

strokePlaneVectorBottom_short = ...
    input_points(index_of_min_angle,:) - ...
    dot(input_points(index_of_min_angle,:), nstrpl) * nstrpl;

angle_to_top_vector = @(p) acosd(dot(p, strokePlaneVectorTop_short)/(norm(p)*norm(strokePlaneVectorTop_short)));
angle_to_top_vector_from_t = @(t) angle_to_top_vector(get_line_point(t));
top_vector_t = fminsearch(angle_to_top_vector_from_t,-10);
strokePlaneVectorTop = get_line_point(top_vector_t);

angle_to_bottom_vector = @(p) acosd(dot(p, strokePlaneVectorBottom_short)/(norm(p)*norm(strokePlaneVectorBottom_short)));
angle_to_bottom_vector_from_t = @(t) angle_to_bottom_vector(get_line_point(t));
bottom_vector_t = fminsearch(angle_to_bottom_vector_from_t,-10);
strokePlaneVectorBottom = get_line_point(bottom_vector_t);

strokLineVector = strokePlaneVectorTop-strokePlaneVectorBottom;
strokePlaneAngleXY = acos(dot(strokLineVector([1 2]), [0 sign(strokLineVector(2))])/norm(strokLineVector([1 2])));
%strokePlaneAngleXZ = (1+sign(strokLineVector(1)))/2*pi-sign(strokLineVector(1))*acos(dot(strokLineVector([1 3]), [-1 0])/norm(strokLineVector([1 3])));
strokePlaneAngleXZ = acos(dot(strokLineVector([1 3]), [-1 0])/norm(strokLineVector([1 3])));
strokePlaneAngleYZ = acos(dot(strokLineVector([2 3]), [sign(strokLineVector(2)) 0])/norm(strokLineVector([2 3])));
strokePlaneAngle = [strokePlaneAngleXY, strokePlaneAngleXZ, strokePlaneAngleYZ];
% if(sign(strokLineVector(1))==1)
% figure; plot3([0 strokLineVector(1)],[0 strokLineVector(2)],[0 strokLineVector(3)]); view(0,0); axis equal
% disp(strokePlaneAngleXZ*57.29578)
% end



