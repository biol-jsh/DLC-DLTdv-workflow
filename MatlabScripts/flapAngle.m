function flapAngle = flapAngle(...
    strokePlane_vector1, ...
    strokePlane_vector2, ...
    y_axis_projected_on_strokePlane, ...
    input_points)
%FLAPANGLE Summary of this function goes here
%   Detailed explanation goes here
normal_of_strokePlane = cross(strokePlane_vector1, strokePlane_vector2);
normal_of_strokePlane = normal_of_strokePlane/norm(normal_of_strokePlane);

ystrpl = y_axis_projected_on_strokePlane;

for i=1:size(input_points,1)
    input_point_projected_on_strokePlane = ...
        input_points(i,:) - ...
        dot(input_points(i,:), normal_of_strokePlane) * normal_of_strokePlane;
    
    
    ipstrpl = input_point_projected_on_strokePlane;    
    
    flapAngle(i) = sign(ipstrpl(3)-ystrpl(3)) * ...
        acosd(dot(ipstrpl, ystrpl)/(norm(ipstrpl)*norm(ystrpl)));
end
end
