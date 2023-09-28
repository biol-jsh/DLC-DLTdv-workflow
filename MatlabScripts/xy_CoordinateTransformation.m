function [xyz_transformed, xy_rotation_angle] = xy_CoordinateTransformation(xyz_raw, direction_vector_or_angle)

vector_or_angle = "not sure";
if(numel(direction_vector_or_angle(:))==size(xyz_raw,1))
    vector_or_angle = "angle";
else
    vector_or_angle = "vector";
end

switch vector_or_angle
    case "angle"
        xy_rotation_angle = direction_vector_or_angle;
    case "vector"
        for i=1:size(xyz_raw,1)
            xy_rotation_angle(i) = sign(direction_vector_or_angle(i,2)) * ...
                acos(dot(direction_vector_or_angle(i,[1 2]), [1 0]) / ...
                (norm(direction_vector_or_angle(i,[1 2]))));
        end
end

for i=1:size(xyz_raw,1)
    
    xy_rotationMatrix = [...
        cos(xy_rotation_angle(i)) sin(xy_rotation_angle(i)) 0; ...
        -sin(xy_rotation_angle(i)) cos(xy_rotation_angle(i)) 0; ...
        0 0 1];
    
    if ndims(xyz_raw) == 3
        for j=1:size(xyz_raw,2)
            xyz_transformed(i,j,:) = reshape( ...
                xy_rotationMatrix * ...
                reshape(xyz_raw(i,j,:), 3, 1),size(xyz_raw(i,j,:)));
        end
    else
        xyz_transformed(i,:) = xy_rotationMatrix*xyz_raw(i,:)';
    end
end