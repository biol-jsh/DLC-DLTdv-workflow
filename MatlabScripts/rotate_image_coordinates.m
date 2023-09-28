function [rotated_x, rotated_y] = rotate_image_coordinates(img_width, img_height, x, y, rotation_angle)
    switch rotation_angle
        case 90
            rotated_x = y;
            rotated_y = img_width - x;
        case 180
            rotated_x = img_width - x;
            rotated_y = img_height - y;
        case 270
            rotated_x = img_height - y;
            rotated_y = x;
        otherwise
            error('Invalid rotation angle. Supported angles are 90, 180, and 270 degrees.');
    end
end