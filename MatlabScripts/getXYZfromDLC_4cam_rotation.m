function [xyzDLC, DLTdvResidual, meanLikelihoodDLC] = ..., [vidWidth1, vidHeight1] , , ,
    getXYZfromDLC_4cam_rotation(easyWandDataFile, camTformsFiles, DLCoutputCSVs, ThruTrackerTranslationFiles, crp, minDLClikelihood, pairs, bodyparts, rotations, vidSizes)
%getXYZfromDLC_4cam_rotation: This function calculates the 3D coordinates
% from 2D DLC outputs from four cameras. The function also computes the 
% DLTdv residuals and, per 3D point prediction, mean likelihood associated
% with the DeepLabCut (DLC) predictions.
%
% This function is similar to getXYZfromDLC_4cam, only difference is that
% this version contains a correction for different sensor rotations used
% during easyWand calibration and DLC analysis in our test case
%
% Author:
%   Jonas Bengt Carina Håkansson, PhD
%   Univeristy of Colorado Colorado Springs (UCCS)
%   Last Modified: 20230908
%   Contact: scientistjonas [at] gmail [dot] com
%   Alt contact: jhaakans [at] uccs [dot] edu
%
% Citation: Please cite (TODO add citation here)
%
% Inputs:
%   easyWandDataFile (string): The path to the file containing the easyWand
%   calibration file
%
%   camTformsFiles (array of strings): The paths to the files containing
%   the transformation matrices for each camera.
%
%   DLCoutputCSVs (array of strings): The paths to the CSV files containing
%   the output from the DeepLabCut analysis.
%
%   ThruTrackerTranslationFiles (cell array of strings): The paths to the
%   files containing translation data from ThruTracker. Replace with [] if
%   no cropping used.
%
%   crp (number): Camera Reduction Penalty, the penalty for using less than
%   all cameras for a 3D prediction. Used to filter out faulty 2D
%   predictions based on the resulting DLTdv residual.
%
%   minDLClikelihood (double): The minimum likelihood threshold for
%   considering a DLC prediction as valid.
%
%   pairs (integer array): A cell array containing pairs of body parts that
%   need to be considered for symmetrical points' corrections.
%
%   bodyparts (integer array): An array containing the label numbers of
%   body parts to be tracked. For paired landmark, only the left landmark 
%   is supplied
%
%   rotations (array): An array specifying the 2D rotations to be applied
%   to the image data from each camera. Valid values are 90, 180, and 270
%
%   vidSizes (cell array of integer arrays): An array specifying the sizes
%   of the videos. Each element is a 1x2 array containint width and height
%
% Outputs:
%   xyzDLC (NxMxL matrix): A matrix containing the 3D coordinates
%   calculated from the 2D DLC outputs. N = frame number, M = landmark
%   number, L = spatial dimnesion (1 = x, 2 = y, 3 = z)
%   
%   DLTdvResidual (N x M matrix): A matrix containing the residuals 
%   associated with the 3D coordinates calculation.
%   
%   meanLikelihoodDLC (N x M matrix): A vector containing the mean likelihoods associated with the DLC predictions for each body part.
%
% Example usage:
%   [xyzDLC, DLTdvResidual, meanLikelihoodDLC] = getXYZfromDLC_4cam_rotation("path/to/easyWandDataFile", ["path/to/cam1", "path/to/cam2", ...], ["path/to/csv1", "path/to/csv2", ...], ["path/to/trans1", "path/to/trans2", ...], 2, 0.3, [1, 15; 2, 13; 3, 14; 4, 12; 5, 10;, 6, 11], [1, 2, 3, 4, 5, 6, 7, 8, 9, 16], [0, 180, 180, 180], {[vidWidth1, vidHeight1], [vidWidth2, vidHeight2], ...});
%
% Note:
%   Ensure all paths are correctly specified.
%   Adjust the minDLClikelihood as necessary based on your data quality.
%   Customize the pairs and bodyparts cell arrays based on the specifics of your study.

    load(easyWandDataFile, "easyWandData");
    coefs = easyWandData.coefs;

    % Define possible camera combinations for 3D predictions (e.g. cam1 and 2, or cam1, 2, 4)
    [possible_combinations1,possible_combinations2,possible_combinations3, possible_combinations4] = ndgrid(0:2); %0,1,2 = none,left,right
    possible_combinations=[possible_combinations1(:), possible_combinations2(:),possible_combinations3(:) possible_combinations4(:)];

    % if more than two cams don't use a point we can't get reconstruction, so get rid of that combination
    possible_combinations((sum(possible_combinations==0,2)>2),:)=[];

    ncams = numel(DLCoutputCSVs); % number of cameras
    [possible_combinations1,possible_combinations2,possible_combinations3, possible_combinations4] = ndgrid(0:1); %0,1 = use, don't use
    possible_bilateral_combinations=[possible_combinations1(:), possible_combinations2(:),possible_combinations3(:) possible_combinations4(:)];

    possible_bilateral_combinations((sum(possible_bilateral_combinations==0,2)>2),:)=[];

midLineBodyparats = bodyparts(~ismember(bodyparts, pairs(:,1)));

if(~isempty(ThruTrackerTranslationFiles))
    [~, confidence_dlc] = getXYfromDLC(DLCoutputCSVs,ThruTrackerTranslationFiles);
else
    [~, confidence_dlc] = getXYfromDLC(DLCoutputCSVs);
end

for cam = 1:ncams
    rotation = rotations(cam);
    vidSize = vidSizes{cam};
    xyCoords{cam} = readtable(DLCoutputCSVs(cam), 'Delimiter', ',');
    xcoordsTEMP = xyCoords{cam}(:,2:3:end);
    xcoordsTEMP = xcoordsTEMP{:,:}; %table to matrix

    ycoordsTEMP = xyCoords{cam}(:,3:3:end);
    ycoordsTEMP = ycoordsTEMP{:,:}; %table to matrix

    if rotation ~= 0 % rotate image coordinates to align with orientation during calibration
        [xcoordsTEMP, ycoordsTEMP] = rotate_image_coordinates(vidSize(1), vidSize(2), xcoordsTEMP, ycoordsTEMP, rotation);
    end

    likelihoodTEMP = xyCoords{cam}(:,4:3:end); %DLC confidences for this camera
    likelihoodTEMP = likelihoodTEMP{:,:}; %table to matrix
    
    if(~isempty(ThruTrackerTranslationFiles))
        load(ThruTrackerTranslationFiles(cam));
    end

    for bodyPoint=1:size(xcoordsTEMP, 2)
        for frameNumber=1:size(xyCoords{cam},1)
            xcoordsDLC(frameNumber,bodyPoint,cam) = xcoordsTEMP(frameNumber,bodyPoint);
            ycoordsDLC(frameNumber,bodyPoint,cam) = ycoordsTEMP(frameNumber,bodyPoint);
            likelihoodDLC(frameNumber,bodyPoint,cam) = likelihoodTEMP(frameNumber,bodyPoint);
        end
        if(~isempty(ThruTrackerTranslationFiles))
            xcoordsDLC(:,bodyPoint,cam) = xcoordsDLC(:,bodyPoint,cam) + XYtop_left(:,1);
            ycoordsDLC(:,bodyPoint,cam) = ycoordsDLC(:,bodyPoint,cam) + XYtop_left(:,2);
        end
    end

    load(camTformsFiles(cam), 'camud');

    image_xvaluesTEMP = xcoordsDLC(:,:,cam);
    image_xvaluesTEMP = image_xvaluesTEMP(:);

    image_yvaluesTEMP = ycoordsDLC(:,:,cam);
    image_yvaluesTEMP = image_yvaluesTEMP(:);

    image_xy_Undistorted = undistort_Tform([image_xvaluesTEMP image_yvaluesTEMP],camud);
    image_xvalues_UndistortedTEMP = image_xy_Undistorted(:,1);
    image_yvalues_UndistortedTEMP = image_xy_Undistorted(:,2);

    image_xvalues_Undistorted(:,:,cam) = reshape(image_xvalues_UndistortedTEMP, size(xcoordsDLC,1), size(xcoordsDLC,2));
    image_yvalues_Undistorted(:,:,cam) = reshape(image_yvalues_UndistortedTEMP, size(xcoordsDLC,1), size(xcoordsDLC,2));
end

nframes = size(confidence_dlc,1);

for frame = 1:nframes
    for bodypart = bodyparts
        if ~ismember(bodypart,midLineBodyparats)
            points = [midLineBodyparats(1) pairs(bodypart,:)]; %we need a point that exists but doesn't make sense, so a midline bodypart works here
            residuals = nan(size(possible_combinations,1),1);
            likely_sides = nan(size(possible_combinations,1),1);
            containsNAN = zeros(size(possible_combinations,1),1);
            containsOtherCamera = zeros(size(possible_combinations,1),1);
            actual_camcombo = nan(size(possible_combinations));

            for camcombo = 1:size(possible_combinations,1)
                current_combination = possible_combinations(camcombo,:);
                current_bodyPoints = points(current_combination+1);
                current_confidences(1) = confidence_dlc(frame,current_bodyPoints(1),1);
                current_confidences(2) = confidence_dlc(frame,current_bodyPoints(2),2);
                current_confidences(3) = confidence_dlc(frame,current_bodyPoints(3),3);
                current_confidences(4) = confidence_dlc(frame,current_bodyPoints(4),4);

                use_or_not = ones(1,4);
                use_or_not(current_bodyPoints==points(1))=nan;
                use_or_not(current_confidences<minDLClikelihood) = nan;
                containsNAN(camcombo)=sum(isnan(use_or_not));

                if(sum(use_or_not,'omitnan')<2)
                    %only one camera left                    
                    continue
                end

                left_likelihood_sum = sum(current_confidences(current_combination == 1 & ~isnan(use_or_not)));
                right_likelihood_sum = sum(current_confidences(current_combination == 2 & ~isnan(use_or_not)));
                if left_likelihood_sum > right_likelihood_sum
                    likely_sides(camcombo) = 1;
                elseif left_likelihood_sum < right_likelihood_sum
                    likely_sides(camcombo) = 2;
                end

                [~, residuals(camcombo)] = dlt_reconstruct(coefs, ...
                    [...
                    use_or_not(1)*image_xvalues_Undistorted(frame,current_bodyPoints(1),1) use_or_not(1)*image_yvalues_Undistorted(frame,current_bodyPoints(1),1) ...
                    use_or_not(2)*image_xvalues_Undistorted(frame,current_bodyPoints(2),2) use_or_not(2)*image_yvalues_Undistorted(frame,current_bodyPoints(2),2) ...
                    use_or_not(3)*image_xvalues_Undistorted(frame,current_bodyPoints(3),3) use_or_not(3)*image_yvalues_Undistorted(frame,current_bodyPoints(3),3) ...
                    use_or_not(4)*image_xvalues_Undistorted(frame,current_bodyPoints(4),4) use_or_not(4)*image_yvalues_Undistorted(frame,current_bodyPoints(4),4) ...
                    ]);
            end

            residuals(find(containsNAN)) = residuals(find(containsNAN))+containsNAN(find(containsNAN))*crp;
            
            if(isnan(min(residuals)))
                likely_side = 0;
            else
            [~, best_combination_index] = min(residuals);
            best_combination = possible_combinations(best_combination_index,:);
            likely_side = likely_sides(best_combination_index);
            end

            impossible_match = zeros(size(possible_combinations,1),1);
            impossible_match = impossible_match + likely_side == likely_sides;
            for i = 1:size(possible_combinations,2)
                if best_combination(i) == 0
                    continue
                end

                impossible_match = impossible_match + (best_combination(i) == actual_camcombo(:,i));
            end

            impossible_match(impossible_match>0)=1;
            impossible_match(isnan(residuals))=1;

            if(sum(impossible_match==0)==0)
                likely_side_partner = 0;
            else
                possible_partner_combinations = possible_combinations;
                partner_residuals = residuals;
                likelys_sides_partner = likely_sides;
                possible_partner_combinations(find(impossible_match),:) = [];
                partner_residuals(find(impossible_match)) = [];
                likelys_sides_partner(find(impossible_match)) = [];
                [~, best_partner_combination_index] = min(partner_residuals);
                best_partner_combination = possible_partner_combinations(best_partner_combination_index,:);
                likely_side_partner = likelys_sides_partner(best_partner_combination_index);
            end

            for lr = 1:2 %main partner
                if lr==1
                    if(likely_side==0)
                        continue
                    end
                    current_bodyPoints = points(best_combination+1);
                    current_outPut_point = points(likely_side+1);
                else
                    if(likely_side_partner==0)
                        continue
                    end
                    current_bodyPoints = points(best_partner_combination+1);
                    current_outPut_point = points(likely_side_partner+1);
                end

                if isempty(current_bodyPoints)
                    continue
                end

                current_confidences(1) = confidence_dlc(frame,current_bodyPoints(1),1);
                current_confidences(2) = confidence_dlc(frame,current_bodyPoints(2),2);
                current_confidences(3) = confidence_dlc(frame,current_bodyPoints(3),3);
                current_confidences(4) = confidence_dlc(frame,current_bodyPoints(4),4);

                use_or_not = ones(1,4);
                use_or_not(current_bodyPoints==points(1))=nan;
                use_or_not(current_confidences<minDLClikelihood) = nan;


                clear output_dlc_confidenceTEMP
                for i=1:ncams
                    if current_bodyPoints(i) == points(1) || isnan(use_or_not(i))
                        output_dlc_confidenceTEMP(i) = nan;
                    else
                        output_dlc_confidenceTEMP(i)=confidence_dlc(frame,current_bodyPoints(i),i);
                    end
                end

                meanLikelihoodDLC(frame,current_outPut_point) = mean(output_dlc_confidenceTEMP,'omitnan');

                [xyzDLC(frame,current_outPut_point,:), DLTdvResidual(frame,current_outPut_point)] = dlt_reconstruct(coefs, ...
                    [...
                    use_or_not(1)*image_xvalues_Undistorted(frame,current_bodyPoints(1),1) use_or_not(1)*image_yvalues_Undistorted(frame,current_bodyPoints(1),1) ...
                    use_or_not(2)*image_xvalues_Undistorted(frame,current_bodyPoints(2),2) use_or_not(2)*image_yvalues_Undistorted(frame,current_bodyPoints(2),2) ...
                    use_or_not(3)*image_xvalues_Undistorted(frame,current_bodyPoints(3),3) use_or_not(3)*image_yvalues_Undistorted(frame,current_bodyPoints(3),3) ...
                    use_or_not(4)*image_xvalues_Undistorted(frame,current_bodyPoints(4),4) use_or_not(4)*image_yvalues_Undistorted(frame,current_bodyPoints(4),4) ...
                    ]);
            end

        else
            current_outPut_point = bodypart;
            points = [pairs(1) bodypart];
            residuals = nan(size(possible_bilateral_combinations,1),1);
            containsNAN = zeros(size(possible_bilateral_combinations,1),1);

            for camcombo = 1:size(possible_bilateral_combinations,1)
                current_combination = possible_bilateral_combinations(camcombo,:);
                current_bodyPoints = points(current_combination+1);
                current_confidences(1) = confidence_dlc(frame,current_bodyPoints(1),1);
                current_confidences(2) = confidence_dlc(frame,current_bodyPoints(2),2);
                current_confidences(3) = confidence_dlc(frame,current_bodyPoints(3),3);
                current_confidences(4) = confidence_dlc(frame,current_bodyPoints(4),4);
                

                use_or_not = [1,1,1,1];
                use_or_not(current_bodyPoints==points(1))=nan;
                use_or_not(current_confidences<minDLClikelihood) = nan;
                containsNAN(camcombo)=sum(isnan(use_or_not))>0;

                [~, residuals(camcombo)] = dlt_reconstruct(coefs, ...
                    [...
                    use_or_not(1)*image_xvalues_Undistorted(frame,current_bodyPoints(1),1) use_or_not(1)*image_yvalues_Undistorted(frame,current_bodyPoints(1),1) ...
                    use_or_not(2)*image_xvalues_Undistorted(frame,current_bodyPoints(2),2) use_or_not(2)*image_yvalues_Undistorted(frame,current_bodyPoints(2),2) ...
                    use_or_not(3)*image_xvalues_Undistorted(frame,current_bodyPoints(3),3) use_or_not(3)*image_yvalues_Undistorted(frame,current_bodyPoints(3),3) ...
                    use_or_not(4)*image_xvalues_Undistorted(frame,current_bodyPoints(4),4) use_or_not(4)*image_yvalues_Undistorted(frame,current_bodyPoints(4),4) ...
                    ]);

            end

            %residuals(find(containsNAN)) = residuals(find(containsNAN))+crp; %penalize using less than all cameras
            residuals(find(containsNAN)) = residuals(find(containsNAN))+containsNAN(find(containsNAN))*crp;

            [~, best_combination_index] = min(residuals);
            best_combination = possible_bilateral_combinations(best_combination_index,:);
            current_bodyPoints = points(best_combination+1);
            use_or_not = [1,1,1,1];
            use_or_not(current_bodyPoints==points(1))=nan;
            use_or_not(current_confidences<minDLClikelihood) = nan;

            clear output_dlc_confidenceTEMP
            for i=1:numel(current_bodyPoints)
                if current_bodyPoints(i) == points(1)
                    output_dlc_confidenceTEMP(i) = nan;
                else
                    output_dlc_confidenceTEMP(i)=confidence_dlc(frame,current_bodyPoints(i),i);
                end
            end

            meanLikelihoodDLC(frame,current_outPut_point) = mean(output_dlc_confidenceTEMP,'omitnan');

            [xyzDLC(frame,current_outPut_point,:), DLTdvResidual(frame,current_outPut_point)] = dlt_reconstruct(coefs, ...
                [...
                use_or_not(1)*image_xvalues_Undistorted(frame,current_bodyPoints(1),1) use_or_not(1)*image_yvalues_Undistorted(frame,current_bodyPoints(1),1) ...
                use_or_not(2)*image_xvalues_Undistorted(frame,current_bodyPoints(2),2) use_or_not(2)*image_yvalues_Undistorted(frame,current_bodyPoints(2),2) ...
                use_or_not(3)*image_xvalues_Undistorted(frame,current_bodyPoints(3),3) use_or_not(3)*image_yvalues_Undistorted(frame,current_bodyPoints(3),3) ...
                use_or_not(4)*image_xvalues_Undistorted(frame,current_bodyPoints(4),4) use_or_not(4)*image_yvalues_Undistorted(frame,current_bodyPoints(4),4) ...
                ]);
            
        end
    end
end
