function [xyzDLC, DLTdvResidual, meanLikelihoodDLC] = ...
    getXYZ_swapFix(easyWandDataFile, camTformsFiles, DLCoutputCSVs, ThruTrackerTranslationFiles, crp, minDLClikelihood)
load(easyWandDataFile);
coefs = easyWandData.coefs;
pairs  = [... % left and right point for...
    1,15; ... % wingtips
    2,13; ... % wrists
    3,14; ... % 5th digits
    4,12; ... % elbows
    5,10; ... % shoulders
    6,11  ... % ankles
    ];

ncams = numel(DLCoutputCSVs);
[possible_combinations2,possible_combinations3,possible_combinations1] = meshgrid(0:2); %0,1,2 = none,left,right
possible_combinations=[possible_combinations1(:), possible_combinations2(:),possible_combinations3(:)];

% if more than one cam doesn't use a point we can't get reconstruction, so get rid of that combination
possible_combinations((sum(possible_combinations==0,2)>1),:)=[];

possible_bilateral_combinations = [...
    1,1,1; ...
    1,1,0; ...
    1,0,1; ...
    0,1,1; ...
    ];

[~, confidence_dlc] = getXYfromDLC(DLCoutputCSVs,ThruTrackerTranslationFiles);

for cam = 1:ncams
    xyCoords{cam} = readtable(DLCoutputCSVs(cam), 'Delimiter', ',');
    xcoordsTEMP = xyCoords{cam}(:,2:3:end);
    xcoordsTEMP = xcoordsTEMP{:,:}; %table to matrix

    ycoordsTEMP = xyCoords{cam}(:,3:3:end);
    ycoordsTEMP = ycoordsTEMP{:,:}; %table to matrix

    likelihoodTEMP = xyCoords{cam}(:,4:3:end);
    likelihoodTEMP = likelihoodTEMP{:,:}; %table to matrix

    load(ThruTrackerTranslationFiles(cam));

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
%         if ismember(bodyPoint, [8])
% 
%             xcoordsDLC(:,bodyPoint,cam) = tybutterNaN(xcoordsDLC(:,bodyPoint,cam), 45, 800, 'low');
%             ycoordsDLC(:,bodyPoint,cam) = tybutterNaN(ycoordsDLC(:,bodyPoint,cam), 45, 800, 'low');
% 
%         end
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
    %for bodypart = [1:9 16]
    for bodypart = [1:9 16]%[1:16]
        if bodypart<7 % symmetrical point
            badPoint = bodypart+1;
            points = [badPoint pairs(bodypart,:)]; %we need a point that exists but doesn't make sense
            residuals = nan(size(possible_combinations,1),1);
            likely_sides = nan(size(possible_combinations,1),1);
            likely_sides(find(sum(possible_combinations==1,2)>1))=1;
            likely_sides(find(sum(possible_combinations==2,2)>1))=2;
            containsNAN = zeros(size(possible_combinations,1),1);

            for camcombo = 1:size(possible_combinations,1)
                current_combination = possible_combinations(camcombo,:);
                current_bodyPoints = points(current_combination+1);
                current_confidences(1) = confidence_dlc(frame,current_bodyPoints(1),1);
                current_confidences(2) = confidence_dlc(frame,current_bodyPoints(2),2);
                current_confidences(3) = confidence_dlc(frame,current_bodyPoints(3),3);
                current_confidences = current_confidences';

                use_or_not = [1,1,1];
                use_or_not(current_bodyPoints==points(1))=nan;
                use_or_not(current_confidences<minDLClikelihood) = nan;
                containsNAN(camcombo)=sum(isnan(use_or_not))>0;

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
                    ]);

            end

            residuals(find(containsNAN)) = residuals(find(containsNAN))+crp; %penalize using less than all cameras

            if(isnan(min(residuals)))
                likely_side = 0;
            else

            [~, best_combination_index] = min(residuals);
            best_combination =  possible_combinations(best_combination_index,:);
            likely_side = likely_sides(best_combination_index);
            end

            impossible_match = zeros(size(possible_combinations,1),1);
            impossible_match = impossible_match + likely_side == likely_sides;
            for i = 1:size(possible_combinations,2)
                if best_combination(i) == 0
                    continue
                end
                impossible_match = impossible_match + (best_combination(i) == possible_combinations(:,i));
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


%             if isnan(likely_side) && isnan(likely_side_partner)
%                 %sides = randperm(2);
%                 sides = [1,2];
%                 main_left_confidence = confidence_dlc(frame,points(2), find(best_combination==1));
%                 main_right_confidence = confidence_dlc(frame,points(3), find(best_combination==2));
% 
%                 main_left_confidence_mean = main_left_confidence + 1 - main_right_confidence;
%                 main_right_confidence_mean = main_right_confidence + 1 - main_left_confidence;
%                 
%                 partner_left_confidence = confidence_dlc(frame,points(2), find(best_partner_combination==1));
%                 partner_right_confidence = confidence_dlc(frame,points(3), find(best_partner_combination==2));
% 
%                 partner_left_confidence_mean = partner_left_confidence + 1 - partner_right_confidence;
%                 partner_right_confidence_mean = partner_right_confidence + 1 - partner_left_confidence;
% 
%                 %[~, I] = max([main_left_confidence_mean, main_right_confidence_mean, partner_left_confidence_mean, partner_right_confidence_mean]);
%                 [~, I] = max([main_left_confidence, main_right_confidence, partner_left_confidence, partner_right_confidence]);
%                 if ismember(I, [1,4])
%                     likely_side = 1;
%                 elseif ismember(I, [2,3])
%                     likely_side_partner = 2;
%                 end
% 
%                 
% 
% %                 current_
% %                 likely_side = sides(1);
% %                 likely_side_partner = sides(2);
% % 
% %                 
% % 
% %                 needRetrack(frame,points(2),:) = 1;
% %                 needRetrack(frame,points(3),:) = 1;
%             end

%             sides = [1,2];
%             if xor(isnan(likely_side), isnan(likely_side_partner))
%                 if isnan(likely_side)
%                     likely_side = 2-(sign(likely_side_partner-1));
%                 else
%                     likely_side_partner = 2-(sign(likely_side-1));
%                 end
%             end

            for lr = 1:2 %left right
                if lr==1
                    if(likely_side==0)
                        %disp("Point not possible to reconstruct")
                        continue
                    end
                    %disp("Point possible to reconstruct")
                    current_bodyPoints = points(best_combination+1);
                    current_outPut_point = points(likely_side+1);
                else
                    if(likely_side_partner==0)
                        %disp("Point not possible to reconstruct")
                        continue
                    end
                    %disp("Point possible to reconstruct")
                    current_bodyPoints = points(best_partner_combination+1);                    
                    current_outPut_point = points(likely_side_partner+1);

                end

                if isempty(current_bodyPoints)
                    disp("Point not possible to reconstruct")
                    continue
                end
                use_or_not = [1,1,1];
                %use_or_not(current_confidences<minDLClikelihood) = nan;
                use_or_not(current_bodyPoints==points(1))=nan;

                current_confidences(1) = confidence_dlc(frame,current_bodyPoints(1),1);
                current_confidences(2) = confidence_dlc(frame,current_bodyPoints(2),2);
                current_confidences(3) = confidence_dlc(frame,current_bodyPoints(3),3);
                current_confidences = current_confidences';

                use_or_not(current_confidences<minDLClikelihood) = nan;

                clear output_dlc_confidenceTEMP
                for i=1:numel(current_bodyPoints)
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
                    ]);
            end

        else % mid body line point (7,8,9, or 16)
            current_outPut_point = bodypart;
            badPoint = bodypart+1;
            if badPoint == 17
                badPoint =15;
            end
            
            points = [badPoint bodypart];
            residuals = nan(size(possible_bilateral_combinations,1),1);
            containsNAN = zeros(size(possible_bilateral_combinations,1),1);

            for camcombo = 1:size(possible_bilateral_combinations,1)
                current_combination = possible_bilateral_combinations(camcombo,:);
                current_bodyPoints = points(current_combination+1);
                current_confidences(1) = confidence_dlc(frame,current_bodyPoints(1),1);
                current_confidences(2) = confidence_dlc(frame,current_bodyPoints(2),2);
                current_confidences(3) = confidence_dlc(frame,current_bodyPoints(3),3);
                current_confidences = current_confidences';

                use_or_not = [1,1,1];
                use_or_not(current_bodyPoints==points(1))=nan;
                use_or_not(current_confidences<minDLClikelihood) = nan;
                containsNAN(camcombo)=sum(isnan(use_or_not))>0;

                [~, residuals(camcombo)] = dlt_reconstruct(coefs, ...
                    [...
                    use_or_not(1)*image_xvalues_Undistorted(frame,current_bodyPoints(1),1) use_or_not(1)*image_yvalues_Undistorted(frame,current_bodyPoints(1),1) ...
                    use_or_not(2)*image_xvalues_Undistorted(frame,current_bodyPoints(2),2) use_or_not(2)*image_yvalues_Undistorted(frame,current_bodyPoints(2),2) ...
                    use_or_not(3)*image_xvalues_Undistorted(frame,current_bodyPoints(3),3) use_or_not(3)*image_yvalues_Undistorted(frame,current_bodyPoints(3),3) ...
                    ]);

            end

            residuals(find(containsNAN)) = residuals(find(containsNAN))+crp; %penalize using less than all cameras

            [~, best_combination_index] = min(residuals);
            best_combination = possible_bilateral_combinations(best_combination_index,:);
            current_bodyPoints = points(best_combination+1);
            use_or_not = [1,1,1];
            use_or_not(current_bodyPoints==points(1))=nan;

            current_confidences(1) = confidence_dlc(frame,current_bodyPoints(1),1);
            current_confidences(2) = confidence_dlc(frame,current_bodyPoints(2),2);
            current_confidences(3) = confidence_dlc(frame,current_bodyPoints(3),3);
            current_confidences = current_confidences';

            use_or_not(current_confidences<minDLClikelihood) = nan;

            clear output_dlc_confidenceTEMP
            for i=1:numel(current_bodyPoints)
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
                ]);
        end
    end
end
