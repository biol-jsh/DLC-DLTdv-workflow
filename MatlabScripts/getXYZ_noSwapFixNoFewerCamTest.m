function [xyzDLC, DLTdvResidual, meanLikelihoodDLC] = ...
    getXYZ_noSwapFixNoFewerCamTest(easyWandDataFile, camTformsFiles, DLCoutputCSVs, ThruTrackerTranslationFiles, crp, minDLClikelihood)
load(easyWandDataFile);
coefs = easyWandData.coefs;
ncams = numel(DLCoutputCSVs);

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
    for bodypart = [1:16]

        current_outPut_point = bodypart;

        current_confidences(1) = confidence_dlc(frame,bodypart,1);
        current_confidences(2) = confidence_dlc(frame,bodypart,2);
        current_confidences(3) = confidence_dlc(frame,bodypart,3);
        current_confidences = current_confidences';
        use_or_not = [1,1,1];
        use_or_not(current_confidences<minDLClikelihood) = nan;


        meanLikelihoodDLC(frame,current_outPut_point) = mean(current_confidences(:).*use_or_not(:),'omitnan');

        [xyzDLC(frame,current_outPut_point,:), DLTdvResidual(frame,current_outPut_point)] = dlt_reconstruct(coefs, ...
            [...
            use_or_not(1)*image_xvalues_Undistorted(frame,bodypart,1) use_or_not(1)*image_yvalues_Undistorted(frame,bodypart,1) ...
            use_or_not(2)*image_xvalues_Undistorted(frame,bodypart,2) use_or_not(2)*image_yvalues_Undistorted(frame,bodypart,2) ...
            use_or_not(3)*image_xvalues_Undistorted(frame,bodypart,3) use_or_not(3)*image_yvalues_Undistorted(frame,bodypart,3) ...
            ]);
    end
end
