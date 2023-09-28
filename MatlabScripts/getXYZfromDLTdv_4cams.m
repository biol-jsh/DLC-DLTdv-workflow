function [xyz,DLTdvResidual] = getXYZfromDLTdv_4cams(DLTdvFile, easyWandDataFile, camTformsFiles, cameraSelectionWeight, ThruTrackerTranslationFiles)
load(DLTdvFile);

NCams = udExport.data.nvid;

xyData = udExport.data.xypts;
xyData = full(xyData);
xyData(xyData==0) = nan;

load(easyWandDataFile);
coefs = easyWandData.coefs;

idx = [];

for cam = 1:NCams
    xcol = cam*2-1; %First column of
    ycol = xcol+1;
    cols = sort([xcol:NCams*2:size(xyData,2) ycol:NCams*2:size(xyData,2)]);
    idxCurrennt = find(sum(xyData(:,cols),2, 'omitnan')>0);
    idx = union(idx,idxCurrennt);
end

for cam = 1:NCams
    xcol = cam*2-1; %First column of
    ycol = xcol+1;
    cols = sort([xcol:NCams*2:size(xyData,2) ycol:NCams*2:size(xyData,2)]);
    if(~isempty(ThruTrackerTranslationFiles))
        load(ThruTrackerTranslationFiles(cam))
        xtranslation = XYtop_left(idx,1);
        ytranslation = XYtop_left(idx,2);
    else
        xtranslation = 0;
        ytranslation = 0;
    end
    xy_curr_cam = xyData(:,cols);

    image_xvalues(:,:,cam) = nan(size(xy_curr_cam(:,1:2:end)));
    image_yvalues(:,:,cam) = nan(size(xy_curr_cam(:,2:2:end)));

    image_xvalues(idx,:,cam) = xy_curr_cam(idx,1:2:end)+xtranslation;
    image_yvalues(idx,:,cam) = xy_curr_cam(idx,2:2:end)+ytranslation;

    load(camTformsFiles(cam), 'camud');

    image_xvaluesTEMP = image_xvalues(:,:,cam);
    image_xvaluesTEMP = image_xvaluesTEMP(:);

    image_yvaluesTEMP = image_yvalues(:,:,cam);
    image_yvaluesTEMP = image_yvaluesTEMP(:);

    image_xy_Undistorted = undistort_Tform([image_xvaluesTEMP image_yvaluesTEMP],camud);
    image_xvalues_UndistortedTEMP = image_xy_Undistorted(:,1);
    image_yvalues_UndistortedTEMP = image_xy_Undistorted(:,2);

    image_xvalues_Undistorted(:,:,cam) = reshape(image_xvalues_UndistortedTEMP, size(image_xvalues,1), size(image_xvalues,2));
    image_yvalues_Undistorted(:,:,cam) = reshape(image_yvalues_UndistortedTEMP, size(image_yvalues,1), size(image_yvalues,2));
end

coordinateDimensions = [size(image_xvalues,1), size(image_xvalues,2), 3];
%
% image_xy_Undistorted = undistort_Tform([image_xvalues(:) image_yvalues(:)],camud);
% image_xvalues_Undistorted = image_xy_Undistorted(:,1);
% image_yvalues_Undistorted = image_xy_Undistorted(:,2);
%
% image_xvalues_Undistorted = reshape(image_xvalues_Undistorted, coordinateDimensions);
% image_yvalues_Undistorted = reshape(image_yvalues_Undistorted, coordinateDimensions);

for bodyPoint=1:coordinateDimensions(2)
    xyInput = [image_xvalues_Undistorted(:,bodyPoint,1) image_yvalues_Undistorted(:,bodyPoint,1) image_xvalues_Undistorted(:,bodyPoint,2) image_yvalues_Undistorted(:,bodyPoint,2) image_xvalues_Undistorted(:,bodyPoint,3) image_yvalues_Undistorted(:,bodyPoint,3) image_xvalues_Undistorted(:,bodyPoint,4) image_yvalues_Undistorted(:,bodyPoint,4)];
%     if(bodyPoint==8)
%         xyInput(:,5:6) = nan;
%     end
    % using all four cameras
    [xyz4(:,bodyPoint,:), res4(:,bodyPoint)] = dlt_reconstruct(coefs, xyInput);

    % with less cameras:
    cams = 1:NCams;
    % 3 cameras:
    permutations = nchoosek(cams,1); %choose one camera to disregard
    for i=1:size(permutations,1)
        camsTemp = cams(permutations(i,:));
        xyInputTEMP = xyInput;
        for currCam = camsTemp
            xyInputTEMP(:,currCam*2) = nan(size(xyInputTEMP(:,currCam*2)));
            xyInputTEMP(:,currCam*2-1) = nan(size(xyInputTEMP(:,currCam*2-1)));
        end
        [xyz3Alternatives(:,bodyPoint,:,i), res3Alternativs(:,bodyPoint,i)] = dlt_reconstruct(coefs, xyInputTEMP);
    end
    for frameNumber = 1:size(xyz3Alternatives,1)
        [res3(frameNumber, bodyPoint), minResIndex] = min(res3Alternativs(frameNumber,bodyPoint,:));
        xyz3(frameNumber,bodyPoint,:) = xyz3Alternatives(frameNumber,bodyPoint,:,minResIndex);
    end

    % 2 cameras:
    permutations = nchoosek(cams,2); %choose two cameras to disregard
    for i=1:size(permutations,1)
        camsTemp = cams(permutations(i,:));
        xyInputTEMP = xyInput;
        for currCam = camsTemp
            xyInputTEMP(:,currCam*2) = nan(size(xyInputTEMP(:,currCam*2)));
            xyInputTEMP(:,currCam*2-1) = nan(size(xyInputTEMP(:,currCam*2-1)));
        end
        [xyz2Alternatives(:,bodyPoint,:,i), res2Alternativs(:,bodyPoint,i)] = dlt_reconstruct(coefs, xyInputTEMP);
    end
    for frameNumber = 1:size(xyz2Alternatives,1)
        [res2(frameNumber, bodyPoint), minResIndex] = min(res2Alternativs(frameNumber,bodyPoint,:));
        xyz2(frameNumber,bodyPoint,:) = xyz2Alternatives(frameNumber,bodyPoint,:,minResIndex);
    end

    %     [xyz12(:,bodyPoint,:),  res12(:,bodyPoint)] =  dlt_reconstruct_fast(coefs, [image_xvalues_Undistorted(:,bodyPoint,1) image_yvalues_Undistorted(:,bodyPoint,1) image_xvalues_Undistorted(:,bodyPoint,2) image_yvalues_Undistorted(:,bodyPoint,2) nan(size(image_xvalues_Undistorted(:,bodyPoint,3))) nan(size(image_yvalues_Undistorted(:,bodyPoint,3)))]);
    %     [xyz13(:,bodyPoint,:),  res13(:,bodyPoint)] =  dlt_reconstruct_fast(coefs, [image_xvalues_Undistorted(:,bodyPoint,1) image_yvalues_Undistorted(:,bodyPoint,1) nan(size(image_xvalues_Undistorted(:,bodyPoint,2))) nan(size(image_yvalues_Undistorted(:,bodyPoint,2))) image_xvalues_Undistorted(:,bodyPoint,3) image_yvalues_Undistorted(:,bodyPoint,3)]);
    %     [xyz23(:,bodyPoint,:),  res23(:,bodyPoint)] =  dlt_reconstruct_fast(coefs, [nan(size(image_xvalues_Undistorted(:,bodyPoint,1))) nan(size(image_yvalues_Undistorted(:,bodyPoint,1))) image_xvalues_Undistorted(:,bodyPoint,2) image_yvalues_Undistorted(:,bodyPoint,2) image_xvalues_Undistorted(:,bodyPoint,3) image_yvalues_Undistorted(:,bodyPoint,3)]);

    for frameNumber = 1:size(xyz4,1)
        [~, I] = min([...
            res4(frameNumber,bodyPoint) ...
            res3(frameNumber,bodyPoint)+cameraSelectionWeight ...
            res2(frameNumber,bodyPoint)+2*cameraSelectionWeight ...
            ]);
        switch I
            case 1
                xyz(frameNumber,bodyPoint,:) = xyz4(frameNumber,bodyPoint,:);
                DLTdvResidual(frameNumber,bodyPoint) = res4(frameNumber,bodyPoint);
            case 2
                xyz(frameNumber,bodyPoint,:) = xyz3(frameNumber,bodyPoint,:);
                DLTdvResidual(frameNumber,bodyPoint) = res3(frameNumber,bodyPoint);
            case 3
                xyz(frameNumber,bodyPoint,:) = xyz2(frameNumber,bodyPoint,:);
                DLTdvResidual(frameNumber,bodyPoint) = res2(frameNumber,bodyPoint);
        end
    end

end

% idx = find(nansum(xyz(:,:,1),1)>0);
% idx
% plot3(xyz(:,1,1),xyz(:,1,2),xyz(:,1,3)); hold on;
%
% coordinateDimensions(1) = size(xyData,1);
% xyzTEMP = nan(coordinateDimensions);
% DLTdvResidualTEMP = nan(coordinateDimensions(1:2));
% xyzTEMP(idx,:,:) = xyz;
% DLTdvResidualTEMP(idx,:,:) = DLTdvResidual;
% xyz = xyzTEMP;
% DLTdvResidual = DLTdvResidualTEMP;
end

