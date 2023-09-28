function [xyz,DLTdvResidual] = getXYZfromDLTdv(DLTdvFile, easyWandDataFile, cam1TformsFile, cameraSelectionWeight, ThruTrackerTranslationFiles)
load(DLTdvFile);

NCams = udExport.data.nvid;

xyData = udExport.data.xypts;
xyData = full(xyData);
xyData(xyData==0) = nan;

load(easyWandDataFile);
coefs = easyWandData.coefs;
load(cam1TformsFile, 'camud');

idx = [];

for cam = 1:NCams
    xcol = cam*2-1;
    ycol = xcol+1;
    cols = sort([xcol:NCams*2:size(xyData,2) ycol:NCams*2:size(xyData,2)]);
    idxCurrennt = find(sum(xyData(:,cols),2, 'omitnan')>0);
    idx = union(idx,idxCurrennt);
end

for cam = 1:NCams
     xcol = cam*2-1;
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
    xy_curr_cam = xyData(idx,cols);
    image_xvalues(:,:,cam) = xy_curr_cam(:,1:2:end)+xtranslation;
    image_yvalues(:,:,cam) = xy_curr_cam(:,2:2:end)+ytranslation; 
end

coordinateDimensions = [size(image_xvalues,1), size(image_xvalues,2), 3];

image_xy_Undistorted = undistort_Tform([image_xvalues(:) image_yvalues(:)],camud);
image_xvalues_Undistorted = image_xy_Undistorted(:,1);
image_yvalues_Undistorted = image_xy_Undistorted(:,2);

image_xvalues_Undistorted = reshape(image_xvalues_Undistorted, coordinateDimensions);
image_yvalues_Undistorted = reshape(image_yvalues_Undistorted, coordinateDimensions);

for bodyPoint=1:coordinateDimensions(2)
    [xyz123(:,bodyPoint,:), res123(:,bodyPoint)] = dlt_reconstruct(coefs, [image_xvalues_Undistorted(:,bodyPoint,1) image_yvalues_Undistorted(:,bodyPoint,1) image_xvalues_Undistorted(:,bodyPoint,2) image_yvalues_Undistorted(:,bodyPoint,2) image_xvalues_Undistorted(:,bodyPoint,3) image_yvalues_Undistorted(:,bodyPoint,3)]);
    [xyz12(:,bodyPoint,:),  res12(:,bodyPoint)] =  dlt_reconstruct(coefs, [image_xvalues_Undistorted(:,bodyPoint,1) image_yvalues_Undistorted(:,bodyPoint,1) image_xvalues_Undistorted(:,bodyPoint,2) image_yvalues_Undistorted(:,bodyPoint,2) nan(size(image_xvalues_Undistorted(:,bodyPoint,3))) nan(size(image_yvalues_Undistorted(:,bodyPoint,3)))]);
    [xyz13(:,bodyPoint,:),  res13(:,bodyPoint)] =  dlt_reconstruct(coefs, [image_xvalues_Undistorted(:,bodyPoint,1) image_yvalues_Undistorted(:,bodyPoint,1) nan(size(image_xvalues_Undistorted(:,bodyPoint,2))) nan(size(image_yvalues_Undistorted(:,bodyPoint,2))) image_xvalues_Undistorted(:,bodyPoint,3) image_yvalues_Undistorted(:,bodyPoint,3)]);
    [xyz23(:,bodyPoint,:),  res23(:,bodyPoint)] =  dlt_reconstruct(coefs, [nan(size(image_xvalues_Undistorted(:,bodyPoint,1))) nan(size(image_yvalues_Undistorted(:,bodyPoint,1))) image_xvalues_Undistorted(:,bodyPoint,2) image_yvalues_Undistorted(:,bodyPoint,2) image_xvalues_Undistorted(:,bodyPoint,3) image_yvalues_Undistorted(:,bodyPoint,3)]);

    for frameNumber = 1:size(xyz123,1)
        [minRes, I] = min([...
            res12(frameNumber,bodyPoint) ...
            res13(frameNumber,bodyPoint) ...
            res23(frameNumber,bodyPoint) ...
            ]);
        if(cameraSelectionWeight+minRes< res123(frameNumber,bodyPoint))
            switch I
                case 1
                    xyz(frameNumber,bodyPoint,:) = xyz12(frameNumber,bodyPoint,:);
                    DLTdvResidual(frameNumber,bodyPoint) = res12(frameNumber,bodyPoint);
                case 2
                    xyz(frameNumber,bodyPoint,:) = xyz13(frameNumber,bodyPoint,:);
                    DLTdvResidual(frameNumber,bodyPoint) = res13(frameNumber,bodyPoint);
                case 3
                    xyz(frameNumber,bodyPoint,:) = xyz23(frameNumber,bodyPoint,:);
                    DLTdvResidual(frameNumber,bodyPoint) = res23(frameNumber,bodyPoint);
            end
        else
            xyz(frameNumber,bodyPoint,:) = xyz123(frameNumber,bodyPoint,:);
            DLTdvResidual(frameNumber,bodyPoint) = res123(frameNumber,bodyPoint);
        end
    end
end

coordinateDimensions(1) = size(xyData,1);
xyzTEMP = nan(coordinateDimensions);
DLTdvResidualTEMP = nan(coordinateDimensions(1:2));
xyzTEMP(idx,:,:) = xyz;
DLTdvResidualTEMP(idx,:,:) = DLTdvResidual;
xyz = xyzTEMP;
DLTdvResidual = DLTdvResidualTEMP;
end

