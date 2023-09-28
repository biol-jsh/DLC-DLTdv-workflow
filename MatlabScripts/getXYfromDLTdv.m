function [xy, digitized_frames] = getXYfromDLTdv(DLTdvFile)
load(DLTdvFile);

NCams = udExport.data.nvid;

xyData = udExport.data.xypts;
xyData = full(xyData);
xyData(xyData==0) = nan;

for cam = 1:NCams
    xcol = cam*2-1; %First column of
    ycol = xcol+1;
    cols = sort([xcol:NCams*2:size(xyData,2) ycol:NCams*2:size(xyData,2)]);
    %digitized_frames = find(nansum(xyData(:,cols),2)>0);
    digitized_frames = find(sum(xyData(:,:),2, 'omitnan')>0);
    xy_curr_cam = xyData(digitized_frames,cols);

    xy(:,:, 1, cam) = xy_curr_cam(:,1:2:end);
    xy(:,:, 2, cam) = xy_curr_cam(:,2:2:end); 
end

