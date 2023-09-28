function [xy, dlc_confidence] = getXYfromDLC(DLCoutputCSVs,ThruTrackerTranslationFiles)

NCams = numel(DLCoutputCSVs);

for cam = 1:NCams
    DLC_data = readtable(DLCoutputCSVs(cam));

    if exist('ThruTrackerTranslationFiles','var')
        tf = load(ThruTrackerTranslationFiles(cam));
        XYtop_left = tf.XYtop_left;
    else
        XYtop_left = zeros(size(DLC_data,1),2);
    end

    numberOfBodyPoints = (size(DLC_data,2)-1)/3;
    xTEMP = DLC_data(:,2:3:end);
    xTEMP = xTEMP{:,:}; %table 2 matrix
    xTEMP = xTEMP + XYtop_left(:,1);
    
    yTEMP = DLC_data(:,3:3:end);
    yTEMP = yTEMP{:,:}; %table 2 matrix
    yTEMP = yTEMP + XYtop_left(:,2);

    xy(:, :, 1, cam) = xTEMP;
    xy(:, :, 2, cam) = yTEMP;    
    
    dlc_confidenceTEMP = DLC_data(:,4:3:end);
    dlc_confidenceTEMP = dlc_confidenceTEMP{:,:}; %table 2 matrix
    
    dlc_confidence(:,:,cam) = dlc_confidenceTEMP;
end