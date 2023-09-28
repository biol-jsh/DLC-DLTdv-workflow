%% 
% This script produces all the figures and statistics related to the wind
% tunnel part.
% TODO: output overlap table in more readable way
% TODO: add clarifications in script
% TODO: maybe get away from all the eval statements. Probably not but would
% be nice. They seemed like a good method at first but as project grew it
% has become rather hard to read so not ideal for people wanting to fully
% understand what's going on here.


%%
% Check if DLTdv available
if exist("dltdv-master") == 7
    addpath(genpath("dltdv-master"))
else
    msgbox({'You need to download DLTdv from e.g.'; ...
        'https://github.com/tlhedrick/dltdv/tree/master'; ...
        'See Command Window for link'}, "Missing DLTdv functions", "Error")
    disp("Web address to DLTdv GitHub: https://github.com/tlhedrick/dltdv/tree/master")
    return
end

% Check if file exchange function shade (filled area plot) exists
if exist("shade") ~= 7
    disp("https://www.mathworks.com/matlabcentral/fileexchange/69652-filled-area-plot")
end


% Check MATLAB dependencies (toolboxes)
if ~any(any(contains(struct2cell(ver), 'Curve Fitting Toolbox')))
    msgbox({'You need to install the Curve Fitting Toolbox.'}, "Missing toolbox", "Error")
    return
end
% Check MATLAB dependencies (toolboxes)
if ~any(any(contains(struct2cell(ver), 'Signal Processing Toolbox')))
    msgbox({'You need to install the Signal Processing Toolbox.'}, "Missing toolbox", "Error")
    return
end
% Check MATLAB dependencies (toolboxes)
if ~any(any(contains(struct2cell(ver), 'Statistics and Machine Learning Toolbox')))
    msgbox({'You need to install the Statistics and Machine Learning Toolbox.'}, "Missing toolbox", "Error")
    return
end

% Check MATLAB dependencies (toolboxes)
if ~any(any(contains(struct2cell(ver), 'Bioinformatics Toolbox')))
    msgbox({'You need to install the Bioinformatics Toolbox.'}, "Missing toolbox", "Error")
    return
end
wingbeats = {}; % struc for saving wingbeats

wt_to_wt_points = [1, 15]; % the two wingtip landmarks
paired_landmarks = [... % left and right point for...
    1,15; ... % wingtips
    2,13; ... % wrists
    3,14; ... % 5th digits
    4,12; ... % elbows
    5,10; ... % shoulders
    6,11  ... % ankles
    ];

land_mark_numbers = [1:9 16];

%% set up filter parameters
speedLimitsFile = "speed_075.mat";
load(speedLimitsFile);
goodSpeeds = maxSpeed*2;
maxSpeeds = maxSpeed*4;
frameRate = 700; % not only for filtering
crp = 5;
smoothspan = 4;
goodDLTdvResidual = 10;
maxDLTdvResidual = 20;
goodDLClikelihood = .85;
minDLClikelihood = .4;
minDLClikelihood_3d = 0.5;

%% 3D calibration files, filmed over two days, so two calibrations
easyWandDataFile_apr1 = fullfile(pwd,"windtunnelCalibration1\extrinsicCalib_20230228_easyWandData.mat");
easyWandDataFile_apr2 = fullfile(pwd,"windtunnelCalibration2\extrinsicCalib_20230228_easyWandData.mat");

% easyWand calibration files for the four individuals
easyWandDataFile(1) = easyWandDataFile_apr1;
easyWandDataFile(2) = easyWandDataFile_apr1;
easyWandDataFile(3) = easyWandDataFile_apr2;
easyWandDataFile(4) = easyWandDataFile_apr2;

x_axis_file_apr1 = "S:\Projects\2023_Methodology_paper\Brown2021_part\calibrations\1 Apr 2021\xaxis20211209_dvProject.mat";
x_axis_file_apr2 = "S:\Projects\2023_Methodology_paper\Brown2021_part\calibrations\2 Apr 2021\xaxis20211209_dvProject.mat";

% x axis files for the four individuals
% x used to rotate the coordinate system such that x axis is along
% direction of flow
x_axis_file(1) = x_axis_file_apr1;
x_axis_file(2) = x_axis_file_apr1;
x_axis_file(3) = x_axis_file_apr2;
x_axis_file(4) = x_axis_file_apr2;

camTformsFile_apr1(1) = fullfile(pwd,"windtunnelCalibration1\extrinsicCalib_20230228_cam1Tforms.mat");
camTformsFile_apr1(2) = fullfile(pwd,"windtunnelCalibration1\extrinsicCalib_20230228_cam2Tforms.mat");
camTformsFile_apr1(3) = fullfile(pwd,"windtunnelCalibration1\extrinsicCalib_20230228_cam3Tforms.mat");
camTformsFile_apr1(4) = fullfile(pwd,"windtunnelCalibration1\extrinsicCalib_20230228_cam4Tforms.mat");

camTformsFile_apr2(1) = fullfile(pwd,"windtunnelCalibration2\extrinsicCalib_20230228_cam1Tforms.mat");
camTformsFile_apr2(2) = fullfile(pwd,"windtunnelCalibration2\extrinsicCalib_20230228_cam2Tforms.mat");
camTformsFile_apr2(3) = fullfile(pwd,"windtunnelCalibration2\extrinsicCalib_20230228_cam3Tforms.mat");
camTformsFile_apr2(4) = fullfile(pwd,"windtunnelCalibration2\extrinsicCalib_20230228_cam4Tforms.mat");

camTforms{1} = camTformsFile_apr1;
camTforms{2} = camTformsFile_apr1;
camTforms{3} = camTformsFile_apr2;
camTforms{4} = camTformsFile_apr2;

%% Analysis of automatically digitized trials
clear speeds % in case we have run it twice, don't want to mix with speeds for manual trials
speeds(1,:) = [3, 2.973, 2.978, 4.477, 4.587, 4.633, 5.96, 5.97, 5.907];
speeds(2,:) = [4.49, 4.547, 4.595, 5.908, 5.92, 5.992, 2.994, 3.013, 3.011];
speeds(3,:) = [6.046, 5.957, 6.079, 4.315, 4.467, 4.537, 2.915, 3.076, 3.092];
speeds(4,:) = [4.366, 4.507, 4.461, 4.614, 2.963, 2.971, 2.965, 5.971, nan];

speeds_nominal(1,:) = [3, 3, 3, 4.5, 4.5, 4.5, 6, 6, 6];
speeds_nominal(2,:) = [4.5, 4.5, 4.5, 6, 6, 6, 3, 3, 3];
speeds_nominal(3,:) = [6, 6, 6, 4.5, 4.5, 4.5, 3, 3, 3];
speeds_nominal(4,:) = [4.5, 4.5, 4.5, 4.5, 3, 3, 3, 6, nan];

batFolders(1) = fullfile(pwd,"windtunnelTrials\April 1\3AE9E\Markerless");
batFolders(2) = fullfile(pwd,"windtunnelTrials\April 1\0790C\Markerless");
batFolders(3) = fullfile(pwd,"windtunnelTrials\April 2\078A3\Markerless");
batFolders(4) = fullfile(pwd,"windtunnelTrials\April 2\FE404\Markerless");

% due to different camera rotations used during calibration vs DLC
% analysis, we need the video dimensions to rotate the coordinates back
load(fullfile(pwd, "windtunnelTrials", "vidSizesWindTunnelAuto.mat"), "vidSizes");

% Based on visual inspection, determine if right or left side is to be
% analuzed. 1 = left wingtip, 15 = right wingtip, nan = that trial number
% is not analyzed

wts_to_analyze(1,:) = [1, 1, 1, 1, 1, 1, 1, 1, 1]; % bat 1 - 3AE9E
wts_to_analyze(2,:) = [1, 1, 15, 1, 1, 15, 15, 1, 15]; % bat 2 - 0790C
wts_to_analyze(3,:) = [1,  1, nan, 1, 1, 15, 1, 15, 1]; % bat 3 - 078A3
wts_to_analyze(4,:) = [nan, 1, 1, 1, 15, 1, 1, 1, nan]; % bat 4 - FE404

Trials_to_use{1} = 1:9; % bat 1 - 3AE9E
Trials_to_use{2} = 1:9; % bat 2 - 0790C
Trials_to_use{3} = [1:2 4:9]; % bat 3 - 078A3
Trials_to_use{4} = 2:8; % bat 4 - FE404

for folderNumber = 1:numel(batFolders)
    currentBatFolder = batFolders(folderNumber);
    current_bat = split(currentBatFolder, "\");
    current_bat = current_bat(end-1); % get name of bat based on folder name

    [x_axis ~] = getXYZfromDLTdv_4cams(x_axis_file(folderNumber), easyWandDataFile(folderNumber), camTforms{folderNumber}, crp, []);
    x_axis = squeeze(x_axis(find(abs(sum(x_axis(:,1,1),2,'omitnan'))>0),:,:));

    % not actually getting velocity, just direction of points along
    % calibration object
    [x_axis,~] = getVelocityAndAcceleration(x_axis, 1, false); % false to not plot
    x_axis = squeeze(x_axis(find(abs(sum(x_axis(:,1,1),2,'omitnan'))>0),:,:));
    x_axis(2) = -x_axis(2); % flip y axis due to bug in easyWand

    for trialNumber = 1:numel(Trials_to_use{folderNumber})% %trial_number_to_check
        trial = strcat("F",num2str(Trials_to_use{folderNumber}(trialNumber)));
        actualTrialNumber = Trials_to_use{folderNumber}(trialNumber);

        wt_to_analyze = wts_to_analyze(folderNumber, actualTrialNumber);
        if(wt_to_analyze == 1)
            points_to_analyze = [1:9];
        elseif(wt_to_analyze == 15)
            points_to_analyze = [15, 13, 14, 12, 10, 11, 7, 8, 9];
        end

        currentFolder = fullfile(currentBatFolder,trial);
        cropFiles = dir(currentFolder);
        cropFiles = string({cropFiles.name})';
        cropFiles = cropFiles(contains(cropFiles, "crop_"));
        cropFiles = sort(cropFiles);
        cropFiles = fullfile(currentFolder, cropFiles);

        CSVs = dir(currentFolder);
        CSVs = string({CSVs.name})';
        CSVs = CSVs(contains(CSVs,"csv"));
        DLCoutputCSVs(1) = CSVs(contains(CSVs, "Bane"));
        DLCoutputCSVs(2) = CSVs(contains(CSVs, "Freeze"));
        DLCoutputCSVs(3) = CSVs(contains(CSVs, "Joker"));
        DLCoutputCSVs(4) = CSVs(contains(CSVs, "Penguin"));
        DLCoutputCSVs = fullfile(currentFolder, DLCoutputCSVs);
        DLCoutputCSVs = DLCoutputCSVs';

        current_vidSizes = vidSizes.(strcat("bat_",current_bat)).(strcat("F", num2str(actualTrialNumber)));

        [xyz_raw, DLTdvResidual, meanLikelihoodDLC] = getXYZfromDLC_4cam_rotation(easyWandDataFile(folderNumber), camTforms{folderNumber}, DLCoutputCSVs, cropFiles, crp,  minDLClikelihood, paired_landmarks, land_mark_numbers,[0, 180, 180, 180], current_vidSizes);

        % this should be put in function
        xyz_raw(xyz_raw==0)=nan;
        DLTdvResidual(DLTdvResidual==0) = inf;
        DLTdvResidual(isnan(DLTdvResidual)) = inf;

        xyz_raw(:,:,2) = -xyz_raw(:,:,2);
        [xyz_raw, xy_rotation_angle] = ...
            xy_CoordinateTransformation(xyz_raw, repmat(x_axis,[size(xyz_raw,1) 1]));

        x_displacement = speeds(folderNumber,trialNumber)/frameRate*(1:size(xyz_raw,1));
        xyz_raw(:,:,1) = xyz_raw(:,:,1) + repmat(x_displacement',[1 size(xyz_raw,2)]);

        xyz_filtered = ...
            filter3Dpredictions(...
            xyz_raw, ...
            DLTdvResidual, ...
            meanLikelihoodDLC, ...
            frameRate, ...
            goodSpeeds, ...
            maxSpeeds, ...
            goodDLTdvResidual, ...
            maxDLTdvResidual, ...
            goodDLClikelihood, ...
            minDLClikelihood_3d);

        xyz_smooth = nan(size(xyz_raw));
        for k=1:size(xyz_raw,2)
            for L=1:3
                try
                    xyz_filtered(find(isoutlier(xyz_filtered(:,k,L),'movmedian',smoothspan*2)),k,L) = nan;
                    firstIndex = find(~isnan(xyz_filtered(:,k,L)),true,'first');
                    lastIndex = find(~isnan(xyz_filtered(:,k,L)),true,'last');
                    xyz_smooth(firstIndex:lastIndex,k,L) = tybutterNaN(smooth(xyz_filtered(firstIndex:lastIndex,k,L),smoothspan, 'rlowess'),45,700,'low');
                end
            end
        end
        xyz_smooth(xyz_smooth==0) = nan;

        wt_to_wt_vector = xyz_smooth(:,wt_to_wt_points(2),:)-xyz_smooth(:,wt_to_wt_points(1),:);
        wt_to_wt_vector = squeeze(wt_to_wt_vector);

        clear wt_to_wt_length wt_to_wt_angle
        for L=1:size(wt_to_wt_vector,1)
            wt_to_wt_length(L) = norm(wt_to_wt_vector(L,:));
        end
        wt_to_wt_length(1:10) = nan;
        [~, midDS_index] = max(wt_to_wt_length(1:end-20));

        nframes = size(xyz_raw,1);

        strn = squeeze(xyz_raw(:,8,:));
        time=(1:size(xyz_raw,1))/frameRate;

        clear strn_smooth
        for k=1:3
            [fitTime, fitData] = prepareCurveData(time', strn(:,k));
            [fitFunction_speed_point, gof] = fit( fitTime, fitData, fittype( 'poly2' ), 'Normalize', 'on');
            strn_smooth(:,k) = fitFunction_speed_point(time);
        end

        [strnvel, ~] = getVelocityAndAcceleration(strn_smooth,frameRate, false);
        hold on;

        jhat = strnvel./vecnorm(strnvel,2,2);

        if false && folderNumber == 2 && actualTrialNumber == 3
            ihatEst = wt_to_wt_vector;
            ihatEst = ihatEst./vecnorm(ihatEst,2,2);
        else
            i_one_frame = wt_to_wt_vector(midDS_index,:);
            i_one_frame = i_one_frame/norm(i_one_frame);
            ihatEst = repmat(i_one_frame, nframes, 1);
        end

        khat = cross(ihatEst,jhat);
        khat = khat./vecnorm(khat,2,2);

        ihat = cross(jhat, khat);

        xyz_BCS_smooth = nan(size(xyz_raw));

        for bodyPoint = 1:16
            xyz_BCS_smooth(:,bodyPoint,:) = changeOfCoordinateSys(squeeze(xyz_smooth(:,bodyPoint,:)), jhat, ihat, khat, strn_smooth);
        end

        newVel = changeOfCoordinateSys(strnvel, jhat, ihat, khat, strn_smooth);

        [begDownStrokeZvalue, begDownStrokeIndexes] = findpeaks(xyz_BCS_smooth(:,points_to_analyze(1),3));
        [begUpStrokeZvalue, begUpStrokeIndexes] = findpeaks(-(xyz_BCS_smooth(:,points_to_analyze(1),3)));

        points = [1 15];

        wt_to_wt_vector = xyz_BCS_smooth(:,wt_to_wt_points(2),:)-xyz_BCS_smooth(:,wt_to_wt_points(1),:);
        wt_to_wt_vector = squeeze(wt_to_wt_vector);

        for L=1:size(wt_to_wt_vector,1)
            wt_to_wt_length(L) = norm(wt_to_wt_vector(L,:));
            wt_to_wt_angle(L) = -sign(wt_to_wt_vector(L,1))*acos(dot(wt_to_wt_vector(L,[1 2]), [0 1])/norm(wt_to_wt_vector(L,[1 2])));
        end

        angle_to_use = wt_to_wt_angle(midDS_index);

        xyz_BCS_smooth2 = xy_CoordinateTransformation(xyz_BCS_smooth,repmat(angle_to_use,1,size(xyz_BCS_smooth,1)));

        wingbeats{end+1}.individual = current_bat;
        currentWB = numel(wingbeats);
        wingbeats{currentWB}.full_wb_BCS_smooth=xyz_BCS_smooth2;
        wingbeats{currentWB}.begDownStrokeIndexes=begDownStrokeIndexes;
        wingbeats{currentWB}.begUpStrokeIndexes=begUpStrokeIndexes;
        wingbeats{currentWB}.torso_fitted_xyz = strn_smooth;
        wingbeats{currentWB}.speed_nominal = speeds_nominal(folderNumber,actualTrialNumber);
        wingbeats{currentWB}.treatment = "auto";
        wingbeats{currentWB}.points_to_analyze = points_to_analyze;
        wingbeats{currentWB}.trial_number = actualTrialNumber;

    end
end

%% Analysis of manually digitized trials

batFolders(1) = "E:\Dropbox\Colorado2020-2023\Papers\Methods Paper\supplement\windtunnelTrials\April 1\3AE9E\Markered";
batFolders(2) = "E:\Dropbox\Colorado2020-2023\Papers\Methods Paper\supplement\windtunnelTrials\April 1\0790C\Markered";
batFolders(3) = "E:\Dropbox\Colorado2020-2023\Papers\Methods Paper\supplement\windtunnelTrials\April 2\078A3\Markered";
batFolders(4) = "E:\Dropbox\Colorado2020-2023\Papers\Methods Paper\supplement\windtunnelTrials\April 2\FE404\Markered";

clear speeds
speeds(1,:) = [2.892, 2.993, 2.928, 4.388, 4.532, 4.382, 6.201, 6.163, 6.103];
speeds(2,:) = [4.497, 4.609, 4.473, 6.032, 6.062, 6.046, 2.898, 2.972, 2.975];
speeds(3,:) = [5.93, 5.902, 6.011, 4.366, 4.496, 4.374, 3.057, 3.022, 2.945];
speeds(4,:) = [4.366, 4.507, 4.461, 3.071, 2.963, 2.971, 5.971, nan, nan];

speeds_nominal(1,:) = [3, 3, 3, 4.5, 4.5, 4.5, 6, 6, 6];
speeds_nominal(2,:) = [4.5, 4.5, 4.5, 6, 6, 6, 3, 3, 3];
speeds_nominal(3,:) = [6, 6, 6, 4.5, 4.5, 4.5, 3, 3, 3];
speeds_nominal(4,:) = [4.5, 4.5, 4.5, 3, 3, 3, 6, nan, nan];


frameLimits(1,:,:) = [... bat 1 - 3AE9E
    423, 506; ... 1
    508, 602; ... 2
    530, 620; ... 3
    426, 510; ... 4
    562, 635; ... 5
    794, 884; ... 6
    528, 677; ... 7
    551, 650; ... 8
    138, 248; ... 9
    ]';

frameLimits(2,:,:) = [... bat 2 - 0790C
    674, 756; ... 1
    155, 252; ... 2
    675, 785; ... 3
    410, 515; ... 4
    125, 215; ... 5
    520, 625; ... 6
    555, 630; ... 7
    520, 605; ... 8
    545, 627; ... 9
    ]';

frameLimits(3,:,:) = [... 078A3 markered
    792, 890; ... F1
    725, 825; ... F2
    785, 875; ... F3
    590, 685; ... F4
    780, 870; ... F5
    825, 925; ... F6
    777, 857; ... F7
    865, 973; ... F8
    820, 910; ... F9
    ]';

frameLimits(4,:,:) = [... bat 4 - FE404
    820, 911; ... F1
    790, 884; ... F2
    810, 935; ... F3
    807, 899; ... F4
    800, 884; ... F5
    850, 935; ... F6
    820, 945; ... F7
    nan, nan; ... F8
    nan, nan; ... F9
    ]';

wts_to_analyze(1,:) = [1, 1, 15, 1, 1, 1, 1, 1, 1]; % bat 1 - 3AE9E
wts_to_analyze(2,:) = [1, 15, 1, 1, 1, 1, 1, 1, 1]; % bat 2 - 0790C
wts_to_analyze(3,:) = [1, 1, 1, 1, 1, 1, 15 ,1 ,1]; % bat 3 - 078A3
wts_to_analyze(4,:) = [1, 1, 1, 1, 1, 1,  1, nan, nan]; % bat 4 - FE404

Trials_to_use{1} = 1:9;
Trials_to_use{2} = 1:9;
Trials_to_use{3} = 1:9;
Trials_to_use{4} = 1:7;

for folderNumber = 1:numel(batFolders)
    currentBatFolder = batFolders(folderNumber);
    current_bat = split(currentBatFolder, "\");
    current_bat = current_bat(end-1); % get name of bat based on folder name

    [x_axis ~] = getXYZfromDLTdv_4cams(x_axis_file(folderNumber), easyWandDataFile(folderNumber), camTforms{folderNumber}, crp, []);
    x_axis = squeeze(x_axis(find(abs(sum(x_axis(:,1,1),2, 'omitnan'))>0),:,:));

    % not actually getting velocity, just direction of points along
    % calibration object
    [x_axis,~] = getVelocityAndAcceleration(x_axis, 1, false); % false to not plot
    x_axis = squeeze(x_axis(find(abs(sum(x_axis(:,1,1),2, 'omitnan'))>0),:,:));
    x_axis(2) = -x_axis(2); % flip y axis due to bug in easyWand

    for trialNumber = 1:numel(Trials_to_use{folderNumber})
        trial = strcat("F",num2str(Trials_to_use{folderNumber}(trialNumber)));
        actualTrialNumber = Trials_to_use{folderNumber}(trialNumber);

        wt_to_analyze = wts_to_analyze(folderNumber, actualTrialNumber);
        if(wt_to_analyze == 1)
            points_to_analyze = [1:9];
        elseif(wt_to_analyze == 15)
            points_to_analyze = [15, 13, 14, 12, 10, 11, 7, 8, 9];
        end

        currentFolder = fullfile(currentBatFolder,trial);
        DLTdv_file = better_dir(currentFolder, 'dvProject.mat');
        DLTdv_file = fullfile(currentFolder, string(DLTdv_file.name));

        [xyz_raw, DLTdvResidual] = getXYZfromDLTdv_4cams(DLTdv_file, easyWandDataFile(folderNumber), camTforms{folderNumber}, crp, []);

        xyz_raw = xyz_raw(frameLimits(folderNumber,1,actualTrialNumber):frameLimits(folderNumber,2,actualTrialNumber),:,:);

        xyz_raw(:,:,2) = -xyz_raw(:,:,2);

        x_displacement = speeds(folderNumber,trialNumber)/frameRate*(1:size(xyz_raw,1));
        xyz_raw(:,:,1) = xyz_raw(:,:,1) + repmat(x_displacement',[1 size(xyz_raw,2)]);


        xyz_smooth = nan(size(xyz_raw));
        for k=1:15
            for L=1:3
                try
                    firstIndex = find(~isnan(xyz_raw(:,k,L)),true,'first');
                    lastIndex = find(~isnan(xyz_raw(:,k,L)),true,'last');
                    xyz_smooth(firstIndex:lastIndex,k,L) = tybutterNaN(smooth(xyz_raw(firstIndex:lastIndex,k,L), smoothspan, 'rlowess'),45,700,'low');
                end
            end
        end
        xyz_smooth(xyz_smooth==0) = nan;

        wt_to_wt_vector = xyz_smooth(:,wt_to_wt_points(2),:)-xyz_smooth(:,wt_to_wt_points(1),:);
        wt_to_wt_vector = squeeze(wt_to_wt_vector);

        clear wt_to_wt_length wt_to_wt_angle
        for L=1:size(wt_to_wt_vector,1)
            wt_to_wt_length(L) = norm(wt_to_wt_vector(L,:));
        end

        [~, midDS_index] = max(wt_to_wt_length);

        nframes = size(xyz_raw,1);

        strn = squeeze(xyz_raw(:,8,:));
        time=(1:size(xyz_raw,1))/frameRate;

        clear strn_smooth
        for k=1:3
            [fitTime, fitData] = prepareCurveData(time, strn(:,k));
            [fitFunction_speed_point, gof] = fit( fitTime, fitData, fittype( 'poly2' ), 'Normalize', 'on');
            strn_smooth(:,k) = fitFunction_speed_point(time);

        end

        [strnvel, ~] = getVelocityAndAcceleration(strn_smooth,frameRate, false);
        jhat = strnvel./vecnorm(strnvel,2,2);


        i_one_frame = wt_to_wt_vector(midDS_index,:);
        i_one_frame = i_one_frame/norm(i_one_frame);
        ihatEst = repmat(i_one_frame, nframes, 1);


        khat = cross(ihatEst,jhat);
        khat = khat./vecnorm(khat,2,2);

        ihat = cross(jhat, khat);

        xyz_BCS_smooth = nan(size(xyz_raw));
        for bodyPoint = 1:15
            xyz_BCS_smooth(:,bodyPoint,:) = changeOfCoordinateSys(squeeze(xyz_smooth(:,bodyPoint,:)), jhat, ihat, khat, strn_smooth);
        end

        [begDownStrokeZvalue, begDownStrokeIndexes] = findpeaks(xyz_BCS_smooth(:,points_to_analyze(1),3));
        [begUpStrokeZvalue, begUpStrokeIndexes] = findpeaks(-(xyz_BCS_smooth(:,points_to_analyze(1),3)));

        points = [1 15];

        wt_to_wt_vector = xyz_BCS_smooth(:,wt_to_wt_points(2),:)-xyz_BCS_smooth(:,wt_to_wt_points(1),:);
        wt_to_wt_vector = squeeze(wt_to_wt_vector);

        for L=1:size(wt_to_wt_vector,1)
            wt_to_wt_length(L) = norm(wt_to_wt_vector(L,:));
            wt_to_wt_angle(L) = -sign(wt_to_wt_vector(L,1))*acos(dot(wt_to_wt_vector(L,[1 2]), [0 1])/norm(wt_to_wt_vector(L,[1 2])));
        end

        angle_to_use = wt_to_wt_angle(midDS_index);

        xyz_BCS_smooth2 = xy_CoordinateTransformation(xyz_BCS_smooth,repmat(angle_to_use,1,size(xyz_BCS_smooth,1)));

        wingbeats{end+1}.individual = current_bat;
        currentWB = numel(wingbeats);
        wingbeats{currentWB}.full_wb_BCS_smooth=xyz_BCS_smooth2;
        wingbeats{currentWB}.begDownStrokeIndexes=begDownStrokeIndexes;
        wingbeats{currentWB}.begUpStrokeIndexes=begUpStrokeIndexes;
        wingbeats{currentWB}.torso_fitted_xyz = strn_smooth;
        wingbeats{currentWB}.speed_nominal = speeds_nominal(folderNumber,actualTrialNumber);
        wingbeats{currentWB}.treatment = "manual";
        wingbeats{currentWB}.points_to_analyze = points_to_analyze;
        wingbeats{currentWB}.trial_number = actualTrialNumber;


    end
end

%% stats and figures
% Check if shade plotting package available
if exist("shade") == 7
    addpath(genpath("shade"))
else
    msgbox({'You need to download shade from e.g.'; ...
        'https://www.mathworks.com/matlabcentral/fileexchange/69652-filled-area-plot'; ...
        'See Command Window for link'}, "Missing Shade function", "Error")
    disp("Web address to Shade on MathWorks File Exchange: https://www.mathworks.com/matlabcentral/fileexchange/69652-filled-area-plot")
    return
end
correlationCoeff = nan(6,9);
correlationP = nan(6,9);

for p = 1:6
    eval(string(['point_', char(num2str(p)), '_flap_amplitudeX.auto_3mps = []']));
    eval(string(['point_', char(num2str(p)), '_flap_amplitudeY.auto_3mps = []']));
    eval(string(['point_', char(num2str(p)), '_flap_amplitudeZ.auto_3mps = []']));

    eval(string(['point_', char(num2str(p)), '_flap_amplitudeX.auto_4p5mps = []']));
    eval(string(['point_', char(num2str(p)), '_flap_amplitudeY.auto_4p5mps = []']));
    eval(string(['point_', char(num2str(p)), '_flap_amplitudeZ.auto_4p5mps = []']));

    eval(string(['point_', char(num2str(p)), '_flap_amplitudeX.auto_6mps = []']));
    eval(string(['point_', char(num2str(p)), '_flap_amplitudeY.auto_6mps = []']));
    eval(string(['point_', char(num2str(p)), '_flap_amplitudeZ.auto_6mps = []']));

    eval(string(['point_', char(num2str(p)), '_flap_amplitudeX.manual_3mps = []']));
    eval(string(['point_', char(num2str(p)), '_flap_amplitudeY.manual_3mps = []']));
    eval(string(['point_', char(num2str(p)), '_flap_amplitudeZ.manual_3mps = []']));

    eval(string(['point_', char(num2str(p)), '_flap_amplitudeX.manual_4p5mps = []']));
    eval(string(['point_', char(num2str(p)), '_flap_amplitudeY.manual_4p5mps = []']));
    eval(string(['point_', char(num2str(p)), '_flap_amplitudeZ.manual_4p5mps = []']));

    eval(string(['point_', char(num2str(p)), '_flap_amplitudeX.manual_6mps = []']));
    eval(string(['point_', char(num2str(p)), '_flap_amplitudeY.manual_6mps = []']));
    eval(string(['point_', char(num2str(p)), '_flap_amplitudeZ.manual_6mps = []']));

    wbinterplength = 100;

    for trial_number = 1:numel(wingbeats)
        current_wingbeat_XYZ = wingbeats{trial_number}.full_wb_BCS_smooth(wingbeats{trial_number}.begDownStrokeIndexes(1):wingbeats{trial_number}.begDownStrokeIndexes(end),:,:);
        points_to_analyze = wingbeats{trial_number}.points_to_analyze;
        current_amplitudeX = current_wingbeat_XYZ(:,points_to_analyze(p),1) - current_wingbeat_XYZ(:,points_to_analyze(5),1)*0;
        current_amplitudeY = current_wingbeat_XYZ(:,points_to_analyze(p),2) - current_wingbeat_XYZ(:,points_to_analyze(5),2)*0;
        current_amplitudeZ = current_wingbeat_XYZ(:,points_to_analyze(p),3) - current_wingbeat_XYZ(:,points_to_analyze(5),3)*0;

        current_amplitudeX = current_amplitudeX*1000; %m to mm
        current_amplitudeY = current_amplitudeY*1000; %m to mm
        current_amplitudeZ = current_amplitudeZ*1000; %m to mm
        speed_and_treatment = strcat(wingbeats{trial_number}.treatment, "_", num2str(wingbeats{trial_number}.speed_nominal),"mps");
        speed_and_treatment = strrep(speed_and_treatment, ".", "p");

        if wingbeats{trial_number}.points_to_analyze(1) == 1
            current_amplitudeY = -current_amplitudeY;
        end

        eval(['point_', char(num2str(p)), '_flap_amplitudeX.','(speed_and_treatment)(end+1,1:wbinterplength) = interp1(linspace(1,wbinterplength,numel(current_amplitudeX)),current_amplitudeX,1:wbinterplength, ''spline'')']);
        eval(['point_', char(num2str(p)), '_flap_amplitudeY.','(speed_and_treatment)(end+1,1:wbinterplength) = interp1(linspace(1,wbinterplength,numel(current_amplitudeY)),current_amplitudeY,1:wbinterplength, ''spline'')']);
        eval(['point_', char(num2str(p)), '_flap_amplitudeZ.','(speed_and_treatment)(end+1,1:wbinterplength) = interp1(linspace(1,wbinterplength,numel(current_amplitudeZ)),current_amplitudeZ,1:wbinterplength, ''spline'')']);
    end
end
%%

DLC_color = [141 47 162]./256;
human_digitizer_color = [0, 0.4470, 0.7410];
lw=2;
fa=0.5;
n_ticks = 5; %define the number of yticks you want

amplitude_fig=figure;

for p =1:6

    switch p
        case 1 % wingtip
            xheight = [-65, 85];
            yheight = [5, 155];
            zheight = [-120, 140];
        case 2 % wrist
            xheight = [-10, 40];
            yheight = [0, 70];
            zheight = [-40, 70];
        case 3 %5th digit
            xheight = [-65, 5];
            yheight = [10, 80];
            zheight = [-55, 75];
        case 4 % Elbow
            xheight = [-30, 0];
            yheight = [0, 40];
            zheight = [-20, 40];
        case 5 % shoulder
            xheight = [-5, 15];
            yheight = [-5, 25];
            zheight = [-5, 35];
        case 6 % ankle
            xheight = [-70, -50];
            yheight = [-10, 20];
            zheight = [-20, 30];
    end
    max_range = max([range(xheight) range(yheight) range(zheight)]);
    xheight = mean(xheight) + [-max_range max_range]/2;
    yheight = mean(yheight) + [-max_range max_range]/2;
    zheight = mean(zheight) + [-max_range max_range]/2;

    panelWidth = 1 / 2;  % Total width divided by the number of columns
    panelHeight = 1 / 3;  % Total height divided by the number of rows

    switch p
        case 1
            outerRow = 1;
            outerCol = 1;
        case 2
            outerRow = 1;
            outerCol = 2;
        case 3
            outerRow = 2;
            outerCol = 1;
        case 4
            outerRow = 2;
            outerCol = 2;
        case 5
            outerRow = 3;
            outerCol = 1;
        case 6
            outerRow = 3;
            outerCol = 2;
    end

    % Calculate the position of each UI panel
    left = (outerCol - 1) * panelWidth;
    bottom = 1 - outerRow * panelHeight;  % The y-coordinate is reversed for UI panels
    position = [left, bottom, panelWidth, panelHeight];

    % Create the UI panel
    panel = uipanel('Position', position);

    eval(['number_of_wingbeats = size(','point_', char(num2str(p)),'_flap_amplitudeX.auto_3mps,1)']);
    subplot(3,3,1, 'Parent',panel); hold on
    subtitle("3 m/s");
    ylabel("x");
    eval(['[trend_auto, x_std, bin_N, bins, norm_time, x_iqr_auto_x_3mps] = paa(point_',char(num2str(p)),'_flap_amplitudeX.auto_3mps, repmat(1:wbinterplength,number_of_wingbeats,1), wbinterplength, 0);']);
    color = DLC_color;
    uppy_auto = trend_auto+1.96*x_std/sqrt(number_of_wingbeats);
    lowy_auto = trend_auto-1.96*x_std/sqrt(number_of_wingbeats);
    h = shade(0:99, uppy_auto, 0:99, lowy_auto, 'FillType',[1 2;2 1], 'FillColor', color, 'FillAlpha',fa);
    lines = h.findobj("Type", "Line");
    for i = 1:numel(lines)
        delete(lines(i));
    end

    eval(['number_of_wingbeats = size(','point_', char(num2str(p)),'_flap_amplitudeX.manual_3mps,1)']);
    eval(['[trend_manual, x_std, bin_N, bins, norm_time, x_iqr_manual_x_3mps] = paa(point_',char(num2str(p)),'_flap_amplitudeX.manual_3mps, repmat(1:wbinterplength,number_of_wingbeats,1), wbinterplength, 0);']);

    [R,P] = corrcoef(trend_auto,trend_manual);
    correlationCoeff(p,1) = R(1,2);
    correlationP(p,1) = P(1,2);

    trend_auto_x_3mps = trend_auto;
    trend_manual_x_3mps = trend_manual;
    trend_diff_x_3mps = abs(trend_manual-trend_auto);

    color = human_digitizer_color;
    uppy_manual = trend_manual+1.96*x_std/sqrt(number_of_wingbeats);
    lowy_manual = trend_manual-1.96*x_std/sqrt(number_of_wingbeats);
    h = shade(0:99, uppy_manual, 0:99, lowy_manual, 'FillType',[1 2;2 1], 'FillColor', color, 'FillAlpha',fa);
    lines = h.findobj("Type", "Line");
    for i = 1:numel(lines)
        delete(lines(i));
    end
    eval(['overlap.point_' char(num2str(p)) '_X_3mps = (uppy_auto >= lowy_manual) & (uppy_manual >= lowy_auto)']);
    eval(['overlap.point_' char(num2str(p)) '_X_3mps_ratio = sum(overlap.point_' char(num2str(p)) '_X_3mps)']);
    plot(0:99, trend_auto, 'color', DLC_color, 'LineWidth',lw)
    plot(0:99, trend_manual, 'color', human_digitizer_color, 'lineWidth', lw)
    ylim(xheight)
    y_limits = ylim(gca); %get the y limits
    yticks(linspace(y_limits(1),y_limits(2),n_ticks)); %set the yticks
    xticklabels([]);

    eval(['number_of_wingbeats = size(','point_', char(num2str(p)),'_flap_amplitudeX.auto_4p5mps,1)']);
    subplot(3,3,2, 'Parent',panel); hold on
    subtitle("4.5 m/s")
    eval(['[trend_auto, x_std, bin_N, bins, norm_time, x_iqr_auto_x_4p5mps] = paa(point_',char(num2str(p)),'_flap_amplitudeX.auto_4p5mps, repmat(1:wbinterplength,number_of_wingbeats,1), wbinterplength, 0);']);
    color = DLC_color;
    uppy_auto = trend_auto+1.96*x_std/sqrt(number_of_wingbeats);
    lowy_auto = trend_auto-1.96*x_std/sqrt(number_of_wingbeats);
    h = shade(0:99, uppy_auto, 0:99, lowy_auto, 'FillType',[1 2;2 1], 'FillColor', color, 'FillAlpha',fa);
    lines = h.findobj("Type", "Line");
    for i = 1:numel(lines)
        delete(lines(i));
    end

    eval(['number_of_wingbeats = size(','point_', char(num2str(p)),'_flap_amplitudeX.manual_4p5mps,1)']);
    eval(['[trend_manual, x_std, bin_N, bins, norm_time, x_iqr_manual_x_4p5mps] = paa(point_',char(num2str(p)),'_flap_amplitudeX.manual_4p5mps, repmat(1:wbinterplength,number_of_wingbeats,1), wbinterplength, 0);']);

    [R,P] = corrcoef(trend_auto,trend_manual);
    correlationCoeff(p,2) = R(1,2);
    correlationP(p,2) = P(1,2);

    trend_diff_x_4p5mps = abs(trend_manual-trend_auto);

    color = human_digitizer_color;
    uppy_manual = trend_manual+1.96*x_std/sqrt(number_of_wingbeats);
    lowy_manual = trend_manual-1.96*x_std/sqrt(number_of_wingbeats);
    h = shade(0:99, uppy_manual, 0:99, lowy_manual, 'FillType',[1 2;2 1], 'FillColor', color, 'FillAlpha',fa);
    lines = h.findobj("Type", "Line");
    for i = 1:numel(lines)
        delete(lines(i));
    end
    eval(['overlap.point_' char(num2str(p)) '_X_4p5mps = (uppy_auto >= lowy_manual) & (uppy_manual >= lowy_auto)']);
    eval(['overlap.point_' char(num2str(p)) '_X_4p5mps_ratio = sum(overlap.point_' char(num2str(p)) '_X_4p5mps)']);
    plot(0:99, trend_auto, 'color', DLC_color, 'LineWidth',lw)
    plot(0:99, trend_manual, 'color', human_digitizer_color, 'lineWidth', lw)
    xticklabels([])

    ylim(xheight)
    y_limits = ylim(gca); %get the y limits
    yticks(linspace(y_limits(1),y_limits(2),n_ticks)); %set the yticks
    yticklabels([])

    eval(['number_of_wingbeats = size(','point_', char(num2str(p)),'_flap_amplitudeX.auto_6mps,1)']);
    subplot(3,3,3, 'Parent',panel); hold on
    subtitle("6 m/s")
    eval(['[trend_auto, x_std, bin_N, bins, norm_time, x_iqr_auto_x_6mps] = paa(point_',char(num2str(p)),'_flap_amplitudeX.auto_6mps, repmat(1:wbinterplength,number_of_wingbeats,1), wbinterplength, 0);']);
    color = DLC_color;
    uppy_auto = trend_auto+1.96*x_std/sqrt(number_of_wingbeats);
    lowy_auto = trend_auto-1.96*x_std/sqrt(number_of_wingbeats);
    h = shade(0:99, uppy_auto, 0:99, lowy_auto, 'FillType',[1 2;2 1], 'FillColor', color, 'FillAlpha',fa);
    lines = h.findobj("Type", "Line");
    for i = 1:numel(lines)
        delete(lines(i));
    end

    eval(['number_of_wingbeats = size(','point_', char(num2str(p)),'_flap_amplitudeX.manual_6mps,1)']);
    eval(['[trend_manual, x_std, bin_N, bins, norm_time, x_iqr_manual_x_6mps] = paa(point_',char(num2str(p)),'_flap_amplitudeX.manual_6mps, repmat(1:wbinterplength,number_of_wingbeats,1), wbinterplength, 0);']);

    [R,P] = corrcoef(trend_auto,trend_manual);
    correlationCoeff(p,3) = R(1,2);
    correlationP(p,3) = P(1,2);

    trend_diff_x_6mps = abs(trend_manual-trend_auto);

    color = human_digitizer_color;
    uppy_manual = trend_manual+1.96*x_std/sqrt(number_of_wingbeats);
    lowy_manual = trend_manual-1.96*x_std/sqrt(number_of_wingbeats);
    h = shade(0:99, uppy_manual, 0:99, lowy_manual, 'FillType',[1 2;2 1], 'FillColor', color, 'FillAlpha',fa);
    lines = h.findobj("Type", "Line");
    for i = 1:numel(lines)
        delete(lines(i));
    end
    eval(['overlap.point_' char(num2str(p)) '_X_6mps = (uppy_auto >= lowy_manual) & (uppy_manual >= lowy_auto)']);
    eval(['overlap.point_' char(num2str(p)) '_X_6mps_ratio = sum(overlap.point_' char(num2str(p)) '_X_6mps)']);
    plot(0:99, trend_auto, 'color', DLC_color, 'LineWidth',lw)
    plot(0:99, trend_manual, 'color', human_digitizer_color, 'lineWidth', lw)
    xticklabels([])

    ylim(xheight)
    y_limits = ylim(gca); %get the y limits
    yticks(linspace(y_limits(1),y_limits(2),n_ticks)); %set the yticks
    yticklabels([])

    eval(['number_of_wingbeats = size(','point_', char(num2str(p)),'_flap_amplitudeY.auto_3mps,1)']);
    subplot(3,3,4, 'Parent',panel); hold on
    ylabel("y")
    eval(['[trend_auto, x_std, bin_N, bins, norm_time, x_iqr_auto_y_3mps] = paa(point_',char(num2str(p)),'_flap_amplitudeY.auto_3mps, repmat(1:wbinterplength,number_of_wingbeats,1), wbinterplength, 0);']);
    color = DLC_color;
    uppy_auto = trend_auto+1.96*x_std/sqrt(number_of_wingbeats);
    lowy_auto = trend_auto-1.96*x_std/sqrt(number_of_wingbeats);
    h = shade(0:99, uppy_auto, 0:99, lowy_auto, 'FillType',[1 2;2 1], 'FillColor', color, 'FillAlpha',fa);
    lines = h.findobj("Type", "Line");
    for i = 1:numel(lines)
        delete(lines(i));
    end

    eval(['number_of_wingbeats = size(','point_', char(num2str(p)),'_flap_amplitudeY.manual_3mps,1)']);
    eval(['[trend_manual, x_std, bin_N, bins, norm_time, x_iqr_manual_y_3mps] = paa(point_',char(num2str(p)),'_flap_amplitudeY.manual_3mps, repmat(1:wbinterplength,number_of_wingbeats,1), wbinterplength, 0);']);

    [R,P] = corrcoef(trend_auto,trend_manual);
    correlationCoeff(p,4) = R(1,2);
    correlationP(p,4) = P(1,2);

    trend_auto_y_3mps = trend_auto;
    trend_manual_y_3mps = trend_manual;
    trend_diff_y_3mps = abs(trend_manual-trend_auto);

    color = human_digitizer_color;
    uppy_manual = trend_manual+1.96*x_std/sqrt(number_of_wingbeats);
    lowy_manual = trend_manual-1.96*x_std/sqrt(number_of_wingbeats);
    h = shade(0:99, uppy_manual, 0:99, lowy_manual, 'FillType',[1 2;2 1], 'FillColor', color, 'FillAlpha',fa);

    lines = h.findobj("Type", "Line");
    for i = 1:numel(lines)
        delete(lines(i));
    end
    eval(['overlap.point_' char(num2str(p)) '_Y_3mps = (uppy_auto >= lowy_manual) & (uppy_manual >= lowy_auto)']);
    eval(['overlap.point_' char(num2str(p)) '_Y_3mps_ratio = sum(overlap.point_' char(num2str(p)) '_Y_3mps)']);
    plot(0:99, trend_auto, 'color', DLC_color, 'LineWidth',lw)
    plot(0:99, trend_manual, 'color', human_digitizer_color, 'lineWidth', lw)
    xticklabels([])
    ylim(yheight)
    y_limits = ylim(gca); %get the y limits
    yticks(linspace(y_limits(1),y_limits(2),n_ticks)); %set the yticks

    eval(['number_of_wingbeats = size(','point_', char(num2str(p)),'_flap_amplitudeY.auto_4p5mps,1)']);
    subplot(3,3,5, 'Parent',panel); hold on
    eval(['[trend_auto, x_std, bin_N, bins, norm_time, x_iqr_auto_y_4p5mps] = paa(point_',char(num2str(p)),'_flap_amplitudeY.auto_4p5mps, repmat(1:wbinterplength,number_of_wingbeats,1), wbinterplength, 0);']);
    color = DLC_color;
    uppy_auto = trend_auto+1.96*x_std/sqrt(number_of_wingbeats);
    lowy_auto = trend_auto-1.96*x_std/sqrt(number_of_wingbeats);
    h = shade(0:99, uppy_auto, 0:99, lowy_auto, 'FillType',[1 2;2 1], 'FillColor', color, 'FillAlpha',fa);
    lines = h.findobj("Type", "Line");
    for i = 1:numel(lines)
        delete(lines(i));
    end

    eval(['number_of_wingbeats = size(','point_', char(num2str(p)),'_flap_amplitudeY.manual_4p5mps,1)']);
    eval(['[trend_manual, x_std, bin_N, bins, norm_time, x_iqr_manual_y_4p5mps] = paa(point_',char(num2str(p)),'_flap_amplitudeY.manual_4p5mps, repmat(1:wbinterplength,number_of_wingbeats,1), wbinterplength, 0);']);

    [R,P] = corrcoef(trend_auto,trend_manual);
    correlationCoeff(p,5) = R(1,2);
    correlationP(p,5) = P(1,2);

    trend_diff_y_4p5mps = abs(trend_manual-trend_auto);

    color = human_digitizer_color;
    uppy_manual = trend_manual+1.96*x_std/sqrt(number_of_wingbeats);
    lowy_manual = trend_manual-1.96*x_std/sqrt(number_of_wingbeats);
    h = shade(0:99, uppy_manual, 0:99, lowy_manual, 'FillType',[1 2;2 1], 'FillColor', color, 'FillAlpha',fa);
    lines = h.findobj("Type", "Line");
    for i = 1:numel(lines)
        delete(lines(i));
    end
    eval(['overlap.point_' char(num2str(p)) '_Y_4p5mps = (uppy_auto >= lowy_manual) & (uppy_manual >= lowy_auto)']);
    eval(['overlap.point_' char(num2str(p)) '_Y_4p5mps_ratio = sum(overlap.point_' char(num2str(p)) '_Y_4p5mps)']);
    plot(0:99, trend_auto, 'color', DLC_color, 'LineWidth',lw)
    plot(0:99, trend_manual, 'color', human_digitizer_color, 'lineWidth', lw)
    xticklabels([])

    ylim(yheight)
    y_limits = ylim(gca); %get the y limits
    yticks(linspace(y_limits(1),y_limits(2),n_ticks)); %set the yticks
    yticklabels([])

    eval(['number_of_wingbeats = size(','point_', char(num2str(p)),'_flap_amplitudeY.auto_6mps,1)']);
    subplot(3,3,6, 'Parent',panel); hold on
    eval(['[trend_auto, x_std, bin_N, bins, norm_time, x_iqr_auto_y_6mps] = paa(point_',char(num2str(p)),'_flap_amplitudeY.auto_6mps, repmat(1:wbinterplength,number_of_wingbeats,1), wbinterplength, 0);']);
    color = DLC_color;
    uppy_auto = trend_auto+1.96*x_std/sqrt(number_of_wingbeats);
    lowy_auto = trend_auto-1.96*x_std/sqrt(number_of_wingbeats);
    h = shade(0:99, uppy_auto, 0:99, lowy_auto, 'FillType',[1 2;2 1], 'FillColor', color, 'FillAlpha',fa);
    lines = h.findobj("Type", "Line");
    for i = 1:numel(lines)
        delete(lines(i));
    end

    eval(['number_of_wingbeats = size(','point_', char(num2str(p)),'_flap_amplitudeY.manual_6mps,1)']);
    eval(['[trend_manual, x_std, bin_N, bins, norm_time, x_iqr_manual_y_6mps] = paa(point_',char(num2str(p)),'_flap_amplitudeY.manual_6mps, repmat(1:wbinterplength,number_of_wingbeats,1), wbinterplength, 0);']);

    [R,P] = corrcoef(trend_auto,trend_manual);
    correlationCoeff(p,6) = R(1,2);
    correlationP(p,6) = P(1,2);

    trend_diff_y_6mps = abs(trend_manual-trend_auto);

    color = human_digitizer_color;
    uppy_manual = trend_manual+1.96*x_std/sqrt(number_of_wingbeats);
    lowy_manual = trend_manual-1.96*x_std/sqrt(number_of_wingbeats);
    h = shade(0:99, uppy_manual, 0:99, lowy_manual, 'FillType',[1 2;2 1], 'FillColor', color, 'FillAlpha',fa);
    lines = h.findobj("Type", "Line");
    for i = 1:numel(lines)
        delete(lines(i));
    end
    eval(['overlap.point_' char(num2str(p)) '_Y_6mps = (uppy_auto >= lowy_manual) & (uppy_manual >= lowy_auto)']);
    eval(['overlap.point_' char(num2str(p)) '_Y_6mps_ratio = sum(overlap.point_' char(num2str(p)) '_Y_6mps)']);
    plot(0:99, trend_auto, 'color', DLC_color, 'LineWidth',lw)
    plot(0:99, trend_manual, 'color', human_digitizer_color, 'lineWidth', lw)
    xticklabels([])

    ylim(yheight)
    y_limits = ylim(gca); %get the y limits
    yticks(linspace(y_limits(1),y_limits(2),n_ticks)); %set the yticks
    yticklabels([])

    eval(['number_of_wingbeats = size(','point_', char(num2str(p)),'_flap_amplitudeZ.auto_3mps,1)']);
    subplot(3,3,7, 'Parent',panel); hold on

    ylabel("z", 'Color','k')
    eval(['[trend_auto, x_std, bin_N, bins, norm_time, x_iqr_auto_z_3mps] = paa(point_',char(num2str(p)),'_flap_amplitudeZ.auto_3mps, repmat(1:wbinterplength,number_of_wingbeats,1), wbinterplength, 0);']);
    color = DLC_color;
    uppy_auto = trend_auto+1.96*x_std/sqrt(number_of_wingbeats);
    lowy_auto = trend_auto-1.96*x_std/sqrt(number_of_wingbeats);
    h = shade(0:99, uppy_auto, 0:99, lowy_auto, 'FillType',[1 2;2 1], 'FillColor', color, 'FillAlpha',fa);
    lines = h.findobj("Type", "Line");
    for i = 1:numel(lines)
        delete(lines(i));
    end

    eval(['number_of_wingbeats = size(','point_', char(num2str(p)),'_flap_amplitudeZ.manual_3mps,1)']);
    eval(['[trend_manual, x_std, bin_N, bins, norm_time, x_iqr_manual_z_3mps] = paa(point_',char(num2str(p)),'_flap_amplitudeZ.manual_3mps, repmat(1:wbinterplength,number_of_wingbeats,1), wbinterplength, 0);']);

    [R,P] = corrcoef(trend_auto,trend_manual);
    correlationCoeff(p,7) = R(1,2);
    correlationP(p,7) = P(1,2);

    trend_auto_z_3mps = trend_auto;
    trend_manual_z_3mps = trend_manual;
    trend_diff_z_3mps = abs(trend_manual-trend_auto);

    color = human_digitizer_color;
    uppy_manual = trend_manual+1.96*x_std/sqrt(number_of_wingbeats);
    lowy_manual = trend_manual-1.96*x_std/sqrt(number_of_wingbeats);
    h = shade(0:99, uppy_manual, 0:99, lowy_manual, 'FillType',[1 2;2 1], 'FillColor', color, 'FillAlpha',fa);
    lines = h.findobj("Type", "Line");
    for i = 1:numel(lines)
        delete(lines(i));
    end
    eval(['overlap.point_' char(num2str(p)) '_Z_3mps = (uppy_auto >= lowy_manual) & (uppy_manual >= lowy_auto)']);
    eval(['overlap.point_' char(num2str(p)) '_Z_3mps_ratio = sum(overlap.point_' char(num2str(p)) '_Z_3mps)']);
    plot(0:99, trend_auto, 'color', DLC_color, 'LineWidth',lw)
    plot(0:99, trend_manual, 'color', human_digitizer_color, 'lineWidth', lw)
    ylim(zheight);
    y_limits = ylim(gca); %get the y limits
    yticks(linspace(y_limits(1),y_limits(2),n_ticks)); %set the yticks

    eval(['number_of_wingbeats = size(','point_', char(num2str(p)),'_flap_amplitudeZ.auto_4p5mps,1)']);
    subplot(3,3,8, 'Parent',panel); hold on
    eval(['[trend_auto, x_std, bin_N, bins, norm_time, x_iqr_auto_z_4p5mps] = paa(point_',char(num2str(p)),'_flap_amplitudeZ.auto_4p5mps, repmat(1:wbinterplength,number_of_wingbeats,1), wbinterplength, 0);']);
    color = DLC_color;
    uppy_auto = trend_auto+1.96*x_std/sqrt(number_of_wingbeats);
    lowy_auto = trend_auto-1.96*x_std/sqrt(number_of_wingbeats);
    h = shade(0:99, uppy_auto, 0:99, lowy_auto, 'FillType',[1 2;2 1], 'FillColor', color, 'FillAlpha',fa);
    lines = h.findobj("Type", "Line");
    for i = 1:numel(lines)
        delete(lines(i));
    end

    eval(['number_of_wingbeats = size(','point_', char(num2str(p)),'_flap_amplitudeZ.manual_4p5mps,1)']);
    eval(['[trend_manual, x_std, bin_N, bins, norm_time, x_iqr_manual_z_4p5mps] = paa(point_',char(num2str(p)),'_flap_amplitudeZ.manual_4p5mps, repmat(1:wbinterplength,number_of_wingbeats,1), wbinterplength, 0);']);

    [R,P] = corrcoef(trend_auto,trend_manual);
    correlationCoeff(p,8) = R(1,2);
    correlationP(p,8) = P(1,2);

    trend_diff_z_4p5mps = abs(trend_manual-trend_auto);

    color = human_digitizer_color;
    uppy_manual = trend_manual+1.96*x_std/sqrt(number_of_wingbeats);
    lowy_manual = trend_manual-1.96*x_std/sqrt(number_of_wingbeats);
    h = shade(0:99, uppy_manual, 0:99, lowy_manual, 'FillType',[1 2;2 1], 'FillColor', color, 'FillAlpha',fa);
    lines = h.findobj("Type", "Line");
    for i = 1:numel(lines)
        delete(lines(i));
    end
    eval(['overlap.point_' char(num2str(p)) '_Z_4p5mps = (uppy_auto >= lowy_manual) & (uppy_manual >= lowy_auto)']);
    eval(['overlap.point_' char(num2str(p)) '_Z_4p5mps_ratio = sum(overlap.point_' char(num2str(p)) '_Z_4p5mps)']);
    plot(0:99, trend_auto, 'color', DLC_color, 'LineWidth',lw)
    plot(0:99, trend_manual, 'color', human_digitizer_color, 'lineWidth', lw)

    ylim(zheight);
    y_limits = ylim(gca); %get the y limits
    yticks(linspace(y_limits(1),y_limits(2),n_ticks)); %set the yticks
    yticklabels([])

    eval(['number_of_wingbeats = size(','point_', char(num2str(p)),'_flap_amplitudeZ.auto_6mps,1)']);
    subplot(3,3,9, 'Parent',panel); hold on
    eval(['[trend_auto, x_std, bin_N, bins, norm_time, x_iqr_auto_z_6mps] = paa(point_',char(num2str(p)),'_flap_amplitudeZ.auto_6mps, repmat(1:wbinterplength,number_of_wingbeats,1), wbinterplength, 0);']);
    color = DLC_color;
    uppy_auto = trend_auto+1.96*x_std/sqrt(number_of_wingbeats);
    lowy_auto = trend_auto-1.96*x_std/sqrt(number_of_wingbeats);
    h = shade(0:99, uppy_auto, 0:99, lowy_auto, 'FillType',[1 2;2 1], 'FillColor', color, 'FillAlpha',fa);
    lines = h.findobj("Type", "Line");
    for i = 1:numel(lines)
        delete(lines(i));
    end

    eval(['number_of_wingbeats = size(','point_', char(num2str(p)),'_flap_amplitudeZ.manual_6mps,1)']);
    eval(['[trend_manual, x_std, bin_N, bins, norm_time, x_iqr_manual_z_6mps] = paa(point_',char(num2str(p)),'_flap_amplitudeZ.manual_6mps, repmat(1:wbinterplength,number_of_wingbeats,1), wbinterplength, 0);']);

    [R,P] = corrcoef(trend_auto,trend_manual);
    correlationCoeff(p,9) = R(1,2);
    correlationP(p,9) = P(1,2);

    trend_diff_z_6mps = abs(trend_manual-trend_auto);

    color = human_digitizer_color;
    uppy_manual = trend_manual+1.96*x_std/sqrt(number_of_wingbeats);
    lowy_manual = trend_manual-1.96*x_std/sqrt(number_of_wingbeats);
    h = shade(0:99, uppy_manual, 0:99, lowy_manual, 'FillType',[1 2;2 1], 'FillColor', color, 'FillAlpha',fa);
    lines = h.findobj("Type", "Line");
    for i = 1:numel(lines)
        delete(lines(i));
    end
    eval(['overlap.point_' char(num2str(p)) '_Z_6mps = (uppy_auto >= lowy_manual) & (uppy_manual >= lowy_auto)']);
    eval(['overlap.point_' char(num2str(p)) '_Z_6mps_ratio = sum(overlap.point_' char(num2str(p)) '_Z_6mps)']);
    plot(0:99, trend_auto, 'color', DLC_color, 'LineWidth',lw)
    plot(0:99, trend_manual, 'color', human_digitizer_color, 'lineWidth', lw)

    ylim(zheight);
    y_limits = ylim(gca); %get the y limits
    yticks(linspace(y_limits(1),y_limits(2),n_ticks)); %set the yticks
    yticklabels([])

    eval(['iqr_auto_x_3mps_point' char(num2str(p)) ' =  x_iqr_auto_x_3mps']);
    eval(['iqr_auto_y_3mps_point' char(num2str(p)) ' =  x_iqr_auto_y_3mps']);
    eval(['iqr_auto_z_3mps_point' char(num2str(p)) ' =  x_iqr_auto_z_3mps']);

    eval(['iqr_auto_x_4p5mps_point' char(num2str(p)) ' =  x_iqr_auto_x_4p5mps']);
    eval(['iqr_auto_y_4p5mps_point' char(num2str(p)) ' =  x_iqr_auto_y_4p5mps']);
    eval(['iqr_auto_z_4p5mps_point' char(num2str(p)) ' =  x_iqr_auto_z_4p5mps']);

    eval(['iqr_auto_x_6mps_point' char(num2str(p)) ' =  x_iqr_auto_x_6mps']);
    eval(['iqr_auto_y_6mps_point' char(num2str(p)) ' =  x_iqr_auto_y_6mps']);
    eval(['iqr_auto_z_6mps_point' char(num2str(p)) ' =  x_iqr_auto_z_6mps']);

    eval(['iqr_manual_x_3mps_point' char(num2str(p)) ' =  x_iqr_manual_x_3mps']);
    eval(['iqr_manual_y_3mps_point' char(num2str(p)) ' =  x_iqr_manual_y_3mps']);
    eval(['iqr_manual_z_3mps_point' char(num2str(p)) ' =  x_iqr_manual_z_3mps']);

    eval(['iqr_manual_x_4p5mps_point' char(num2str(p)) ' =  x_iqr_manual_x_4p5mps']);
    eval(['iqr_manual_y_4p5mps_point' char(num2str(p)) ' =  x_iqr_manual_y_4p5mps']);
    eval(['iqr_manual_z_4p5mps_point' char(num2str(p)) ' =  x_iqr_manual_z_4p5mps']);

    eval(['iqr_manual_x_6mps_point' char(num2str(p)) ' =  x_iqr_manual_x_6mps']);
    eval(['iqr_manual_y_6mps_point' char(num2str(p)) ' =  x_iqr_manual_y_6mps']);
    eval(['iqr_manual_z_6mps_point' char(num2str(p)) ' =  x_iqr_manual_z_6mps']);
end

%%
frameRate = 700;
% Define speeds and treatments
speeds = [3, 4.5, 6];
treatments = {'automatic', 'manual'};

% Initialize empty cell arrays for each parameter
wingbeat_frequency = cell(length(speeds), length(treatments));
stroke_amplitude = cell(length(speeds), length(treatments));
stroke_plane_angle_ds = cell(length(speeds), length(treatments));
stroke_plane_angle_us = cell(length(speeds), length(treatments));
mid_downstroke_speed = cell(length(speeds), length(treatments));
average_wingtip_speed = cell(length(speeds), length(treatments));
angle_of_attack = cell(length(speeds), length(treatments));
wing_area_midDS = cell(length(speeds), length(treatments));
individual = cell(length(speeds), length(treatments));
num_frames = 0;

% figure; hold on; axis equal

for trial_number = 1:numel(wingbeats)
    current_wingbeat_XYZ = wingbeats{trial_number}.full_wb_BCS_smooth(wingbeats{trial_number}.begDownStrokeIndexes(1):wingbeats{trial_number}.begDownStrokeIndexes(end),:,:);
    % whole_wingbeats{trial_number} = current_wingbeat_XYZ;

    switch wingbeats{trial_number}.treatment
        case "auto"
            treatment_I = 1;
            % plot_color = 'r'
        case "manual"
            treatment_I = 2;
            % plot_color = 'b'
    end

    switch wingbeats{trial_number}.speed_nominal
        case 3
            speed_I = 1;
        case 4.5
            speed_I = 2;
        case 6
            speed_I = 3;
    end

    points_to_analyze = wingbeats{trial_number}.points_to_analyze;
    current_amplitudeZ = current_wingbeat_XYZ(:,points_to_analyze(1),3); % wingtip amplitude
    current_amplitudeZ = current_amplitudeZ*1000; %m to mm

    [max_amplitude, max_I] = max(current_amplitudeZ, [],'omitnan');
    [min_amplitude, min_I] = min(current_amplitudeZ, [],'omitnan');

    stroke_amplitude_current = max_amplitude-min_amplitude;
    stroke_amplitude{speed_I,treatment_I} = [stroke_amplitude{speed_I,treatment_I}, stroke_amplitude_current];

    wingbeat_frequency_current = 1/(size(current_wingbeat_XYZ,1)*1/frameRate);
    wingbeat_frequency{speed_I,treatment_I} = [wingbeat_frequency{speed_I,treatment_I}, wingbeat_frequency_current];
    if wingbeats{trial_number}.treatment == "auto"
        num_frames = num_frames + size(current_wingbeat_XYZ,1);
    end
    current_downstroke = current_wingbeat_XYZ(1:min_I,:,:);
    current_upstroke = current_wingbeat_XYZ(min_I+1:end,:,:);

    [~, ~, ~, stroke_plane_angle_ds_current] = strokePlaneVectors((current_downstroke(:,points_to_analyze(1),:)));
    [~, ~, ~, stroke_plane_angle_us_current] = strokePlaneVectors((current_upstroke(:,points_to_analyze(1),:)));
    stroke_plane_angle_ds_current = stroke_plane_angle_ds_current(2)*57.29578;
    stroke_plane_angle_us_current = stroke_plane_angle_us_current(2)*57.29578;

    stroke_plane_angle_ds{speed_I,treatment_I} = [stroke_plane_angle_ds{speed_I,treatment_I}, stroke_plane_angle_ds_current];
    stroke_plane_angle_us{speed_I,treatment_I} = [stroke_plane_angle_us{speed_I,treatment_I}, stroke_plane_angle_us_current];

    wt_to_wt_vector = squeeze(current_wingbeat_XYZ(:,15,:)-current_wingbeat_XYZ(:,1,:));

    if(sum(isnan(wt_to_wt_vector),'all') > 0)
        n_nanpoints = char(num2str(sum(isnan(wt_to_wt_vector),'all')/3));
        for i=1:3
            wt_to_wt_vector(:,i) = fillmissing(wt_to_wt_vector(:,i),"linear");
        end
        warning(['Input points contained ', n_nanpoints, ' nan points. These have been filled in using linear interpolation of neighboring, nonmissing values.'])
        disp(['Trial number: ' char(num2str(trial_number))])
    end

    wt_to_wt_length = vecnorm(wt_to_wt_vector',2);

    [~, midDS_index] = max(wt_to_wt_length);

    angle_of_attack_current = ...
        atand((current_wingbeat_XYZ(midDS_index,points_to_analyze(2),3)-current_wingbeat_XYZ(midDS_index,points_to_analyze(3),3))...
        /...
        (current_wingbeat_XYZ(midDS_index,points_to_analyze(2),1)-current_wingbeat_XYZ(midDS_index,points_to_analyze(3),1)));
    if(false && angle_of_attack_current<-10)
        figure
        better_plot3(current_wingbeat_XYZ); hold on; axis equal
        line([current_wingbeat_XYZ(midDS_index,points_to_analyze(2),1) current_wingbeat_XYZ(midDS_index,points_to_analyze(3),1)],[current_wingbeat_XYZ(midDS_index,points_to_analyze(2),2) current_wingbeat_XYZ(midDS_index,points_to_analyze(3),2)],[current_wingbeat_XYZ(midDS_index,points_to_analyze(2),3) current_wingbeat_XYZ(midDS_index,points_to_analyze(3),3)], 'Color', 'k')
        view(0,0)
        subtitle(trial_number)
        disp(angle_of_attack_current)
    end

    angle_of_attack{speed_I,treatment_I} = [angle_of_attack{speed_I,treatment_I}, angle_of_attack_current];

    [average_wingtip_speed_current, ~] = getVelocityAndAcceleration(current_wingbeat_XYZ(:,points_to_analyze(1),:),frameRate,false);
    average_wingtip_speed_current = mean(vecnorm(average_wingtip_speed_current',2),'omitnan');
    average_wingtip_speed{speed_I,treatment_I} = [average_wingtip_speed{speed_I,treatment_I}, average_wingtip_speed_current];

    current_hand_area = ...
        triangle_area(...
        current_wingbeat_XYZ(midDS_index,points_to_analyze(1), :), ...
        current_wingbeat_XYZ(midDS_index,points_to_analyze(2), :), ...
        current_wingbeat_XYZ(midDS_index,points_to_analyze(3), :));
    current_arm_area = ...
        triangle_area(...
        current_wingbeat_XYZ(midDS_index,points_to_analyze(2), :), ...
        current_wingbeat_XYZ(midDS_index,points_to_analyze(3), :), ...
        current_wingbeat_XYZ(midDS_index,points_to_analyze(4), :)) + ...
        triangle_area(...
        current_wingbeat_XYZ(midDS_index,points_to_analyze(2), :), ...
        current_wingbeat_XYZ(midDS_index,points_to_analyze(4), :), ...
        current_wingbeat_XYZ(midDS_index,points_to_analyze(5), :)) + ...
        triangle_area(...
        current_wingbeat_XYZ(midDS_index,points_to_analyze(3), :), ...
        current_wingbeat_XYZ(midDS_index,points_to_analyze(4), :), ...
        current_wingbeat_XYZ(midDS_index,points_to_analyze(6), :)) + ...
        triangle_area(...
        current_wingbeat_XYZ(midDS_index,points_to_analyze(4), :), ...
        current_wingbeat_XYZ(midDS_index,points_to_analyze(5), :), ...
        current_wingbeat_XYZ(midDS_index,points_to_analyze(6), :));
    current_wing_area = 2*(current_hand_area+current_arm_area);

    wing_area_midDS{speed_I,treatment_I} = [wing_area_midDS{speed_I,treatment_I}, current_wing_area*10000];

    individual{speed_I,treatment_I} = [individual{speed_I,treatment_I}, wingbeats{trial_number}.individual];
end

%plot stroke amplitude
% Create a tiled layout
kinparamfig = figure;
layout = tiledlayout(3, 2); % Adjust the numbers according to your desired layout (rows x columns)

% Call the plot_wing_kinematic function for each parameter
plot_wing_kinematic(stroke_amplitude, speeds, DLC_color, human_digitizer_color, 'Stroke Amplitude (mm)', nexttile(layout));
xlabel([]);
xticks()
xticklabels([]);
legend off
plot_wing_kinematic(wingbeat_frequency, speeds, DLC_color, human_digitizer_color, 'Wingbeat frequency (Hz)',nexttile(layout));
xlabel([]);
xticklabels([]);
legend off
plot_wing_kinematic(stroke_plane_angle_ds, speeds, DLC_color, human_digitizer_color, ['Stroke plane angle DS (' char(176) ')'],nexttile(layout));
xlabel([]);
xticklabels([]);
ylims_sp = get(gca,'ylim');
legend off
plot_wing_kinematic(stroke_plane_angle_us, speeds, DLC_color, human_digitizer_color, ['Stroke plane angle US (' char(176) ')'], nexttile(layout));
xlabel([]);
xticklabels([]);
ylim(ylims_sp)
legend off
plot_wing_kinematic(angle_of_attack, speeds, DLC_color, human_digitizer_color, ['Angle of attack mid DS (' char(176) ')'] ,nexttile(layout));
legend off
xlabel([]);
plot_wing_kinematic(wing_area_midDS, speeds, DLC_color, human_digitizer_color, 'Wing area mid DS (cm^3)',nexttile(layout));
xlabel([]);

% Customize the layout
layout.TileSpacing = 'compact';
layout.Padding = 'compact';
xlabel(layout,"Speed (m/s)")

%%
fnames = fieldnames(overlap);
fnames(contains(fnames,"ratio"))=[];
allOverlap = [];
overlapRatio = {};
for i = 1:numel(fnames)
    fname = fnames(i);
    allOverlap = [allOverlap overlap.(string(fname))];
end
overlapRatioTotal = sum(allOverlap)/numel(allOverlap)

fnames = fieldnames(overlap);
fnames(~contains(fnames,"ratio"))=[];
overlapRatio = {};
for i = 1:numel(fnames)
    fname = fnames(i);
    overlapRatio.(string(fname)) = overlap.(string(fname));
end

%%
diff_iqr_x_3mps_point1 = mean(iqr_manual_x_3mps_point1-iqr_auto_x_3mps_point1);
diff_iqr_y_3mps_point1 = mean(iqr_manual_y_3mps_point1-iqr_auto_y_3mps_point1);
diff_iqr_z_3mps_point1 = mean(iqr_manual_z_3mps_point1-iqr_auto_z_3mps_point1);

diff_iqr_x_4p5mps_point1 = mean(iqr_manual_x_4p5mps_point1-iqr_auto_x_4p5mps_point1);
diff_iqr_y_4p5mps_point1 = mean(iqr_manual_y_4p5mps_point1-iqr_auto_y_4p5mps_point1);
diff_iqr_z_4p5mps_point1 = mean(iqr_manual_z_4p5mps_point1-iqr_auto_z_4p5mps_point1);

diff_iqr_x_6mps_point1 = mean(iqr_manual_x_6mps_point1-iqr_auto_x_6mps_point1);
diff_iqr_y_6mps_point1 = mean(iqr_manual_y_6mps_point1-iqr_auto_y_6mps_point1);
diff_iqr_z_6mps_point1 = mean(iqr_manual_z_6mps_point1-iqr_auto_z_6mps_point1);

diff_iqr_x_3mps_point2 = mean(iqr_manual_x_3mps_point2-iqr_auto_x_3mps_point2);
diff_iqr_y_3mps_point2 = mean(iqr_manual_y_3mps_point2-iqr_auto_y_3mps_point2);
diff_iqr_z_3mps_point2 = mean(iqr_manual_z_3mps_point2-iqr_auto_z_3mps_point2);

diff_iqr_x_4p5mps_point2 = mean(iqr_manual_x_4p5mps_point2-iqr_auto_x_4p5mps_point2);
diff_iqr_y_4p5mps_point2 = mean(iqr_manual_y_4p5mps_point2-iqr_auto_y_4p5mps_point2);
diff_iqr_z_4p5mps_point2 = mean(iqr_manual_z_4p5mps_point2-iqr_auto_z_4p5mps_point2);

diff_iqr_x_6mps_point2 = mean(iqr_manual_x_6mps_point2-iqr_auto_x_6mps_point2);
diff_iqr_y_6mps_point2 = mean(iqr_manual_y_6mps_point2-iqr_auto_y_6mps_point2);
diff_iqr_z_6mps_point2 = mean(iqr_manual_z_6mps_point2-iqr_auto_z_6mps_point2);

diff_iqr_x_3mps_point3 = mean(iqr_manual_x_3mps_point3-iqr_auto_x_3mps_point3);
diff_iqr_y_3mps_point3 = mean(iqr_manual_y_3mps_point3-iqr_auto_y_3mps_point3);
diff_iqr_z_3mps_point3 = mean(iqr_manual_z_3mps_point3-iqr_auto_z_3mps_point3);

diff_iqr_x_4p5mps_point3 = mean(iqr_manual_x_4p5mps_point3-iqr_auto_x_4p5mps_point3);
diff_iqr_y_4p5mps_point3 = mean(iqr_manual_y_4p5mps_point3-iqr_auto_y_4p5mps_point3);
diff_iqr_z_4p5mps_point3 = mean(iqr_manual_z_4p5mps_point3-iqr_auto_z_4p5mps_point3);

diff_iqr_x_6mps_point3 = mean(iqr_manual_x_6mps_point3-iqr_auto_x_6mps_point3);
diff_iqr_y_6mps_point3 = mean(iqr_manual_y_6mps_point3-iqr_auto_y_6mps_point3);
diff_iqr_z_6mps_point3 = mean(iqr_manual_z_6mps_point3-iqr_auto_z_6mps_point3);

diff_iqr_x_3mps_point4 = mean(iqr_manual_x_3mps_point4-iqr_auto_x_3mps_point4);
diff_iqr_y_3mps_point4 = mean(iqr_manual_y_3mps_point4-iqr_auto_y_3mps_point4);
diff_iqr_z_3mps_point4 = mean(iqr_manual_z_3mps_point4-iqr_auto_z_3mps_point4);

diff_iqr_x_4p5mps_point4 = mean(iqr_manual_x_4p5mps_point4-iqr_auto_x_4p5mps_point4);
diff_iqr_y_4p5mps_point4 = mean(iqr_manual_y_4p5mps_point4-iqr_auto_y_4p5mps_point4);
diff_iqr_z_4p5mps_point4 = mean(iqr_manual_z_4p5mps_point4-iqr_auto_z_4p5mps_point4);

diff_iqr_x_6mps_point4 = mean(iqr_manual_x_6mps_point4-iqr_auto_x_6mps_point4);
diff_iqr_y_6mps_point4 = mean(iqr_manual_y_6mps_point4-iqr_auto_y_6mps_point4);
diff_iqr_z_6mps_point4 = mean(iqr_manual_z_6mps_point4-iqr_auto_z_6mps_point4);

diff_iqr_x_3mps_point5 = mean(iqr_manual_x_3mps_point5-iqr_auto_x_3mps_point5);
diff_iqr_y_3mps_point5 = mean(iqr_manual_y_3mps_point5-iqr_auto_y_3mps_point5);
diff_iqr_z_3mps_point5 = mean(iqr_manual_z_3mps_point5-iqr_auto_z_3mps_point5);

diff_iqr_x_4p5mps_point5 = mean(iqr_manual_x_4p5mps_point5-iqr_auto_x_4p5mps_point5);
diff_iqr_y_4p5mps_point5 = mean(iqr_manual_y_4p5mps_point5-iqr_auto_y_4p5mps_point5);
diff_iqr_z_4p5mps_point5 = mean(iqr_manual_z_4p5mps_point5-iqr_auto_z_4p5mps_point5);

diff_iqr_x_6mps_point5 = mean(iqr_manual_x_6mps_point5-iqr_auto_x_6mps_point5);
diff_iqr_y_6mps_point5 = mean(iqr_manual_y_6mps_point5-iqr_auto_y_6mps_point5);
diff_iqr_z_6mps_point5 = mean(iqr_manual_z_6mps_point5-iqr_auto_z_6mps_point5);

diff_iqr_x_3mps_point6 = mean(iqr_manual_x_3mps_point6-iqr_auto_x_3mps_point6);
diff_iqr_y_3mps_point6 = mean(iqr_manual_y_3mps_point6-iqr_auto_y_3mps_point6);
diff_iqr_z_3mps_point6 = mean(iqr_manual_z_3mps_point6-iqr_auto_z_3mps_point6);

diff_iqr_x_4p5mps_point6 = mean(iqr_manual_x_4p5mps_point6-iqr_auto_x_4p5mps_point6);
diff_iqr_y_4p5mps_point6 = mean(iqr_manual_y_4p5mps_point6-iqr_auto_y_4p5mps_point6);
diff_iqr_z_4p5mps_point6 = mean(iqr_manual_z_4p5mps_point6-iqr_auto_z_4p5mps_point6);

diff_iqr_x_6mps_point6 = mean(iqr_manual_x_6mps_point6-iqr_auto_x_6mps_point6);
diff_iqr_y_6mps_point6 = mean(iqr_manual_y_6mps_point6-iqr_auto_y_6mps_point6);
diff_iqr_z_6mps_point6 = mean(iqr_manual_z_6mps_point6-iqr_auto_z_6mps_point6);

CI_diff_iqr_x_3mps_point1 = [diff_iqr_x_3mps_point1-1.96*std(iqr_manual_x_3mps_point1-iqr_auto_x_3mps_point1)/sqrt(numel(iqr_manual_x_3mps_point1)), ...
    diff_iqr_x_3mps_point1+1.96*std(iqr_manual_x_3mps_point1-iqr_auto_x_3mps_point1)/sqrt(numel(iqr_manual_x_3mps_point1))];
CI_diff_iqr_y_3mps_point1 = [diff_iqr_y_3mps_point1-1.96*std(iqr_manual_y_3mps_point1-iqr_auto_y_3mps_point1)/sqrt(numel(iqr_manual_y_3mps_point1)), ...
    diff_iqr_y_3mps_point1+1.96*std(iqr_manual_y_3mps_point1-iqr_auto_y_3mps_point1)/sqrt(numel(iqr_manual_y_3mps_point1))];
CI_diff_iqr_z_3mps_point1 = [diff_iqr_z_3mps_point1-1.96*std(iqr_manual_z_3mps_point1-iqr_auto_z_3mps_point1)/sqrt(numel(iqr_manual_z_3mps_point1)), ...
    diff_iqr_z_3mps_point1+1.96*std(iqr_manual_z_3mps_point1-iqr_auto_z_3mps_point1)/sqrt(numel(iqr_manual_z_3mps_point1))];

CI_diff_iqr_x_4p5mps_point1 = [diff_iqr_x_4p5mps_point1-1.96*std(iqr_manual_x_4p5mps_point1-iqr_auto_x_4p5mps_point1)/sqrt(numel(iqr_manual_x_4p5mps_point1)), ...
    diff_iqr_x_4p5mps_point1+1.96*std(iqr_manual_x_4p5mps_point1-iqr_auto_x_4p5mps_point1)/sqrt(numel(iqr_manual_x_4p5mps_point1))];
CI_diff_iqr_y_4p5mps_point1 = [diff_iqr_y_4p5mps_point1-1.96*std(iqr_manual_y_4p5mps_point1-iqr_auto_y_4p5mps_point1)/sqrt(numel(iqr_manual_y_4p5mps_point1)), ...
    diff_iqr_y_4p5mps_point1+1.96*std(iqr_manual_y_4p5mps_point1-iqr_auto_y_4p5mps_point1)/sqrt(numel(iqr_manual_y_4p5mps_point1))];
CI_diff_iqr_z_4p5mps_point1 = [diff_iqr_z_4p5mps_point1-1.96*std(iqr_manual_z_4p5mps_point1-iqr_auto_z_4p5mps_point1)/sqrt(numel(iqr_manual_z_4p5mps_point1)), ...
    diff_iqr_z_4p5mps_point1+1.96*std(iqr_manual_z_4p5mps_point1-iqr_auto_z_4p5mps_point1)/sqrt(numel(iqr_manual_z_4p5mps_point1))];

CI_diff_iqr_x_6mps_point1 = [diff_iqr_x_6mps_point1-1.96*std(iqr_manual_x_6mps_point1-iqr_auto_x_6mps_point1)/sqrt(numel(iqr_manual_x_6mps_point1)), ...
    diff_iqr_x_6mps_point1+1.96*std(iqr_manual_x_6mps_point1-iqr_auto_x_6mps_point1)/sqrt(numel(iqr_manual_x_6mps_point1))];
CI_diff_iqr_y_6mps_point1 = [diff_iqr_y_6mps_point1-1.96*std(iqr_manual_y_6mps_point1-iqr_auto_y_6mps_point1)/sqrt(numel(iqr_manual_y_6mps_point1)), ...
    diff_iqr_y_6mps_point1+1.96*std(iqr_manual_y_6mps_point1-iqr_auto_y_6mps_point1)/sqrt(numel(iqr_manual_y_6mps_point1))];
CI_diff_iqr_z_6mps_point1 = [diff_iqr_z_6mps_point1-1.96*std(iqr_manual_z_6mps_point1-iqr_auto_z_6mps_point1)/sqrt(numel(iqr_manual_z_6mps_point1)), ...
    diff_iqr_z_6mps_point1+1.96*std(iqr_manual_z_6mps_point1-iqr_auto_z_6mps_point1)/sqrt(numel(iqr_manual_z_6mps_point1))];

CI_diff_iqr_x_3mps_point2 = [diff_iqr_x_3mps_point2-1.96*std(iqr_manual_x_3mps_point2-iqr_auto_x_3mps_point2)/sqrt(numel(iqr_manual_x_3mps_point2)), ...
    diff_iqr_x_3mps_point2+1.96*std(iqr_manual_x_3mps_point2-iqr_auto_x_3mps_point2)/sqrt(numel(iqr_manual_x_3mps_point2))];
CI_diff_iqr_y_3mps_point2 = [diff_iqr_y_3mps_point2-1.96*std(iqr_manual_y_3mps_point2-iqr_auto_y_3mps_point2)/sqrt(numel(iqr_manual_y_3mps_point2)), ...
    diff_iqr_y_3mps_point2+1.96*std(iqr_manual_y_3mps_point2-iqr_auto_y_3mps_point2)/sqrt(numel(iqr_manual_y_3mps_point2))];
CI_diff_iqr_z_3mps_point2 = [diff_iqr_z_3mps_point2-1.96*std(iqr_manual_z_3mps_point2-iqr_auto_z_3mps_point2)/sqrt(numel(iqr_manual_z_3mps_point2)), ...
    diff_iqr_z_3mps_point2+1.96*std(iqr_manual_z_3mps_point2-iqr_auto_z_3mps_point2)/sqrt(numel(iqr_manual_z_3mps_point2))];

CI_diff_iqr_x_4p5mps_point2 = [diff_iqr_x_4p5mps_point2-1.96*std(iqr_manual_x_4p5mps_point2-iqr_auto_x_4p5mps_point2)/sqrt(numel(iqr_manual_x_4p5mps_point2)), ...
    diff_iqr_x_4p5mps_point2+1.96*std(iqr_manual_x_4p5mps_point2-iqr_auto_x_4p5mps_point2)/sqrt(numel(iqr_manual_x_4p5mps_point2))];
CI_diff_iqr_y_4p5mps_point2 = [diff_iqr_y_4p5mps_point2-1.96*std(iqr_manual_y_4p5mps_point2-iqr_auto_y_4p5mps_point2)/sqrt(numel(iqr_manual_y_4p5mps_point2)), ...
    diff_iqr_y_4p5mps_point2+1.96*std(iqr_manual_y_4p5mps_point2-iqr_auto_y_4p5mps_point2)/sqrt(numel(iqr_manual_y_4p5mps_point2))];
CI_diff_iqr_z_4p5mps_point2 = [diff_iqr_z_4p5mps_point2-1.96*std(iqr_manual_z_4p5mps_point2-iqr_auto_z_4p5mps_point2)/sqrt(numel(iqr_manual_z_4p5mps_point2)), ...
    diff_iqr_z_4p5mps_point2+1.96*std(iqr_manual_z_4p5mps_point2-iqr_auto_z_4p5mps_point2)/sqrt(numel(iqr_manual_z_4p5mps_point2))];

CI_diff_iqr_x_6mps_point2 = [diff_iqr_x_6mps_point2-1.96*std(iqr_manual_x_6mps_point2-iqr_auto_x_6mps_point2)/sqrt(numel(iqr_manual_x_6mps_point2)), ...
    diff_iqr_x_6mps_point2+1.96*std(iqr_manual_x_6mps_point2-iqr_auto_x_6mps_point2)/sqrt(numel(iqr_manual_x_6mps_point2))];
CI_diff_iqr_y_6mps_point2 = [diff_iqr_y_6mps_point2-1.96*std(iqr_manual_y_6mps_point2-iqr_auto_y_6mps_point2)/sqrt(numel(iqr_manual_y_6mps_point2)), ...
    diff_iqr_y_6mps_point2+1.96*std(iqr_manual_y_6mps_point2-iqr_auto_y_6mps_point2)/sqrt(numel(iqr_manual_y_6mps_point2))];
CI_diff_iqr_z_6mps_point2 = [diff_iqr_z_6mps_point2-1.96*std(iqr_manual_z_6mps_point2-iqr_auto_z_6mps_point2)/sqrt(numel(iqr_manual_z_6mps_point2)), ...
    diff_iqr_z_6mps_point2+1.96*std(iqr_manual_z_6mps_point2-iqr_auto_z_6mps_point2)/sqrt(numel(iqr_manual_z_6mps_point2))];

CI_diff_iqr_x_3mps_point3 = [diff_iqr_x_3mps_point3-1.96*std(iqr_manual_x_3mps_point3-iqr_auto_x_3mps_point3)/sqrt(numel(iqr_manual_x_3mps_point3)), ...
    diff_iqr_x_3mps_point3+1.96*std(iqr_manual_x_3mps_point3-iqr_auto_x_3mps_point3)/sqrt(numel(iqr_manual_x_3mps_point3))];
CI_diff_iqr_y_3mps_point3 = [diff_iqr_y_3mps_point3-1.96*std(iqr_manual_y_3mps_point3-iqr_auto_y_3mps_point3)/sqrt(numel(iqr_manual_y_3mps_point3)), ...
    diff_iqr_y_3mps_point3+1.96*std(iqr_manual_y_3mps_point3-iqr_auto_y_3mps_point3)/sqrt(numel(iqr_manual_y_3mps_point3))];
CI_diff_iqr_z_3mps_point3 = [diff_iqr_z_3mps_point3-1.96*std(iqr_manual_z_3mps_point3-iqr_auto_z_3mps_point3)/sqrt(numel(iqr_manual_z_3mps_point3)), ...
    diff_iqr_z_3mps_point3+1.96*std(iqr_manual_z_3mps_point3-iqr_auto_z_3mps_point3)/sqrt(numel(iqr_manual_z_3mps_point3))];

CI_diff_iqr_x_4p5mps_point3 = [diff_iqr_x_4p5mps_point3-1.96*std(iqr_manual_x_4p5mps_point3-iqr_auto_x_4p5mps_point3)/sqrt(numel(iqr_manual_x_4p5mps_point3)), ...
    diff_iqr_x_4p5mps_point3+1.96*std(iqr_manual_x_4p5mps_point3-iqr_auto_x_4p5mps_point3)/sqrt(numel(iqr_manual_x_4p5mps_point3))];
CI_diff_iqr_y_4p5mps_point3 = [diff_iqr_y_4p5mps_point3-1.96*std(iqr_manual_y_4p5mps_point3-iqr_auto_y_4p5mps_point3)/sqrt(numel(iqr_manual_y_4p5mps_point3)), ...
    diff_iqr_y_4p5mps_point3+1.96*std(iqr_manual_y_4p5mps_point3-iqr_auto_y_4p5mps_point3)/sqrt(numel(iqr_manual_y_4p5mps_point3))];
CI_diff_iqr_z_4p5mps_point3 = [diff_iqr_z_4p5mps_point3-1.96*std(iqr_manual_z_4p5mps_point3-iqr_auto_z_4p5mps_point3)/sqrt(numel(iqr_manual_z_4p5mps_point3)), ...
    diff_iqr_z_4p5mps_point3+1.96*std(iqr_manual_z_4p5mps_point3-iqr_auto_z_4p5mps_point3)/sqrt(numel(iqr_manual_z_4p5mps_point3))];

CI_diff_iqr_x_6mps_point3 = [diff_iqr_x_6mps_point3-1.96*std(iqr_manual_x_6mps_point3-iqr_auto_x_6mps_point3)/sqrt(numel(iqr_manual_x_6mps_point3)), ...
    diff_iqr_x_6mps_point3+1.96*std(iqr_manual_x_6mps_point3-iqr_auto_x_6mps_point3)/sqrt(numel(iqr_manual_x_6mps_point3))];
CI_diff_iqr_y_6mps_point3 = [diff_iqr_y_6mps_point3-1.96*std(iqr_manual_y_6mps_point3-iqr_auto_y_6mps_point3)/sqrt(numel(iqr_manual_y_6mps_point3)), ...
    diff_iqr_y_6mps_point3+1.96*std(iqr_manual_y_6mps_point3-iqr_auto_y_6mps_point3)/sqrt(numel(iqr_manual_y_6mps_point3))];
CI_diff_iqr_z_6mps_point3 = [diff_iqr_z_6mps_point3-1.96*std(iqr_manual_z_6mps_point3-iqr_auto_z_6mps_point3)/sqrt(numel(iqr_manual_z_6mps_point3)), ...
    diff_iqr_z_6mps_point3+1.96*std(iqr_manual_z_6mps_point3-iqr_auto_z_6mps_point3)/sqrt(numel(iqr_manual_z_6mps_point3))];

CI_diff_iqr_x_3mps_point4 = [diff_iqr_x_3mps_point4-1.96*std(iqr_manual_x_3mps_point4-iqr_auto_x_3mps_point4)/sqrt(numel(iqr_manual_x_3mps_point4)), ...
    diff_iqr_x_3mps_point4+1.96*std(iqr_manual_x_3mps_point4-iqr_auto_x_3mps_point4)/sqrt(numel(iqr_manual_x_3mps_point4))];
CI_diff_iqr_y_3mps_point4 = [diff_iqr_y_3mps_point4-1.96*std(iqr_manual_y_3mps_point4-iqr_auto_y_3mps_point4)/sqrt(numel(iqr_manual_y_3mps_point4)), ...
    diff_iqr_y_3mps_point4+1.96*std(iqr_manual_y_3mps_point4-iqr_auto_y_3mps_point4)/sqrt(numel(iqr_manual_y_3mps_point4))];
CI_diff_iqr_z_3mps_point4 = [diff_iqr_z_3mps_point4-1.96*std(iqr_manual_z_3mps_point4-iqr_auto_z_3mps_point4)/sqrt(numel(iqr_manual_z_3mps_point4)), ...
    diff_iqr_z_3mps_point4+1.96*std(iqr_manual_z_3mps_point4-iqr_auto_z_3mps_point4)/sqrt(numel(iqr_manual_z_3mps_point4))];

CI_diff_iqr_x_4p5mps_point4 = [diff_iqr_x_4p5mps_point4-1.96*std(iqr_manual_x_4p5mps_point4-iqr_auto_x_4p5mps_point4)/sqrt(numel(iqr_manual_x_4p5mps_point4)), ...
    diff_iqr_x_4p5mps_point4+1.96*std(iqr_manual_x_4p5mps_point4-iqr_auto_x_4p5mps_point4)/sqrt(numel(iqr_manual_x_4p5mps_point4))];
CI_diff_iqr_y_4p5mps_point4 = [diff_iqr_y_4p5mps_point4-1.96*std(iqr_manual_y_4p5mps_point4-iqr_auto_y_4p5mps_point4)/sqrt(numel(iqr_manual_y_4p5mps_point4)), ...
    diff_iqr_y_4p5mps_point4+1.96*std(iqr_manual_y_4p5mps_point4-iqr_auto_y_4p5mps_point4)/sqrt(numel(iqr_manual_y_4p5mps_point4))];
CI_diff_iqr_z_4p5mps_point4 = [diff_iqr_z_4p5mps_point4-1.96*std(iqr_manual_z_4p5mps_point4-iqr_auto_z_4p5mps_point4)/sqrt(numel(iqr_manual_z_4p5mps_point4)), ...
    diff_iqr_z_4p5mps_point4+1.96*std(iqr_manual_z_4p5mps_point4-iqr_auto_z_4p5mps_point4)/sqrt(numel(iqr_manual_z_4p5mps_point4))];

CI_diff_iqr_x_6mps_point4 = [diff_iqr_x_6mps_point4-1.96*std(iqr_manual_x_6mps_point4-iqr_auto_x_6mps_point4)/sqrt(numel(iqr_manual_x_6mps_point4)), ...
    diff_iqr_x_6mps_point4+1.96*std(iqr_manual_x_6mps_point4-iqr_auto_x_6mps_point4)/sqrt(numel(iqr_manual_x_6mps_point4))];
CI_diff_iqr_y_6mps_point4 = [diff_iqr_y_6mps_point4-1.96*std(iqr_manual_y_6mps_point4-iqr_auto_y_6mps_point4)/sqrt(numel(iqr_manual_y_6mps_point4)), ...
    diff_iqr_y_6mps_point4+1.96*std(iqr_manual_y_6mps_point4-iqr_auto_y_6mps_point4)/sqrt(numel(iqr_manual_y_6mps_point4))];
CI_diff_iqr_z_6mps_point4 = [diff_iqr_z_6mps_point4-1.96*std(iqr_manual_z_6mps_point4-iqr_auto_z_6mps_point4)/sqrt(numel(iqr_manual_z_6mps_point4)), ...
    diff_iqr_z_6mps_point4+1.96*std(iqr_manual_z_6mps_point4-iqr_auto_z_6mps_point4)/sqrt(numel(iqr_manual_z_6mps_point4))];

CI_diff_iqr_x_3mps_point5 = [diff_iqr_x_3mps_point5-1.96*std(iqr_manual_x_3mps_point5-iqr_auto_x_3mps_point5)/sqrt(numel(iqr_manual_x_3mps_point5)), ...
    diff_iqr_x_3mps_point5+1.96*std(iqr_manual_x_3mps_point5-iqr_auto_x_3mps_point5)/sqrt(numel(iqr_manual_x_3mps_point5))];
CI_diff_iqr_y_3mps_point5 = [diff_iqr_y_3mps_point5-1.96*std(iqr_manual_y_3mps_point5-iqr_auto_y_3mps_point5)/sqrt(numel(iqr_manual_y_3mps_point5)), ...
    diff_iqr_y_3mps_point5+1.96*std(iqr_manual_y_3mps_point5-iqr_auto_y_3mps_point5)/sqrt(numel(iqr_manual_y_3mps_point5))];
CI_diff_iqr_z_3mps_point5 = [diff_iqr_z_3mps_point5-1.96*std(iqr_manual_z_3mps_point5-iqr_auto_z_3mps_point5)/sqrt(numel(iqr_manual_z_3mps_point5)), ...
    diff_iqr_z_3mps_point5+1.96*std(iqr_manual_z_3mps_point5-iqr_auto_z_3mps_point5)/sqrt(numel(iqr_manual_z_3mps_point5))];

CI_diff_iqr_x_4p5mps_point5 = [diff_iqr_x_4p5mps_point5-1.96*std(iqr_manual_x_4p5mps_point5-iqr_auto_x_4p5mps_point5)/sqrt(numel(iqr_manual_x_4p5mps_point5)), ...
    diff_iqr_x_4p5mps_point5+1.96*std(iqr_manual_x_4p5mps_point5-iqr_auto_x_4p5mps_point5)/sqrt(numel(iqr_manual_x_4p5mps_point5))];
CI_diff_iqr_y_4p5mps_point5 = [diff_iqr_y_4p5mps_point5-1.96*std(iqr_manual_y_4p5mps_point5-iqr_auto_y_4p5mps_point5)/sqrt(numel(iqr_manual_y_4p5mps_point5)), ...
    diff_iqr_y_4p5mps_point5+1.96*std(iqr_manual_y_4p5mps_point5-iqr_auto_y_4p5mps_point5)/sqrt(numel(iqr_manual_y_4p5mps_point5))];
CI_diff_iqr_z_4p5mps_point5 = [diff_iqr_z_4p5mps_point5-1.96*std(iqr_manual_z_4p5mps_point5-iqr_auto_z_4p5mps_point5)/sqrt(numel(iqr_manual_z_4p5mps_point5)), ...
    diff_iqr_z_4p5mps_point5+1.96*std(iqr_manual_z_4p5mps_point5-iqr_auto_z_4p5mps_point5)/sqrt(numel(iqr_manual_z_4p5mps_point5))];

CI_diff_iqr_x_6mps_point5 = [diff_iqr_x_6mps_point5-1.96*std(iqr_manual_x_6mps_point5-iqr_auto_x_6mps_point5)/sqrt(numel(iqr_manual_x_6mps_point5)), ...
    diff_iqr_x_6mps_point5+1.96*std(iqr_manual_x_6mps_point5-iqr_auto_x_6mps_point5)/sqrt(numel(iqr_manual_x_6mps_point5))];
CI_diff_iqr_y_6mps_point5 = [diff_iqr_y_6mps_point5-1.96*std(iqr_manual_y_6mps_point5-iqr_auto_y_6mps_point5)/sqrt(numel(iqr_manual_y_6mps_point5)), ...
    diff_iqr_y_6mps_point5+1.96*std(iqr_manual_y_6mps_point5-iqr_auto_y_6mps_point5)/sqrt(numel(iqr_manual_y_6mps_point5))];
CI_diff_iqr_z_6mps_point5 = [diff_iqr_z_6mps_point5-1.96*std(iqr_manual_z_6mps_point5-iqr_auto_z_6mps_point5)/sqrt(numel(iqr_manual_z_6mps_point5)), ...
    diff_iqr_z_6mps_point5+1.96*std(iqr_manual_z_6mps_point5-iqr_auto_z_6mps_point5)/sqrt(numel(iqr_manual_z_6mps_point5))];

CI_diff_iqr_x_3mps_point6 = [diff_iqr_x_3mps_point6-1.96*std(iqr_manual_x_3mps_point6-iqr_auto_x_3mps_point6)/sqrt(numel(iqr_manual_x_3mps_point6)), ...
    diff_iqr_x_3mps_point6+1.96*std(iqr_manual_x_3mps_point6-iqr_auto_x_3mps_point6)/sqrt(numel(iqr_manual_x_3mps_point6))];
CI_diff_iqr_y_3mps_point6 = [diff_iqr_y_3mps_point6-1.96*std(iqr_manual_y_3mps_point6-iqr_auto_y_3mps_point6)/sqrt(numel(iqr_manual_y_3mps_point6)), ...
    diff_iqr_y_3mps_point6+1.96*std(iqr_manual_y_3mps_point6-iqr_auto_y_3mps_point6)/sqrt(numel(iqr_manual_y_3mps_point6))];
CI_diff_iqr_z_3mps_point6 = [diff_iqr_z_3mps_point6-1.96*std(iqr_manual_z_3mps_point6-iqr_auto_z_3mps_point6)/sqrt(numel(iqr_manual_z_3mps_point6)), ...
    diff_iqr_z_3mps_point6+1.96*std(iqr_manual_z_3mps_point6-iqr_auto_z_3mps_point6)/sqrt(numel(iqr_manual_z_3mps_point6))];

CI_diff_iqr_x_4p5mps_point6 = [diff_iqr_x_4p5mps_point6-1.96*std(iqr_manual_x_4p5mps_point6-iqr_auto_x_4p5mps_point6)/sqrt(numel(iqr_manual_x_4p5mps_point6)), ...
    diff_iqr_x_4p5mps_point6+1.96*std(iqr_manual_x_4p5mps_point6-iqr_auto_x_4p5mps_point6)/sqrt(numel(iqr_manual_x_4p5mps_point6))];
CI_diff_iqr_y_4p5mps_point6 = [diff_iqr_y_4p5mps_point6-1.96*std(iqr_manual_y_4p5mps_point6-iqr_auto_y_4p5mps_point6)/sqrt(numel(iqr_manual_y_4p5mps_point6)), ...
    diff_iqr_y_4p5mps_point6+1.96*std(iqr_manual_y_4p5mps_point6-iqr_auto_y_4p5mps_point6)/sqrt(numel(iqr_manual_y_4p5mps_point6))];
CI_diff_iqr_z_4p5mps_point6 = [diff_iqr_z_4p5mps_point6-1.96*std(iqr_manual_z_4p5mps_point6-iqr_auto_z_4p5mps_point6)/sqrt(numel(iqr_manual_z_4p5mps_point6)), ...
    diff_iqr_z_4p5mps_point6+1.96*std(iqr_manual_z_4p5mps_point6-iqr_auto_z_4p5mps_point6)/sqrt(numel(iqr_manual_z_4p5mps_point6))];

CI_diff_iqr_x_6mps_point6 = [diff_iqr_x_6mps_point6-1.96*std(iqr_manual_x_6mps_point6-iqr_auto_x_6mps_point6)/sqrt(numel(iqr_manual_x_6mps_point6)), ...
    diff_iqr_x_6mps_point6+1.96*std(iqr_manual_x_6mps_point6-iqr_auto_x_6mps_point6)/sqrt(numel(iqr_manual_x_6mps_point6))];
CI_diff_iqr_y_6mps_point6 = [diff_iqr_y_6mps_point6-1.96*std(iqr_manual_y_6mps_point6-iqr_auto_y_6mps_point6)/sqrt(numel(iqr_manual_y_6mps_point6)), ...
    diff_iqr_y_6mps_point6+1.96*std(iqr_manual_y_6mps_point6-iqr_auto_y_6mps_point6)/sqrt(numel(iqr_manual_y_6mps_point6))];
CI_diff_iqr_z_6mps_point6 = [diff_iqr_z_6mps_point6-1.96*std(iqr_manual_z_6mps_point6-iqr_auto_z_6mps_point6)/sqrt(numel(iqr_manual_z_6mps_point6)), ...
    diff_iqr_z_6mps_point6+1.96*std(iqr_manual_z_6mps_point6-iqr_auto_z_6mps_point6)/sqrt(numel(iqr_manual_z_6mps_point6))];

pre_table = [[diff_iqr_x_3mps_point1,diff_iqr_y_3mps_point1,diff_iqr_z_3mps_point1,diff_iqr_x_4p5mps_point1,diff_iqr_y_4p5mps_point1,diff_iqr_z_4p5mps_point1,diff_iqr_x_6mps_point1,diff_iqr_y_6mps_point1,diff_iqr_z_6mps_point1];...
    [diff_iqr_x_3mps_point2,diff_iqr_y_3mps_point2,diff_iqr_z_3mps_point2,diff_iqr_x_4p5mps_point2,diff_iqr_y_4p5mps_point2,diff_iqr_z_4p5mps_point2,diff_iqr_x_6mps_point2,diff_iqr_y_6mps_point2,diff_iqr_z_6mps_point2];...
    [diff_iqr_x_3mps_point3,diff_iqr_y_3mps_point3,diff_iqr_z_3mps_point3,diff_iqr_x_4p5mps_point3,diff_iqr_y_4p5mps_point3,diff_iqr_z_4p5mps_point3,diff_iqr_x_6mps_point3,diff_iqr_y_6mps_point3,diff_iqr_z_6mps_point3];...
    [diff_iqr_x_3mps_point4,diff_iqr_y_3mps_point4,diff_iqr_z_3mps_point4,diff_iqr_x_4p5mps_point4,diff_iqr_y_4p5mps_point4,diff_iqr_z_4p5mps_point4,diff_iqr_x_6mps_point4,diff_iqr_y_6mps_point4,diff_iqr_z_6mps_point4];...
    [diff_iqr_x_3mps_point5,diff_iqr_y_3mps_point5,diff_iqr_z_3mps_point5,diff_iqr_x_4p5mps_point5,diff_iqr_y_4p5mps_point5,diff_iqr_z_4p5mps_point5,diff_iqr_x_6mps_point5,diff_iqr_y_6mps_point5,diff_iqr_z_6mps_point5];...
    [diff_iqr_x_3mps_point6,diff_iqr_y_3mps_point6,diff_iqr_z_3mps_point6,diff_iqr_x_4p5mps_point6,diff_iqr_y_4p5mps_point6,diff_iqr_z_4p5mps_point6,diff_iqr_x_6mps_point6,diff_iqr_y_6mps_point6,diff_iqr_z_6mps_point6]];

pre_table_CI = [[CI_diff_iqr_x_3mps_point1,CI_diff_iqr_y_3mps_point1,CI_diff_iqr_z_3mps_point1,CI_diff_iqr_x_4p5mps_point1,CI_diff_iqr_y_4p5mps_point1,CI_diff_iqr_z_4p5mps_point1,CI_diff_iqr_x_6mps_point1,CI_diff_iqr_y_6mps_point1,CI_diff_iqr_z_6mps_point1];...
    [CI_diff_iqr_x_3mps_point2,CI_diff_iqr_y_3mps_point2,CI_diff_iqr_z_3mps_point2,CI_diff_iqr_x_4p5mps_point2,CI_diff_iqr_y_4p5mps_point2,CI_diff_iqr_z_4p5mps_point2,CI_diff_iqr_x_6mps_point2,CI_diff_iqr_y_6mps_point2,CI_diff_iqr_z_6mps_point2];...
    [CI_diff_iqr_x_3mps_point3,CI_diff_iqr_y_3mps_point3,CI_diff_iqr_z_3mps_point3,CI_diff_iqr_x_4p5mps_point3,CI_diff_iqr_y_4p5mps_point3,CI_diff_iqr_z_4p5mps_point3,CI_diff_iqr_x_6mps_point3,CI_diff_iqr_y_6mps_point3,CI_diff_iqr_z_6mps_point3];...
    [CI_diff_iqr_x_3mps_point4,CI_diff_iqr_y_3mps_point4,CI_diff_iqr_z_3mps_point4,CI_diff_iqr_x_4p5mps_point4,CI_diff_iqr_y_4p5mps_point4,CI_diff_iqr_z_4p5mps_point4,CI_diff_iqr_x_6mps_point4,CI_diff_iqr_y_6mps_point4,CI_diff_iqr_z_6mps_point4];...
    [CI_diff_iqr_x_3mps_point5,CI_diff_iqr_y_3mps_point5,CI_diff_iqr_z_3mps_point5,CI_diff_iqr_x_4p5mps_point5,CI_diff_iqr_y_4p5mps_point5,CI_diff_iqr_z_4p5mps_point5,CI_diff_iqr_x_6mps_point5,CI_diff_iqr_y_6mps_point5,CI_diff_iqr_z_6mps_point5];...
    [CI_diff_iqr_x_3mps_point6,CI_diff_iqr_y_3mps_point6,CI_diff_iqr_z_3mps_point6,CI_diff_iqr_x_4p5mps_point6,CI_diff_iqr_y_4p5mps_point6,CI_diff_iqr_z_4p5mps_point6,CI_diff_iqr_x_6mps_point6,CI_diff_iqr_y_6mps_point6,CI_diff_iqr_z_6mps_point6]];

pre_table_CI = [...
    [strjoin(string(CI_diff_iqr_x_3mps_point1), " to "),strjoin(string(CI_diff_iqr_y_3mps_point1), " to "),strjoin(string(CI_diff_iqr_z_3mps_point1), " to "),strjoin(string(CI_diff_iqr_x_4p5mps_point1), " to "),strjoin(string(CI_diff_iqr_y_4p5mps_point1), " to "),strjoin(string(CI_diff_iqr_z_4p5mps_point1), " to "),strjoin(string(CI_diff_iqr_x_6mps_point1), " to "),strjoin(string(CI_diff_iqr_y_6mps_point1), " to "),strjoin(string(CI_diff_iqr_z_6mps_point1), " to ")];...
    [strjoin(string(CI_diff_iqr_x_3mps_point2), " to "),strjoin(string(CI_diff_iqr_y_3mps_point2), " to "),strjoin(string(CI_diff_iqr_z_3mps_point2), " to "),strjoin(string(CI_diff_iqr_x_4p5mps_point2), " to "),strjoin(string(CI_diff_iqr_y_4p5mps_point2), " to "),strjoin(string(CI_diff_iqr_z_4p5mps_point2), " to "),strjoin(string(CI_diff_iqr_x_6mps_point2), " to "),strjoin(string(CI_diff_iqr_y_6mps_point2), " to "),strjoin(string(CI_diff_iqr_z_6mps_point2), " to ")];...
    [strjoin(string(CI_diff_iqr_x_3mps_point3), " to "),strjoin(string(CI_diff_iqr_y_3mps_point3), " to "),strjoin(string(CI_diff_iqr_z_3mps_point3), " to "),strjoin(string(CI_diff_iqr_x_4p5mps_point3), " to "),strjoin(string(CI_diff_iqr_y_4p5mps_point3), " to "),strjoin(string(CI_diff_iqr_z_4p5mps_point3), " to "),strjoin(string(CI_diff_iqr_x_6mps_point3), " to "),strjoin(string(CI_diff_iqr_y_6mps_point3), " to "),strjoin(string(CI_diff_iqr_z_6mps_point3), " to ")];...
    [strjoin(string(CI_diff_iqr_x_3mps_point4), " to "),strjoin(string(CI_diff_iqr_y_3mps_point4), " to "),strjoin(string(CI_diff_iqr_z_3mps_point4), " to "),strjoin(string(CI_diff_iqr_x_4p5mps_point4), " to "),strjoin(string(CI_diff_iqr_y_4p5mps_point4), " to "),strjoin(string(CI_diff_iqr_z_4p5mps_point4), " to "),strjoin(string(CI_diff_iqr_x_6mps_point4), " to "),strjoin(string(CI_diff_iqr_y_6mps_point4), " to "),strjoin(string(CI_diff_iqr_z_6mps_point4), " to ")];...
    [strjoin(string(CI_diff_iqr_x_3mps_point5), " to "),strjoin(string(CI_diff_iqr_y_3mps_point5), " to "),strjoin(string(CI_diff_iqr_z_3mps_point5), " to "),strjoin(string(CI_diff_iqr_x_4p5mps_point5), " to "),strjoin(string(CI_diff_iqr_y_4p5mps_point5), " to "),strjoin(string(CI_diff_iqr_z_4p5mps_point5), " to "),strjoin(string(CI_diff_iqr_x_6mps_point5), " to "),strjoin(string(CI_diff_iqr_y_6mps_point5), " to "),strjoin(string(CI_diff_iqr_z_6mps_point5), " to ")];...
    [strjoin(string(CI_diff_iqr_x_3mps_point6), " to "),strjoin(string(CI_diff_iqr_y_3mps_point6), " to "),strjoin(string(CI_diff_iqr_z_3mps_point6), " to "),strjoin(string(CI_diff_iqr_x_4p5mps_point6), " to "),strjoin(string(CI_diff_iqr_y_4p5mps_point6), " to "),strjoin(string(CI_diff_iqr_z_4p5mps_point6), " to "),strjoin(string(CI_diff_iqr_x_6mps_point6), " to "),strjoin(string(CI_diff_iqr_y_6mps_point6), " to "),strjoin(string(CI_diff_iqr_z_6mps_point6), " to ")]...
    ];

string([mean(pre_table(:)), mean(pre_table(:))-1.96*std(pre_table(:))/sqrt(numel(pre_table(:))), ...
    mean(pre_table(:))+1.96*std(pre_table(:))/sqrt(numel(pre_table(:)))]);

diff_iqr_table = array2table(pre_table,'RowNames',["1", "2", "3", "4", "5", "6"], 'VariableNames',["3 m/s x","3 m/s y","3 m/s z","4.5 m/s x","4.5 y","4.5 z", "6 m/s x", "6 m/s y", "6 m/s z"])
diff_iqr_table_CI = array2table(pre_table_CI,'RowNames',["1", "2", "3", "4", "5", "6"], 'VariableNames',["3 m/s x","3 m/s y","3 m/s z","4.5 m/s x","4.5 y","4.5 z", "6 m/s x", "6 m/s y", "6 m/s z"])

figure(amplitude_fig);
fontsize(amplitude_fig,8,"points");

amplitude_fig.Units = 'centimeters';
amplitude_fig.OuterPosition = [.1 .1 18.1 18.1];
% exportgraphics(amplitude_fig,'amplitude_fig.eps','Resolution',600)

figure(kinparamfig);
fontsize(kinparamfig,8,"points")

kinparamfig.Units = 'centimeters';
kinparamfig.OuterPosition = [.1 .1 8.1 15.1];
% exportgraphics(kinparamfig,'kinparamfig_fig.eps','Resolution',600)

lme_stroke_amplitude = formatForLME(stroke_amplitude,speeds,treatments,individual,"stroke_amplitude");
lme_wingbeat_frequency = formatForLME(wingbeat_frequency,speeds,treatments,individual,"wingbeat_frequency");
lme_stroke_plane_angle_ds = formatForLME(stroke_plane_angle_ds,speeds,treatments,individual,"stroke_plane_angle_ds");
lme_stroke_plane_angle_us = formatForLME(stroke_plane_angle_us,speeds,treatments,individual,"stroke_plane_angle_us");
lme_angle_of_attack = formatForLME(angle_of_attack,speeds,treatments,individual,"angle_of_attack");
lme_wing_area_midDS = formatForLME(wing_area_midDS,speeds,treatments,individual,"wing_area_midDS");

varNames = ["stroke_amplitude","wingbeat_frequency","stroke_plane_angle_ds","stroke_plane_angle_us","AoA","wing_area_midDS"];

tests = {lme_stroke_amplitude,lme_wingbeat_frequency,lme_stroke_plane_angle_ds,lme_stroke_plane_angle_us,lme_angle_of_attack,lme_wing_area_midDS};

pValues = []; % speed, treatment, speed, treatment, etc

for i = 1:numel(varNames)
    test = tests{i}.anova;
    terms = string(test.Term);
    speed_I = find(matches(terms,"speed"));
    treatment_I = find(matches(terms,"treatment"));
    pValues = [pValues; test.pValue([speed_I treatment_I])];
end

adjusted_p_values = mafdr(pValues, 'BHFDR', true, 'showplot', true);

for i=1:numel(adjusted_p_values)
    if adjusted_p_values(i) < 0.001
        p_valueColumn(i) = "<0.001";
    else
        p_valueColumn(i) = string(num2str(adjusted_p_values(i), "%.3f"));
    end
end

pValues_speed = p_valueColumn(1:2:end);
pValues_treatment = p_valueColumn(2:2:end);

kinparamtable = table(pValues_speed',pValues_treatment', 'RowNames',["stroke amplitude", "wingbeat frequency", "stroke plane angle (ds)", "stroke plane angle (us)", "AoA", "wing area (mid ds)"], 'variableNames', ["p-value speed", "p-value treatment"]); % adjusted
%writetable(kinparamtable, "pValuesKinParams20230616.csv")

parameters = {stroke_amplitude,wingbeat_frequency,stroke_plane_angle_ds,stroke_plane_angle_us,angle_of_attack,wing_area_midDS};
preTable=string(nan(6,2));

units = ["mm", string(char(176)),string(char(176)),string(char(176)),string(char(176)),strcat("cm",string(char(178)))];

for i=1:6
    parameter = parameters{i};
    parameter_auto = cell2mat({parameter{1:3,1}});
    parameter_manual = cell2mat({parameter{1:3,2}});
    parameter_average_auto =  mean(parameter_auto);
    parameter_average_manual =  mean(parameter_manual);
    parameter_95CI_auto = parameter_average_auto+[-1,1]*(1.96*std(parameter_auto)/sqrt(length(parameter_auto)));
    parameter_95CI_manual = parameter_average_manual+[-1,1]*(1.96*std(parameter_manual)/sqrt(length(parameter_manual)));
    parameter_95CI_auto = sort(parameter_95CI_auto);
    parameter_95CI_manual = sort(parameter_95CI_manual);
    preTable(i,1) = strcat(num2str(parameter_average_auto, "%.2f"), " (", [num2str(parameter_95CI_auto(1), "%.2f"), ' to ', num2str(parameter_95CI_auto(2), "%.2f")], ")");
    preTable(i,2) = strcat(num2str(parameter_average_manual, "%.2f"), " (", [num2str(parameter_95CI_manual(1), "%.2f"), ' to ', num2str(parameter_95CI_manual(2), "%.2f")], ")");
end

% array2table(preTable,'VariableNames', ["Auto", "Manual"], "RowNames", ...
%     ["stroke amplitude (mm)", ...
%     "wingbeat frequency (Hz)", ...
%     strcat("stroke plane angle, ds (", units(3), ")"), ...
%     strcat("stroke plane angle, us (", units(4), ")"), ...
%     strcat("angle of attack (", units(5), ")"), ...
%     strcat("wing area mid ds (", units(6), ")")])


kinParamTable=array2table([preTable,pValues_speed',pValues_treatment'],'VariableNames', ["Auto", "Manual","p-value speed", "p-value treatment"], "RowNames", ...
    ["stroke amplitude (mm)", ...
    "wingbeat frequency (Hz)", ...
    strcat("stroke plane angle, ds (", units(3), ")"), ...
    strcat("stroke plane angle, us (", units(4), ")"), ...
    strcat("angle of attack (", units(5), ")"), ...
    strcat("wing area mid ds (", units(6), ")")])

estimate = [];
tstat = [];
df = [];
pvalue = [];
for i = 1:numel(tests)
    test = tests{i};
    estimate = [estimate; table2array(dataset2table(test.Coefficients(:,2)))];
    tstat = [tstat; table2array(dataset2table(test.Coefficients(:,4)))];
    df = [df; table2array(dataset2table(test.Coefficients(:,5)))];
    pvalue = [pvalue; table2array(dataset2table(test.Coefficients(:,6)))];
end

fullstat_table = table([estimate tstat df pvalue])
%writetable(fulltable, 'fullkinparstat.xls')