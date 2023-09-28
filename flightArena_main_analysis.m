%%
% This script analyzes the DLC output from our trained network. It
% calculates the human-human variation from paired human digitization 
% and then compares that to the human-auto variation.
% This is the code used to produce figures 6-7 and table 2.

% TODO: add check for dependencies here like in the windtunnel script
% TODO: improve commenting for better readability.


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

% filter parameters
speedLimitsFile = "speed_075.mat"
load(speedLimitsFile);
goodSpeeds = maxSpeed*3;
maxSpeeds = maxSpeed*9;
smoothspan = 4;
goodDLTdvResidual = 4;
maxDLTdvResidual = 16;
goodDLClikelihood = 0.55;
minDLClikelihood_2d = .475;
minDLClikelihood_3d = .525;
frameRate = 800;
crp=2;

easyWandDataFile = fullfile(pwd, "flightArenaTrials\calibration09_1maxRes_easyWandData.mat");
load(easyWandDataFile);

cam1TformsFile = fullfile(pwd, "flightArenaTrials/calibration09_1maxRes_cam1Tforms.mat");

dlc_network_folder = fullfile(pwd, "flightArenaTrials/LRSWAP_20230407_shuffle4_snapshot22")';

dlc_network_folder_noswap = fullfile(pwd, "flightArenaTrials/NOLRSWAP_20230516_shuffle3_snapshot22")';


dlc_predicted_bouts = dir(dlc_network_folder);
dlc_predicted_bouts = fullfile(string(dlc_predicted_bouts(1).folder),string({dlc_predicted_bouts.name}))';
dlc_predicted_bouts = dlc_predicted_bouts(3:end);

dlc_predicted_bouts_noswap = dir(dlc_network_folder_noswap);
dlc_predicted_bouts_noswap = fullfile(string(dlc_predicted_bouts_noswap(1).folder),string({dlc_predicted_bouts_noswap.name}))';
dlc_predicted_bouts_noswap = dlc_predicted_bouts_noswap(3:end);

translations_base_folder = fullfile(pwd, "flightArenaTrials/crop_files/");
translation_files_folders = dir(translations_base_folder);
translation_files_folders = fullfile(translations_base_folder, string({translation_files_folders.name}))';
translation_files_folders = translation_files_folders(3:end);


validations1_base_folder = fullfile(pwd,"flightArenaTrials/validation_files_1/");
validations2_base_folder = fullfile(pwd,"flightArenaTrials/validation_files_2/");

validation_files1 = dir(validations1_base_folder);
validation_files1 = fullfile(validations1_base_folder, string({validation_files1.name}))';
validation_files1 = validation_files1(3:end);

validation_files2 = dir(validations2_base_folder);
validation_files2 = fullfile(validations2_base_folder, string({validation_files2.name}))';
validation_files2 = validation_files2(3:end);

point_distance_humans = [];
point_distance_swapfix = [];
point_distance_noswapfix = [];
point_distance_noswapaug = [];

humanToUse =  [1 2 1 2 2 1];

for trial_number = 1:(numel(dlc_predicted_bouts))
    current_bout = dlc_predicted_bouts(trial_number);
    current_bout_noswap = dlc_predicted_bouts_noswap(trial_number);

    h1_dlt = validation_files1(trial_number);
    h2_dlt = validation_files2(trial_number);

    [~, digitized_frames1] = getXYfromDLTdv(h1_dlt);
    [~, digitized_frames2] = getXYfromDLTdv(h2_dlt);

    % the frames that are digitized by humans
    digitized_frames1(mod((digitized_frames1-25),50)~=0,:,:,:) = [];
    digitized_frames2(mod((digitized_frames2-25),50)~=0,:,:,:) = [];

    % get the DLC output CSVs (one per camera)
    CSVs = dir(current_bout);
    CSVs = string({CSVs.name})';
    DLCoutputCSVs(1) = CSVs(contains(CSVs, "Cam1") & contains(CSVs, "csv"));
    DLCoutputCSVs(2) = CSVs(contains(CSVs, "Cam2") & contains(CSVs, "csv"));
    DLCoutputCSVs(3) = CSVs(contains(CSVs, "Cam3") & contains(CSVs, "csv"));
    DLCoutputCSVs = fullfile(current_bout, DLCoutputCSVs);
    DLCoutputCSVs = DLCoutputCSVs';

    % get the DLC output CSVs noswap (one per camera)
    CSVs = dir(current_bout_noswap);
    CSVs = string({CSVs.name})';
    DLCoutputCSVs_noswap(1) = CSVs(contains(CSVs, "Cam1") & contains(CSVs, "csv"));
    DLCoutputCSVs_noswap(2) = CSVs(contains(CSVs, "Cam2") & contains(CSVs, "csv"));
    DLCoutputCSVs_noswap(3) = CSVs(contains(CSVs, "Cam3") & contains(CSVs, "csv"));
    DLCoutputCSVs_noswap = fullfile(current_bout_noswap, DLCoutputCSVs_noswap);
    DLCoutputCSVs_noswap = DLCoutputCSVs_noswap';

    % get the thrutracker translation files (one per camera)
    mats = dir(translation_files_folders(trial_number));
    mats = string({mats.name});
    ThruTrackerTranslationFiles(1) = mats(contains(mats, "Cam1"));
    ThruTrackerTranslationFiles(2) = mats(contains(mats, "Cam2"));
    ThruTrackerTranslationFiles(3) = mats(contains(mats, "Cam3"));
    ThruTrackerTranslationFiles = fullfile(translation_files_folders(trial_number), ThruTrackerTranslationFiles);
    ThruTrackerTranslationFiles = ThruTrackerTranslationFiles';

    % get the XYZ from human digitizers
    XYZ_h1 = getXYZfromDLTdv(h1_dlt, easyWandDataFile, cam1TformsFile, crp*1, ThruTrackerTranslationFiles);
    XYZ_h2 = getXYZfromDLTdv(h2_dlt, easyWandDataFile, cam1TformsFile, crp*1, ThruTrackerTranslationFiles);

    % get the XYZ from dlc with swap fix
    [XYZ_DLC_swapfix, DLTdvResidual_dlc_swapfix, meanLikelihoodDLC_swapfix] = ...
        getXYZ_swapFix(easyWandDataFile, [cam1TformsFile, cam1TformsFile, cam1TformsFile], DLCoutputCSVs, ThruTrackerTranslationFiles, crp, minDLClikelihood_2d);
    
    XYZ_DLC_swapfix(XYZ_DLC_swapfix(:)==0) = nan;
    DLTdvResidual_dlc_swapfix(DLTdvResidual_dlc_swapfix==0) = inf;

    [XYZ_DLC_noswapaug, DLTdvResidual_dlc_noswapaug, meanLikelihoodDLC_noswapaug] = ...
        getXYZ_noSwapFixNoFewerCamTest(easyWandDataFile, [cam1TformsFile, cam1TformsFile, cam1TformsFile], DLCoutputCSVs_noswap, ThruTrackerTranslationFiles, crp, minDLClikelihood_2d);
    XYZ_DLC_noswapaug(XYZ_DLC_noswapaug(:)==0) = nan;
    DLTdvResidual_dlc_noswapaug(DLTdvResidual_dlc_noswapaug==0) = inf;

    [XYZ_DLC_noswapfix, DLTdvResidual_dlc_noswapfix, meanLikelihoodDLC_noswapfix] = ...
        getXYZ_noSwapFixNoFewerCamTest(easyWandDataFile, [cam1TformsFile, cam1TformsFile, cam1TformsFile], DLCoutputCSVs, ThruTrackerTranslationFiles, crp, minDLClikelihood_2d);

    XYZ_DLC_noswapfix(XYZ_DLC_noswapfix(:)==0) = nan;
    DLTdvResidual_dlc_noswapfix(DLTdvResidual_dlc_noswapfix==0) = inf;

    % filter out bad predictions
    XYZ_DLC_filtered_CCS_swapfix = ...
        filter3Dpredictions(...
        XYZ_DLC_swapfix, ...
        DLTdvResidual_dlc_swapfix, ...
        meanLikelihoodDLC_swapfix, ...
        frameRate, ...
        goodSpeeds, ...
        maxSpeeds, ...
        goodDLTdvResidual, ...
        maxDLTdvResidual, ...
        goodDLClikelihood, ...
        minDLClikelihood_3d);

    XYZ_DLC_smoothed_noswapaug = nan(size(XYZ_DLC_swapfix));
    XYZ_DLC_smoothed_swapfix = nan(size(XYZ_DLC_swapfix));
    XYZ_DLC_smoothed_noswapfix = nan(size(XYZ_DLC_swapfix));

    minLastIndex_swapfix = inf;

    for k=1:size(XYZ_DLC_filtered_CCS_swapfix,2)
        for j=1:3
            XYZ_DLC_filtered_CCS_swapfix(find(isoutlier(XYZ_DLC_filtered_CCS_swapfix(:,k,j),'movmedian',smoothspan*2)),k,j) = nan;
            firstIndex = find(~isnan(XYZ_DLC_filtered_CCS_swapfix(:,k,j)),true,'first');
            lastIndex = find(~isnan(XYZ_DLC_filtered_CCS_swapfix(:,k,j)),true,'last');
            minLastIndex_swapfix = min(minLastIndex_swapfix, lastIndex);
            XYZ_DLC_smoothed_swapfix(firstIndex:lastIndex,k,j) = tybutterNaN(smooth(XYZ_DLC_filtered_CCS_swapfix(firstIndex:lastIndex,k,j),smoothspan,'rlowess'), 45, frameRate, 'low');

            XYZ_DLC_noswapfix(find(isoutlier(XYZ_DLC_noswapfix(:,k,j),'movmedian',smoothspan*2)),k,j) = nan; %maybe remove...
            firstIndex = find(~isnan(XYZ_DLC_noswapfix(:,k,j)),true,'first');
            lastIndex = find(~isnan(XYZ_DLC_noswapfix(:,k,j)),true,'last');
            XYZ_DLC_smoothed_noswapfix(firstIndex:lastIndex,k,j) = tybutterNaN(smooth(XYZ_DLC_noswapfix(firstIndex:lastIndex,k,j),smoothspan,'rlowess'), 45, frameRate, 'low');

            XYZ_DLC_noswapaug(find(isoutlier(XYZ_DLC_noswapaug(:,k,j),'movmedian',smoothspan*2)),k,j) = nan; %maybe remove...
            firstIndex = find(~isnan(XYZ_DLC_noswapaug(:,k,j)),true,'first');
            lastIndex = find(~isnan(XYZ_DLC_noswapaug(:,k,j)),true,'last');
            XYZ_DLC_smoothed_noswapaug(firstIndex:lastIndex,k,j) = tybutterNaN(smooth(XYZ_DLC_noswapaug(firstIndex:lastIndex,k,j),smoothspan,'rlowess'), 45, frameRate, 'low');
        end
    end
    XYZ_DLC_smoothed_swapfix(XYZ_DLC_smoothed_swapfix==0)=nan;
    XYZ_DLC_noswapfix(XYZ_DLC_noswapfix==0)=nan;
    XYZ_DLC_noswapaug(XYZ_DLC_noswapaug==0)=nan;

colors = [107,11,225;... fixed so that tail is not red
    87,61,233; ...
    76,94,226; ...
    32,146,216; ...
    15,190,215; ...
    57,209,173; ...
    94,232,185; ...
    124,242,200; ...
    163,238,197; ...
    196,226,174; ...
    222,193,130; ...
    225,158,105; ...
    219,109,82; ...
    232,53,55; ...
    222,12,31; ...
    205,218,158; ...
    ]/255;

    figure(trial_number); hold on; axis equal    
    frameEnd = minLastIndex_swapfix;%size(XYZ_DLC_smoothed,1)-0;
    plot_3D_bat(XYZ_DLC_smoothed_swapfix,frameEnd)
    for k=1:size(XYZ_DLC_filtered_CCS_swapfix,2) 
        plot3(XYZ_DLC_smoothed_swapfix(1:frameEnd,k,1),XYZ_DLC_smoothed_swapfix(1:frameEnd,k,2),XYZ_DLC_smoothed_swapfix(1:frameEnd,k,3),'LineWidth',1.5,'Color', colors(k,:));
    end
    grid on
    if(trial_number==6)
        xlim([-.6 .8]); xticks(-.6:.2:.8); xticklabels({})
        ylim([-.7 .5]); yticks(-.7:.2:.5); yticklabels({})
        zlim([-.3 .1]); zticks(-.3:.2:.1); zticklabels({})
        %plot3([-0.2, 0],[.1 .1],[-.3, -.3],'k', 'LineWidth',3)
        plot3([0.8, 0.6],[-.7 -.7],[.1, .1],'k', 'LineWidth',3)
        plot3([0.8, 0.8],[-.7 -.5],[.1, .1],'k', 'LineWidth',3)
        plot3([0.8, 0.8],[-.7 -.7],[.1, -.1],'k', 'LineWidth',3)
        view(-140,20)
        currfig = gcf;

    currfig.Units = 'centimeters';
    currfig.OuterPosition = [.1 .1 18.1 12.1];
    end

     I = humanToUse(trial_number);
%     if I == 1
%         plot3(XYZ_h1(digitized_frames1,:,1),XYZ_h1(digitized_frames1,:,2),XYZ_h1(digitized_frames1,:,3),'kx');
%     else
%         plot3(XYZ_h2(digitized_frames1,:,1),XYZ_h2(digitized_frames1,:,2),XYZ_h2(digitized_frames1,:,3), 'kx');
%     end

    current_points = [1:16];
    XYZ_h1_test = squeeze(XYZ_h1(digitized_frames1,current_points, :));
    XYZ_h2_test = squeeze(XYZ_h2(digitized_frames1,current_points, :));

    XYZ_DLC_swapfix_reduced = squeeze(XYZ_DLC_smoothed_swapfix(digitized_frames1,current_points, :));
    XYZ_DLC_noswapfix_reduced = squeeze(XYZ_DLC_smoothed_noswapfix(digitized_frames1,current_points, :));
    XYZ_DLC_noswapaug_reduced = squeeze(XYZ_DLC_smoothed_noswapaug(digitized_frames1,current_points, :));

    diff3D_h = XYZ_h1_test - XYZ_h2_test;
    
    if I == 1
        XYZ_human = XYZ_h1_test;
    else
        XYZ_human = XYZ_h2_test;
    end

    diff3D_swapfix = XYZ_human - XYZ_DLC_swapfix_reduced;
    diff3D_noswapfix = XYZ_human - XYZ_DLC_noswapfix_reduced;
    diff3D_noswapaug = XYZ_human - XYZ_DLC_noswapaug_reduced;

    clear   point_distance_h_TEMP ...
            point_distance_swapfix_TEMP ...
            point_distance_noswapfix_TEMP ...
            point_distance_noswapaug_TEMP 

    for k = 1:size(diff3D_h,1)
        for j = 1:size(diff3D_h,2)
            point_distance_h_TEMP(k,j) = norm(squeeze(diff3D_h(k,j,:)));
            point_distance_swapfix_TEMP(k,j) = norm(squeeze(diff3D_swapfix(k,j,:)));
            point_distance_noswapfix_TEMP(k,j) = norm(squeeze(diff3D_noswapfix(k,j,:)));
            point_distance_noswapaug_TEMP(k,j) = norm(squeeze(diff3D_noswapaug(k,j,:)));
        end
    end

    point_distance_humans = [point_distance_humans; point_distance_h_TEMP];
    point_distance_swapfix = [point_distance_swapfix; point_distance_swapfix_TEMP];
    point_distance_noswapfix = [point_distance_noswapfix; point_distance_noswapfix_TEMP];
    point_distance_noswapaug = [point_distance_noswapaug; point_distance_noswapaug_TEMP];
end

if trial_number ~=6
    return
end

bat_image = fullfile(pwd, "flightArenaTrials/n_humeralis morphology.png");

wingspan_in_m = norm(squeeze(XYZ_h1_test(6,1,:)-XYZ_h1_test(6,15,:)))/2;
resultsCircleFig = figure; imshow(bat_image); hold on;

wingspan_bar_position = [195 344; 188 188];

wingspan_in_iu = range(wingspan_bar_position(1,:)); % image units

m_to_iu = wingspan_in_iu/wingspan_in_m;

line_color = 'w'; %standard MATLAB red
DLC_color = [141 47 162]./256;
human_digitizer_color = [0, 0.4470, 0.7410];

line(wingspan_bar_position(1,:), wingspan_bar_position(2,:), 'Color', line_color, 'LineWidth', 2);
plot(wingspan_bar_position(1,1), wingspan_bar_position(2,1), 'MarkerEdgeColor', line_color, 'MarkerFaceColor', line_color, 'Marker', '|', 'LineWidth',2);
plot(wingspan_bar_position(1,2), wingspan_bar_position(2,2), 'MarkerEdgeColor', line_color, 'MarkerFaceColor', line_color, 'Marker', '|', 'LineWidth',2);

text_position = [wingspan_bar_position(1,1)*1.3, wingspan_bar_position(2,1)*0.94];

wingspan_string = strcat(num2str(wingspan_in_m, 2), " m");

point_x_positions = [...
    46, ...     % t3L
    121, ...    % wstL
    108, ...    % t5L
    159, ...    % elbL
    179, ...    % shdL
    172, ...    % ankL
    202, ...    % nl
    198, ...    % str
    193, ...    % lmb
    215, ...    % shdR
    215, ...    % ankR
    230, ...    % elbR
    273, ...    % wstR
    279, ...    % t5R
    343, ...    % t3R
    194, ...    % tail
    ];

point_y_positions = [...
    98, ...     % t3L
    100, ...    % wstL
    145, ...    % t5L
    112, ...    % elbL
    107, ...    % shdL
    159, ...    % ankL
    93, ...     % nl
    116, ...    % str
    148, ...    % lmb
    106, ...    % shdR
    159, ...    % ankR
    111, ...    % elbR
    101, ...    % wstR
    145, ...    % t5R
    101, ...    % t3R
    180, ...    % tail
    ];

ref_x_pos = [311.5 317 324 333 343];
ref_y_pos = [165 165 165 165 165];

circle_alpha = 0.75;

meansMachine = [];
meansHuman = [];
seMachine = [];
seHuman = [];
nMachine = [];
nHuman = [];
p_values = [];
CIs_diff = nan(11,2);
CIs_meanMachine = nan(11,2);
CIs_meanHuman = nan(11,2);
for i = [7 8 9 16]
    r = mean(point_distance_swapfix(:,i), 'omitnan')*m_to_iu/2;
    if(i==16)
        nMachine(10) = numel(point_distance_swapfix(:,10))-sum(isnan(point_distance_swapfix(:,16)));
        nHuman(10) = numel(point_distance_humans(:,10))-sum(isnan(point_distance_humans(:,16)));
        meansMachine(10) = mean(point_distance_swapfix(:,16),'omitnan');
        seMachine(10) = std(point_distance_swapfix(:,16),'omitnan')/sqrt(nMachine(10));
        meansHuman(10) = mean(point_distance_humans(:,16),'omitnan');
        seHuman(10) = std(point_distance_humans(:,16),'omitnan')/sqrt(nHuman(10));

        CIs_meanHuman(10,1) = meansHuman(10) - 1.96*seHuman(10);
        CIs_meanHuman(10,2) = meansHuman(10) + 1.96*seHuman(10);
        CIs_meanMachine(10,1) = meansMachine(10) - 1.96*seMachine(10);
        CIs_meanMachine(10,2) = meansMachine(10) + 1.96*seMachine(10);

        [h,p_values(10),CIs_diff(10,:)] = ttest2(point_distance_humans(:,16)*1000, point_distance_swapfix(:,16)*1000,'Vartype','unequal');
    else
        nMachine(i) = numel(point_distance_swapfix(:,i))-sum(isnan(point_distance_swapfix(:,i)));
        nHuman(i) = numel(point_distance_humans(:,i))-sum(isnan(point_distance_humans(:,i)));

        meansMachine(i) = mean(point_distance_swapfix(:,i), 'omitnan');
        seMachine(i) = std(point_distance_swapfix(:,i),'omitnan')/sqrt(nMachine(i));
        meansHuman(i) = mean(point_distance_humans(:,i),'omitnan');
        seHuman(i) = std(point_distance_humans(:,i),'omitnan')/sqrt(nHuman(i));

        CIs_meanHuman(i,1) = meansHuman(i) - 1.96*seHuman(i);
        CIs_meanHuman(i,2) = meansHuman(i) + 1.96*seHuman(i);
        CIs_meanMachine(i,1) = meansMachine(i) - 1.96*seMachine(i);
        CIs_meanMachine(i,2) = meansMachine(i) + 1.96*seMachine(i);

        [h,p_values(i),CIs_diff(i,:)] = ttest2(point_distance_humans(:,i)*1000, point_distance_swapfix(:,i)*1000,'Vartype','unequal');
        
    end
    th = linspace( pi/2, -pi/2, 100);
    semicrc = r.*[cos(th); sin(th)];
    x = semicrc(1,:) + point_x_positions(i);
    y = semicrc(2,:) + point_y_positions(i);
    fill(x, y, DLC_color, 'EdgeAlpha', 0, 'FaceAlpha', circle_alpha)

    r = nanmean(point_distance_humans(:,i))*m_to_iu/2;
    th = linspace( -pi/2, -pi*1.5, 100);
    semicrc = r.*[cos(th); sin(th)];
    x = semicrc(1,:) + point_x_positions(i);
    y = semicrc(2,:) + point_y_positions(i);
    fill(x, y, human_digitizer_color, 'EdgeAlpha', 0, 'FaceAlpha', circle_alpha)
end

pair_points = [15 13 14 12 10 11];
for i = [1 2 3 4 5 6]
    hdiff = point_distance_humans(:,[i pair_points(i)]);
    hdiff = hdiff(:);

    mdiff = point_distance_swapfix(:,[i pair_points(i)]);
    mdiff = mdiff(:);

    nMachine(i) = numel(point_distance_swapfix(:,[i pair_points(i)]))-sum(isnan(point_distance_swapfix(:,[i pair_points(i)])),'all');
    nHuman(i) = numel(point_distance_humans(:,[i pair_points(i)]))-sum(isnan(point_distance_humans(:,[i pair_points(i)])),'all');

    meansMachine(i) = mean(point_distance_swapfix(:,[i pair_points(i)]), 'all', 'omitnan');
    seMachine(i) = std(point_distance_swapfix(:,[i pair_points(i)]), [],'all', 'omitnan')/sqrt(nMachine(i));    
    meansHuman(i) = mean(point_distance_humans(:,[i pair_points(i)]), 'all', 'omitnan');
    seHuman(i) = std(point_distance_humans(:,[i pair_points(i)]), [],'all', 'omitnan')/sqrt(nHuman(i));

    CIs_meanHuman(i,1) = meansHuman(i) - 1.96*seHuman(i);
    CIs_meanHuman(i,2) = meansHuman(i) + 1.96*seHuman(i);
    CIs_meanMachine(i,1) = meansMachine(i) - 1.96*seMachine(i);
    CIs_meanMachine(i,2) = meansMachine(i) + 1.96*seMachine(i);

    [h,p_values(i),CIs_diff(i,:)] = ttest2(hdiff*1000, mdiff*1000,'Vartype','unequal');

    r = nanmean(point_distance_swapfix(:,[i pair_points(i)]), 'all')*m_to_iu/2;
    th = linspace( pi/2, -pi/2, 100);
    semicrc = r.*[cos(th); sin(th)];
    x = semicrc(1,:) + point_x_positions(pair_points(i));
    y = semicrc(2,:) + point_y_positions(pair_points(i));
    fill(x, y, DLC_color, 'EdgeAlpha', 0, 'FaceAlpha', circle_alpha)

    r = nanmean(point_distance_humans(:,[i pair_points(i)]), 'all')*m_to_iu/2;
    th = linspace( -pi/2, -pi*1.5, 100);
    semicrc = r.*[cos(th); sin(th)];
    x = semicrc(1,:) + point_x_positions(pair_points(i));
    y = semicrc(2,:) + point_y_positions(pair_points(i));
    fill(x, y, human_digitizer_color, 'EdgeAlpha', 0, 'FaceAlpha', circle_alpha)
end

for i = [1 2 3 4 5]
    r = i/1000*m_to_iu/2;
    th = linspace(pi/2, -pi*1.5, 100);
    semicrc = r.*[cos(th); sin(th)];
    x = semicrc(1,:) + ref_x_pos(i);
    y = semicrc(2,:) + ref_y_pos(i);
    fill(x, y, "w", 'EdgeAlpha', 0, 'FaceAlpha', circle_alpha)

end

nMachine(11) = numel(point_distance_swapfix(:))-sum(isnan(point_distance_swapfix(:)),'all');
nHuman(11) = numel(point_distance_humans(:))-sum(isnan(point_distance_humans(:)),'all');

meansMachine(end+1) = mean(point_distance_swapfix(:), 'all', 'omitnan');
seMachine(end+1) = std(point_distance_swapfix(:), [],'all', 'omitnan')/sqrt(nMachine(11));
meansHuman(end+1) = mean(point_distance_humans(:), 'all', 'omitnan');
seHuman(end+1) = std(point_distance_humans(:), [],'all', 'omitnan')/sqrt(nHuman(11));

CIs_meanHuman(11,1) = meansHuman(11) - 1.96*seHuman(11);
CIs_meanHuman(11,2) = meansHuman(11) + 1.96*seHuman(11);
CIs_meanMachine(11,1) = meansMachine(11) - 1.96*seMachine(11);
CIs_meanMachine(11,2) = meansMachine(11) + 1.96*seMachine(11);

[h,p_values(11),CIs_diff(11,:)] = ttest2(point_distance_humans(:)*1000, point_distance_swapfix(:)*1000,'Vartype','unequal');

%% make stat table
adjusted_p_values = mafdr(p_values, 'BHFDR', true, 'showplot', false);
clear differenceInMeanColumn
for i=1:numel(meansMachine)
    humanColumn(i) = strcat(num2str(meansHuman(i)*1000,"%.2f")," (", num2str(CIs_meanHuman(i,1)*1000,"%.2f"), " to ", num2str(CIs_meanHuman(i,2)*1000,"%.2f"), ")");
    
    machineColumn(i) = strcat(num2str(meansMachine(i)*1000,"%.2f")," (", num2str(CIs_meanMachine(i,1)*1000,"%.2f"), " to ", num2str(CIs_meanMachine(i,2)*1000,"%.2f"), ")");

    differenceInMeanColumn(i) = strcat(num2str((meansHuman(i)-meansMachine(i))*1000,"%.2f")," (", num2str(CIs_diff(i,1),"%.2f"), " to ", num2str(CIs_diff(i,2),"%.2f"), ")");

    if adjusted_p_values(i) < 0.001
        p_valueColumn(i) = "<0.001";
    else
        p_valueColumn(i) = string(num2str(adjusted_p_values(i), "%.3f"));
    end

    if adjusted_p_values(i) < 0.001
        p_valueColumn(i) = strcat(p_valueColumn(i), "***");
    elseif adjusted_p_values(i) < 0.01
        p_valueColumn(i) = strcat(p_valueColumn(i), "**");
    elseif adjusted_p_values(i) < 0.05
        p_valueColumn(i) = strcat(p_valueColumn(i), "*");
    end
end

output_table = table(humanColumn', machineColumn', differenceInMeanColumn', p_valueColumn','variableNames', ["Human-human difference (mm)","Human-Auto difference (mm)", "Difference in means (mm)", "p-value"], 'rowNames', ["Wingtip","Wrist","5th digit","Elbow","Shoulder","Ankle","Nose", "Sternum", "Base of tail", "Tip of tail", "All"])

hdiff = point_distance_humans(:,[1:6,10:15]);
hdiff = hdiff(:);

mdiff = point_distance_swapfix(:,[1:6,10:15]);
mdiff = mdiff(:);

disp("Humans")
mean(point_distance_humans(:), 'omitnan')
disp("swapfix")
mean(point_distance_swapfix(:), 'omitnan')
disp("noswapfix")
mean(point_distance_noswapfix(:), 'omitnan')
disp("noswapaug")
mean(point_distance_noswapaug(:), 'omitnan')

meanHumanAll = mean(point_distance_humans(:), 'omitnan')*1000;
nHumanAll = numel(point_distance_humans(:))-sum(isnan(point_distance_humans(:)));
CIRangeHumanAll = 1.96*std(point_distance_humans(:), [],'all', 'omitnan')/sqrt(nHumanAll)*1000;

meanHumanWingtips = mean(point_distance_humans(:,[1,15]), 'all', 'omitnan')*1000;
nHumanWingtips = numel(point_distance_humans(:,[1,15]))-sum(isnan(point_distance_humans(:,[1,15])),'all');
CIRangeHumanWingtips = 1.96*std(point_distance_humans(:,[1,15]), [],'all', 'omitnan')/sqrt(nHumanWingtips)*1000;

meanSwapFixAll = mean(point_distance_swapfix(:), 'omitnan')*1000;
nSwapFixAll = numel(point_distance_swapfix(:))-sum(isnan(point_distance_swapfix(:)));
CIRangeSwapFixAll = 1.96*std(point_distance_swapfix(:), [],'all', 'omitnan')/sqrt(nSwapFixAll)*1000;

meanSwapFixWingtips = mean(point_distance_swapfix(:,[1,15]), 'all', 'omitnan')*1000;
nSwapFixWingtips = numel(point_distance_swapfix(:,[1,15]))-sum(isnan(point_distance_swapfix(:,[1,15])),'all');
CIRangeSwapFixWingtips = 1.96*std(point_distance_swapfix(:,[1,15]), [],'all', 'omitnan')/sqrt(nSwapFixWingtips)*1000;

meanNoSwapFixAll = mean(point_distance_noswapfix(:), 'omitnan')*1000;
nNoSwapFixAll = numel(point_distance_noswapfix(:))-sum(isnan(point_distance_noswapfix(:)));
CIRangeNoSwapFixAll = 1.96*std(point_distance_noswapfix(:), [],'all', 'omitnan')/sqrt(nNoSwapFixAll)*1000;

meanNoSwapFixWingtips = mean(point_distance_noswapfix(:,[1,15]), 'all', 'omitnan')*1000;
nNoSwapFixWingtips = numel(point_distance_noswapfix(:,[1,15]))-sum(isnan(point_distance_noswapfix(:,[1,15])),'all');
CIRangeNoSwapFixWingtips = 1.96*std(point_distance_noswapfix(:,[1,15]), [],'all', 'omitnan')/sqrt(nNoSwapFixWingtips)*1000;

meanNoSwapAugAll = mean(point_distance_noswapaug(:), 'omitnan')*1000;
nNoSwapAugAll = numel(point_distance_noswapaug(:))-sum(isnan(point_distance_noswapaug(:)));
CIRangeNoSwapAugAll = 1.96*std(point_distance_noswapaug(:), [],'all', 'omitnan')/sqrt(nNoSwapAugAll)*1000;

meanNoSwapAugWingtips = mean(point_distance_noswapaug(:,[1,15]), 'all', 'omitnan')*1000;
nNoSwapAugWingtips = numel(point_distance_noswapaug(:,[1,15]))-sum(isnan(point_distance_noswapaug(:,[1,15])),'all');
CIRangeNoSwapAugWingtips = 1.96*std(point_distance_noswapaug(:,[1,15]), [],'all', 'omitnan')/sqrt(nNoSwapAugWingtips)*1000;

improvement_means_all = [meanNoSwapAugAll, meanNoSwapFixAll, meanSwapFixAll, meanHumanAll];
improvement_means_wts = [meanNoSwapAugWingtips, meanNoSwapFixWingtips, meanSwapFixWingtips, meanHumanWingtips];

improvement_CIs_all = [CIRangeNoSwapAugAll, CIRangeNoSwapFixAll, CIRangeSwapFixAll, CIRangeHumanAll];
improvement_CIs_wts = [CIRangeNoSwapAugWingtips, CIRangeNoSwapFixWingtips, CIRangeSwapFixWingtips, CIRangeHumanWingtips];

improvements_x = 1:4;
improvementsfigure = figure; hold on
plot(improvements_x-0.05, improvement_means_all, 'Color',[mean(colors)*0.6, 1]);
plot(improvements_x+0.05, improvement_means_wts, 'Color',[colors(15,:), 1]);
errorbar(improvements_x-0.05, improvement_means_all,improvement_CIs_all, 'LineStyle','none', 'Color',mean(colors)*0.6);
errorbar(improvements_x+0.05, improvement_means_wts,improvement_CIs_wts','LineStyle','none', 'Color',colors(15,:));

xlim([0.75,4.25]);
xticks([1:4]);

xlim([0.75,4.25]);
ylim([0,70]);
yticks([0:5:70]);
yticklabel_array = string([0:5:70]);
yticklabel_array(2:2:end)="";

yticklabels(yticklabel_array);


xticklabels(["b", "fl", "sf ", "h-h"]);
set(gca, 'YGrid', 'on', 'XGrid', 'off');
box off;
ylabel("Error (mm)");
legend(["Average of all landmarks", "Average of wingtips only"]);

ax=gca;
ax.FontSize=9;
ax.FontName='SansSerif';
improvementsfigure.Units = 'centimeters';
improvementsfigure.OuterPosition = [.1 .1 8.1 10.1];

figure(resultsCircleFig);
set(gca,'xlim', [165, 360]);
set(gca,'ylim', [75, 200]);
resultsCircleFig.Units = 'centimeters';
resultsCircleFig.OuterPosition(3:4) = [18 14];