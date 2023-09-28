function xyzFiltered = filter3Dpredictions(xyzInputCoords, DLTdvResidual, meanLikelihoodDLC,frameRate, acceptableSpeeds, maxSpeeds, acceptableDLTdvResidual, maxDLTdvResidual, acceptableDLClikelihood, minDLClikelihood)

% get speed of input points (xyzInputCoords)
diffMeanxyzCoord = diff(xyzInputCoords, [], 1);
diffMeanxyzCoord(2:end+1,:,:) = diffMeanxyzCoord(1:end,:,:);
diffMeanxyzCoord(1,:,:) = nan;

for i=1:size(diffMeanxyzCoord,1)
    for j = 1:size(xyzInputCoords,2)
        speed(i,j) = norm(reshape(diffMeanxyzCoord(i, j, :),1,3));
    end
end

speed = speed * frameRate;

acceptableSpeed = (speed < acceptableSpeeds);
badSpeed = speed > maxSpeeds;
acceptableResid = DLTdvResidual < acceptableDLTdvResidual;
badResid = DLTdvResidual > maxDLTdvResidual;
acceptableLikel = meanLikelihoodDLC > acceptableDLClikelihood;
badLikel = meanLikelihoodDLC < minDLClikelihood;


goodPoints = zeros(size(speed));

% check that no measures of quality are bad and that at least two are
% acceptable

goodPoints(...
    ~badLikel & ...
    ~badResid & ...
    ~badSpeed & ...
    ((acceptableSpeed & acceptableResid) | ...
    (acceptableSpeed & acceptableLikel) | ...
    (acceptableResid & acceptableLikel)) ...
    ) = true;

% reform to match sice of positonal data (spatial data)
goodPoints = logical(goodPoints);
goodPoints = cat(3,goodPoints, goodPoints, goodPoints);

% good points remain, bad points are now nans
xyzFiltered = nan(size(xyzInputCoords));
xyzFiltered(goodPoints) = xyzInputCoords(goodPoints);
