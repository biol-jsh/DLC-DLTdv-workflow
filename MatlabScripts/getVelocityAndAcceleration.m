function [velocity,acceleration] = getVelocityAndAcceleration(trajetory, frameRate, plot_on)

if size(trajetory,2) == 1
    trajetory = reshape(trajetory, size(trajetory,1), size(trajetory,3));
end

velocity = diff(trajetory, [], 1);
velocity(2:end+1,:) = velocity(1:end,:);
velocity(1,:) = nan;
for i = 1:3
    velocity(:,i) = fillmissing(velocity(:,i), 'spline');
end
velocity = velocity * frameRate;

acceleration = diff(velocity, [], 1);
acceleration(2:end+1,:) = acceleration(1:end,:);
acceleration(1,:) = nan;
acceleration = acceleration * frameRate;

if(plot_on)
    for i=1:size(acceleration,1)
        normVel(i)= norm(velocity(i,:));
        normAcc(i)= norm(acceleration(i,:));
    end
    quiverResolution = 15;
    figure; hold on; axis equal
    plot3(trajetory(:,1),trajetory(:,2),trajetory(:,3));
    quiver3(trajetory(1:quiverResolution:end,1),trajetory(1:quiverResolution:end,2), trajetory(1:quiverResolution:end,3), velocity(1:quiverResolution:end,1), velocity(1:quiverResolution:end,2), velocity(1:quiverResolution:end,3));
    quiver3(trajetory(1:quiverResolution:end,1),trajetory(1:quiverResolution:end,2), trajetory(1:quiverResolution:end,3), acceleration(1:quiverResolution:end,1), acceleration(1:quiverResolution:end,2), acceleration(1:quiverResolution:end,3));        
end
end

