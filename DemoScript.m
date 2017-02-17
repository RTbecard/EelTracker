%% Set tracking parameters
startFrame = 150;
endFrame = 600;
video = 'ExampleVideo/sprint_20160722_024359_PM_20160722_025355_PM.avi';

%% Start tracking
buffer = 400;  % load multiple video frames at a time (increases processing speed)
refOffset = 150; % The number of frames it takes for the fish to leave the bounding box
% How many frames ahead (or behind) the reference image will be grabbed from
% The reference frames will be __ +/- the current frame
midlineResolution = 3; % units (in pixels) between each midline point
loessSmooth = 0.25; % Intensity of smoothing as a percent of total points
maxMidlineLength = 150; %% maximum length of straight midline (used for object detection)
frameInterval = 3; % analyse every __ frame
thresh = 0.02; % sensitivity for object detection (between 0 and 1).  Lower to increase sensitivity
minChunkSize = 5; % The smnallest size (in pixels) a detected object can have

% File name paramaters
Individual = '1';
Treatment = 'A';

writeVideo = 1;  %% save video output, 1 = yes, 0 = no

Tracking(video,startFrame,endFrame,refOffset,midlineResolution,loessSmooth,...
    maxMidlineLength,buffer,frameInterval,thresh,Individual,...
    Treatment,writeVideo,minChunkSize);

%% Write tracking results to a new video file
% This function can server as a good example on how to use the tracking
% data on a frame-by-frame basis.
file = 'TrackingResuts1A_17-Feb-2017.csv';
figure(1) % Manually resize figure before creating video
writeEelPreview(file,10)