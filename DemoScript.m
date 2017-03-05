%% Small fish 1
video = 'ExampleVideo/sprint_20140630_121833_PM_20140630_012519_PM.avi';
framerate = 10;
startFrame = 5*framerate;
endFrame = 26*framerate;
loessSmooth = 9; % Intensity of smoothing (points on either side to use in smoothing)
buffer = 400;  % load multiple video frames at a time (increases processing speed)
refOffset = 5*framerate; % The number of frames it takes for the fish to leave the bounding box
% How many frames ahead (or behind) the reference image will be grabbed from
% The reference frames will be __ +/- the current frame
midlineResolution = 30; % units (in pixels) between each midline point
maxMidlineLength = 80; %% maximum length of straight midline (used for object detection)
frameInterval = 3; % analyse every __ frame
thresh = 0.1; % sensitivity for object detection (between 0 and 1).  Lower to increase sensitivity
minChunkSize = 5; % The smallest size (in pixels) a detected object can have

% File name paramaters
Individual = '1';
Treatment = 'A';

%% Large fish 1
video = 'ExampleVideo/sprint_20150521_012239_PM_20150521_020457_PM.avi';
framerate = 10;
startFrame = 15*framerate;
endFrame = 22*framerate;
loessSmooth = 9; % Intensity of smoothing (points on either side to use in smoothing)
buffer = 400;  % load multiple video frames at a time (increases processing speed)
refOffset = 15*framerate; % The number of frames it takes for the fish to leave the bounding box
% How many frames ahead (or behind) the reference image will be grabbed from
% The reference frames will be __ +/- the current frame
midlineResolution = 40; % units (in pixels) between each midline point
maxMidlineLength = 320; %% maximum length of straight midline (used for object detection)
frameInterval = 3; % analyse every __ frame
thresh = 0.28; % sensitivity for object detection (between 0 and 1).  Lower to increase sensitivity
minChunkSize = 50; % The smallest size (in pixels) a detected object can have
% File name paramaters
Individual = '1';
Treatment = 'B';

%% Eel
video = 'ExampleVideo/sprint_20160722_024359_PM_20160722_025355_PM.avi';
framerate = 10;
startFrame = 31*framerate;
endFrame = 62*framerate;
loessSmooth = 9; % Intensity of smoothing (points on either side to use in smoothing)
buffer = 1;  % load multiple video frames at a time (increases processing speed)
refOffset = 15*framerate; % The number of frames it takes for the fish to leave the bounding box
% How many frames ahead (or behind) the reference image will be grabbed from
% The reference frames will be __ +/- the current frame
midlineResolution = 40; % units (in pixels) between each midline point
frameInterval = 3; % analyse every __ frame
thresh = 0.04; % sensitivity for object detection (between 0 and 1).  Lower to increase sensitivity
% File name paramaters
Individual = '1';
Treatment = 'C';


%% Run Tracking script
writeVideo = 0;  %% save video output, 1 = yes, 0 = no
Tracking(video,startFrame,endFrame,refOffset,midlineResolution,loessSmooth,...
    buffer,frameInterval,thresh,Individual,...
    Treatment,writeVideo);

%% Write tracking results to a new video file
% This function can server as a good example on how to use the tracking
% data on a frame-by-frame basis.

%file = 'TrackingResuts1A_17-Feb-2017.csv';
%figure(1) % Manually resize figure before creating video
%writeEelPreview(file,10)