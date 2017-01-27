%% Set tracking aprameters
startFrame = 150;
video = 'ExampleVideo/sprint_20160722_024359_PM_20160722_025355_PM.avi';

%% Start tracking
Tracking(video,startFrame)

%% Write tracking results to a new video file
% This function can server as a good example on how to use the tracking
% data on a frame-by-frame basis.
file = 'TrackingResuts27-Jan-2017.csv';
writeEelPreview(file,1600, 500, 30)
 
%% Show video in VLC (only works on linux systems with vlc installed)
system('vlc EelPreview.avi')