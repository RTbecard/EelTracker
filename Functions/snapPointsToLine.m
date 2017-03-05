function [xOut,yOut] = snapPointsToLine(x,y,m,b)
    % Snaps points to a line
    % Points will be snapped along the shortest distance (perpendicular to the
    % line)
    %
    % x: x coorinates of the points to snap
    % y: coordinates of the points to snap
    % m: slope fo the line to snap to
    % b: y-intercept of the line to snap to

    transposedPointsTemp = (x + y*m - m*b) / (1 + m^2);
    transposedPointsTemp  = [transposedPointsTemp (m*transposedPointsTemp  + b)];
    xOut = transposedPointsTemp(:,1);
    yOut = transposedPointsTemp(:,2);
    end
