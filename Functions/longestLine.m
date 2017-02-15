function [longestline,D] = longestLine(XY)
    maxDist = 0;
    
    % loop through all rows
    for i = 1:size(XY,1)
        % Compare this point to each poin and measure distance
        lines = sqrt((XY(i,1) - XY(:,1)).^2 + (XY(i,2) - XY(:,2)).^2);
        [D,I] = max(lines);
        if D > maxDist
           maxDist = D;
           longestline = [XY(i,:);XY(I,:)];
        end
    end
end

