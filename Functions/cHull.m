function cv = cHull(bw)
%% Draw convec hull around bw image. return points
    cv = bwboundaries(bw);
    cv = cell2mat(cv);
    CH = convhull(cv(:,1),cv(:,2));
    cv = cv(CH,[2 1]);
end