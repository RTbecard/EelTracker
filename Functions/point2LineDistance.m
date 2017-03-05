function [distance] = point2LineDistance(Q1,Q2,P)
    % Calculates the distance of a point to a line
    % Q1: The x and y coordinates of a point along the line
    % Q2: The x and y coordinates of a second point along the line
    % P: The point which you'd like to calculate the distance from

    % Absolute distance
    distance= abs( det([P-Q1;Q2-Q1]) )/norm(Q2-Q1);

    % Find slope and intercept of defined line
    m = (Q1(2) - Q2(2))/(Q1(1) - Q2(1));
    b = Q1(2) - m*Q1(1);
    % Correct negative values
    if find(P(2) > (m*P(1) + b));
        distance = distance*-1;
    end
end