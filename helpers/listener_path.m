function [position, heading] = listener_path(t, turn, route)
% turn = 1 for left, -1 for right, and 0 for straight.
t = min(t, route.walkStop);
distance = route.speed * max(0, t - route.walkStart);
angle = min(1, max(0, (t - route.turnStart) / route.turnDuration)) * pi/2;
radius = route.speed * route.turnDuration / (pi/2);
entryLength = route.speed * (route.turnStart - route.walkStart);
exitLength = route.speed * max(0, t - route.turnStart - route.turnDuration);

x = entryLength + radius * sin(angle) + exitLength .* cos(angle);
y = turn * (radius * (1 - cos(angle)) + exitLength .* sin(angle));
x(t <= route.turnStart) = distance(t <= route.turnStart);
heading = turn * angle;

if turn == 0
    x = distance;
    y = zeros(size(t));
end
position = [x, y];
end
