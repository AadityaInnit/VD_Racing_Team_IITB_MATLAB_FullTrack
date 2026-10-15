function track = make_full_track(xy,closed,nStations)
%MAKE_FULL_TRACK Smooth a prescribed driving line and calculate signed curvature.
% xy: ordered N-by-2 XY points [m]; closed: one lap or open course.
% Do not pass sharp corners: the spline defines a smooth driving line.
arguments
    xy {mustBeNumeric,mustBeReal,mustBeFinite}
    closed (1,1) logical = true
    nStations (1,1) {mustBeNumeric,mustBeReal,mustBeInteger,mustBeGreaterThan(nStations,20)} = 600
end
if size(xy,2) ~= 2 || size(xy,1) < 5
    error('full_track:Track','Provide at least five ordered XY points.');
end
if closed && norm(xy(end,:)-xy(1,:)) < 1e-8, xy(end,:) = []; end
if any(vecnorm(diff(xy),2,2) < 1e-8)
    error('full_track:Track','Consecutive track points must be distinct.');
end
if closed
    xy = [xy;xy(1,:)];
end
q = [0;cumsum(vecnorm(diff(xy),2,2))];
if closed
    % Three repeated laps give a smooth seam without Spline Toolbox.
    L = q(end); qi = [q(1:end-1)-L;q;q(2:end)+L];
    yi = [xy(1:end-1,:);xy;xy(2:end,:)];
else
    qi = q; yi = xy;
end
px = spline(qi,yi(:,1)); py = spline(qi,yi(:,2));
u = linspace(0,q(end),max(4000,8*nStations)).';
x = ppval(px,u); y = ppval(py,u);
arc = [0;cumsum(hypot(diff(x),diff(y)))];
s = linspace(0,arc(end),nStations+1).';
u = interp1(arc,u,s);
xp = ppval(derivative(px),u); yp = ppval(derivative(py),u);
xpp = ppval(derivative(derivative(px)),u);
ypp = ppval(derivative(derivative(py)),u);
if any(hypot(xp,yp) < 1e-6)
    error('full_track:Track','Spline has a stationary point; revise the driving line.');
end
track.s = s; track.X = ppval(px,u); track.Y = ppval(py,u);
track.heading = unwrap(atan2(yp,xp));
track.curvature = (xp.*ypp-yp.*xpp)./(xp.^2+yp.^2).^(3/2);
track.length = s(end); track.closed = closed;
track.sourceXY = xy;
if closed
    % Make the finite repeated-spline seam exactly consistent for interpolation.
    track.X(end) = track.X(1); track.Y(end) = track.Y(1);
    track.curvature(end) = track.curvature(1);
    turn = round((track.heading(end)-track.heading(1))/(2*pi));
    track.heading(end) = track.heading(1)+2*pi*turn;
end
end

function dp = derivative(p)
[b,a,n,k,d] = unmkpp(p);
if k <= 1
    dp = mkpp(b,zeros(n*d,1),d);
else
    dp = mkpp(b,a(:,1:k-1).*(k-1:-1:1),d);
end
end
