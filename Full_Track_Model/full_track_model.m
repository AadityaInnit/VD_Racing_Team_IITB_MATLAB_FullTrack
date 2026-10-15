function out = full_track_model(c,track,options)
%FULL_TRACK_MODEL Four-wheel planar RWD simulation and matching GGV physics.
% States [Vx Vy r s ey epsi delta]; positive y, steer and yaw point left.
% Wheel order FL FR RL RR. All planner/plant/trim/GGV physics is local here.
% Build GGV: L = full_track_model(c,[],BuildGGVOnly=true).
% Reuse GGV: o = full_track_model(c,track,GGV=L).
% Probe plant: d = full_track_model(c,[],PlantProbe=struct('state',[u;v;r;delta],'force',F)).
arguments
    c (1,1) struct
    track = []
    options.GGV (1,1) struct = struct()
    options.BuildGGVOnly (1,1) logical = false
    options.PlantProbe (1,1) struct = struct()
end
validateConfig(c);
q = capacities(c);
if isfield(options.PlantProbe,'state')
    z = options.PlantProbe.state;
    validateattributes(z,{'numeric'},{'real','finite','vector','numel',4});
    validateattributes(options.PlantProbe.force,{'numeric'},{'real','finite','scalar'});
    out = plant(z(1),z(2),z(3),z(4),options.PlantProbe.force,c,q);
    return
end
L = options.GGV;
if ~isfield(L,'ax_acc'), L = buildGGV(c,q); end
if ~isequaln(L.signature,signature(c))
    error('full_track:GGVMismatch','Vehicle/tyre/steering/drivetrain settings changed: rebuild the GGV.');
end
if options.BuildGGVOnly, out = L; return; end
validateTrack(track,c);
plan = speedPlan(c,track,L);
trim = trimMap(c,track,plan,q,L);
k0 = track.curvature(1);
if c.startMode == "flying"
    V0 = plan.speed(1); beta0 = trim.beta(1); delta0 = trim.delta(1);
else
    V0 = min(c.entrySpeed,plan.speed(1)); beta0 = 0;
    delta0 = trim.delta(1);
end
beta0 = beta0+c.initialSideslipOffset;
x0 = [V0*cos(beta0);V0*sin(beta0);V0*k0+c.initialYawRateOffset; ...
    0;c.initialLateralError;-beta0;delta0];
opts = odeset('RelTol',c.relTol,'AbsTol',c.absTol,'MaxStep',c.maxStep, ...
    'Events',@(t,x) stopEvents(t,x,c,track));
[t,x,~,~,events] = ode45(@(t,x) dynamics(t,x,c,track,plan,trim,q),[0 c.maxTime],x0,opts);
n = numel(t); D = zeros(n,12);
Fz = zeros(n,4); alpha = Fz; Fx = Fz; Fy = Fz; steer = Fz;
for j = 1:n
    [dx,d] = dynamics(t(j),x(j,:).',c,track,plan,trim,q);
    V = hypot(x(j,1),x(j,2)); bd = 0; ka = NaN;
    if V > 1e-6, bd = (x(j,1)*dx(2)-x(j,2)*dx(1))/V^2; end
    if V > 0.25, ka = (x(j,3)+bd)/V; end
    D(j,:) = [dx(1:3).' bd ka d.FX d.FY d.Mz d.targetSpeed ...
        d.utilization d.loadResidual d.steerCommand];
    Fz(j,:) = d.Fz.'; alpha(j,:) = d.alpha.';
    Fx(j,:) = d.Fx.'; Fy(j,:) = d.Fy.'; steer(j,:) = d.wheelSteer.';
end
out.config = c; out.track = track; out.plan = plan; out.trim = trim; out.ggv = L;
out.t = t; out.Vx = x(:,1); out.Vy = x(:,2); out.r = x(:,3);
out.s = x(:,4); out.e_y = x(:,5); out.e_psi = x(:,6); out.delta = x(:,7);
out.speed = hypot(out.Vx,out.Vy); out.beta = atan2(out.Vy,out.Vx);
out.VxDot = D(:,1); out.VyDot = D(:,2); out.rDot = D(:,3);
out.betaDot = D(:,4); out.curvature = D(:,5);
out.FX = D(:,6); out.FY = D(:,7); out.Mz = D(:,8);
out.targetSpeed = D(:,9); out.tireUtilization = D(:,10);
out.loadResidual = D(:,11); out.steerCommand = D(:,12);
out.axBody = out.FX/c.vehicle.m; out.ayBody = out.FY/c.vehicle.m;
out.Fz = Fz; out.slipAngles = alpha; out.FxWheel = Fx; out.FyWheel = Fy;
out.wheelSteer = steer;
sr = min(max(out.s,0),track.length);
heading = interp1(track.s,track.heading,sr);
out.X = interp1(track.s,track.X,sr,'spline')-out.e_y.*sin(heading);
out.Y = interp1(track.s,track.Y,sr,'spline')+out.e_y.*cos(heading);
out.psi = heading+out.e_psi; out.course = out.psi+out.beta;
out.referenceCurvature = interp1(track.s,track.curvature,sr);
out.tireFitRangeExceeded = any(Fz < 200 | Fz > 1100,2) ...
    | any(abs(alpha) > deg2rad(30),2);
out.wheelLiftPredicted = any(Fz < 1,2);
out.completed = any(events == 1); out.lapTime = NaN;
if out.completed
    out.lapTime = t(end); out.status = "completed";
elseif any(events == 2), out.status = "failed: track boundary exceeded";
elseif any(events == 3), out.status = "failed: sideslip limit exceeded";
elseif any(events == 4), out.status = "failed: reversed body motion";
elseif any(events == 5), out.status = "failed: Frenet chart singularity";
else, out.status = "failed: time limit reached";
end
if any(out.loadResidual > 1e-3)
    out.completed = false; out.lapTime = NaN;
    out.status = "failed: load-transfer iteration did not converge";
end
if c.makePlots, plotResult(out); end
end

function validateConfig(c)
p = c.vehicle;
positive = [p.m p.g p.wb p.a p.b p.Iz p.track_f p.track_r p.rho ...
    p.tire.P c.powerLimit c.maxTime c.maxStep c.relTol c.absTol ...
    c.trackHalfWidth c.steerLimit c.steerRateLimit c.steerTimeConstant ...
    c.maxSideslip c.lowSpeedBlend c.ggvAccelerationBound];
validateattributes(positive,{'numeric'},{'real','finite','positive'});
validateattributes([p.h_cg p.ClA p.CdA c.entrySpeed],{'numeric'},{'real','finite','nonnegative'});
validateattributes([p.rwb p.aero_bal_f c.frontRollTransferFraction ...
    c.ackermannFraction c.brakeBiasFront],{'numeric'},{'real','finite','>=',0,'<=',1});
validateattributes(c.longitudinalMargin,{'numeric'},{'real','finite','scalar','positive','<=',1});
validateattributes(c.speedScale,{'numeric'},{'real','finite','scalar','positive','<=',1});
validateattributes([c.initialSideslipOffset c.initialYawRateOffset c.initialLateralError ...
    c.courseGain c.lateralGain c.yawGain c.speedGain p.tire.gamma_f p.tire.gamma_r], ...
    {'numeric'},{'real','finite'});
validateattributes(c.ggvSpeeds,{'numeric'},{'real','finite','vector','positive','increasing'});
validateattributes(c.ggvLateral,{'numeric'},{'real','finite','vector','increasing'});
if numel(c.ggvSpeeds) < 3 || numel(c.ggvLateral) < 3 || ~any(c.ggvLateral == 0)
    error('full_track:Grid','GGV needs at least three speeds and three lateral samples, including zero.');
end
if c.ggvLateral(1) ~= -1 || c.ggvLateral(end) ~= 1
    error('full_track:Grid','GGV must cover both left and right turns.');
end
if abs(p.a+p.b-p.wb)>1e-9 || abs(p.a/p.wb-p.rwb)>1e-9
    error('full_track:Geometry','Keep a, b, wheelbase and rear weight fraction consistent.');
end
if c.driveType ~= "RWD" || ~any(c.differential == ["open" "independent"])
    error('full_track:Drive','Use RWD with open or independent rear force allocation.');
end
if ~any(c.startMode == ["standing" "flying"])
    error('full_track:Start','Use standing or flying start.');
end
if abs(c.initialLateralError) >= c.trackHalfWidth
    error('full_track:Initial','Initial position must lie inside the corridor.');
end
end

function validateTrack(t,c)
if ~isstruct(t) || ~all(isfield(t,{'s','X','Y','heading','curvature','length','closed'}))
    error('full_track:Track','Create a track with make_full_track.');
end
validateattributes(t.s,{'numeric'},{'real','finite','vector','increasing'});
if numel(t.s) < 21 || t.s(1) ~= 0 || abs(t.s(end)-t.length)>1e-8
    error('full_track:Track','Track stations must run from zero to track length.');
end
for name = ["X" "Y" "heading" "curvature"]
    validateattributes(t.(name),{'numeric'},{'real','finite','vector','numel',numel(t.s)});
end
if any(abs(t.curvature)*c.trackHalfWidth >= 0.8)
    error('full_track:Track','Corridor is too wide for the local curvature and Frenet chart.');
end
end

function s = signature(c)
% Numerical grids are included so coarse caches cannot silently replace finer ones.
s = struct('vehicle',c.vehicle,'power',c.powerLimit,'drive',c.driveType, ...
    'differential',c.differential,'roll',c.frontRollTransferFraction, ...
    'ackermann',c.ackermannFraction,'brakeBias',c.brakeBiasFront, ...
    'steerLimit',c.steerLimit,'lowSpeed',c.lowSpeedBlend, ...
    'speeds',c.ggvSpeeds,'lateral',c.ggvLateral,'axBound',c.ggvAccelerationBound);
end

function q = capacities(c)
p = c.vehicle.tire;
z = linspace(200,2000,300).'; angles = deg2rad(linspace(-30,30,401));
kappas = linspace(-0.198,0.13,401);
f = max(abs(pacejka_MF_core(angles,z,p.gamma_f,p.P)),[],2);
r = max(abs(pacejka_MF_core(angles,z,p.gamma_r,p.P)),[],2);
[~,~,fx] = pacejka_MF_core(0,z,p.gamma_f,p.P,kappas);
[~,~,rx] = pacejka_MF_core(0,z,p.gamma_r,p.P,kappas);
q.FzMin = z(1); q.FzMax = z(end);
q.Fyf = griddedInterpolant(z,f,'linear','nearest');
q.Fyr = griddedInterpolant(z,r,'linear','nearest');
q.Fxdf = griddedInterpolant(z,max(fx,[],2),'linear','nearest');
q.Fxdr = griddedInterpolant(z,max(rx,[],2),'linear','nearest');
q.Fxbf = griddedInterpolant(z,abs(min(fx,[],2)),'linear','nearest');
q.Fxbr = griddedInterpolant(z,abs(min(rx,[],2)),'linear','nearest');
end

function d = plant(u,w,r,delta,force,c,q)
p = c.vehicle; V = hypot(u,w);
xw = [p.a;p.a;-p.b;-p.b]; yw = [p.track_f;-p.track_f;p.track_r;-p.track_r]/2;
wheelSteer = zeros(4,1);
% Ackermann: cot(delta_i)=cot(delta_centre)-y_i/wheelbase.
ideal = atan2(tan(delta),1-yw(1:2)*tan(delta)/p.wb);
wheelSteer(1:2) = (1-c.ackermannFraction)*delta+c.ackermannFraction*ideal;
wheelSteer = min(max(wheelSteer,-c.steerLimit),c.steerLimit);
uw = u-r*yw; vw = w+r*xw;
cs = cos(wheelSteer); sn = sin(wheelSteer);
uwheel = uw.*cs+vw.*sn; vwheel = -uw.*sn+vw.*cs;
alpha = atan2(vwheel,max(uwheel,0.25));
blend = V/sqrt(V^2+c.lowSpeedBlend^2);
gamma = [p.tire.gamma_f;p.tire.gamma_f;p.tire.gamma_r;p.tire.gamma_r];
down = 0.5*p.rho*p.ClA*V^2;
axle = [p.m*p.g*(1-p.rwb)+p.aero_bal_f*down; ...
    p.m*p.g*p.rwb+(1-p.aero_bal_f)*down];
Fz = [axle(1);axle(1);axle(2);axle(2)]/2;
drag = 0.5*p.rho*p.CdA*V*[u;w];
context = struct('axle',axle,'alpha',alpha,'force',force,'uwheel',uwheel, ...
    'blend',blend,'gamma',gamma,'cs',cs,'sn',sn);
for j = 1:100
    target = loadTarget(Fz,context,c,q);
    residual = max(abs(target-Fz));
    if residual < 1e-5, break; end
    Fz = 0.8*Fz+0.2*target;
end
if residual > 1e-4
    % A bounded root solve handles stiff load/combined-force feedback.
    scale = sum(axle);
    opts = optimoptions('lsqnonlin','Display','off','MaxIterations',80, ...
        'FiniteDifferenceStepSize',1e-5,'FiniteDifferenceType','central', ...
        'FunctionTolerance',1e-14,'OptimalityTolerance',1e-11);
    z = lsqnonlin(@(z) (loadTarget(z*scale,context,c,q)-z*scale)/scale, ...
        Fz/scale,zeros(4,1),ones(4,1),opts);
    Fz = z*scale;
end
residual = max(abs(loadTarget(Fz,context,c,q)-Fz));
% Re-evaluate at final loads so reported forces and normal loads match exactly.
[Fx,Fy,util] = tireForces(Fz,alpha,force,uwheel,blend,gamma,c,q);
xb = Fx.*cs-Fy.*sn; yb = Fx.*sn+Fy.*cs;
d.FX = sum(xb)-drag(1); d.FY = sum(yb)-drag(2);
d.Mz = sum(xw.*yb-yw.*xb);
d.Fz = Fz; d.alpha = alpha; d.Fx = Fx; d.Fy = Fy;
d.xBody = xb; d.yBody = yb; d.wheelSteer = wheelSteer;
d.wheelLongitudinalSpeed = uwheel;
d.utilization = max(util); d.loadResidual = residual;
d.loadConverged = residual < 1e-3;
end

function target = loadTarget(Fz,d,c,q)
p = c.vehicle;
[Fx,Fy] = tireForces(Fz,d.alpha,d.force,d.uwheel,d.blend,d.gamma,c,q);
xb=Fx.*d.cs-Fy.*d.sn; yb=Fx.*d.sn+Fy.*d.cs;
pitch=sum(xb)*p.h_cg/p.wb;
rollMoment=sum(yb)*p.h_cg;
frontRoll=c.frontRollTransferFraction*rollMoment/p.track_f;
rearRoll=(1-c.frontRollTransferFraction)*rollMoment/p.track_r;
ft=min(max(d.axle(1)-pitch,0),sum(d.axle)); rt=sum(d.axle)-ft;
target=[ft/2-frontRoll;ft/2+frontRoll;rt/2-rearRoll;rt/2+rearRoll];
target(1)=min(max(target(1),0),ft); target(2)=ft-target(1);
target(3)=min(max(target(3),0),rt); target(4)=rt-target(3);
end

function [Fx,Fy,util] = tireForces(Fz,alpha,force,uwheel,blend,gamma,c,q)
z = min(max(Fz,q.FzMin),q.FzMax); scale = min(1,Fz/q.FzMin);
ycap = [q.Fyf(z(1:2));q.Fyr(z(3:4))].*scale;
if force >= 0
    xcap = [q.Fxdf(z(1:2));q.Fxdr(z(3:4))].*scale;
    force = min(force,c.powerLimit/max(0.25,mean(uwheel(3:4))));
    Fx = zeros(4,1);
    if c.differential == "open"
        Fx(3:4) = min(force/2,min(xcap(3:4)));
    else
        % Ideal independent rear force allocation, proportional to capacities.
        Fx(3:4) = min(force,sum(xcap(3:4)))*xcap(3:4)/max(sum(xcap(3:4)),1e-6);
        actualPower = sum(Fx.*max(uwheel,0));
        Fx = Fx*min(1,c.powerLimit/max(actualPower,1e-6));
    end
else
    xcap = [q.Fxbf(z(1:2));q.Fxbr(z(3:4))].*scale;
    weights = [c.brakeBiasFront;c.brakeBiasFront;1-c.brakeBiasFront;1-c.brakeBiasFront]/2;
    Fx = -min(abs(force)*weights,xcap);
end
% Clamp load for model evaluation and scale unloaded tyres down to zero.
pure = pacejka_MF_core(alpha,z,gamma,c.vehicle.tire.P).*scale*blend;
pure = min(max(pure,-ycap),ycap);
Fy = pure.*sqrt(max(0,1-(Fx./max(xcap,1e-6)).^2));
util = (Fx./max(xcap,1e-6)).^2+(Fy./max(ycap,1e-6)).^2;
end

function L = buildGGV(c,q)
% All physics comes from the SAME plant used for transient integration.
% Local quasi-steady GGV uses constant beta and zero yaw moment.
v = c.ggvSpeeds(:); ay = c.ggvLateral(:);
L.v = v; L.fractions = ay; L.ay = nan(numel(ay),numel(v)); L.ax_acc = nan(numel(ay),numel(v)); L.ax_brk = L.ax_acc;
L.ay_left = nan(size(v)); L.ay_right = L.ay_left;
L.residual_acc = L.ax_acc; L.residual_brk = L.ax_acc;
L.state_acc = nan(numel(ay),numel(v),4); L.state_brk = L.state_acc;
for iv = 1:numel(v)
    fprintf('Four-wheel GGV: %.2f m/s (%d/%d)\n',v(iv),iv,numel(v));
    for turn = [-1 1]
        [limit,residual] = lateralBoundary(v(iv),turn,c,q);
        fprintf(' lateral sign %+d: ay %.3f, residual %.3g\n',turn,limit,residual);
        if turn>0, L.ay_left(iv)=limit; else, L.ay_right(iv)=limit; end
    end
    % Start with zero lateral demand; continuation moves outward in each direction.
    [~,order] = sort(abs(ay)); seeds = [0;0;0;0];
    for mode = [1 -1]
        previous = repmat(seeds,1,2);
        for ia = order.'
            turn = sign(ay(ia)); if turn == 0, turn = 1; end
            side = 1+(turn<0);
            limit = L.ay_left(iv); if turn < 0, limit = L.ay_right(iv); end
            if ~isfinite(limit), continue; end
            demand = ay(ia)*0.95*limit; L.ay(ia,iv) = demand;
            guess = previous(:,side);
            [z,flag,residual] = longitudinalBoundary(v(iv),demand,mode,guess,c,q);
            ci = -1;
            if flag > 0 && residual < 1e-4 && max(ci) < 1e-3
                previous(:,side) = z;
                if mode > 0
                    L.ax_acc(ia,iv) = z(4)*c.vehicle.g; L.residual_acc(ia,iv) = residual;
                    L.state_acc(ia,iv,:) = z.*[1;1;c.vehicle.m*c.vehicle.g;c.vehicle.g];
                else
                    L.ax_brk(ia,iv) = z(4)*c.vehicle.g; L.residual_brk(ia,iv) = residual;
                    L.state_brk(ia,iv,:) = z.*[1;1;c.vehicle.m*c.vehicle.g;c.vehicle.g];
                end
            end
        end
    end
end
L.signature = signature(c); L.speed_frame = "path"; L.acceleration_frame = "path";
L.description = "Four independent wheel kinematics, same force-controlled RWD plant";
if any(~isfinite(L.ay_left)) || any(~isfinite(L.ay_right))
    error('full_track:GGV','A steady lateral envelope solve failed; inspect grid or vehicle settings.');
end
end

function [best,flag,residual] = longitudinalBoundary(V,ay,mode,guess,c,q)
% Continue force demand along the feasible low-slip equilibrium branch.
% Retain the greatest signed tangential acceleration; stop at saturation.
fun=@(z) envelopeConstraint(z,V,ay,1,c,q,false);
lb=[-c.maxSideslip;-c.steerLimit;-c.ggvAccelerationBound/c.vehicle.g];
ub=[c.maxSideslip;c.steerLimit;c.ggvAccelerationBound/c.vehicle.g];
opts=optimoptions('lsqnonlin','Display','off','MaxIterations',120, ...
    'FiniteDifferenceStepSize',1e-4,'FiniteDifferenceType','central', ...
    'FunctionTolerance',1e-12,'OptimalityTolerance',1e-9);
guess(1)=atan(c.vehicle.b*ay/V^2);
guess(2)=atan(c.vehicle.wb*ay/V^2);
guess(2)=min(max(guess(2),-c.steerLimit),c.steerLimit);
y=guess([1 2 4]); best=nan(4,1); flag=-1; residual=Inf;
lo=0; hi=5; last=y; bestScore=-Inf; flat=0;
for force=0:0.15:5
    [y,~,res]=lsqnonlin(@(y) fixedForceResidual(y,mode*force,fun),last,lb,ub,opts);
    candidate=[y(1:2);mode*force;y(3)]; [ci,~]=fun(candidate);
    if max(abs(res))>1e-5 || max(ci)>1e-5
        hi=force; break;
    end
    score=mode*y(3);
    if score>bestScore
        improvement=score-bestScore; bestScore=score; best=candidate;
        residual=max(abs(res)); flag=1;
    else
        improvement=0;
    end
    last=y; lo=force;
    if improvement<1e-5, flat=flat+1; else, flat=0; end
    if flat>=2, return; end
end
if flag<0, return; end
for j=1:16
    force=(lo+hi)/2;
    [y,~,res]=lsqnonlin(@(y) fixedForceResidual(y,mode*force,fun),last,lb,ub,opts);
    candidate=[y(1:2);mode*force;y(3)]; [ci,~]=fun(candidate);
    if max(abs(res))<1e-5 && max(ci)<1e-5
        lo=force; last=y;
        if mode*y(3)>bestScore
            bestScore=mode*y(3); best=candidate; residual=max(abs(res));
        end
    else
        hi=force;
    end
    if hi-lo<1e-4, break; end
end
end

function e = fixedForceResidual(y,force,fun)
[~,e] = fun([y(1:2);force;y(3)]);
end

function [limit,residual] = lateralBoundary(V,turn,c,q)
% Continue the low-slip steady branch and bracket its first feasible boundary.
% This avoids jumping to disconnected post-peak tyre equilibria.
F=5; lb=[-c.maxSideslip;-c.steerLimit;0]; ub=[c.maxSideslip;c.steerLimit;F];
opts=optimoptions('lsqnonlin','Display','off','MaxIterations',150, ...
    'FiniteDifferenceStepSize',1e-4,'FiniteDifferenceType','central', ...
    'FunctionTolerance',1e-12,'OptimalityTolerance',1e-9);
last=[0;0;0.02]; lo=0; hi=3.5*c.vehicle.g;
step=max(0.01,min(0.15*c.vehicle.g,0.1*V^2/c.vehicle.wb));
residual=Inf;
for target=0:step:hi
    fun=@(y) fixedLateralResidual(y,target,turn,V,c,q);
    [z,~,res]=lsqnonlin(fun,last,lb,ub,opts);
    [ci,~]=envelopeConstraint([z;target/c.vehicle.g],V,0,turn,c,q,true);
    if max(abs(res))>1e-5 || max(ci)>1e-5
        hi=target; break;
    end
    lo=target; last=z; residual=max(abs(res));
end
if ~isfinite(residual), error('full_track:GGV','Zero-lateral trim failed.'); end
for j=1:22
    target=(lo+hi)/2;
    [z,~,res]=lsqnonlin(@(y) fixedLateralResidual(y,target,turn,V,c,q),last,lb,ub,opts);
    [ci,~]=envelopeConstraint([z;target/c.vehicle.g],V,0,turn,c,q,true);
    if max(abs(res))<1e-5 && max(ci)<1e-5
        lo=target; last=z; residual=max(abs(res));
    else
        hi=target;
    end
    if hi-lo<1e-4, break; end
end
limit=lo;
end

function e = fixedLateralResidual(y,target,turn,V,c,q)
[~,e] = envelopeConstraint([y;target/c.vehicle.g],V,0,turn,c,q,true);
end

function [ci,ce] = envelopeConstraint(z,V,ay,turn,c,q,lateral)
if lateral, ay = turn*z(4)*c.vehicle.g; ax = 0; else, ax = z(4)*c.vehicle.g; end
beta = z(1); u = V*cos(beta); w = V*sin(beta); r = ay/V;
d = plant(u,w,r,z(2),z(3)*c.vehicle.m*c.vehicle.g,c,q); p = c.vehicle;
at = (d.FX*cos(beta)+d.FY*sin(beta))/p.m;
an = (-d.FX*sin(beta)+d.FY*cos(beta))/p.m;
ce = [(at-ax)/p.g;(an-ay)/p.g;d.Mz/(p.m*p.g*p.wb)];
ci = [d.loadResidual-1e-3;max(abs(d.alpha))-deg2rad(30)];
end

function plan = speedPlan(c,t,L)
% Preserve corner cap -> forward -> reverse -> minimum -> trapezoidal time.
s = t.s(:); k = t.curvature(:); ds = diff(s); v = L.v;
cap = zeros(size(s));
for j = 1:numel(s)
    if abs(k(j)) < 1e-8, cap(j) = v(end); continue; end
    ay = L.ay_left; if k(j)<0, ay = L.ay_right; end
    f = v.^2*abs(k(j))-ay;
    first = find(f >= 0,1);
    if isempty(first), cap(j) = v(end);
    elseif first == 1, error('full_track:Corner','Track requires speed below the GGV minimum.');
    else
        cap(j) = interp1(f(first-1:first),v(first-1:first),0);
    end
end
cap = cap*c.speedScale;
forward = cap; reverse = cap;
if c.startMode == "standing", forward(1) = min(c.entrySpeed,cap(1)); end
% Closed flying laps require cyclic forward/backward constraints at the seam.
for pass = 1:30
    old = min(forward,reverse);
    if t.closed && c.startMode == "flying"
        forward(1) = min([forward(1),forward(end),reverse(1),reverse(end)]);
    end
    for j = 1:numel(ds)
        a = c.longitudinalMargin*max(0,ggvLimit(L,max(forward(j),v(1)),forward(j)^2*k(j),1));
        forward(j+1) = min([forward(j+1),reverse(j+1),sqrt(max(0,forward(j)^2+2*a*ds(j)))]);
    end
    if t.closed && c.startMode == "flying"
        reverse(end) = min([reverse(end),forward(end),forward(1)]);
    end
    for j = numel(ds):-1:1
        ceiling = min(forward(j),reverse(j));
        if ceiling <= reverse(j+1), reverse(j) = ceiling; continue; end
        lo = reverse(j+1); hi = ceiling;
        for b = 1:35
            vm = (lo+hi)/2; a = c.longitudinalMargin*min(0,ggvLimit(L,max(vm,v(1)),vm^2*k(j),-1));
            if vm^2+2*a*ds(j)>reverse(j+1)^2, hi = vm; else, lo = vm; end
        end
        reverse(j) = (lo+hi)/2;
    end
    if max(abs(min(forward,reverse)-old)) < 1e-6, break; end
end
speed = min([cap forward reverse],[],2);
at = diff(speed.^2)./(2*ds);
plan.s = s; plan.speed = speed; plan.cap = cap; plan.length = t.length;
plan.acceleration = [at;at(end)]; plan.lapTime = sum(2*ds./max(speed(1:end-1)+speed(2:end),1e-9));
plan.cyclicResidual = abs(speed(end)-speed(1));
end

function a = ggvLimit(L,V,ay,mode)
% Interpolate only within contiguous feasible rows; never nearest-fill holes.
if mode>0, A = L.ax_acc; else, A = L.ax_brk; end
V = min(max(V,L.v(1)),L.v(end));
iv = find(L.v<=V,1,'last'); iv = min(iv,numel(L.v)-1);
envelope = L.ay_left; if ay<0, envelope = L.ay_right; end
localLimit = interp1(L.v,envelope,V);
fraction = ay/max(localLimit,1e-9);
a0 = slice(L.ay(:,iv),A(:,iv),fraction*envelope(iv));
a1 = slice(L.ay(:,iv+1),A(:,iv+1),fraction*envelope(iv+1));
a = a0+(a1-a0)*(V-L.v(iv))/(L.v(iv+1)-L.v(iv));
end

function a = slice(y,z,query)
if query < y(1)-1e-9 || query > y(end)+1e-9
    error('full_track:GGVRange','Requested lateral demand exceeds the sampled envelope; reduce speedScale.');
end
query = min(max(query,y(1)),y(end));
exact = find(abs(y-query)<1e-10,1);
if ~isempty(exact) && isfinite(z(exact)), a=z(exact); return; end
j = find(y<=query,1,'last');
if isempty(j), j=1; end
j = min(j,numel(y)-1);
if ~all(isfinite(z(j:j+1)))
    error('full_track:GGVHole','Requested acceleration lies outside sampled feasible GGV; refine the grid or reduce speedScale.');
end
a = interp1(y(j:j+1),z(j:j+1),query);
end

function trim = trimMap(c,t,plan,q,L)
% Solve plant steady trim along the track; do not reset dynamic states.
% Sample up to 100 stations for affordable steering feedforward.
idx = unique(round(linspace(1,numel(t.s),min(100,numel(t.s)))));
ss = t.s(idx); beta = zeros(size(ss)); delta = beta; residual = beta;
opts = optimoptions('lsqnonlin','Display','off','MaxIterations',100, ...
    'FunctionTolerance',1e-10,'OptimalityTolerance',1e-8, ...
    'FiniteDifferenceStepSize',1e-4,'FiniteDifferenceType','central');
F = 5;
guess = [0;0;0];
for j = 1:numel(idx)
    V = max(plan.speed(idx(j)),0.5); k = t.curvature(idx(j));
    [~,iv]=min(abs(L.v-V));
    [~,ia]=min(abs(L.ay(:,iv)-V^2*k));
    if plan.acceleration(idx(j))>=0, states=L.state_acc; else, states=L.state_brk; end
    seed=reshape(states(ia,iv,:),4,1);
    if all(isfinite(seed)), guess=[seed(1:2);seed(3)/(c.vehicle.m*c.vehicle.g)]; end
    at = plan.acceleration(idx(j));
    % Commanded force in a saturated GGV state can lie on a flat plateau.
    % Start inside the actuator range so force sensitivity is available.
    guess(3) = (at*cos(guess(1))-V^2*k*sin(guess(1)))/c.vehicle.g ...
        +0.5*c.vehicle.rho*c.vehicle.CdA*V^2/(c.vehicle.m*c.vehicle.g) ...
        +V^2*k*c.vehicle.b/c.vehicle.wb*tan(guess(2))/c.vehicle.g;
    fun = @(z) trimResidual(z,V,k,at,c,q);
    [z,~,r] = lsqnonlin(fun,guess,[-c.maxSideslip;-c.steerLimit;-F], ...
        [c.maxSideslip;c.steerLimit;F],opts);
    residual(j) = max(abs(r));
    if residual(j)>1e-3, fprintf('Trim fail s %.2f V %.2f k %.4f at %.2f residual %.3g\n',ss(j),V,k,at,residual(j)); end
    beta(j) = z(1); delta(j) = z(2); guess = z;
end
if any(residual > 1e-3)
    error('full_track:Trim','Steady feedforward trim failed; lower speedScale or inspect GGV grid.');
end
trim.s = ss; trim.beta = beta; trim.delta = delta; trim.residual = residual;
end

function e = trimResidual(z,V,k,at,c,q)
p = c.vehicle; u = V*cos(z(1)); w = V*sin(z(1)); r = V*k;
d = plant(u,w,r,z(2),z(3)*p.m*p.g,c,q);
e = [((d.FX*cos(z(1))+d.FY*sin(z(1)))/p.m-at)/p.g; ...
    ((-d.FX*sin(z(1))+d.FY*cos(z(1)))/p.m-V^2*k)/p.g; ...
    d.Mz/(p.m*p.g*p.wb)];
end

function [dx,d] = dynamics(~,x,c,t,plan,trim,q)
u = x(1); w = x(2); r = x(3); s = x(4); ey = x(5); epsi = x(6); delta = x(7);
p = c.vehicle; V = hypot(u,w); beta = atan2(w,u);
sr = min(max(s,0),t.length); k = interp1(t.s,t.curvature,sr);
vtarget = interp1(plan.s,plan.speed,sr); atref = interp1(plan.s,plan.acceleration,sr);
ff = interp1(trim.s,trim.delta,sr);
courseError = atan2(sin(epsi+beta),cos(epsi+beta));
steer = ff-c.courseGain*courseError-c.lateralGain*ey-c.yawGain*(r-V*k);
steer = min(max(steer,-c.steerLimit),c.steerLimit);
deltaDot = min(max((steer-delta)/c.steerTimeConstant,-c.steerRateLimit),c.steerRateLimit);
at = atref+c.speedGain*(vtarget-V); an = V^2*k;
force = p.m*(at*cos(beta)-an*sin(beta))+0.5*p.rho*p.CdA*V^2*cos(beta) ...
    +p.m*an*p.b/p.wb*tan(delta);
d = plant(u,w,r,delta,force,c,q);
sDot = (u*cos(epsi)-w*sin(epsi))/max(0.05,1-k*ey);
eyDot = u*sin(epsi)+w*cos(epsi);
dx = [d.FX/p.m+r*w;d.FY/p.m-r*u;d.Mz/p.Iz;sDot;eyDot;r-k*sDot;deltaDot];
d.targetSpeed = vtarget; d.steerCommand = steer;
end

function [value,terminal,direction] = stopEvents(~,x,c,t)
betaMargin = c.maxSideslip-abs(atan2(x(2),x(1)));
if hypot(x(1),x(2))<0.5, betaMargin = 1; end
k = interp1(t.s,t.curvature,min(max(x(4),0),t.length));
value = [t.length-x(4);c.trackHalfWidth-abs(x(5));betaMargin;x(1)+0.1;1-k*x(5)-0.1];
terminal = ones(5,1); direction = -ones(5,1);
end

function plotResult(o)
fig=figure('Name','Four-wheel full track','Color','w','Position',[100 100 1200 850]); tiledlayout(3,2);
nexttile; plot(o.track.X,o.track.Y,'--',o.X,o.Y,'LineWidth',1.2);
axis equal; grid on; xlabel('X [m]'); ylabel('Y [m]'); legend('Reference','Actual');
nexttile; plot(o.t,[o.speed o.targetSpeed]); grid on;
xlabel('Time [s]'); ylabel('Speed [m/s]'); legend('Actual','Target');
nexttile; plot(o.t,rad2deg([o.beta o.delta])); grid on;
xlabel('Time [s]'); ylabel('Angle [deg]'); legend('Sideslip','Centre steer');
nexttile; plot(o.t,o.Fz); grid on; xlabel('Time [s]'); ylabel('Normal load [N]');
legend('FL','FR','RL','RR');
nexttile; plot(o.t,[o.r o.speed.*o.referenceCurvature]); grid on;
xlabel('Time [s]'); ylabel('Yaw rate [rad/s]'); legend('Actual','Reference');
nexttile; plot(o.t,o.e_y); grid on; xlabel('Time [s]'); ylabel('Path error [m]');
sgtitle("Four-wheel RWD: "+o.status);
set(findall(fig,'Type','axes'),'Color','w','XColor',[0.15 0.15 0.15], ...
    'YColor',[0.15 0.15 0.15],'GridColor',[0.65 0.65 0.65],'FontSize',10);
set(findall(fig,'Type','text'),'Color',[0.15 0.15 0.15]);
set(findall(fig,'Type','legend'),'Color','w','TextColor',[0.15 0.15 0.15]);
end
