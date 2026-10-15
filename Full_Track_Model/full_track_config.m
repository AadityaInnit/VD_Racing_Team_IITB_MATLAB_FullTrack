function c = full_track_config()
%FULL_TRACK_CONFIG Four-wheel RWD configuration; assumptions are listed below.
d = load(fullfile(fileparts(mfilename('fullpath')),'data','vehicle_snapshot.mat'),'ggvLUT');
c.vehicle = d.ggvLUT.vehicle;
c.vehicle.Iz = c.vehicle.m*c.vehicle.a*c.vehicle.b; % UNKNOWN: estimate, kg m^2
c.powerLimit = 70000; % Usable wheel power [W], inherited nominal value.
c.driveType = "RWD";
c.differential = "open"; % UNKNOWN: equal rear force, limited by weaker tyre.
c.frontRollTransferFraction = 0.50; % UNKNOWN: front share of lateral roll moment.
c.ackermannFraction = 1.0; % UNKNOWN: 0 parallel, 1 geometric Ackermann.
c.brakeBiasFront = 0.60; % UNKNOWN: front brake force fraction.
c.steerTimeConstant = 0.05; % UNKNOWN: steering actuator time constant [s].
c.steerRateLimit = deg2rad(180); % UNKNOWN: steering rate [rad/s].
c.steerLimit = deg2rad(35); % Virtual centre steer and individual wheel limit [rad].
c.courseGain = 1.8; % Driver tuning: rad/rad.
c.lateralGain = 0.18; % Driver tuning: rad/m.
c.yawGain = 0.12; % Driver tuning: seconds.
c.speedGain = 2.0; % Driver tuning: 1/s.
c.lowSpeedBlend = 0.20; % Numerical tyre regularization speed [m/s].
c.speedScale = 0.70; % Initial transient tracking margin on GGV corner caps.
c.longitudinalMargin = 0.70; % Driver reserve on GGV acceleration/braking limits.
c.startMode = "standing";
c.entrySpeed = 0;
c.initialSideslipOffset = 0; % rad
c.initialYawRateOffset = 0; % rad/s
c.initialLateralError = 0; % m, positive left
c.trackHalfWidth = 2; % Reference-path corridor half-width [m].
c.maxSideslip = deg2rad(35);
c.maxTime = 120;
c.maxStep = 0.02;
c.relTol = 1e-6;
c.absTol = 1e-8;
c.makePlots = true;
c.ggvSpeeds = [0.5 3 6 10 15 20 25 30]; % m/s, increase resolution as needed.
c.ggvLateral = linspace(-1,1,15); % Signed fractions of each speed's lateral envelope.
c.ggvAccelerationBound = 4*c.vehicle.g;
c.assumptions = table( ...
    ["vehicle.Iz";"frontRollTransferFraction";"differential"; ...
     "ackermannFraction";"brakeBiasFront";"steerTimeConstant"; ...
     "steerRateLimit";"steerLimit";"vehicle.h_cg";"vehicle.track_f";"vehicle.track_r"], ...
    ["Estimate m*a*b; replace with measured yaw inertia"; ...
     "Assumed 50% front roll moment; replace using suspension data"; ...
     "Assumed open differential; confirm rear differential type"; ...
     "Assumed ideal Ackermann; replace using measured steering geometry"; ...
     "Assumed 60% front braking; replace using brake system data"; ...
     "Assumed actuator response; replace using steering response data"; ...
     "Assumed rate limit; replace using measured steering rate"; ...
     "Assumed angle limit; replace using measured steering travel"; ...
     "Inherited nominal 0.30 m; verify CG height"; ...
     "Inherited nominal 1.20 m; verify front track width"; ...
     "Inherited nominal 1.25 m; verify rear track width"], ...
    'VariableNames',{'Parameter','UpdateLater'});
end
