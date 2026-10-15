# Four-wheel full-track model

An independent experimental planar RWD vehicle model. The original skidpad and `Transient_Skidpad` packages are unchanged. MATLAB plus Optimization Toolbox (`lsqnonlin`) are required.

## Run

Open `run_full_track.m` and run it. The example is a smooth ellipse with varying curvature. The script builds/caches a matching four-wheel GGV, constructs the spatial speed plan, integrates the vehicle and saves the results and overview figure.

```matlab
addpath('/Users/aaditya/Desktop/VD_Racing_Team/Full_Track_Model');
c = full_track_config();
disp(c.assumptions);                 % Parameters requiring measurements/confirmation.
xy = readmatrix('my_driving_line.csv'); % Ordered [X Y] points, metres.
track = make_full_track(xy,true,600); % true: closed lap; false: open course.
ggv = full_track_model(c,[],BuildGGVOnly=true);
out = full_track_model(c,track,GGV=ggv);
```

Coordinates describe the desired CG driving line, not just track boundaries. Input points must describe a smooth route without repeated consecutive points; spline smoothing changes the line between points. Inspect the resulting shape/curvature before using measured track data. Sharp corners and self-intersecting driving lines require care; this is not a racing-line optimizer. The corridor is a constant-width tolerance around that line, not a physical boundary map.

`make_full_track` constructs a cubic spline, resamples approximately uniformly in arc length and calculates signed curvature from its derivatives. Positive curvature turns left; negative turns right. Closed tracks use three repeated laps to smooth the seam. Increase point/station resolution and check convergence for tight or complex routes.

## Files and packaging

- `full_track_config.m`: vehicle, controls, GGV grids and labelled assumptions.
- `make_full_track.m`: prescribed path geometry and signed curvature.
- `full_track_model.m`: **all tyre, load-transfer, steering, drivetrain, GGV equilibrium, trim, planner and simulation physics**, as local functions.
- `run_full_track.m`: editable track example, GGV cache and saved outputs.
- `test_full_track_model.m`: physical balances and integration checks.
- `private/pacejka_MF_core.m`: independent copy of the existing tyre fit.
- `data/vehicle_snapshot.mat`: compact snapshot of nominal vehicle parameters; no old bicycle GGV is used for planning.
- `data/four_wheel_ggv.mat`: generated envelope for this four-wheel plant.

A cached GGV must match the vehicle, tyre, differential, Ackermann, roll-transfer, brake-bias, steering-limit, low-speed regularization and numerical grid settings. A mismatch is rejected. Set `rebuildGGV=true` in the run script after changing those fields. Driver gains and yaw inertia do not determine a steady zero-yaw-moment envelope, although the conservative signature also includes inertia.

## Physics

States: `[Vx, Vy, r, s, ey, epsi, delta]`. Wheel order is **FL, FR, RL, RR**. Vehicle axes are x forward, y left, positive yaw counterclockwise. Angles are radians; lengths metres; time seconds; tyre pressure kPa.

For wheel location `(xi,yi)`, its body-frame contact velocity is `(Vx-r*yi, Vy+r*xi)`. Rotate by its individual steering angle into wheel coordinates before calculating slip angle. The tyre fit's existing slip/force sign convention is retained. Each wheel has its own slip angle and normal load; front steering follows a configurable fraction of ideal Ackermann.

The integrated body equations are:

```
VxDot = FX/m + r*Vy
VyDot = FY/m - r*Vx
rDot  = Mz/Iz
Mz    = sum(xi*Fyi - yi*Fxi)
```

Tyre forces are rotated back into the body frame before summation. Aerodynamic drag acts opposite CG velocity. Aero downforce has a prescribed front/rear split. The assumed aerodynamic application heights equal CG height for pitch-load-transfer purposes; additional aero pitching moments are omitted.

Vertical loads combine static axle weights, downforce, ground-force longitudinal transfer and lateral roll moment. `frontRollTransferFraction` splits total lateral roll moment between axles; axle track widths convert moment to inside/outside load transfer. This is a lumped replacement for unknown roll stiffness/roll-centre data, **not a suspension model**. An unloaded tyre's load is clipped to zero and the axle's remaining vertical load transferred to its partner. Loads and forces iterate to a consistent fixed point, and convergence residuals are recorded.

Pure lateral Pacejka forces are reduced using a combined-force ellipse after longitudinal allocation. RWD acceleration acts only at the rear. With equal assumed rear rolling radii, the open differential requests equal rear forces and limits both to the weaker wheel's capacity. `differential="independent"` models ideal independently allocated rear forces; it is not a locked or limited-slip differential. Power is limited at the wheels. Braking uses fixed front/rear bias and individual tyre caps. No wheel rotation, longitudinal-slip states, ABS or traction-control dynamics are included.

Low-speed lateral tyre forces blend toward zero. Capacity lookup loads are bounded at 200–2000 N and forces scale to zero below the lower boundary. This is an approximation outside measured tyre data; the stated longitudinal fit range is 200–1100 N. Outputs flag loads outside that range and slip angles beyond 30 degrees. Lateral/longitudinal fits and the combined ellipse need experimental validation.

## GGV and lap algorithm

The GGV generator calls the **same four-wheel plant** as the simulation. It solves body tangential/normal force balance and zero yaw moment with sideslip and centre steering as unknowns. It follows the low-slip steady branch to its first feasibility boundary separately for left/right turns, then follows drive/brake force demand along the feasible branch to construct acceleration/braking envelopes on signed fractions of those limits. `ggvLateral` is dimensionless, from -1 to +1; samples use 95% of the solved limit for numerical margin. Tables store actual SI acceleration coordinates in `ggv.ay`.

The lateral boundary is a conservative continuation result rather than a certified global maximum over all tyre branches; longitudinal boundaries similarly follow force demand to feasibility or saturation and retain the greatest signed acceleration found. This is local quasi-steady GGV physics: beta is constant and yaw moment is zero. It does not represent yaw acceleration during changing curvature, steering-rate limits or tyre relaxation. Those discrepancies are allowed to evolve during transient execution and can cause tracking failures.

The lap planner retains the existing sequence: curvature corner-speed caps, forward acceleration pass, backward braking pass, minimum profile and trapezoidal segment time. Flying closed laps use repeated passes to propagate seam constraints; standing starts retain an unconstrained finishing speed, as in the skidpad planner. Lookup interpolation is performed at equal fractions of the local lateral envelope across speed columns. Failed GGV samples are not nearest-filled: requests touching them fail clearly. Increase grid resolution or lower `speedScale` if the sampled envelope cannot support a requested plan.

`speedScale=0.7` is an initial transient tracking margin on corner caps, including the straight-line speed ceiling. It is configurable and is not an optimized lap-time result. The driver uses acceleration-aware equilibrium trim feedforward plus course, lateral-position and yaw feedback; it follows the plan rather than enforcing its speed or curvature. `longitudinalMargin=0.7` retains 30% reserve on the GGV acceleration/braking limits while preserving the same planning passes.

Exact Frenet kinematics include Vy:

```
sDot    = (Vx*cos(epsi)-Vy*sin(epsi))/(1-kappa*ey)
eyDot   = Vx*sin(epsi)+Vy*cos(epsi)
epsiDot = r-kappa*sDot
```

Actual CG curvature is `(r+betaDot)/hypot(Vx,Vy)` at adequate speed. The solver stops on completion, excessive path error/sideslip, reversed motion, a Frenet singularity or the time limit. Failed runs have `lapTime=NaN`; they never masquerade as completed laps.

## Parameters to update later

`c.assumptions` is the editable measurement checklist. In particular, confirm the rear differential type, yaw inertia, roll-transfer distribution, Ackermann fraction, brake bias and steering actuator limits/response. CG height and track widths are inherited nominal values needing confirmation. Power, aero coefficients/balance and tyre pressure/camber also retain nominal project values and should be checked against the actual setup. Controller gains are simulation tuning parameters.

The first version assumes a flat track, rigid planar body and instantaneous tyre forces. Suspension roll/pitch/heave dynamics, tyre relaxation, camber gain, compliance, wheel rotation, differential torque dynamics, rolling resistance, road banking/elevation and racing-line optimization are omitted intentionally from this initial scope.

## Validation

```matlab
results = runtests('Full_Track_Model/test_full_track_model.m');
assert(all([results.Passed]));
```

These tests check balances, individual wheel kinematics, RWD power allocation, matched GGV equilibrium, changing-curvature motion, clockwise laps, open left/right courses, failure reporting and time-step convergence. They establish numerical/implementation consistency, not measured vehicle accuracy.

Latest default demonstration: 181.407 m ellipse, 11.016 s simulated versus 11.004 s planned, maximum path error 0.080 m. All 12 tests passed and MATLAB code analysis reported no issues for the package entry points/tests. At least one tyre exceeds the fit's stated load range in 100% of recorded demo samples; the fit extrapolation/clamping assumptions require better tyre data before claiming measured vehicle accuracy. The initial driver margins are deliberately conservative, and more aggressive settings can fail tracking or load convergence.
