function [Fy, C, Fx, Cx] = pacejka_MF_core(alpha, Fz, gamma, P, kappa)
% PACEJKA_MF_CORE  Nonlinear Pacejka Magic-Formula tire model.
%
% Backward-compatible pure-lateral usage:
%   Fy        = pacejka_MF_core(alpha, Fz)
%   Fy        = pacejka_MF_core(alpha, Fz, gamma, P)
%   [Fy, C]   = pacejka_MF_core(alpha, Fz, gamma, P)
%
% Optional pure-longitudinal usage (added from Fx_70_16_in.m):
%   [Fy, C, Fx]     = pacejka_MF_core(alpha, Fz, gamma, P, kappa)
%   [Fy, C, Fx, Cx] = pacejka_MF_core(alpha, Fz, gamma, P, kappa)
%
% INPUTS (scalars or arrays of matching/broadcastable size):
%   alpha - lateral slip angle [rad]
%   Fz    - vertical (normal) load per tire [N], positive
%   gamma - inclination/camber angle [rad] (default: 0)
%   P     - inflation pressure [kPa] (default: 70)
%   kappa - longitudinal slip ratio [-] (optional; if omitted Fx = [])
%
% OUTPUTS:
%   Fy - pure-slip lateral force [N]
%   C  - lateral diagnostic coefficient struct
%   Fx - pure-slip longitudinal force [N] if kappa is supplied, else []
%   Cx - longitudinal diagnostic coefficient struct if kappa is supplied,
%        else an empty struct
%
% NOTES:
%   * The lateral model is the original MF 6.1 fit already used by this
%     project.
%   * The longitudinal branch is taken directly from Fx_70_16_in.m.
%   * Fx_70_16_in.m was fitted at P = 70 kPa and IA = 0 deg. Therefore the
%     longitudinal model currently has no pressure correction and no camber
%     dependence. gamma and P are retained for interface consistency only.
%   * Lateral and longitudinal forces are PURE-SLIP outputs. This function
%     does NOT yet implement Pacejka combined-slip weighting functions Gx/Gy.
%
% Longitudinal fit validity stated in Fx_70_16_in.m:
%   kappa approximately -0.198 to 0.13
%   Fz    approximately 200 to 1100 N
%   P = 70 kPa, IA = 0 deg

if nargin < 3 || isempty(gamma), gamma = 0; end
if nargin < 4 || isempty(P),     P = 70;    end
if nargin < 5,                   kappa = []; end

%% ========================================================================
% LATERAL MAGIC FORMULA COEFFICIENTS (existing project fit)
% ========================================================================
PCY1 =  1.5004;
PDY1 =  2.3101;
PDY2 = -0.4792;
PDY3 =  5.7663;
PEY1 = -0.0227;
PEY2 = -0.5330;
PEY3 = -0.0070;
PEY4 = -19.5150;
PEY5 = 289.5528;
PKY1 = -38.3895;
PKY2 =  0.8174;
PKY3 =  0.4086;
PKY4 =  1.0500;
PKY5 = 14.5710;
PKY6 =  4.3507;
PKY7 =  1.8557;
PHY1 = -8.7429e-04;
PHY2 = -5.3980e-04;
PVY1 = -0.0646;
PVY2 = -0.0124;
PVY3 = -0.5581;
PVY4 =  0.9574;
PPY1 =  0.5268;
PPY2 =  1.3435;
PPY3 = -0.2115;
PPY4 = -0.0992;
PPY5 = -0.6321;

Fz0_y = 1080;
P0    = 70;

lambda_Cy       = 1;
lambda_Ky       = 1;
lambda_Ey       = 1;
lambda_muy      = 1;
lambda_Ky_gamma = 1;
lambda_Hy       = 1;
lambda_Vy       = 1;

zeta0 = 1;
zeta2 = 1;
zeta3 = 1;
zeta4 = 1;

%% ========================================================================
% NUMERICAL GUARD
% ========================================================================
% Prevent zero/negative vertical load from making the MF coefficients
% singular during a nonlinear solver iterate.
Fz_safe = max(Fz, 1);

%% ========================================================================
% LATERAL PURE-SLIP MODEL
% ========================================================================
dfz_y = (Fz_safe - Fz0_y) ./ Fz0_y;
dpi   = (P - P0) ./ P0;

Ky_alpha = PKY1 .* Fz0_y .* (1 + PPY1 .* dpi) .* ...
    sin(PKY4 .* atan(Fz_safe ./ ...
    ((PKY2 + PKY5 .* gamma.^2) .* Fz0_y .* (1 + PPY2 .* dpi)))) .* ...
    (1 - PKY3 .* abs(gamma)) .* lambda_Ky .* zeta3;

Ky_gamma_0 = (PKY6 + PKY7 .* dfz_y) .* Fz_safe .* ...
    lambda_Ky_gamma .* (1 + PPY5 .* dpi);

SVy_0 = Fz_safe .* (PVY1 + PVY2 .* dfz_y) .* ...
    lambda_Vy .* lambda_muy;
SVy_gamma = Fz_safe .* (PVY3 + PVY4 .* dfz_y) .* gamma .* ...
    lambda_Ky .* lambda_muy .* zeta2;
SVy = SVy_0 .* zeta2 + SVy_gamma;

SHy_0 = (PHY1 + PHY2 .* dfz_y) .* lambda_Hy;
SHy_gamma = (Ky_gamma_0 .* gamma - SVy_gamma) ./ Ky_alpha .* zeta0 + zeta4 - 1;
SHy = SHy_0 + SHy_gamma;

alpha_y = alpha + SHy;

muy = (PDY1 + PDY2 .* dfz_y) .* (1 - PDY3 .* gamma.^2) .* ...
      (1 + PPY3 .* dpi + PPY4 .* dpi.^2) .* lambda_muy;

Cy = PCY1 .* lambda_Cy .* ones(size(Fz_safe));
Dy = muy .* Fz_safe .* zeta2;
By = Ky_alpha ./ (Cy .* Dy + eps);

Ey = (PEY1 + PEY2 .* dfz_y) .* ...
     (1 + PEY5 .* gamma.^2 - ...
     (PEY3 + PEY4 .* gamma) .* sign(alpha_y)) .* lambda_Ey;
Ey = min(Ey, 1);

ByAlpha = By .* alpha_y;
Fy0 = Dy .* sin(Cy .* atan(ByAlpha - Ey .* (ByAlpha - atan(ByAlpha)))) + SVy;
%Fy  = 0.65 .* Fy0; 
Fy = Fy0;

if nargout > 1
    C = struct( ...
        'Ky_alpha', Ky_alpha, ...
        'Ky_gamma_0', Ky_gamma_0, ...
        'SVy_0', SVy_0, ...
        'SVy_gamma', SVy_gamma, ...
        'SVy', SVy, ...
        'SHy_0', SHy_0, ...
        'SHy_gamma', SHy_gamma, ...
        'SHy', SHy, ...
        'alpha_y', alpha_y, ...
        'muy', muy, ...
        'Cy', Cy, ...
        'Dy', Dy, ...
        'By', By, ...
        'Ey', Ey, ...
        'dfz', dfz_y, ...
        'dpi', dpi, ...
        'Fz0', Fz0_y, ...
        'Fy0', Fy0, ...
        'scale', 1); %was 0.65
end

%% ========================================================================
% LONGITUDINAL PURE-SLIP MODEL FROM Fx_70_16_in.m
% ========================================================================
if isempty(kappa)
    Fx = [];
    Cx = struct();
    return
end

% Optimized coefficients copied from Fx_70_16_in.m
pCx1   = 1.549294;
pDx1   = 2.468610;
pDx2   = -0.558514;
pDx3   = 0.000000;
pEx1   = -0.083271;
pEx2   = -0.343884;
pEx3   = -0.323038;
pEx4   = -4.996932;
pKx1   = 36.690891;
pKx2   = -0.000042;
pKx3   = -0.274972;
pHx1   = 0.030364;
pHx2   = 0.037008;
pVx1   = -0.020686;

% IMPORTANT: this fit uses its own reference load.
Fz0_x = 1080;
dfz_x = (Fz_safe - Fz0_x) ./ Fz0_x;

% Shape factor
Cx_mf = pCx1 .* ones(size(Fz_safe));

% Peak factor.  The source fit explicitly used gamma = 0; pDx3 is zero,
% so there is no camber effect in this longitudinal fit.
mux = (pDx1 + pDx2 .* dfz_x) .* (1 - pDx3 .* 0);
Dx = mux .* Fz_safe;

% Longitudinal slip stiffness
Kx = Fz_safe .* (pKx1 + pKx2 .* dfz_x) .* exp(pKx3 .* dfz_x);

% Stiffness factor
Bx = Kx ./ (Cx_mf .* Dx + eps);

% Horizontal/vertical shifts
SHx = pHx1 + pHx2 .* dfz_x;
SVx = Fz_safe .* pVx1;
kappa_x = kappa + SHx;

% Curvature factor
Ex = (pEx1 + pEx2 .* dfz_x + pEx3 .* dfz_x.^2) .* ...
     (1 - pEx4 .* sign(kappa_x));
Ex = min(Ex, 0.999);

% Magic Formula
BxKappa = Bx .* kappa_x;
Fx0 = Dx .* sin(Cx_mf .* atan(BxKappa - Ex .* (BxKappa - atan(BxKappa)))) + SVx;
%Fx  = 0.65 .* Fx0;
Fx = Fx0;

if nargout > 3
    Cx = struct( ...
        'mux', mux, ...
        'Cx', Cx_mf, ...
        'Dx', Dx, ...
        'Kx', Kx, ...
        'Bx', Bx, ...
        'Ex', Ex, ...
        'SHx', SHx, ...
        'SVx', SVx, ...
        'kappa_x', kappa_x, ...
        'dfz', dfz_x, ...
        'Fz0', Fz0_x, ...
        'Fx0', Fx0, ...
        'scale', 1, ... %was 0.65
        'fit_pressure_kPa', 70, ...
        'fit_camber_rad', 0, ...
        'valid_kappa_min', -0.198, ...
        'valid_kappa_max', 0.13, ...
        'valid_Fz_min_N', 200, ...
        'valid_Fz_max_N', 2000);
end

end
