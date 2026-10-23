% hit_stats.m
% Single-point and spectral statistics of a MHIT36 velocity field
% (staggered C-grid, u(i) on the left face of cell i, triply periodic box of side 2*pi).
%
% Reports: u_rms per component, k, dissipation (spectral and with the code's own
% 2nd-order finite differences), eta, lambda, Re_lambda, integral scale, time scales,
% kmax*eta and dx/eta, derivative skewness, discrete divergence.
% Writes the shell-averaged energy spectrum to spectrum.txt and plots it.
%
% Optional: physical scales and particle Stokes numbers for a rescaled box (L, nuP).
clear

%% Input
dir_in = 'mhit36_fields';
files  = {'u.dat', 'v.dat', 'w.dat'};
prec   = 'double';     % MHIT36 writes double precision, x fastest
nu     = 0.015;        % MHIT36 kinematic viscosity, mu/rho from input.inp

% optional rescaling to a physical box (L = [] to skip)
L      = 0.01;         % physical box size [m]
nuP    = 1e-5;         % physical kinematic viscosity [m^2/s]
rhoP   = 1000;         % particle density [kg/m^3]
rhoG   = 1;            % gas density [kg/m^3]
dp     = [1:10, 1.5:1:9.5]*1e-6;   % diameters [m]

%% Read
U = cell(1, 3);
for c = 1:3
    fid = fopen(fullfile(dir_in, files{c}), 'r');
    a = fread(fid, inf, prec);
    fclose(fid);
    n = round(numel(a)^(1/3));
    if n^3 ~= numel(a), error('%s: %d values is not a cube', files{c}, numel(a)); end
    U{c} = reshape(a, [n n n]);
end
clear a
dx = 2*pi/n;
fprintf('Grid %d^3, dx = %.5g, nu = %g\n', n, dx, nu);

%% Discrete divergence (MHIT36 stencil) and FD velocity gradients
D = zeros(n, n, n);  gfd = 0;  dudx = cell(1, 3);
for c = 1:3
    for d = 1:3
        G = (circshift(U{c}, -1, d) - U{c})/dx;   % staggered difference
        gfd = gfd + mean(G(:).^2);
        if c == d, dudx{c} = G; D = D + G; end
    end
end
g2 = mean(cellfun(@(x) mean(x(:).^2), dudx));
fprintf('Discrete divergence: rms(div)/rms(du_i/dx_i) = %.2e\n', sqrt(mean(D(:).^2)/g2));
clear D G

%% Spectral energy (staggering only changes phases: no effect on these statistics)
k1 = [0:n/2-1, -n/2:-1];
[KX, KY, KZ] = ndgrid(k1, k1, k1);
K2 = KX.^2 + KY.^2 + KZ.^2;
clear KX KY KZ
E3 = zeros(n, n, n);
for c = 1:3
    E3 = E3 + abs(fftn(U{c})/n^3).^2;
end

%% Single-point statistics
urms_c = cellfun(@(u) sqrt(mean((u(:) - mean(u(:))).^2)), U);
urms   = sqrt(mean(urms_c.^2));
kin    = 1.5*urms^2;
epsS   = nu*sum(K2(:).*E3(:));        % spectral: nu <du_i/dx_j du_i/dx_j>
epsFD  = nu*gfd;                       % code's 2nd-order FD gradients
eta    = (nu^3/epsS)^(1/4);
tauEta = sqrt(nu/epsS);
lam    = sqrt(15*nu*urms^2/epsS);
Relam  = urms*lam/nu;

lamc = zeros(1, 3);  Sk = zeros(1, 3);
for c = 1:3
    dd = dudx{c}(:);
    lamc(c) = sqrt(urms_c(c)^2/mean(dd.^2));   % longitudinal Taylor microscale
    Sk(c)   = mean(dd.^3)/mean(dd.^2)^1.5;     % longitudinal derivative skewness
end
clear dudx

%% Spectrum E(k), integral scale
kk = round(sqrt(K2));
nk = n/2;
Ek = accumarray(kk(:) + 1, E3(:)/2);
Ek = Ek(1:nk + 1);
ks = (0:nk)';
clear K2 kk E3
Lint = pi/(2*urms^2)*sum(Ek(2:end)./ks(2:end));
TL   = Lint/urms;
Te   = kin/epsS;

fid = fopen('spectrum.txt', 'w');
fprintf(fid, '# k  E(k)   (n = %d, nu = %g, eps = %g, eta = %g)\n', n, nu, epsS, eta);
fprintf(fid, '%d %.8e\n', [ks(2:end)'; Ek(2:end)']);
fclose(fid);

%% Report
fprintf('\nVelocity\n');
fprintf('  u_rms (x,y,z)        = %.4g %.4g %.4g   (max/min = %.3f)\n', urms_c, max(urms_c)/min(urms_c));
fprintf('  u'' = u_rms           = %.4g\n', urms);
fprintf('  k = 3/2 u''^2         = %.4g   (from spectrum: %.4g)\n', kin, sum(Ek(2:end)));
fprintf('Dissipation\n');
fprintf('  eps (spectral)       = %.4g\n', epsS);
fprintf('  eps (2nd-order FD)   = %.4g   (FD/spectral = %.3f)\n', epsFD, epsFD/epsS);
fprintf('Scales\n');
fprintf('  eta                  = %.4g   (2*pi/eta = %.1f)\n', eta, 2*pi/eta);
fprintf('  lambda               = %.4g   (longitudinal x,y,z: %.4g %.4g %.4g)\n', lam, lamc);
fprintf('  L_int                = %.4g   (L_int/eta = %.1f)\n', Lint, Lint/eta);
fprintf('  tau_eta              = %.4g\n', tauEta);
fprintf('  T_L = L_int/u''       = %.4g,   T_e = k/eps = %.4g\n', TL, Te);
fprintf('Reynolds numbers\n');
fprintf('  Re_lambda            = %.1f\n', Relam);
fprintf('  Re_L = u'' L_int/nu   = %.1f\n', urms*Lint/nu);
fprintf('Resolution\n');
fprintf('  dx/eta               = %.2f\n', dx/eta);
fprintf('  kmax*eta (kmax=n/2)  = %.2f\n', nk*eta);
fprintf('  derivative skewness  = %.3f %.3f %.3f   (resolved HIT: about -0.5)\n', Sk);

%% Optional: physical box
if ~isempty(L)
    ell = L/(2*pi);  V = (nuP/nu)/ell;
    etaP = eta*ell;  teP = etaP^2/nuP;
    fprintf('\nRescaled to L = %g m, nu = %g m^2/s (same Re_lambda)\n', L, nuP);
    fprintf('  u'' = %.4g m/s, eta = %.4g m, tau_eta = %.4g s, lambda = %.4g m, L_int = %.4g m, eps = %.4g m^2/s^3\n', ...
        urms*V, etaP, teP, lam*ell, Lint*ell, nuP^3/etaP^4);
    tauP = rhoP*dp.^2/(18*nuP*rhoG);
    fprintf('  d [um]   St_eta\n');
    fprintf('  %6.2f   %.4f\n', [dp*1e6; tauP/teP]);
    fprintf('  box size giving St_eta = 1 at d = 10 um: L = %.4g m\n', ...
        L*sqrt(rhoP*(10e-6)^2/(18*nuP*rhoG)/teP));
end

%% Plot
figure('Name', 'Energy spectrum');
loglog(ks(2:end)*eta, Ek(2:end)/(epsS*nu^5)^(1/4), 'k-', 'LineWidth', 1.5); hold on
kr = logspace(log10(2*eta), log10(0.3), 20);
loglog(kr, 1.6*kr.^(-5/3), 'r--');
xlabel('k\eta'); ylabel('E(k)/(\epsilon\nu^5)^{1/4}');
legend('MHIT36', '1.6 (k\eta)^{-5/3}', 'Location', 'southwest');
title(sprintf('Re_\\lambda = %.1f, k_{max}\\eta = %.2f', Relam, nk*eta));
grid on
