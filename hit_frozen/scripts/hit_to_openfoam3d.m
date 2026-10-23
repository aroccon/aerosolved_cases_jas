% hit_to_openfoam3d.m
% Map a MHIT36 velocity field (staggered Arakawa C-grid, triply periodic 2*pi box)
% onto the cell centres of a uniform OpenFOAM mesh, rescaled to physical units.
%
%  - Fourier interpolation: exact for the resolved modes, sharp spectral cut-off
%    when the OpenFOAM mesh is coarser than the MHIT36 grid.
%  - Staggering handled exactly with phase shifts (face positions -> cell centres).
%  - Projection with the modified wavenumber sin(k*dx)/dx of OpenFOAM's linear
%    face interpolation, so that div(phi) is zero to round-off on the target mesh
%    (uniform mesh, cyclic patches, uniform density). No extra OpenFOAM step needed.
%  - Rescaling preserves all Reynolds numbers: lengths x L/(2*pi), velocities x V,
%    with V = (nuP/nuD)/(L/(2*pi)).
%
% Output: <outdir>/U (volVectorField, cells ordered x fastest, then y, then z,
% as written by blockMesh for a single hex block). Keep 0/p uniform.
clear

%% MHIT36 input
nd     = 256;                         % MHIT36 grid points per direction
nuD    = 0.015;                       % MHIT36 kinematic viscosity (code units)
dir_in = 'mhit36_fields';              % folder with the MHIT36 fields
files  = {'u.dat', 'v.dat', 'w.dat'}; % binary velocity files
prec   = 'double';                    % precision of the binary files
endian = 'ieee-le';
face   = 'auto';  % staggering: 0 -> u(i) on the face at x = (i-1)*dx (left face of cell i)
                  %             1 -> u(i) on the face at x = i*dx     (right face of cell i)
                  %        'auto' -> chosen as the one with zero discrete divergence

%% OpenFOAM output
no     = 150;      % OpenFOAM cells per direction (blockMesh: (no no no))
L      = 0.01;     % physical box size [m] (blockMesh: cube of side L)
nuP    = 1e-5;     % physical kinematic viscosity [m^2/s]
rhoG   = 1;        % gas density [kg/m^3] (only for the Stokes numbers below)
rhoP   = 1000;     % particle density [kg/m^3]
dp     = [1:10, 1.5:1:9.5]*1e-6;  % diameters for the Stokes-number table [m]
outdir = '0';

if isempty(nuD), error('Set nuD, the MHIT36 viscosity in code units.'); end

%% Read MHIT36 velocity (array index 1 = x, 2 = y, 3 = z)
U = cell(1, 3);
for c = 1:3
    fname = fullfile(dir_in, files{c});
    fid = fopen(fname, 'r', endian);
    if fid < 0, error('Cannot open %s (current folder: %s)', fname, pwd); end
    U{c} = reshape(double(fread(fid, nd^3, prec)), [nd nd nd]);
    fclose(fid);
end
dD = 2*pi/nd;

%% Staggering: face offsets (in units of dx) of each component, cell centres at +1/2
if strcmp(face, 'auto')
    r = zeros(1, 2);
    for f = 0:1
        r(f+1) = rms_div(U, f, dD);
    end
    [~, i] = min(r);
    face = i - 1;
    fprintf('Staggering: relative discrete divergence %.2e (face=0), %.2e (face=1) -> face = %d\n', ...
        r(1), r(2), face);
end
s = 0.5*ones(3);              % s(c,d): offset of component c along direction d
for c = 1:3, s(c, c) = face; end

%% Fourier coefficients at the true sample positions (x = (n + s)*dx, n = 0..nd-1)
k1 = [0:nd/2-1, -nd/2:-1];
[KX, KY, KZ] = ndgrid(k1, k1, k1);
C = cell(1, 3);
for c = 1:3
    C{c} = fftn(U{c})/nd^3 .* exp(-1i*dD*(KX*s(c,1) + KY*s(c,2) + KZ*s(c,3)));
    C{c}(abs(KX) == nd/2 | abs(KY) == nd/2 | abs(KZ) == nd/2) = 0;   % drop Nyquist
end
U = []; %#ok<NASGU>

st = flow_stats(C, KX, KY, KZ, nuD);
fprintf('MHIT36: u_rms = %.4g, eps = %.4g, eta = %.4g, lambda = %.4g, Re_lambda = %.1f, kmax*eta = %.2f\n', ...
    st.urms, st.eps, st.eta, st.lambda, st.Relam, (nd/2)*st.eta);
clear KX KY KZ

%% Keep the modes resolved by both grids, shift to the OpenFOAM cell centres
kmax = floor((min(nd, no) - 1)/2);
ks   = -kmax:kmax;
src  = mod(ks, nd) + 1;
dst  = mod(ks, no) + 1;
dO   = 2*pi/no;
[KX, KY, KZ] = ndgrid(ks, ks, ks);
ph   = exp(1i*dO/2*(KX + KY + KZ));          % sample at (m + 1/2)*dO
for c = 1:3
    C{c} = C{c}(src, src, src) .* ph;
end
clear ph

% projection with OpenFOAM's modified wavenumber (central difference over 2*dO)
TX = sin(KX*dO)/dO;  TY = sin(KY*dO)/dO;  TZ = sin(KZ*dO)/dO;
T2 = TX.^2 + TY.^2 + TZ.^2;  T2(T2 == 0) = 1;
P  = (TX.*C{1} + TY.*C{2} + TZ.*C{3})./T2;
e0 = sum(abs(C{1}(:)).^2 + abs(C{2}(:)).^2 + abs(C{3}(:)).^2);
C{1} = C{1} - TX.*P;  C{2} = C{2} - TY.*P;  C{3} = C{3} - TZ.*P;
e1 = sum(abs(C{1}(:)).^2 + abs(C{2}(:)).^2 + abs(C{3}(:)).^2);
fprintf('Projection removed %.3f%% of the kinetic energy\n', 100*(1 - e1/e0));
clear TX TY TZ T2 P KX KY KZ

%% Back to physical space on the OpenFOAM mesh, rescaled
ell = L/(2*pi);
V   = (nuP/nuD)/ell;
Uo  = cell(1, 3);
for c = 1:3
    G = complex(zeros(no, no, no));
    G(dst, dst, dst) = C{c};
    Uo{c} = V*real(ifftn(G))*no^3;
end
clear G C

r = rms_div_collocated(Uo, L/no);
fprintf('OpenFOAM mesh %d^3: discrete div(phi) / |grad U| = %.2e\n', no, r);

%% Physical scales
etaP = st.eta*ell;  tauEta = etaP^2/nuP;  uP = st.urms*V;
fprintf(['Physical: L = %.4g m, u_rms = %.4g m/s, eta = %.4g m, tau_eta = %.4g s, ' ...
         'lambda = %.4g m, dx/eta = %.2f\n'], L, uP, etaP, tauEta, st.lambda*ell, (L/no)/etaP);
tauP = rhoP*dp.^2/(18*nuP*rhoG);
fprintf('  d [um]   St_eta\n');
fprintf('  %6.2f   %.4f\n', [dp*1e6; tauP/tauEta]);

%% Write U
if ~exist(outdir, 'dir'), mkdir(outdir); end
fid = fopen(fullfile(outdir, 'U'), 'w');
fprintf(fid, ['FoamFile\n{\n    version     2.0;\n    format      ascii;\n' ...
              '    class       volVectorField;\n    object      U;\n}\n\n' ...
              'dimensions      [0 1 -1 0 0 0 0];\n\n' ...
              'internalField   nonuniform List<vector>\n%d\n(\n'], no^3);
fprintf(fid, '(%.8e %.8e %.8e)\n', [Uo{1}(:)'; Uo{2}(:)'; Uo{3}(:)']);
fprintf(fid, ');\n\nboundaryField\n{\n    ".*"\n    {\n        type cyclic;\n    }\n}\n');
fclose(fid);
fprintf('Written %s\n', fullfile(outdir, 'U'));


%% ---------------------------------------------------------------------------
function r = rms_div(U, face, d)
% relative rms of the staggered-grid divergence for a given face convention
D = zeros(size(U{1}));  g = 0;
for c = 1:3
    if face == 0
        dU = (circshift(U{c}, -1, c) - U{c})/d;   % faces at the cell's left
    else
        dU = (U{c} - circshift(U{c}, 1, c))/d;    % faces at the cell's right
    end
    D = D + dU;  g = g + mean(dU(:).^2);
end
r = sqrt(mean(D(:).^2)/g);
end

function r = rms_div_collocated(U, d)
% div(phi) with linear face interpolation on a uniform mesh = central difference over 2d
D = zeros(size(U{1}));  g = 0;
for c = 1:3
    dU = (circshift(U{c}, -1, c) - circshift(U{c}, 1, c))/(2*d);
    D = D + dU;  g = g + mean(dU(:).^2);
end
r = sqrt(mean(D(:).^2)/g);
end

function st = flow_stats(C, KX, KY, KZ, nu)
% single-point statistics from the Fourier coefficients (incompressible HIT)
K2 = KX.^2 + KY.^2 + KZ.^2;
E  = abs(C{1}).^2 + abs(C{2}).^2 + abs(C{3}).^2;
st.urms   = sqrt(sum(E(:))/3);
st.eps    = nu*sum(K2(:).*E(:));
st.eta    = (nu^3/st.eps)^(1/4);
st.lambda = sqrt(15*nu*st.urms^2/st.eps);
st.Relam  = st.urms*st.lambda/nu;
end
