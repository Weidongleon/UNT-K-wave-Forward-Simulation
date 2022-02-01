clear
%%
% create the computational grid
scale = 1;

PML_size = 10;                  % size of the PML in grid points
Nx = 200 * scale - 2 * PML_size; % number of grid points in the x direction
Ny = 200 * scale - 2 * PML_size; % number of grid points in the y direction
Nz = 200 * scale - 2 * PML_size; % number of grid points in the z direction
dx = 1e-3 / scale;            % grid point spacing in the x direction [m]
dy = 1e-3 / scale;            % grid point spacing in the y direction [m]
dz = 1e-3 / scale;            % grid point spacing in the z direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);


% define the properties of the propagation medium
medium.sound_speed = 1500; %[m/s]

% create initial pressure distribution using makeDisc
disc_radius = 3; % [grid points]
disc_magnitude = 10; %[Pa]
p0 = disc_magnitude * makeBall(Nx, Ny, Nz, Nx/2, Ny/2, Nz/2, disc_radius);

% smooth the initial pressure distribution and restore the magnitude
source.p0 = smooth(p0, true);

%% sensor mask creation 
sensor.mask = zeros(kgrid.Nx,kgrid.Ny, kgrid.Nz);
sensor.mask (:,:,141:180) = 1;
%% visualisation
voxelPlot(double(p0 | sensor.mask));

%% run simulation 
% create the time array
kgrid.makeTime(medium.sound_speed);

% input arguments **inclue PML
input_args = {'PMLSize', PML_size, 'PMLInside', false, ...
    'PlotPML', true, 'Smooth', true, 'DataCast', 'single'};
% run simulation 
sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:});

%% preparing for reconstruction
% creation of reshape sensor array for FFT
sensor_data_rs_pre = zeros(180,30,180);
sensor_data_rs_pre(:,:,180) = 1;

sensor_data_rs = reshape(a, Nx, Ny, kgrid.Nt);

%% reconstruction simulation in FFT
% reconstruct the initial pressure
p_xyz = kspacePlaneRecon(sensor_data_rs, kgrid.dx, kgrid.dy, kgrid.dt, ...
    medium.sound_speed, 'DataOrder', 'yzt', 'PosCond', true, 'Plot', true);

%% pre-plots
% define a k-space grid using the dimensions of p_xyz
[Nx_recon, Ny_recon, Nz_recon] = size(p_xyz);
kgrid_recon = kWaveGrid(Nx_recon, kgrid.dt * medium.sound_speed, Ny_recon, kgrid.dy, Nz_recon, kgrid.dz);

% define a k-space grid with the same z-spacing as p0
[Nx_p0, Ny_p0, Nz_p0] = size(source.p0);
kgrid_interp = kWaveGrid(Nx_p0, kgrid.dx, Ny_p0, kgrid.dy, Nz_p0, kgrid.dz);

% resample the p_xyz to be the same size as p0; for a matrix indexed as 
% [M, N, P], the axis variables passed to interp3 must be given in the 
% order N, M, P
p_xyz_rs = interp3(kgrid_recon.y - min(kgrid_recon.y(:)), ...
                   kgrid_recon.x - min(kgrid_recon.x(:)), ...
                   kgrid_recon.z - min(kgrid_recon.z(:)), ...
                   p_xyz, ...
                   kgrid_interp.y - min(kgrid_interp.y(:)), ...
                   kgrid_interp.x - min(kgrid_interp.x(:)), ...
                   kgrid_interp.z - min(kgrid_interp.z(:)));
%% plots
% plot the initial pressure
voxelPlot(double(p0 | sensor.mask));
set(gca, 'Projection', 'perspective');
view([0, 99]);

% plot the initial pressure
figure;
plot_scale = [-10, 10];
subplot(2, 2, 1);
imagesc(kgrid_interp.y_vec * 1e3, kgrid_interp.x_vec * 1e3, squeeze(source.p0(:, :, Nz/2)), plot_scale);
title('x-y plane');
axis image;

subplot(2, 2, 2);
imagesc(kgrid_interp.z_vec * 1e3, kgrid_interp.x_vec * 1e3, squeeze(source.p0(:, Ny/2, :)), plot_scale);
title('x-z plane');
axis image;
xlabel('(All axes in mm)');

subplot(2, 2, 3);
imagesc(kgrid_interp.z_vec * 1e3, kgrid_interp.y_vec * 1e3, squeeze(source.p0(Nx/2, :, :)), plot_scale);
title('y-z plane');
axis image;
colormap(getColorMap);
sgtitle('Initial pressure source');

% plot the reconstructed initial pressure
figure;
subplot(2, 2, 1);
imagesc(kgrid_interp.y_vec * 1e3, kgrid_interp.x_vec * 1e3, squeeze(p_xyz_rs(:, :, Nz/2)), plot_scale);
title('x-y plane');
axis image;

subplot(2, 2, 2);
imagesc(kgrid_interp.z_vec * 1e3, kgrid_interp.x_vec * 1e3, squeeze(p_xyz_rs(:, Ny/2, :)), plot_scale);
title('x-z plane');
axis image;
xlabel('(All axes in mm)');

subplot(2, 2, 3);
imagesc(kgrid_interp.z_vec * 1e3, kgrid_interp.y_vec * 1e3, squeeze(p_xyz_rs(Nx/2, :, :)), plot_scale);
title('y-z plane');
axis image;
colormap(getColorMap);
sgtitle('Reconstructed Pressure Source');

%% view the reconstruction slice by slice
flyThrough(p_xyz_rs);









