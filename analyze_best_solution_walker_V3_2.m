% analyze_best_solution_walker_V3_2 is an example script for analyzing the 
% revisit time of an optimized Walker-style satellite constellation using 
% previously computed solutions or custom parameters.
%
% Specific to this implementation:
% - The ground points are assumed to be at a flying altitude of 6 km above 
%   the Earth's surface, representing typical aircraft altitudes.
% - Each aircraft is assumed to transmit a signal within a conical region 
%   defined by a maximum elevation angle on the aircraft's local horizon. 
%   This affects the calculation of visibility from the satellite but is 
%   negligible for most practical FOVs.
% - The script computes worst-case, average, and best revisit times for 
%   all grid points globally, with specific attention to latitudes 
%   relevant for air traffic (|lat| ≤ 75°).
%
%--- Copyright notice ---%
% Copyright (C) 2026 Alessandro Zamboni
%
% Written for educational and research purposes
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the 
% Free Software Foundation, either version 3 of the License, or (at your 
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but 
% WITHOUT ANY WARRANTY; without even the implied warranty of 
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General 
% Public License for more details.
%
% You should have received a copy of the GNU General Public License along 
% with this program. If not, see <http://www.gnu.org/licenses/>.
%

clear; clc;
tic



%% ============================================================
% LOADS THE LAST OPTIMIZATION RESULT OR THE CHOSEN ONE
%% ============================================================

fprintf('[ANALYSIS] Loading Walker best solution...\n');

%UNCOMMENT TO LOAD THE CHOSEN ONE
%load('Optimizes_orbit_V3_2_#1.mat');

%UNCOMMENT TO LOAD THE LAST OPTIMIZATION RESULT
files = dir('Optimized_constellation_V3_2_#*.mat');

if isempty(files)
    error('No matching files found.')
end

nums = zeros(length(files),1);

for i = 1:length(files)
    token = regexp(files(i).name,'#(\d+)\.mat','tokens');
    nums(i) = str2double(token{1});
end

[~,idx] = max(nums);
latest_file = files(idx).name;

load(latest_file)
fprintf('[ANALYSIS] Loaded solution: %s \n',latest_file);

%% ============================================================
% CONSTANTS
%% ============================================================

mu = 398600.4418;          % km^3/s^2
Re = 6378.137;             % km
omegaE = 7.2921159e-5;     % rad/s
fa=6;                      %flying altitude

SIM_DT = 50;               % seconds

%% ============================================================
% DECODE WALKER SOLUTION
% ============================================================

% %PARAMETERS FROM THE LOADED FILE:
% Np = round(xbest(1));         %number of orbital planes
% Ns = round(xbest(2));         %satellite per plane
% F  = round(xbest(3));         %phasing
% a  = xbest(4);                %semiaxis
% inc  = xbest(5);              %inclination


% CHOOSE THE PARAMETERS YOURSELF:
Np = 5;                %number of orbital planes
Ns = 4;                %satellite per plane
F  = 1;                %phasing
a  = Re + 863;         %semiaxis
inc  = deg2rad(80.53); %inclination

Nt = Np * Ns;

fprintf('[INFO] Planes: %d\n',Np);
fprintf('[INFO] Satellites per plane: %d\n',Ns);
fprintf('[INFO] Total satellites: %d\n',Nt);

% ============================================================
% PRINT CONCISE ORBIT CHARACTERISTICS
% ============================================================
n = sqrt(mu / a^3);                % mean motion [rad/s]
T_orbit = 2*pi / n;                % orbital period [s]
altitude = a - Re;                 % circular orbit altitude [km]
incl_deg = rad2deg(inc);             % inclination [deg]
raan_spacing_deg = 360 / Np;       % RAAN spacing between planes [deg]
phasing_deg = 360 * F / Nt;        % relative phasing per satellite [deg]

fprintf('\n[ORBIT SUMMARY]\n');
fprintf(' Semi-major axis (a): %.3f km\n', a);
fprintf(' Altitude above mean Earth radius: %.3f km\n', altitude);
fprintf(' Inclination: %.3f deg\n', incl_deg);
fprintf(' Orbital period: %.3f min\n', T_orbit/60);
fprintf(' Mean motion: %.6f rev/day\n', n*86400/(2*pi));
fprintf(' Number of planes: %d\n', Np);
fprintf(' Satellites/plane: %d\n', Ns);
fprintf(' Total satellites: %d\n', Nt);
fprintf(' Walker phasing (F): %d -> phase offset per sat: %.3f deg\n', F, phasing_deg);
fprintf(' RAAN spacing between planes: %.3f deg\n\n', raan_spacing_deg);

%% ============================================================
% SIMULATION TIME (approx groundtrack repeat)
%% ============================================================

T_orbit = 2*pi*sqrt(a^3/mu);
repeat_factor = lcm(Np,Ns);

SIM_TIME = repeat_factor * T_orbit*10;

times = 0:SIM_DT:SIM_TIME;
Nt_steps = length(times);

fprintf('[INFO] Sim time: %.2f hours\n',SIM_TIME/3600);

%% ============================================================
% SENSOR PARAMETERS
%% ============================================================

max_elv=70; %max elevation angle on the local horizon: the plane only send the signal the satellite below this elevation angle

x=asind(((Re+fa)/(altitude+Re+fa))*sind(70+90))+70;
fov_internal=x-70;
phi=asind((Re+fa)/(Re+fa+altitude)); %halfbeam width to the satellite horizon
phi_central=90-phi; %central angle of external fov
phi_null=90-x; %central angle of null signal

FOV_HALF_ANGLE = deg2rad(phi);      % sensor half-angle
ELEV_MASK      = deg2rad(0);       % minimum elevation

% Sperical caps of FOV and NULL_FOV
A_fov=2*pi*(Re+fa)^2*(1-cosd(phi_central)); %area on ground seen by the antenna
A_null=2*pi*(Re+fa)^2*(1-cosd(phi_null));   %area on ground right below the satellite where ...
% the VHF signal cannot reach the satellite---> referred to as blindspot

sat_angle_ratio=phi/(x-70);

null_ratio=(A_null/A_fov)*100;

%With 978.25km of altitude The null ratio becomes phi/fov_internal=60.127/17.25=3.48 ...
% and the null ratio is 0.8654. However, the central angles phi_central=28.97 and phi_null=2.7480 ...
% both increase with altitude, and the magnitude of the increse is respectively about 10 deg ...
% and 1.5 deg. While with 400km of altitude phi/fov_internal=70.22/18.7751=3.74 and the null_ratio ...
% is 0.3876. because the fov_internal slightly decreases with altitude while phi decreases ...
% significantly with altitude. Considering that the area contained in the spherical caps defined ...
% by these central angles increases in a more-than-linear way for small central angles and then in ...
% an almost-linear way for central angles~>=40deg, with the increase of the central angles then the ...
% null_ratio is bound to increase in our case since phi_central increases linearly and phi_null ...
% increases more than linearly.
%IN CONCLUSION: it makes sense that the null_ratio increases with altitude, but it is still very ...
%low so in this preliminary analisys we can ignore the we can ignore this blindspot.

cos_fov  = cos(FOV_HALF_ANGLE);
sin_mask = sin(ELEV_MASK);

%% ============================================================
% EQUAL AREA SPHERE GRID (Carnegie Mellon method)
% ============================================================

%SOURCE: How to generate equidistributed points on the surface of a sphere
%IMPORTANT: The algorithm here described distributes a lot of point along ...
% each longitude and few along each latitudide. So it was decided to adapt it...
% including a latitudinalization_factor whic decreases the latitude step.

% If this section is uncommented also the end that closes the "else" must be uncommented
% grid_file = 'EarthGrid.h5';
%
% if isfile(grid_file)
%     fprintf('[INFO] Found existing grid file: %s. Loading...\n', grid_file);
%     lat_vec = h5read(grid_file, '/lat');
%     lon_vec = h5read(grid_file, '/lon');
%     Npoints = length(lat_vec);
%     fprintf('[INFO] Loaded %d grid points from file.\n', Npoints);
% else
%    fprintf('[INFO] No grid file found. Building ground grid...\n');

A_earth=4*pi*(Re+fa)^2;
N_mingridpoints=ceil(A_earth/A_fov)*500;

fprintf('[INFO] Generating equispaced Earth point grid ...\n');

N = N_mingridpoints;

r = 1; %unitary radius of the sphere
Ncount = 0; %initial point counter

area_cell_tentativo = 4*pi*r^2 / N; %first tentative area of a cell= Spehere area divided by number of points
d = sqrt(area_cell_tentativo);
lat_factor=1;
Mtheta = round(lat_factor*pi / d);
dtheta = pi / Mtheta;
dphi = area_cell_tentativo / dtheta;

lat_vec = zeros( max(N,1), 1 );
lon_vec = zeros( max(N,1), 1 );

pts = zeros( max(N,1), 3 ); % preallocate approximate size

for m = 0:(Mtheta-1)
    theta = pi*(m + 0.5) / Mtheta;           % ϑ
    lat = pi/2 - theta;                      %because longitude is defined in [-pi, pi] with the zero at the equator
    Mphi = round( 2*pi * sin(theta) / dphi );% number of longitudes at this latitude

    for n = 0:(Mphi-1)
        phi = 2*pi * n / Mphi;               % ϕ

        lon_vec(Ncount+1,1) = phi - pi;         %because longitude is defined in [-pi, pi] with the zero at the Greenwich Meridian
        lat_vec(Ncount+1,1) = lat;

        Ncount = Ncount + 1;

    end
end
lon_vec = lon_vec(1:Ncount,1); % trim unused preallocated rows
lat_vec = lat_vec(1:Ncount,1);

Npoints = Ncount;
fprintf('[INFO] Generated %d grid points.\n', Npoints);

%Cell Area Heatmap:
flag_plots=0; %if zero it does not plot anything, it just calculates
[avg_cell_area, max_cell_area, min_cell_area]=cell_area_heatmap(lat_vec,lon_vec,Re, fa,flag_plots);
%fprintf('[INFO] Grid cell areas [km^2]: avg = %.3e, min = %.3e, max = %.3e\n', avg_cell_area, min_cell_area, max_cell_area);
fov_to_cell_area_ratio = A_fov/avg_cell_area;

% %calculates the intervals in latitude "dlat" and in longitude "dlon"
% dlat=zeros( max(Ncount,1), 1 );
% temp_count=0;
% for i=1:Ncount-1
%     temp=lat_vec(i)-lat_vec(i+1);
%     if temp == 0
%         temp_count=temp_count+1;
%         continue
%     end
%     dlat(i-temp_count)=temp;
% end
% dlat=dlat(1:Ncount-temp_count-1,1);
% dlat=rad2deg(dlat);
% fprintf('[INFO] LatGrid resolution: %.3f deg\n',dlat);
%
% dlon=zeros( max(Ncount,1), 1 );
% temp_count=0;
% for i=1:Ncount-1
%     if lat_vec(i+1)-lat_vec(i) ~= 0
%         continue
%     end
%     temp_count=temp_count+1;
%     dlon(temp_count) = lon_vec(i+1)-lon_vec(i);
%
% end
% dlon=dlon(1:temp_count,1);
% dlon=rad2deg(dlon);
% %fprintf('[INFO] LonGrid resolution: %.3f deg\n',dlon);

% % Plot sphere surface and points
% figure;
% % sphere surface
% [Xs, Ys, Zs] = sphere(60);
% surf(r*Xs, r*Ys, r*Zs, 'FaceAlpha', 1, 'EdgeColor', 'none', 'FaceColor', [0.8 0.8 1]);
% hold on;
% % points
% x = r * cos(lat_vec) .* cos(lon_vec);
% y = r * cos(lat_vec) .* sin(lon_vec);
% z = r * sin(lat_vec);
% scatter3(x, y, z, 20, 'r', 'filled');
% axis equal;
% xlabel('X'); ylabel('Y'); zlabel('Z');
% title(sprintf('Nearly Equal-Area Points on Sphere (N requested = %d, generated = %d)', N, Ncount));
% view(30,20);
% light; lighting gouraud;


% % Save to HDF5
% fprintf('[INFO] Saving grid points to %s...\n', grid_file);
% if isfile(grid_file), delete(grid_file); end
% h5create(grid_file,'/lat',size(lat_vec),'Datatype','double');
% h5create(grid_file,'/lon',size(lon_vec),'Datatype','double');
% h5write(grid_file,'/lat',lat_vec);
% h5write(grid_file,'/lon',lon_vec);
% fprintf('[INFO] Grid saved.\n');

%end %uncomment if the "grid_file = 'EarthGrid.h5';" function is uncommented


%% Ground positions (ECEF)


rg = (Re+fa) * [cos(lat_vec).*cos(lon_vec), ...
    cos(lat_vec).*sin(lon_vec), ...
    sin(lat_vec)];

rg = single(rg);

zenith = rg ./ vecnorm(rg,2,2);

%% ============================================================
% VISIBILITY MATRIX
%% ============================================================

vis_matrix = false(Npoints,Nt_steps);

%% ============================================================
% PRECOMPUTE ORBIT QUANTITIES
%% ============================================================

n = sqrt(mu/a^3);
M_all = n*times;

theta = omegaE * times;
c = cos(theta);
s = sin(theta);

%% ============================================================
% CONSTELLATION PROPAGATION
%% ============================================================

%iteration global counter:
total_iter = Np * Ns * Nt_steps;
iter = 0;
n_updates = 20; % how many prints you want
update_steps = round(linspace(1,total_iter,n_updates));


fprintf('[INFO] Propagating constellation...\n');

for p = 0:Np-1

    RAAN = 2*pi*p/Np;

    R3 = [cos(RAAN) -sin(RAAN) 0;
        sin(RAAN)  cos(RAAN) 0;
        0 0 1];

    R1 = [1 0 0;
        0 cos(inc) -sin(inc);
        0 sin(inc)  cos(inc)];

    Q = R3*R1;

    for sidx = 0:Ns-1

        M0 = 2*pi*sidx/Ns + 2*pi*F*p/Nt;

        for k = 1:Nt_steps

            M = M0 + M_all(k);

            r_orb = [a*cos(M);
                a*sin(M);
                0];

            r_eci = Q * r_orb;

            %% ECI -> ECEF

            R = [ c(k)  s(k) 0;
                -s(k)  c(k) 0;
                0     0    1 ];

            rs = R * r_eci;
            rs_vec = rs';

            %% ==================================================
            % ELEVATION CHECK
            %% ==================================================

            ground_to_sat = rs_vec - rg;
            %size(ground_to_sat)
            dist = vecnorm(ground_to_sat,2,2);
            ground_to_sat_unit = ground_to_sat ./ dist;

            el_sin = sum(ground_to_sat_unit .* zenith,2);
            visible_elev = el_sin > sin_mask;

            %% ==================================================
            % SENSOR FOV CHECK (nadir pointing) - FIXED
            %% ==================================================

            sat_to_ground = rg - rs_vec;
            sat_to_ground_norm = vecnorm(sat_to_ground,2,2);

            % Correct dot product with satellite radius
            cos_view = sum(sat_to_ground .* (-rs_vec),2) ./ (sat_to_ground_norm .* a);

            visible_fov = cos_view >= cos_fov;

            n_vis_elev = sum(visible_elev);
            n_vis_fov  = sum(visible_fov);
            %fprintf('Example step: visible by elevation = %d, by fov = %d\n', n_vis_elev, n_vis_fov);

            %% ==================================================
            % FINAL VISIBILITY
            %% ==================================================

            vis = visible_elev & visible_fov;
            vis_matrix(:,k) = vis_matrix(:,k) | vis;

            % --- Progress bar update ---
            iter = iter + 1;

            if ismember(iter,update_steps)
                fprintf('[PROGRESS] %.0f%% complete\n',100*iter/total_iter);
            end

        end
    end
end

fprintf('[INFO] Coverage propagation finished\n');

%% ============================================================
% REVISIT CALCULATION
%% ============================================================

fprintf('[INFO] Computing revisit metrics...\n');

worst_revisit = NaN(Npoints,1);
avg_revisit   = NaN(Npoints,1);
best_revisit  = NaN(Npoints,1);

valid_mask = false(Npoints,1);

for g = 1:Npoints

    visibility = vis_matrix(g,:);

    if any(visibility)

        valid_mask(g) = true;

        zero_runs = diff([0 visibility==0 0]);

        start_idx = find(zero_runs==1);
        end_idx   = find(zero_runs==-1)-1;

        run_lengths = end_idx - start_idx + 1;

        run_times = run_lengths * SIM_DT;

        if isempty(run_times) %when a point is alwauys visible we observe visibility=[1,1,1,...,1,1] ...
            % so run_times is empty and causes an error. This only happens in rare cases when we have ...
            % polar orbits + large FOV or the simulation time is too short

            worst_revisit(g) = 0;
            avg_revisit(g)   = 0;
            best_revisit(g)  = 0;

        else

            worst_revisit(g) = max(run_times)/60;
            avg_revisit(g)   = mean(run_times)/60;
            best_revisit(g)  = min(run_times)/60;

        end

    end
end

%% ============================================================
% COVERAGE METRIC
%% ============================================================

threshold = 15; %[min]

good_mask = (worst_revisit<=threshold) & valid_mask;

covered_area = sum(good_mask);
total_area   = sum(valid_mask);

surface_percent = covered_area/total_area;

avg_revisit_min   = mean(avg_revisit(valid_mask));
worst_revisit_min = max(worst_revisit(valid_mask));
revisit_variance  = mean(worst_revisit(valid_mask)) - mean(avg_revisit(valid_mask));

% 80% interested area max 15 min revisit constraint
LAT_THRESHOLD=deg2rad(75);
interesting_area_mask = (abs(lat_vec) <=  LAT_THRESHOLD);

ordered_AVG_interesting_rev_time   = sort(avg_revisit(interesting_area_mask), 'ascend');
avg_revisit80 = ordered_AVG_interesting_rev_time(1:round(0.8*length(ordered_AVG_interesting_rev_time))); %we consider the best 80%
max_avg_revisit80 = avg_revisit80(end);
fprintf("The max AVG revisit time interval on the best 80 %% of the sampled points in the area interested by air traffic is %3.5f \n", max_avg_revisit80);

ordered_WORST_interesting_rev_time = sort(worst_revisit(interesting_area_mask), 'ascend'); %we consider the worst revisit time of each location of interest, and we order them from smallest to biggest
worst_revisit80 = ordered_WORST_interesting_rev_time(1:round(0.8*length(ordered_WORST_interesting_rev_time))); %we consider the best 80%
max_worst_revisit80 = worst_revisit80(end);
fprintf("The worts revisit time interval on the best 80 %% of the sampled points in the area interested by air traffic is %3.5f \n", max_worst_revisit80);
toc

%% ============================================================
% PLOTS
%% ============================================================

lon_deg = rad2deg(lon_vec);
lat_deg = rad2deg(lat_vec);

figure
scatter(lon_deg,lat_deg,30,worst_revisit,'filled')
colorbar
title('Worst Revisit Time [min]')
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
axis equal
xlim([-180 180])
ylim([-90 90])

figure
scatter(lon_deg,lat_deg,30,avg_revisit,'filled')
colorbar
title('Average Revisit Time [min]')
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
axis equal
xlim([-180 180])
ylim([-90 90])

figure
scatter(lon_deg,lat_deg,30,best_revisit,'filled')
colorbar
title('Best Revisit Time [min]')
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
axis equal
xlim([-180 180])
ylim([-90 90])

%% ============================================================
% PLOTS (INTERESTING AREA ONLY)
%% ============================================================

figure
scatter(lon_deg(interesting_area_mask), ...
        lat_deg(interesting_area_mask), ...
        30, ...
        worst_revisit(interesting_area_mask), ...
        'filled')
colorbar
title('Worst Revisit Time [min] (|Lat| ≤ 75°)')
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
axis equal
xlim([-180 180])
ylim([-90 90])

figure
scatter(lon_deg(interesting_area_mask), ...
        lat_deg(interesting_area_mask), ...
        30, ...
        avg_revisit(interesting_area_mask), ...
        'filled')
colorbar
title('Average Revisit Time [min] (|Lat| ≤ 75°)')
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
axis equal
xlim([-180 180])
ylim([-90 90])

figure
scatter(lon_deg(interesting_area_mask), ...
        lat_deg(interesting_area_mask), ...
        30, ...
        best_revisit(interesting_area_mask), ...
        'filled')
colorbar
title('Best Revisit Time [min] (|Lat| ≤ 75°)')
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')
axis equal
xlim([-180 180])
ylim([-90 90])

fprintf('[ANALYSIS COMPLETE]\n');

%% ============================================================
% SAVE MAP DATA
% ============================================================

fprintf('[INFO] Saving heatmap data...\n');
any vi
save('BEST_SOLUTION_YET_map_data.mat', ...
    'worst_revisit', ...
    'avg_revisit', ...
    'best_revisit', ...
    'lat_vec', ...
    'lon_vec', ...
    'Np','Ns','F','a','inc','SIM_DT');

fprintf('[INFO] Heatmap data saved\n');

fprintf('Npoints = %d\n', Npoints);
fprintf('Any visibility at all? %d\n', any(vis_matrix(:)));
fprintf('Number of grid points ever visible: %d\n', sum(any(vis_matrix,2)));
