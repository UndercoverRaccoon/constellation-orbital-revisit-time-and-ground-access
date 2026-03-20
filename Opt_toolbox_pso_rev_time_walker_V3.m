% Opt_toolbox_pso_rev_time_walker_V3 is an example script for optimizing Walker-style 
% satellite constellations for global revisit time using a Particle Swarm 
% Optimization (PSO) approach. The script evaluates constellation 
% configurations over an equal-area Earth grid and computes worst-case and 
% average revisit times for each ground point.
%
% Specific to this implementation:
% - The grid points are assumed to be at a height of 6 km above the Earth's 
%   surface, representing typical aircraft flight altitudes.
% - Each aircraft transmits a signal within a local conical coverage defined 
%   by a maximum elevation angle relative to the aircraft's local horizon. 
%   This feature is included for completeness, but it has negligible effect 
%   on the results for typical satellite FOVs.
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

fprintf('\n[INIT] Starting Walker Circular PSO (equal-area grid version)\n');

%% ============================================================
% CONSTANTS
%% ============================================================

mu = 398600.4418;
Re = 6378.137;
omegaE = 7.2921159e-5;

MIN_ALT = 900;
MAX_ALT = 1100;
MAX_SATS = 20;

SIM_DT = 120;
h_air = 6;

lowest_inclination = 50;
highest_inclination = 70;

ELEV_MASK = deg2rad(0);


%% ============================================================
% OPTIMIZATION VARIABLES
%% ============================================================

nvars = 5;

lb = [2 4 1 Re+MIN_ALT deg2rad(lowest_inclination)];
ub = [4 6 3 Re+MAX_ALT deg2rad(highest_inclination)];

%% ============================================================
% FITNESS HANDLE
%% ============================================================

fitness = @(x) walkerFitness_vec( ...
    x,mu,Re,omegaE, ...
    MIN_ALT,MAX_ALT,MAX_SATS, ...
    SIM_DT,ELEV_MASK, ...
    h_air,highest_inclination);

%% ============================================================
% PSO OPTIONS
%% ============================================================

opts = optimoptions('particleswarm', ...
    'SwarmSize',60,...
    'MaxIterations',30,...
    'Display','iter',...
    'UseParallel',true);

if isempty(gcp('nocreate'))
    parpool;
end

[xbest,fbest] = particleswarm(fitness,nvars,lb,ub,opts);



fprintf('\n[FINAL] Best fitness: %.4f\n',fbest);
disp(xbest);

%% ============================================================
% SAVING THE SOLUTION
%% ============================================================

baseName = 'Optimized_constellation_V3_2_#';
k = 1;

while true
    filename = sprintf('%s%d.mat', baseName, k);
    
    if ~exist(filename, 'file')
        break
    end
    
    k = k + 1;
end

save(filename,'xbest','fbest');

%%-------------------------------------------------------------------------
function fitness = walkerFitness_vec(x,mu,Re,omegaE, ...
    MIN_ALT,MAX_ALT,MAX_SATS, ...
    SIM_DT,ELEV_MASK, ...
    h_air,highest_inclination)

Np = round(x(1));
Ns = round(x(2));
F  = round(x(3));
a  = x(4);
inc = x(5);

Nt = Np*Ns;
alt = a-Re;

if Nt > MAX_SATS
    fitness = 1e6;
    return
end

if alt < MIN_ALT || alt > MAX_ALT
    fitness = 1e6;
    return
end

%% ============================================================
% SENSOR GEOMETRY (same as validated script)
%% ============================================================

%x = asind(((Re+h_air)/(alt+Re+h_air))*sind(70+90))+70;
phi = asind((Re+h_air)/(Re+h_air+alt));

FOV_HALF_ANGLE = deg2rad(phi);
cos_fov = cos(FOV_HALF_ANGLE);

%% ============================================================
% EQUAL AREA GRID (validated method)
%% ============================================================

%fprintf('[INFO] Generating equal-area Earth grid...\n')

A_earth = 4*pi*(Re+h_air)^2;
A_fov   = 2*pi*(Re+h_air)^2*(1-cosd(90-phi));

N_mingridpoints = ceil(A_earth/A_fov)*10;

r = 1;
N = N_mingridpoints;

area_cell = 4*pi*r^2/N;
d = sqrt(area_cell);

Mtheta = round(pi/d);
dtheta = pi/Mtheta;
dphi = area_cell/dtheta;

lat_vec = zeros(N,1);
lon_vec = zeros(N,1);

Ncount = 0;

for m = 0:(Mtheta-1)

    theta = pi*(m+0.5)/Mtheta;
    lat = pi/2 - theta;

    Mphi = round(2*pi*sin(theta)/dphi);

    for n = 0:(Mphi-1)

        phi2 = 2*pi*n/Mphi;

        Ncount = Ncount + 1;

        lon_vec(Ncount) = phi2 - pi;
        lat_vec(Ncount) = lat;

    end
end

lon_vec = lon_vec(1:Ncount);
lat_vec = lat_vec(1:Ncount);

Npoints = Ncount;

%fprintf('[INFO] Grid points: %d\n',Npoints);

%% ============================================================
% GROUND POSITIONS (ECEF)
%% ============================================================

rg = (Re+h_air) * [cos(lat_vec).*cos(lon_vec), ...
    cos(lat_vec).*sin(lon_vec), ...
    sin(lat_vec)]; %rg=[x1,y1,z1;
                   %    x2,y2,z2;
                   %    x3,...] one row for each point

zenith = rg ./ vecnorm(rg,2,2);

%% ============================================================
% Simulation time (repeat approximation)
%% ============================================================

Torb = 2*pi*sqrt(a^3/mu);
repeat_factor = lcm(Np,Ns);

SIM_TIME = repeat_factor*Torb;

times = 0:SIM_DT:SIM_TIME;
Nt_steps = numel(times);

Npoints = size(rg,1);

vis_matrix = false(Npoints,Nt_steps);

%% ============================================================
% Orbit propagation
%% ============================================================

n = sqrt(mu/a^3);

for p = 0:Np-1

    RAAN = 2*pi*p/Np;

    R3 = [cos(RAAN) -sin(RAAN) 0;
        sin(RAAN)  cos(RAAN) 0;
        0 0 1];

    R1 = [1 0 0;
        0 cos(inc) -sin(inc);
        0 sin(inc) cos(inc)];

    Q = R3*R1;

    for sidx = 0:Ns-1

        M0 = 2*pi*sidx/Ns + 2*pi*F*p/Nt;

        for k = 1:Nt_steps

            M = M0 + n*times(k);

            r_orb = [a*cos(M); a*sin(M); 0];

            r_eci = Q*r_orb;

            theta = omegaE*times(k);

            R = [ cos(theta)  sin(theta) 0;
                -sin(theta)  cos(theta) 0;
                0 0 1];

            rs = R*r_eci;
            rs_vec = rs';

            %% Elevation
            %we check at one instant the elevation this satellite is seen at by every gs
            ground_to_sat = rs_vec - rg;
            dist = vecnorm(ground_to_sat,2,2);

            ground_to_sat_unit = ground_to_sat./dist;

            el_sin = sum(ground_to_sat_unit.*zenith,2);

            visible_elev = el_sin > sin(ELEV_MASK);

            %% FOV check

            sat_to_ground = rg - rs_vec;
            sat_norm = vecnorm(sat_to_ground,2,2);

            cos_view = sum(sat_to_ground.*(-rs_vec),2) ./ (sat_norm.*a);

            visible_fov = cos_view >= cos_fov; %logical: (example) giorgio = 1 > 2 --> false --> giorgio=0 

            vis = visible_elev & visible_fov;

            vis_matrix(:,k) = vis_matrix(:,k) | vis; %if anyone of vis_matrix or vis = true then also ...
            % vis_matrix is true. The OR is necessary to not overwrite the vis_matrix across different satellites. 
            % Each row of the vis_matrix is relative to a point, the columns are for successive time steps

        end
    end
end

%% ============================================================
% REVISIT TIME (YOUR VALIDATED METHOD)
%% ============================================================

worst_revisit = NaN(Npoints,1);
avg_revisit   = NaN(Npoints,1);

valid_mask = false(Npoints,1);

for g = 1:Npoints

    visibility = vis_matrix(g,:); % visibility from any satellite of a specific point identified by "g" along the propagation time

    if any(visibility)

        valid_mask(g) = true; %tell me if that specific point identified by g is seen at all or not

        zero_runs = diff([0 visibility==0 0]); %at the beginning of a no-visibility-period ...
        % it places a 1, at the end places a -1. all the other values are turned to zeros. ...
        % If visibility ends with a zero then dim(zero_runs)=1+dim(visibility) and the last ...
        % element added is a -1 to reflect the last period of not being seen

        start_idx = find(zero_runs==1);
        end_idx   = find(zero_runs==-1)-1;

        run_lengths = end_idx-start_idx+1;

        run_times = run_lengths*SIM_DT;

        if isempty(run_times) %when a point is alwauys visible we observe visibility=[1,1,1,...,1,1] ...
            % so run_times is empty and causes an error. This only happens in rare cases when we have ...
            % polar orbits + large FOV or the simulation time is too short

            worst_revisit(g) = 0;
            avg_revisit(g)   = 0;

        else

            worst_revisit(g) = max(run_times);
            avg_revisit(g)   = mean(run_times);

        end


    end
end


if ~any(valid_mask) % tests along the first dimension of the array (column). ...
    % If not even one of the gridpoints is seen at least one time gives heavy penalty in the fitness function
    fitness = 1e6;
    return
end

%% ============================================================
% COVERAGE METRIC
%% ============================================================

threshold = 15*60; %seconds

good_mask = (worst_revisit<=threshold) & valid_mask;

covered_area = sum(good_mask);
total_area   = sum(valid_mask);

surface_percent = covered_area/total_area;

avg_revisit_min   = mean(avg_revisit(valid_mask))/60;
worst_revisit_min = max(worst_revisit(valid_mask))/60;
revisit_variance=mean(worst_revisit(valid_mask)) - mean(avg_revisit(valid_mask));

% 80% interested area max 15 min revisit constraint
LAT_THRESHOLD=deg2rad(75);
interesting_area_mask = (abs(lat_vec) <=  LAT_THRESHOLD);

ordered_AVG_interesting_rev_time   = sort(avg_revisit(interesting_area_mask)/60, 'ascend');
avg_revisit80_min = ordered_AVG_interesting_rev_time(1:round(0.8*length(ordered_AVG_interesting_rev_time))); %we consider the best 80%

ordered_WORST_interesting_rev_time = sort(worst_revisit(interesting_area_mask)/60, 'ascend'); %we consider the worst revisit time of each location of interest, and we order them from smallest to biggest
worst_revisit80_min = ordered_WORST_interesting_rev_time(1:round(0.8*length(ordered_WORST_interesting_rev_time))); %we consider the best 80%


%% ============================================================
% FITNESS
%% ============================================================

fitness = 5*(1-surface_percent) + ...
    3*(inc/highest_inclination) + ...
    2*revisit_variance + ...
    0.5*(Nt/MAX_SATS) + ...
    0.4*(alt) + ...
    0.1*(max(avg_revisit80_min)/15) + ...
    0.05*(max(worst_revisit80_min)/15);

end
