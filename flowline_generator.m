%% HEADER - flowline_generator_v1

% Description: This code makes shapefiles of flowlines given velocity data,
% and user input of a cross section to make starting points from.

% Instructions:
% 1. use the fastice subsetting tool here
%   (https://github.com/fastice/GrIMPNotebooks) to subset velocity data in
%   time / space.
% 2. Set the vel_fi path to that file.
% 3. Use the ncdisp function to see the attributes of time1 and time2 in 
%   that file. Set t1_ref and t2_ref times as those.
% 4. Set t_target as the time you want to use.
% 5. Set xSection_ends_fi as '' to pick new points. Or set as a .mat file
%   where you've saved the cursor_info variable with two points.
% 6. Set margin_detection type and include_margins (see description)
% 6. Set vel_ratio_cut as the fraction of maximum velocity to use as the
%   margins. (e.g., glacier margins are when velocity is 5% of the max
%   velocity in the center)
% 7. Set num_flows as the number of flowlines to be created
% 8. Set shp_out as the filename of the shapefile that will hold the
%   flowlines

% Author: Ben Reynolds (Aug 7, 2024)

close all; clear; clc

%% USER INPUTS

% velocity data inputs
vel_fi = './sermeq_kujalleq_2018/GrIMPSubset.NSIDC-0731.nc'; % file path
t1_ref = datetime(2018,1,1, 12,0,0); % reference time for start dates (see time1 attributes with ncdisp)
t2_ref = datetime(2018,1,31, 12,0,0); % reference time for end dates (see time2 attributes with ncdisp)

% target time
t_target = datetime(2018,11,15, 12,0,0); % date of velocity to use (year, month, day, hour, min, sec) 

% saved cross section points - to manually pick new points, set as empty
% character array (''). To save cross section points, save the cursor_info
% variable that gets made by exporting dataTips.
xSection_ends_fi = '';

% flowline creation paramaters
margin_detection = 3; % options: 1 = percentage of max value
                               % 2 = percentage of max-min per side
                               % 3 = maximum slope     
include_margins = 1; % 0 = git rid of margin points (e.g., for 5 flowlines make 7 points margin to margin and keep inner 5)
                     % 1 = include margin points
vel_ratio_cut = .05; % fraction of maximum velocity used to as margins of glacier (if margin_detection = 1 or 2)
num_flows = 5; % number of flowlines to make

% output file name
shp_out = 'test.shp';

%% GET VELOCITY DATA

fprintf('Getting velocity data.\n')

% get coordinates
x = double(ncread(vel_fi, 'x'));
y = double(ncread(vel_fi, 'y'));
[x_mesh, y_mesh] = meshgrid(x, y); % put into meshgrid

% get velocity data and corresponding time1 (start of avg) and t2 (end of
% avg)
vels = double(ncread(vel_fi, 'VelocitySeries'));
vels = pagetranspose(vels); % this neat function transposes the 1st two
                            % dimensions (each page of a 4d book)
time1 = double(ncread(vel_fi, 'time1'));
time2 = double(ncread(vel_fi, 'time2'));
band = (ncread(vel_fi, 'band'));

% put dates into excel format which refs of off 1900 in days
t1_ref_days = convertTo(t1_ref, 'excel');
time1_excel = time1 + t1_ref_days;
t2_ref_days = convertTo(t2_ref, 'excel');
time2_excel = time2 + t2_ref_days;

% center date for each velocity mosaic
vel_dates = (time1_excel + time2_excel) ./ 2; 

% find the velocity mosaic date closest to target date
t_target_excel = convertTo(t_target, 'excel');
t_misses = abs(vel_dates - t_target_excel);
[min_t_miss, t_use_ind] = min(t_misses);
used_date = vel_dates(t_use_ind);
used_date = datetime(used_date, 'ConvertFrom', 'excel');

% grab the velocity mosaic at the nearest date
vx = vels(:,:,1,t_use_ind);
vy = vels(:,:,2,t_use_ind);
v_mag = sqrt(vx.^2 + vy.^2);

% plot the velocity data so that user can pick points for cross section
figure(1)
s = pcolor(x,y,v_mag);
set(s, 'edgeColor', 'none')
c = colorbar;
c.Label.String = 'velocity (m yr^-^1)';
clim([0 inf])
hold on
ax1 = gca;
ax1.DataAspectRatio = [1,1,1];
title(sprintf(['Velocity ' char(string(used_date, 'yyyy-MMM-dd'))]))

if isempty(xSection_ends_fi)
    % instruct user
    fprintf('\tSelect two points (dataTips) to be endpoints of the cross section.\n')
    fprintf('\tThen right click and use ''export cursor data to workspace'' with default name.\n')
    fprintf('\tPress any key with cursor in command window once points selected\n\t\t')
    pause % pauses script for any input into the command line
    fprintf('\n')
    
    fprintf('Making Flowlines.\n')
else
    load(xSection_ends_fi)
    fprintf('\tUsing saved cross section end points from file.\n')
end


% fill in cross section struct from cursor info struct
cross_section = struct('x', {}, 'y', {});
cross_section(1).x = cursor_info(1).Position(1,1);
cross_section(1).y = cursor_info(1).Position(1,2);
cross_section(2).x = cursor_info(2).Position(1,1);
cross_section(2).y = cursor_info(2).Position(1,2);

% plot the cross section on the map view
plot([cross_section(1:end).x], [cross_section(1:end).y], 'k', 'LineWidth', 2);
plot(cross_section(2).x, cross_section(2).y, 'g.')
plot(cross_section(1).x, cross_section(1).y, 'r.')

% make dense set of points as cross section, then interp velocity data
cross_x = linspace(cross_section(1).x, cross_section(2).x, 5000);
cross_y = linspace(cross_section(1).y, cross_section(2).y, 5000);
cross_vel = interp2(x_mesh, y_mesh, v_mag, cross_x, cross_y);

% get the distance of the cross section from x and y poins
cross_x_del = diff(cross_x);
cross_y_del = diff(cross_y);
cross_diffs = sqrt(cross_x_del.^2 + cross_y_del.^2);
cross_dist = cumsum(cross_diffs); % a running sum
cross_dist = [0, cross_dist]; % add back the first point (distance = 0)

% plot the velocity on the cross section
figure(2)
plot(cross_dist, cross_vel)
ylabel('Velocity Magnitude (m yr^-^1)')
xlabel('Cross Section Distance (m)')
title('Cross Section Velocity Magnitude')
grid on; grid minor
hold on

% determine the margins of the glacier based on user-selected method
[max_vel, max_vel_ind] = max(cross_vel);
if margin_detection == 1 % ratio of max vel
    vel_margin = vel_ratio_cut * max_vel;
    glac_start_ind = find(cross_vel>vel_margin, 1, 'first'); 
    glac_end_ind = find(cross_vel>vel_margin, 1, 'last');
elseif margin_detection == 2 % ratio of per side max vel minus min vel
    vel_min1 = min(cross_vel(1:max_vel_ind));
    vel_min2 = min(cross_vel(max_vel_ind:end));
    vel_margin1 = (max_vel - vel_min1) * vel_ratio_cut + vel_min1;
    vel_margin2 = (max_vel - vel_min2) * vel_ratio_cut + vel_min2;
    glac_start_ind = find(cross_vel>vel_margin1, 1, 'first'); 
    glac_end_ind = find(cross_vel>vel_margin2, 1, 'last');
elseif margin_detection == 3 % highest slope
    grad_cross_vel = gradient(cross_vel, (cross_dist(2)-cross_dist(1)));
    [~, glac_start_ind] = max(grad_cross_vel);
    [~, glac_end_ind] = min(grad_cross_vel);
else
    fprintf('\Error: margin_detection must be set as 1, 2, or 3.\n')
end

% put the glacier margin points on the cross section velocity plot
plot(cross_dist(glac_start_ind), cross_vel(glac_start_ind), 'rx', 'MarkerSize', 10)
plot(cross_dist(glac_end_ind), cross_vel(glac_end_ind), 'rx', 'MarkerSize', 10)

if include_margins == 0
    flow_x = linspace(cross_x(glac_start_ind), cross_x(glac_end_ind), (num_flows+2));
    flow_y = linspace(cross_y(glac_start_ind), cross_y(glac_end_ind), (num_flows+2));
    flow_dist = linspace(cross_dist(glac_start_ind), cross_dist(glac_end_ind), (num_flows+2));
    
    % get rid of the end points (points on the margin)
    flow_x = flow_x(2:num_flows+1);
    flow_y = flow_y(2:num_flows+1);
    flow_dist = flow_dist(2:num_flows+1);
elseif include_margins == 1
    flow_x = linspace(cross_x(glac_start_ind), cross_x(glac_end_ind), num_flows);
    flow_y = linspace(cross_y(glac_start_ind), cross_y(glac_end_ind), num_flows);
    flow_dist = linspace(cross_dist(glac_start_ind), cross_dist(glac_end_ind), num_flows);
else
    fprintf('\tError: select 0 or 1 for include_margins.\n')
end

% interpolate velocity onto points for plotting
flow_vel = interp1(cross_dist, cross_vel, flow_dist);

% plot the flowline points on the cross section plot
plot(flow_dist, flow_vel, 'bs', 'MarkerSize', 10)

% highlight margin to margin in red on the map view plot
figure(1)
plot([cross_x(glac_start_ind), cross_x(glac_end_ind)],...
    [cross_y(glac_start_ind), cross_y(glac_end_ind)], 'r', 'LineWidth', 2)


% plot for the flowlines
plot(flow_x, flow_y, 'bs')

% now loop through each of the points, create the flowline, plot it, and
% save it to the struct that will be used to write the shape file
S = struct();
for ii = 1:num_flows

    % get vertices of streamlines (as cell arrays with one cell that holds
    % two column vectors that are the x and y points)
    verts_fwd = stream2(x_mesh, y_mesh, vx, vy, flow_x(ii), flow_y(ii));
    verts_bkwd = stream2(x_mesh, y_mesh, -vx, -vy, flow_x(ii), flow_y(ii));
    
    % pull em out of the cell arrays
    verts_fwd = verts_fwd{1};
    verts_bkwd = verts_bkwd{1};
    
    % flip the one going forward towards the terminus (such that the first
    % point is at the terminus)
    verts_fwd = flipud(verts_fwd(2:end, :));

    % split into individual x and y column vectors
    verts_fwd_x = verts_fwd(:,1);
    verts_fwd_y = verts_fwd(:,2);
    verts_bkwd_x = verts_bkwd(:,1);
    verts_bkwd_y = verts_bkwd(:,2);

    % get rid of all the nans
    verts_fwd_x = verts_fwd_x(~isnan(verts_fwd_x));
    verts_fwd_y = verts_fwd_y(~isnan(verts_fwd_y));
    verts_bkwd_x = verts_bkwd_x(~isnan(verts_bkwd_x));
    verts_bkwd_y = verts_bkwd_y(~isnan(verts_bkwd_x));
    
    % put the upstream and downstream flowline segments together
    verts_x = [verts_fwd_x; verts_bkwd_x];
    verts_y = [verts_fwd_y; verts_bkwd_y];
    
    % plot the flowlines
    plot(verts_x, verts_y)

    % fill out the struct with the necessary fields to write a shapefile
    S(ii).Geometry = 'Line';
    S(ii).BoundingBox = [min(verts_x), min(verts_y), max(verts_x), max(verts_y)];
    S(ii).X = verts_x;
    S(ii).Y = verts_y;
    S(ii).flowline = char(string(ii));
end

% get rid of the data tips to see cross section and flowlines
figure(1)
dt = findobj(ax1,'Type','DataTip');
delete(dt)

% make the shapefile
fprintf('Writing shape file.\n')
shapewrite(S, shp_out)


