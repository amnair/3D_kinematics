function plot_traj(seq_path)


%% Prompt for seq directory

if nargin < 1
   seq_path = uigetdir(pwd,'Choose sequence to plot');
end

data_path = [seq_path filesep 'coord_data.mat'];

if ~isempty(dir(data_path))
    
    % Load 'd'
    load(data_path)
      
else
    error('There must be a coord_data.mat file in the selected directory')  
    
end

%% Load other data

cd(seq_path)

load(['..' filesep 'cal_data.mat'])


% Locate index of frames with complete data

idx = ~isnan(d.top.leftEye(:,1)) & ~isnan(d.side.leftEye(:,1));


%% Calculate coordinates in global system

d3.t = (d.frames(idx)-min(d.frames(idx)))./d.frame_rate;

d3.x = cal.camT.const .* (d.top.leftEye(idx,1)-cal.camT.pnt(1));
d3.y = cal.camT.const .* (d.top.leftEye(idx,2)-cal.camT.pnt(2));
d3.z = cal.camS.const .* (d.side.leftEye(idx,2)-cal.camS.pnt(2));


%% Plot
figure;

yrange = [0 max([range(d3.y) range(d3.z) range(d3.x)])];

subplot(3,2,1)
plot(d3.t,d3.x-d3.x(1),'r')
xlabel('time (s)')
ylabel('x (mm)')
ylim(yrange)

subplot(3,2,3)
plot(d3.t,d3.y-d3.y(1),'r')
xlabel('time (s)')
ylabel('y (mm)')
ylim(yrange)

subplot(3,2,5)
plot(d3.t,d3.z-d3.z(1),'r')
xlabel('time (s)')
ylabel('z (mm)')
set(gca,'YDir','reverse')
ylim(yrange)

subplot(3,2,[2 4 6])
plot3(d3.x,d3.y,d3.z,'r',d3.x(1),d3.y(1),d3.z(1),'ro')

xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')
set(gca,'ZDir','reverse')
axis equal
grid on

ttt= 3