function ana_all(root_path)
% Analyzes the trajectory data for all sequences 


%% Parameters


% Number of divisions in each quadrant for summary polar plot
num_div = 3;


% Number of points for drawing the full circle
num_pts = 1000;

% Time to evaluate direction (ms)
eval_time = 19.6;


%% Prompt for seq directory

if nargin < 1
   root_path = uigetdir(pwd,'Choose root directory');
end

a = dir(root_path);


%% Check for pool_data file

if ~isempty(dir([root_path filesep 'pool_data.mat']))

    load([root_path filesep 'pool_data.mat'])
    data_pooled = 1;
    
else
    data_pooled = 0;
end


%% Take inventory of directories

if ~data_pooled
    k = 1;
    for i = 3:length(a)
        if a(i).isdir && strcmp(a(i).name(1:9),'Zebrafish') && ...
                ~isempty(dir([root_path filesep a(i).name filesep 'cal_data.mat']))
            
            b = dir([root_path filesep a(i).name]);
            
            for j = 3:length(b)
                if strcmp(b(j).name(1),'0') && length(b(j).name)==3 && ...
                        ~isempty(dir([root_path filesep a(i).name filesep b(j).name filesep 'coord_data.mat']))
                    
                    seq(k).path = [root_path filesep a(i).name filesep b(j).name filesep 'coord_data.mat'];
                    seq(k).cal_path = [root_path filesep a(i).name filesep 'cal_data.mat'];
                    seq(k).indiv = a(i).name;
                    seq(k).seq_name = b(j).name;
                    
                    k = k+1;
                end
            end
            
            clear b
        end
    end
    
    if k == 1
        error('There are no directories that start with "Zebrafish" in the selected directory')
    end
    
    clear a k i
    
end


%% Analyze data for each sequence
if ~data_pooled
    for i = 1:length(seq)
        
        if ~data_pooled
            % Load calibration data in 'cal'
            load(seq(i).cal_path)
            
            % Load coordinate data in 'd'
            load(seq(i).path)
            
            % Locate index of frames with complete data
            
            idx = ~isnan(d.top.leftEye(:,1)) & ~isnan(d.side.leftEye(:,1));
            
            
            % Calculate coordinates in global system
            
            d3(i).t = (d.frames(idx)-min(d.frames(idx)))./d.frame_rate;
            
            d3(i).x = cal.camT.const .* (d.top.leftEye(idx,1)-cal.camT.pnt(1));
            d3(i).y = cal.camT.const .* (d.top.leftEye(idx,2)-cal.camT.pnt(2));
            d3(i).z = cal.camS.const .* (d.side.leftEye(idx,2)-cal.camS.pnt(2));
        end
        
        
        
        
        clear cal d
    end
end


% Calculate angle traveled at end of trajectory
for i = 1:length(d3)
    
    tmp = d3(i).t-eval_time./1000;
    time_idx = find(tmp==min(abs(tmp)),1,'first');
    
    d3(i).time_idx = time_idx;
    
    
    %atan2(d3(i).z(end),d3(i).x(end))
    plr.theta_xz(i,1)  = atan2(d3(i).z(time_idx)-d3(i).z(1),...
                               d3(i).x(time_idx)-d3(i).x(1));
                           
    plr.theta_yz(i,1)  = atan2(d3(i).z(time_idx)-d3(i).z(1),...
                                d3(i).y(time_idx)-d3(i).y(1));
end


%% Save data for all sequences

%save(['/Users/mmchenry/Dropbox/Matlab/mueller_3d' filesep 'pool_data.mat'],'d3','plr')



%% Plot all

msize = 2;

figure

for i = 1:length(d3)
    subplot(1,2,1)
    plot(d3(i).x-d3(i).x(1),d3(i).z-d3(i).z(1),'k-')
    hold on
    h = plot(d3(i).x(d3(i).time_idx)-d3(i).x(1),...
         d3(i).z(d3(i).time_idx)-d3(i).z(1),'ro');
    set(h,'MarkerSize',msize)
    set(h,'MarkerFaceColor','r')
    axis equal
    set(gca,'YDir','reverse')
    xlabel('x (mm)')
    ylabel('z (mm)')
    
    subplot(1,2,2)
    plot(d3(i).y-d3(i).y(1),d3(i).z-d3(i).z(1),'k-')
     hold on
    h = plot(d3(i).y(d3(i).time_idx)-d3(i).y(1),...
         d3(i).z(d3(i).time_idx)-d3(i).z(1),'ro');
    set(h,'MarkerSize',msize)
    set(h,'MarkerFaceColor','r')
    xlabel('y (mm)')
    ylabel('z (mm)')
    axis equal
    set(gca,'YDir','reverse')
end

subplot(1,2,1)
plot(0,0,'ro')

subplot(1,2,2)
plot(0,0,'ro')



%% Polar plot

figure
subplot(1,2,1)
polar(-plr.theta_xz,ones(length(plr.theta_xz),1),'o')
title('xz plane')

subplot(1,2,2)
polar(-plr.theta_yz,ones(length(plr.theta_yz),1),'o')
title('yz plane')


%% Summary polar plot

figure

subplot(1,2,1)
plot_sum_bar(plr.theta_xz,num_div,num_pts)
ylims = xlim;
text(0,1.1*ylims(2),'The xz plane')

subplot(1,2,2)
plot_sum_bar(plr.theta_yz,num_div,num_pts)
ylims = xlim;
text(0,1.1*ylims(2),'The yz plane')


function plot_sum_bar(theta_xz,num_div,num_pts)

% Increment for angular bins
ang_div = pi/2/num_div;

% Vector of angle groups
div_vect = [0:ang_div:2*pi]-ang_div/2;

% Step through each angular bin
for i = 2:length(div_vect)
    
    % Define angle values that span the bin
    theta = linspace(div_vect(i-1),div_vect(i),round(num_pts/num_div))';
    
    % Index for values in current bin
    idx = (theta_xz > div_vect(i-1)) & (theta_xz <= div_vect(i)); 
    
    % Find mean and 95 CI for the bin
    [phat,pci] = binofit(length(theta_xz(idx)),length(theta_xz));
    
    % Fine x & z values for that angular span
    x = phat.*cos(theta);
    z = -phat.*sin(theta);
    
    % Plot x & z values for mean of bin
    h = fill([0; x; 0],[0; z; 0],[1 0 0]);
    set(h,'EdgeColor','w')
    set(h,'FaceColor',[1 .4 .2])
    hold on
    
    clear x y
    
    % Define error flag for 95% CI of bin
    x_ci = pci.*cos(mean(theta));
    z_ci = -pci.*sin(mean(theta));
    
    % Plot error flags
    h = plot(x_ci,z_ci,'k');
    
    set(h,'Color',.5.*[1 1 1]);
    axis equal
    
    % Clear for next loop
    clear mu s muci sci
end


% Modify appearence
 set(gca,'Visible','off')
 
 xtick = get(gca,'XTick');
 xtick = xtick(xtick>0);
 xtick = round(xtick*1000)/1000;
 theta = linspace(0,2*pi,1000);
 offset = range(xtick)/10;

for i = 1:length(xtick) 
    h = plot(xtick(i).*cos(theta),xtick(i).*sin(theta),'k');
    set(h,'Color',.7.*[1 1 1])
    text(xtick(i)+offset,-offset/2,num2str(xtick(i)));
end
 
 
