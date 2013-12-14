function vis_3Dkine


%% Define paths, load data

% Root path
if isdir('/Users/mmchenry/Dropbox/Projects/3D kinematics/')
    root = '/Users/mmchenry/Dropbox/Projects/3D kinematics/';
elseif isdir('/Volumes/Flow HD/Dropbox/Projects/3D kinematics')
    root = '/Volumes/Flow HD/Dropbox/Projects/3D kinematics';
else
    error('You need to specify a root directory')
end

% Path to save animation files
ani_path = [root filesep 'new animation'];

% Load kinematic data ('dorsalTailCoord', 'headCoord', 'ventralTailCoord')
load([root filesep 'kinematicData.mat'])

% Load morphological data ('m')
load([root filesep '6_02_L19_metrics.mat']);
    

%% Parameters

% Number of points around the periphery of the prey body & tail fin
numPts_circ = 50;

% Color of surface of prey body
preyColor = .5.*[1 1 1];

% Execute rendering of prey in global FOR
render_global = 1;

% Execute rendering of prey in local FOR
render_local = 1;

% Constant by which to expand the limits of axes
lim_const = 0.1;

% Magnfier of figure window size
mag_fig = 2;


%% Determine moprhological parameters

% Step thru time
for i = 1:length(headCoord{1})
    
    % Extract trunk coordinates for current frame
    Leye    = headCoord{1}(:,i)'./1000;
    Reye    = headCoord{2}(:,i)'./1000;
    SB      = headCoord{3}(:,i)'./1000;
    
    % Extract tail cooridnate for current frame
    Dtail = dorsalTailCoord{i}'./1000;
    Vtail = ventralTailCoord{i}'./1000;
    
    % Make origin between eyes
    origin = mean([Leye;Reye]);
    
    % Define transformation matrix 
    S = localSystem(Leye,origin,SB);
    
    % Transform tail data into local FOR
    Dtail_L = globalToLocal(Dtail,origin,S);
    Vtail_L = globalToLocal(Vtail,origin,S);
    
    % Distance to swim bladder
    s_SB(i) = sqrt((SB(1)-origin(1)).^2 + ...
                   (SB(2)-origin(2)).^2 + ...
                   (SB(3)-origin(3)).^2);
            
    % Calculate arclength positions of tail margins
    s_D = cumsum(sqrt(diff(Dtail_L(:,1)).^2 + diff(Dtail_L(:,2)).^2));
    s_V = cumsum(sqrt(diff(Dtail_L(:,1)).^2 + diff(Dtail_L(:,2)).^2));
          
    b_length(i) = min([(Dtail_L(1,1) + max(s_D));(Vtail_L(1,1) + max(s_V))]);
    tail_length(i) = min([ max(s_D);max(s_V)]);
end

% Recalculate, considering all values and convert to meters
%b_length = max(b_length);
SB_pos     = mean(s_SB);

% Clear unneeded for next 
clear Leye Reye s_D s_V Dtail Vtail origin S Dtail_L Vtail_L s_SB


%% Calculate polygons

xMax = []; yMax = []; zMax = [];
xMin = []; yMin = []; zMin = [];

% Step thru time
for i = 1:length(headCoord{1})
    
    % Extract trunk coordinate for current frame (in m)
    Leye    = headCoord{1}(:,i)'./1000;
    Reye    = headCoord{2}(:,i)'./1000;
    SB      = headCoord{3}(:,i)'./1000;
    
    % Extract tail cooridnate for current frame (in m)
    Dtail = dorsalTailCoord{i}'./1000;
    Vtail = ventralTailCoord{i}'./1000;
    
    % Make origin between eyes
    origin = mean([Leye;Reye]);
    
    % Define transformation matrix 
    S = localSystem(Reye,origin,SB);
    
    % Transform tail data into local FOR
    Dtail_L = globalToLocal(Dtail,origin,S);
    Vtail_L = globalToLocal(Vtail,origin,S);
    
    % Define 3d data for prey in local FOR
    [pX,pY,pZ] = prey_surf(m,Dtail_L,Vtail_L,numPts_circ,b_length(i),SB_pos);
    
    % Returns coordinates for prey body in global FOR
    [pXg,pYg,pZg] = prey_global(pX,pY,pZ,origin,S);
       
    % Store results
    R(i).L.x       = pX;
    R(i).L.y       = pY;
    R(i).L.z       = pZ;
    R(i).G.x       = pXg;
    R(i).G.y       = pYg;
    R(i).G.z       = pZg;
    R(i).S         = S;
    R(i).origin    = origin;
    
    % Set limits of global domain
    xMax    = max([xMax max(pXg(:))]);
    yMax    = max([yMax max(pYg(:))]);
    zMax    = max([zMax max(pZg(:))]);
    xMin    = min([xMin min(pXg(:))]);
    yMin    = min([yMin min(pYg(:))]);
    zMin    = min([zMin min(pZg(:))]);
 
    % Clear values for next iteration
    clear Leye Reye SB Dtail Vtail origin S Dtail_L Vtail_L
end
    

%% Render in global FOR

if render_global

    % Set figure window
    hF = figure;
    set(hF,'DoubleBuffer','on');
    set(hF,'WindowStyle','modal') 
    
     % Adjust position
    f_pos = get(hF,'Position');
    set(hF,'Position',[f_pos(1) f_pos(2) mag_fig*f_pos(3) mag_fig*f_pos(4)])
    
    % Set lights
    hL(1) = light('position',[0 0 20]);
    hL(2) = light('position',[0 0 -20]);

    for i = 1:length(R)

        % Render the prey 
        h = patch(R(i).G.x,R(i).G.y,R(i).G.z,R(i).G.z*0);
        
        % Set properties
        set(h,'FaceLighting','gouraud',...
            'LineStyle','none',...
            'BackFaceLighting','reverselit',...
            'FaceColor',preyColor,...
            'AmbientStrength',.5);
        hold on
        
        % Set axes
        axis square
        xlim([xMin*(1-lim_const) xMax*(1+lim_const)])
        ylim([yMin*(1-lim_const) yMax*(1+lim_const)])
        zlim([zMin*(1-lim_const) zMax*(1+lim_const)])
        
        % Set axes & titles
        xlabel('X');ylabel('Y');zlabel('Z');
        set(gca,'XColor','w','YColor','w','ZColor','w')
        title(['Frame ' num2str(i)])
        view([162 21])
        
        % Capture frame
        figure(hF)
        F = getframe(hF);
        pause(0.01);
        
        % Save image
        fnum = ['00' num2str(i)];
        imname = ['Global frame ' fnum(end-2:end) '.tif'];      
        imwrite(F.cdata,colormap,[ani_path filesep imname],'tif'); 
        
        % Delete for next iteration
        delete(h)
    end 
    
    close(hF)
end


%% Render in local FOR

if render_local

    % Set figure window
    hF = figure;
    set(hF,'DoubleBuffer','on');
    set(hF,'WindowStyle','modal')
    
    % Adjust position
    f_pos = get(hF,'Position');
    set(hF,'Position',[f_pos(1) f_pos(2) mag_fig*f_pos(3) mag_fig*f_pos(4)])
    
    ax_range = 4.5e-3;
    
    % Set lights
    subplot(2,1,1)
    hL(1) = light('position',[0 0 20]);
    hL(2) = light('position',[0 0 -20]);
    
    subplot(2,1,2)
    hL(1) = light('position',[0 0 20]);
    hL(2) = light('position',[0 0 -20]);
    
    for i = 1:length(R)

        % Render the prey 
        subplot(2,1,1)
        h(1) = patch(R(i).L.x,R(i).L.y,R(i).L.z,R(i).L.z*0);
        view([0 90])
        ha(1) = gca;
            
        subplot(2,1,2)
        h(2) = patch(R(i).L.x,R(i).L.y,R(i).L.z,R(i).L.z*0);
        view([0 0])
        ha(2) = gca;
          
        for j = 1:2
            
            % Set properties
            set(h(j),'FaceLighting','gouraud',...
                'LineStyle','none',...
                'BackFaceLighting','reverselit',...
                'FaceColor',preyColor,...
                'AmbientStrength',.5);
            hold on
            
            % Set axes
            axes(ha(j))
            axis square
            xlim([0 ax_range])
            ylim([-ax_range/2 ax_range/2])
            zlim([-ax_range/2 ax_range/2])
            
            % Set axes & titles
            xlabel('X');ylabel('Y');zlabel('Z');
            set(gca,'XColor','w','YColor','w','ZColor','w')
            
        end
        
        axes(ha(1))
        title(['Frame ' num2str(i)])
        
        % Capture frame
        figure(hF)
        F = getframe(hF);
        pause(0.01);
        
        % Save image
        fnum = ['00' num2str(i)];
        imname = ['Local frame ' fnum(end-2:end) '.tif'];      
        imwrite(F.cdata,colormap,[ani_path filesep imname],'tif'); 
        
        % Delete for next iteration
        delete(h)
    end 
    
    close(hF)
end




return
     
    
   
    

function [X,Y,Z]= prey_surf(m,tailD,tailV,numPts,b_length,s_SB)
% Provides 3D coordinates of the surface of the prey body

eye_pos = 0.32e-3;

preyColor = .5.*[1 1 1];


s_SB = s_SB + eye_pos;

% Define radial positions along vector
theta = linspace(0,2*pi,numPts)';

% Define empty vectors for coordinates
x=[];y=[];z=[];

% Extract and scale morph data
h = m.h;
w = m.w;
s = m.s;
c = -m.c;

% Clear for next
clear m 

% Index of tail start
iTailStart = find(s >= s_SB,1,'first');

% Offset all centers wrt point of tail start 
c = c - c(iTailStart);

% Matching arclength positions of dorsal & ventral tail coord _____________

% Number of fin points should match number of body points posterior to the
% start of the fin
numFin = sum(s > s(iTailStart));

% X-position of kinematic data to start considering for interpolaton
%v_post = tailV(1,1) > tailD(1,1);
xStart = max([tailD(1,1) tailV(1,1)]);
i_tailD = find(tailD(:,1)>=xStart,1,'first'):length(tailD(:,1));
i_tailV = find(tailV(:,1)>=xStart,1,'first'):length(tailV(:,1));

% Calculate arclength positions of tail margins
s_D = cumsum(sqrt(diff(tailD(i_tailD,1)).^2 + diff(tailD(i_tailD,2)).^2));
s_V = cumsum(sqrt(diff(tailV(i_tailV,1)).^2 + diff(tailV(i_tailV,2)).^2));

% Scaling constant 
const = (max(s)-s(iTailStart))./min([s_D(end) s_V(end)]);

% Scale and position the front end of the kinematic data
tailV(:,1) = const.* (tailV(:,1) - xStart) + s(iTailStart);
tailV(:,2) = const .* (tailV(:,2) - tailV(1,2));
tailV(:,3) = const .* tailV(:,3);

tailD(:,1) = const.* (tailD(:,1) - xStart) + s(iTailStart);
tailD(:,2) = const .* (tailD(:,2) - tailD(1,2));
tailD(:,3) = const .* tailD(:,3);

% Measured arclength positions of the two margins
sBackD = [0; cumsum(sqrt(diff(tailD(i_tailD,1)).^2 + ...
                         diff(tailD(i_tailD,2)).^2 + ...
                         diff(tailD(i_tailD,3)).^2))];
sBackV = [0; cumsum(sqrt(diff(tailV(i_tailV,1)).^2 + ...
                         diff(tailV(i_tailV,2)).^2 + ...
                         diff(tailV(i_tailV,3)).^2))];  

% Arclength positions to interpolate tail
s_tail = linspace(0, min([max(sBackD) max(sBackV)]),numFin);

% Run interpolation
iTailV(:,1) = interp1(sBackV,tailV(i_tailV,1),s_tail);
iTailV(:,2) = interp1(sBackV,tailV(i_tailV,2),s_tail);
iTailV(:,3) = interp1(sBackV,tailV(i_tailV,3),s_tail);

iTailD(:,1) = interp1(sBackD,tailD(i_tailD,1),s_tail);
iTailD(:,2) = interp1(sBackD,tailD(i_tailD,2),s_tail);
iTailD(:,3) = interp1(sBackD,tailD(i_tailD,3),s_tail);

clear tailV tailD numFin minZ

% Overwrite old tail margin coordinates
tailV  = iTailV; 
tailD  = iTailD; 

% Adjust z-positions
midZ = mean([tailV(1,3) tailD(1,3)]);
tailV(:,3) = tailV(:,3) - midZ;
tailD(:,3) = tailD(:,3) - midZ;

% Define tail midline
tail(:,1) = mean([tailV(:,1) tailD(:,1)],2);
tail(:,2) = mean([tailV(:,2) tailD(:,2)],2);
tail(:,3) = mean([tailV(:,3) tailD(:,3)],2);

% Clear for next
clear iTailD iTailV sBackD sBackV idxD idxV

% Visualize raw data
if 0
%     %figure
     subplot(1,2,1)
     plot(s,w/2,'k',s,-w/2,'k',tailD(:,1),tailD(:,2),'r',...
          tailV(:,1),tailV(:,2),'b',tail(:,1),tail(:,2),'g--')
     axis equal
     xlabel('X');ylabel('Y')
    grid on
    
    subplot(1,2,2)
    plot(s,c+h/2,'k',s,c-h/2,'k',tailD(:,1),tailD(:,3),'r',...
         tailV(:,1),tailV(:,3),'b',tail(:,1),tail(:,3),'g--')
    axis equal
    xlabel('X');ylabel('Z')
    grid on
    pause(.3)
end


% Make mouth cap  _________________________________________________________
n = numPts/10;
phi = linspace(0,.75*pi/2,n)';
ds = .02.*range(s); %2*s(2)-s(1);
%sC = linspace(s(1)-ds,s(1),n);
hC = h(1) .* sin(phi)./max(sin(phi));
wC = w(1) .* sin(phi)./max(sin(phi));
sC = -(ds.*cos(phi)-ds.*cos(phi(end)));

% Loop down the body length
for i=1:length(sC)-1  
    
    % Draw first ellipse   
    xTemp1 = sC(i)*ones(size(theta));
    yTemp1 = (wC(i)/2) .* cos(theta);
    zTemp1 = (hC(i)/2) .* sin(theta) + c(1);
    
    % Draw second ellipse  
    xTemp2 = sC(i+1)*ones(size(theta));
    yTemp2 = (wC(i+1)/2) .* cos(theta);
    zTemp2 = (hC(i+1)/2) .* sin(theta) + c(1);
    
    % Combine data (works with 'patch')
    x	= [x [xTemp1(1:end-1)';... 
              xTemp2(1:end-1)';... 
              xTemp2(2:end)';... 
              xTemp1(2:end)']];
                      
    y   = [y [yTemp1(1:end-1)';... 
              yTemp2(1:end-1)';...
              yTemp2(2:end)';...
              yTemp1(2:end)']];
                      
    z   = [z [zTemp1(1:end-1)';...
              zTemp2(1:end-1)';...
              zTemp2(2:end)';...
              zTemp1(2:end)']];
end 

clear xTemp1 yTemp1 zTemp1 xTemp2 yTemp2 zTemp2
clear n phi ds sC hC wC


% Make rigid trunk coordinates __________________________________

% Loop down the body length
for i=1:length(s)-1  
    
    % Stop rendering where tail fin data begins
    if s(i) >= s(iTailStart)
        i_last = i;
        break
    end
    
  % Draw first ellipse  
    xTemp1      = s(i)*ones(size(theta));
    yTemp1      = (w(i)/2) .* cos(theta);
    zTemp1      = (h(i)/2) .* sin(theta) + c(i);
    
  % Draw second ellipse    
    xTemp2      = s(i+1)*ones(size(theta));
    yTemp2      = (w(i+1)/2) .* cos(theta);
    zTemp2      = (h(i+1)/2) .* sin(theta) + c(i+1);
    
  % Combine data (works with 'patch')
    x	= [x [xTemp1(1:end-1)';... 
              xTemp2(1:end-1)';... 
              xTemp2(2:end)';... 
              xTemp1(2:end)']];
                      
    y   = [y [yTemp1(1:end-1)';... 
              yTemp2(1:end-1)';...
              yTemp2(2:end)';...
              yTemp1(2:end)']];
                      
    z   = [z [zTemp1(1:end-1)';...
              zTemp2(1:end-1)';...
              zTemp2(2:end)';...
              zTemp1(2:end)']];
    if 0
        hand = patch(x,y,z,z.*0);
        title(['i =' num2str(i)])
         % Set properties
%         set(hand,'FaceLighting','gouraud',...
%             'LineStyle','none',...
%             'BackFaceLighting','reverselit',...
%             'FaceColor',preyColor,...
%             'AmbientStrength',.5);
        hold on
    end
end 

% This links the trunk to tail coordinates    
E_last = [xTemp2 yTemp2 zTemp2];

% Make cellular tail coordinates __________________________________

% Initiate index for tail coordinates
j = 1;

% Loop down the body length
for i = (i_last+1):length(s)-1  
    
    % Define local system for tail segment
    S = localSystemZ(tail(j+1,:),tail(j,:),tailD(j,:));
    %S = localSystemZ(tail(j+1,:),tail(j,:),[0 0 1]);
    
    % Draw first ellipse
    if 0 %i == (i_last+1)
        E1(:,1) = zeros(size(theta));
        E1(:,2) = (w(i)/2) .* cos(theta);
        E1(:,3) = (h(i)/2) .* sin(theta);
        origin  = [tail(j,1:2) tail(j,3)];
        E1      = localToGlobal(E1,origin,S);
    else
        E1 = E_last;
    end
   
    % Draw second ellipse    
    E2(:,1) = zeros(size(theta));
    E2(:,2) = (w(i+1)/2) .* cos(theta);
    E2(:,3) = (h(i+1)/2) .* sin(theta);
    origin  = [tail(j+1,1:2) tail(j+1,3)];
    E2      = localToGlobal(E2,origin,S);
    
    % Combine data (works with 'patch')
    x	= [x [E1(1:end-1,1)';... 
              E2(1:end-1,1)';... 
              E2(2:end,1)';... 
              E1(2:end,1)']];
                      
    y   = [y [E1(1:end-1,2)';... 
              E2(1:end-1,2)';...
              E2(2:end,2)';...
              E1(2:end,2)']];
                      
    z   = [z [E1(1:end-1,3)';...
              E2(1:end-1,3)';...
              E2(2:end,3)';...
              E1(2:end,3)']];
          
    E_last = E2;
    
    %clear E1 E2 S
    
    j = j + 1;
    
    if 0
    hand = patch(x,y,z,z.*0);
    title(['i =' num2str(i)])
    % Set properties
    %         set(hand,'FaceLighting','gouraud',...
    %             'LineStyle','none',...
    %             'BackFaceLighting','reverselit',...
    %             'FaceColor',preyColor,...
    %             'AmbientStrength',.5);
    axis equal
    view([180 0])
    tttt=2;
    %pause
    
    end
end 



clear E_last i j


% % Make tail cap  _______________________________________
% n = numPts/10;
% phi = linspace(0,0.75*pi/2,n)';
% ds = .02.*range(s); %2*s(2)-s(1);
% %sC = linspace(s(1)-ds,s(1),n);
% hC = h(end) .* sin(phi)./max(sin(phi));
% wC = w(end) .* sin(phi)./max(sin(phi));
% sC = s(end) + ds.*cos(phi)-+ ds.*cos(phi(end));
% 
% % Loop down the body length
% for i=1:length(sC)-1  
%     
%   % Draw first ellipse   
%     xTemp1 = sC(i)*ones(size(theta));
%     yTemp1 = (wC(i)/2) .* cos(theta);
%     zTemp1 = (hC(i)/2) .* sin(theta) + c(end);
%     
%   % Draw second ellipse  
%     xTemp2 = sC(i+1)*ones(size(theta));
%     yTemp2 = (wC(i+1)/2) .* cos(theta);
%     zTemp2 = (hC(i+1)/2) .* sin(theta) + c(end);
%     
%   % Combine data (works with 'patch')
%     x	= [x [xTemp1(1:end-1)';... 
%               xTemp2(1:end-1)';... 
%               xTemp2(2:end)';... 
%               xTemp1(2:end)']];
%                       
%     y   = [y [yTemp1(1:end-1)';... 
%               yTemp2(1:end-1)';...
%               yTemp2(2:end)';...
%               yTemp1(2:end)']];
%                       
%     z   = [z [zTemp1(1:end-1)';...
%               zTemp2(1:end-1)';...
%               zTemp2(2:end)';...
%               zTemp1(2:end)']];
% end 

clear xTemp1 yTemp1 zTemp1 xTemp2 yTemp2 zTemp2
clear n phi ds sC hC wC

% Transform for units and wrt rostrum
% X  = (y-min(y(:)));
% Y  = x;
% Z  = -(z-z(1)+max(h(:))/2);
X = x; Y = y; Z = z;


function ptsT = local_to_global(rost,com,tail,pts)
% Transforms coordinates (pts, in n x 3) from global to local coordinate
% system, assuming the y-axis of larvae runs perpendicular to gravity

% Check dimensions
if size(rost,1)~=1 || size(rost,2)~=3 || size(com,1)~=1 || size(com,2)~=3 ...
   size(tail,1)~=1 || size(tail,2)~=3 || size(pts,2)~=3 
    error('inputs have incorrect dimensions')
end

% Retrieve local x axis to determine coordinate system
xaxis(1,1) = tail(1) - rost(1);
xaxis(1,2) = tail(2) - rost(2);
xaxis(1,3) = tail(3) - rost(3);

% Normalize to create a unit vector
xaxis = xaxis./norm(xaxis);

%Determine local y axis
%Short hand of cross product of inertial z axis and local x axis
yaxis = [-xaxis(2) xaxis(1) 0];

% Normalize to create a unit vector
yaxis = yaxis./norm(yaxis);

%Determine local z axis
zaxis = cross(xaxis,yaxis);

% Normalize to create a unit vector
zaxis = zaxis./norm(zaxis);

%Create rotation matrix (from inertial axes to local axes)
R = [xaxis' yaxis' zaxis'];

% Rotate points
ptsT = [R * pts']';

% Translate global coordinates wrt rostrum
ptsT(:,1) = ptsT(:,1) + rost(1);
ptsT(:,2) = ptsT(:,2) + rost(2);
ptsT(:,3) = ptsT(:,3) + rost(3);

% Visualize to test
if 0
    
    blength = norm([tail(1)-rost(1) tail(2)-rost(2) tail(3)-rost(3)]);    
    
    figure
    
    subplot(2,2,[1 3])
    plot3([tail(1) rost(1)],[tail(2) rost(2)],[tail(3) rost(3)],'b',...
          rost(1),rost(2),rost(3),'bo');
    hold on
    plot3(ptsT(:,1),ptsT(:,2),ptsT(:,3),'ro')
    xlabel('x'); ylabel('y'); zlabel('z')
    hold off
    grid on;axis equal
    view(3)
    title('global')
    
    subplot(2,2,2)
    plot([0 blength],[0 0],'b',0,0,'ob',pts(:,1),pts(:,2),'ro')
    xlabel('x');ylabel('y')
    grid on; axis equal
    title('local')
    
    subplot(2,2,4)
    plot([0 blength],[0 0],'b',0,0,'ob',pts(:,1),pts(:,3),'ro')
    xlabel('x');ylabel('z')
    grid on; axis equal
end



function S = localSystem(P1,P2,P3)
% Defines a transformation matrix for a local coordinate system in an
% inertial frame of reference.  Uses P1 as the yaxis and P2 as the origin, and 
% P3 as the x-axis. Coordinates must be (1x3) vectors. Note: if theses axes 
% are not orthogonal, the x-axis direction is assumed to be more accurate
% than the z-axis and the y-axis direction is adjusted to make the coordinates 
% orthoganal.
 
% Check dimensions of inputs
if size(P1,1)~=1 || size(P1,2)~=3 ||...
   size(P2,1)~=1 || size(P2,2)~=3 ||...
   size(P3,1)~=1 || size(P3,2)~=3
    error('Coordinates must be 1x3 vectors');
end
 
% Define units vectors for x and y axes
xAxis   = (P3-P2)./norm(P3-P2);
yAxis   = (P1-P2)./norm(P1-P2);
 
% Define yaxis from the cross product
zAxis   = cross(xAxis,yAxis);
zAxis   = zAxis./norm(zAxis);
 
% Redefine the yaxis, so all axes are orthoganal
yAxis   = cross(zAxis,xAxis);
 
% Define transformation matrix
S       = [xAxis' yAxis' zAxis'];
 
 
function S = localSystemZ(P1,P2,P3)
% Defines a transformation matrix for a local coordinate system in an
% inertial frame of reference.  Uses P1 as the xaxis and P2 as the origin, and 
% P3 as the z-axis. Coordinates must be (1x3) vectors. Note: if theses axes 
% are not orthogonal, the z-axis direction is assumed to be more accurate
% than the x-axis and the x-axis direction is adjusted to make the coordinates 
% orthoganal.
 
% Check dimensions of inputs
if size(P1,1)~=1 || size(P1,2)~=3 ||...
   size(P2,1)~=1 || size(P2,2)~=3 ||...
   size(P3,1)~=1 || size(P3,2)~=3
    error('Coordinates must be 1x3 vectors');
end
 
% Define units vectors for x and y axes
xAxis   = (P1-P2)./norm(P1-P2);
zAxis   = (P3-P2)./norm(P3-P2);
 
% Define yaxis from the cross product
yAxis   = cross(zAxis,xAxis);
yAxis   = yAxis./norm(yAxis);
 
% Redefine the xaxis, so all axes are orthoganal
xAxis   = cross(yAxis,zAxis);
 
% Define transformation matrix
S       = [xAxis' yAxis' zAxis'];
 

function pts_out = localToGlobal(pts,origin,S)
% Transforms coordinates from the local coordinate system to the global
% system. pts must be a n x 3 matrix with x, y & z coordinates

pts     = [inv(S)'*pts']';

pts_out(:,1) = pts(:,1) + origin(1);
pts_out(:,2) = pts(:,2) + origin(2);
pts_out(:,3) = pts(:,3) + origin(3);

 
function pts_out = globalToLocal(pts,origin,S)
% Transforms coordinates from the global coordinate system to the local
% system. pts must be a n x 3 matrix with x,y & z coordinates

pts(:,1)    = pts(:,1)-origin(1);
pts(:,2)    = pts(:,2)-origin(2);
pts(:,3)    = pts(:,3)-origin(3);

pts_out     = [S'*pts']';


function [pXg,pYg,pZg] = prey_global(pX,pY,pZ,origin,S)
% Transforms prey points in global FOR

tmp1 = localToGlobal([pX(1,:)' pY(1,:)' pZ(1,:)'],origin,S);
tmp2 = localToGlobal([pX(2,:)' pY(2,:)' pZ(2,:)'],origin,S);
tmp3 = localToGlobal([pX(3,:)' pY(3,:)' pZ(3,:)'],origin,S);
tmp4 = localToGlobal([pX(4,:)' pY(4,:)' pZ(4,:)'],origin,S);

pXg = [tmp1(:,1)'; tmp2(:,1)'; tmp3(:,1)'; tmp4(:,1)'];
pYg = [tmp1(:,2)'; tmp2(:,2)'; tmp3(:,2)'; tmp4(:,2)'];
pZg = [tmp1(:,3)'; tmp2(:,3)'; tmp3(:,3)'; tmp4(:,3)'];


 

