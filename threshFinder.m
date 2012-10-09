function varargout = threshFinder(varargin)
% THRESHFINDER M-file for threshFinder.fig
%      THRESHFINDER, is a gui that interactively selects a threhols level
%      for a tiff image
%
%      H = THRESHFINDER returns the handle to a new THRESHFINDER or the handle to
%      the existing singleton*.
%
%      THRESHFINDER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in THRESHFINDER.M with the given input arguments.
%
%      THRESHFINDER('Property','Value',...) creates a new THRESHFINDER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before threshFinder_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to threshFinder_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help threshFinder

% Last Modified by GUIDE v2.5 30-Oct-2006 13:19:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @threshFinder_OpeningFcn, ...
                   'gui_OutputFcn',  @threshFinder_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before threshFinder is made visible.
function threshFinder_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
handles.imData = varargin{1};

d = varargin{2};

if nargin < 4
    error('threshFinder needs 4 inputs')
end

if strcmp(varargin{3},'top')
    
    handles.tVal = d.T.tVal(d.curr_frame_idx);
    handles.view = 'top';
    
elseif strcmp(varargin{3},'side')
    
    handles.tVal = d.S.tVal(d.curr_frame_idx);
    handles.view = 'side';
end

% Store d in handles structure
handles.d = d;

% Store sequence path
handles.seq_path = varargin{4};

% Update gui
set(handles.edit1,'String',num2str(handles.tVal));
set(handles.slider1,'Value',handles.tVal);

makeImage(handles)

% Update handles structure
guidata(hObject, handles);



% --- Outputs from this function are returned to the command line.
function varargout = threshFinder_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
%varargout{1} = handles.p.tVal;
%varargout{2} = str2num(get(handles.edit1,'String'));


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
tVal    = get(handles.slider1,'Value');
        set(handles.edit1,'String',num2str(tVal));
makeImage(handles)

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit1_Callback(hObject, eventdata, handles)
tVal 	= str2num(get(handles.edit1,'String'));
set(handles.slider1,'Value',tVal);
makeImage(handles)

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)

tVal = get(handles.slider1,'Value');
close

% Retrieve 'd'
d = handles.d;

% Define idx
idx = d.curr_frame_idx;

% Frames after current (logical)
mod_remain = [[1:length(d.frames)] > idx]';

if strcmp(handles.view,'top')
    
    % Find next modified tVal
    mod_next = min([find(mod_remain & d.T.mod_tVal,1,'first') ...
                length(d.frames)]);
    
    % Ajust tVal for current and future frames
    d.T.tVal(idx:mod_next) = tVal.*ones(length(idx:mod_next),1);
    
    % Record tVal mod
    d.T.mod_tVal(idx) = 1; 
   
elseif strcmp(handles.view,'side')
       
    % Find next modified tVal
    mod_next = min([find(mod_remain & d.S.mod_tVal,1,'first') ...
                length(d.frames)]);
    
    % Ajust tVal for current and future frames
    d.S.tVal(idx:mod_next) = tVal.*ones(length(idx:mod_next),1);
    
    % Record tVal mod
    d.S.mod_tVal(idx) = 1;

end

% Save updated d
save([handles.seq_path filesep 'coord_data.mat'],'d')


% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function makeImage(handles)
warning off all
tVal    = get(handles.slider1,'Value');
im  = handles.imData;
cmap = [linspace(0,1,256)' linspace(0,1,256)' linspace(0,1,256)'];
imshow(im,cmap);
hold on;
%fillColor           = [.43 .49 1];
fillColor           = [1 0 0];
imBW                = ~im2bw(im,cmap,tVal);
im(find(imBW))    = 244.*ones(length((find(imBW))),1);
%cMap                = [[0:1./255:1]' [0:1./255:1]' [0:1./255:1]'];
cmap(245,:)         = fillColor;
imshow(im,cmap)
hold off;
warning on all

function yes = isfile(fname,dirname)
a = dir(dirname);
yes = 0;
for i = 3:length(a)
    if strcmp(a(i).name,fname)
        yes = 1;
        return
    end
end
