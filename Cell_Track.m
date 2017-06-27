function varargout = Cell_Track(varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A graphical user interface for tracking bacterial motion from bright    %
% field microscopy time series (in TIFF format.) Interface allows loading %
% of a portion of the image stack, so memory requirements are minimal.    %
% Once loaded, the stack can be looped in the software, and parameters    %
% for the tracking can be set. The tracking is performed on the entire    %
% image sequence (not just the part being displayed) and the resultant    %
% trajectories are shown in the GUI. The data can then be exported to a   %
% text file for further processing.
% Images and processes it to determine the number of objects as well as 
% create a list of positions according to the centroid location.
% Last edited March 19th, 2012. 
% GUI implementation by Corey Kelly - coreyjkelly@gmail.com
% Tracking and image processing by Robin White - robint.white90@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Cell_Track_OpeningFcn, ...
                   'gui_OutputFcn',  @Cell_Track_OutputFcn, ...
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

% --- Executes just before Cell_Track is made visible.
function Cell_Track_OpeningFcn(hObject, eventdata, handles, varargin)

% Declare and set default values for global variables
addpath(pwd);
handles.output = hObject;
handles.ImageFile = '';
handles.tiffInfo = [];
handles.filepath = '' ;
handles.fileSize = 0;
handles.imagestack = [];
handles.currimage = [];
handles.dims = [];
handles.maxint = [];
handles.dt = 0;
handles.mag = 0;
handles.trajectories = [];
clear handles.rois;

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = Cell_Track_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

function Cell_Track_CloseRequestFcn(hObject, eventdata, handles)
delete(hObject);
close all;

function select_button_Callback(hObject, eventdata, handles)

if strcmp(handles.filepath,'')
    [ImageFile,FilePath] = uigetfile('*.tif','Select the image file');
else
    [ImageFile,FilePath] = uigetfile('*.tif',...
        'Select the image file',handles.filepath);
end
if ImageFile ~= 0
    handles.ImageFile = ImageFile;
    handles.filepath = FilePath;
    addpath(handles.filepath);
    
    tiffInfo = imfinfo(handles.ImageFile);  % TIFF file information
    handles.tiffInfo = tiffInfo;
    width  = tiffInfo.Width;
    height = tiffInfo.Height;
    n_frames = numel(tiffInfo);    % Number of frames in the file
    handles.dims = [height width n_frames];
    size = tiffInfo.FileSize;
    size = round(size/1048576);
        
    %Writes file path/name in static text box
    set(handles.filename_text,'String',[FilePath,ImageFile]);
    
    %Fill image dimensions text box
    set(handles.dim_text,'String',['Resolution: ',...
        num2str(width),'x', num2str(height),...
        '  Frames: ',num2str(n_frames),'  File Size: ',...
        num2str(size),' MB']);
    
    set(handles.openbutton,'Enable','on');
end

guidata(hObject,handles);

function openbutton_Callback(hObject, eventdata, handles)

n = get(handles.n_text,'Value');

switch get(handles.load_menu,'Value')
    case 1
        loadframes = 1:n;
    case 2
        loadframes = 1:n:handles.dims(3);
    case 3
        loadframes = 1:handles.dims(3);
end

wb = waitbar(0,'Loading image stack...');

% Pre-allocate array
handles.imagestack = zeros([handles.dims(1:2) length(loadframes)]);

for i = loadframes   % Load each frame, adding it to handles.imagestack()
    waitbar((i-1)/length(loadframes));
    handles.imagestack(:,:,find(loadframes==i))...
        = squeeze(imread(handles.ImageFile,'Index',i,'Info',handles.tiffInfo));
end

close(wb);

% Grabs dimensions.
handles.dims = size(handles.imagestack);
handles.maxint = max(max(max(handles.imagestack)));

%Resizes the sliders based on the properties of the image stack
slider_step(1) = 1/(handles.dims(3));
slider_step(2) = 1/(handles.dims(3));
set(handles.frame_slider1,...
    'sliderstep',slider_step,...
    'max',handles.dims(3),...
    'min',1,...
    'Value',1);
set(handles.frame_text,'String',['Frame: 1/'...
    num2str(handles.dims(3))]);
set(handles.time_text,'String','Time: 0s');


%displays the image, sets callback for clicking on image
handles.currimage = imshow(handles.imagestack(:,:,1));
axis image;
set(handles.axes1,'CLim',[0, handles.maxint]);
set(handles.currimage,'ButtonDownFcn',@axes1_ButtonDownFcn,...
    'HitTest','on');

%Prompts for timestep input
handles.dt = str2double(inputdlg('Input timestep, dt in seconds:',...
    'Timestep input',[1 35],{'0'}));
set(handles.dt_text,'String',num2str(handles.dt));

%Prompts for magnification input
[magchoice,~] = listdlg('ListString',{'150x (.0430 micron/pixel)'...
    '60x (.1075 micron/pixel)'},'SelectionMode','single','Name',...
    'Objective Selection','PromptString','Select an objective:',...
    'ListSize',[160 160]);
switch magchoice
    case 1
        handles.mag = .0430;
        set(handles.mag_text,'String','150x');
    case 2
        handles.mag = .1075;
        set(handles.mag_text,'String','60x');
end
guidata(hObject, handles);

function frame_slider1_Callback(hObject, eventdata, handles)
%Updates which image frame is displayed based on the slider position.
set(handles.currimage,'CData',...
    handles.imagestack(:,:,round(get(hObject,'Value'))));
set(handles.frame_text,'String',['Frame: '...
    num2str(round(get(hObject,'Value'))) '/' num2str(handles.dims(3))]);
set(handles.time_text,'String',['Time: '...
    num2str((round(get(hObject,'Value'))-1)*handles.dt) 's']);
guidata(hObject, handles);

function play_button_Callback(hObject, eventdata, handles)
%Loops through all frames of the time series.
i = round(get(handles.frame_slider1,'Value'));
while get(hObject,'Value')
    set(handles.currimage,'CData',handles.imagestack(:,:,i));
    % Set frame and relative time values
    set(handles.frame_text,'String',['Frame: ' num2str(i)...
        '/' num2str(handles.dims(3))]);
    set(handles.time_text,'String',['Time: '...
    num2str((i-1)*handles.dt) 's']);
    set(handles.frame_slider1,'Value',i);
    %Determines playback speed
    pause(1/get(handles.speed_slider,'Value'));
    if i == handles.dims(3)
        i = 1;
    else i = i+1;
    end
end

function track_button_Callback(hObject, eventdata, handles)
% pass variables (handles) to tracking code
handles.trajectories = image_analyse(handles);

% write number of unique trajectories to text box
set(handles.num_traj_text,'String',...
    ['Unique Trajectories Found: '...
    num2str(length(unique(handles.trajectories(:,4))))]);
guidata(hObject, handles);

function dt_text_Callback(hObject, eventdata, handles)
% Gets time step from text box and stores it in handles.dt
handles.dt = round(str2double(get(hObject,'String')));
guidata(hObject, handles);

function export_button_Callback(hObject, eventdata, handles)
% Exports data from table to a text file, with the same name as the image 
% file. 
exp_traj = handles.trajectories;
exp_traj(:,1:2) = exp_traj(:,1:2)*handles.mag;
fn = strrep(get(handles.filename_text,'String'),'.tif','.txt');
fmt = '%4d %4d %d %d\r\n';
fid = fopen(fn, 'w');
fprintf(fid,fmt,exp_traj');
fclose(fid);
msgbox(['Data saved to ' fn],'Save Successful!');

function reset_button_Callback(hObject, eventdata, handles)
if strcmp(questdlg('This will delete current data. Continue?',...
        'Reset?','Yes','No','No'),'Yes')
    closeGUI = handles.figure1;
    guiName = get(handles.figure1,'Name'); %get the name of the GUI
    close(closeGUI); %close the old GUI
    eval(guiName) %call the GUI again
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRACKING FUNCTION CODE BEGINS HERE %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function trajectories = image_analyse(handles)
%
%   This funtion imports a image and processes it to determine
%   the number of objects as well as create a list of positions according
%   to the centroid position
%   Last edit March 12, 2012 Robin White

tiffInfo = handles.tiffInfo;  % Get the TIFF file information
width  = tiffInfo.Width;
height = tiffInfo.Height;
n_frames = numel(tiffInfo);

a = zeros(height,width); % Pre-allocate array

t = handles.dt; %time step
%tmax = t*n_frames;
count = 0;

wb = waitbar(0,'Computing centroids...');
for j = 1:2:n_frames
    
    waitbar(j/n_frames);
    a = double(imread(handles.ImageFile,'Index',j,'Info',tiffInfo));
    max_a = max(max(max(a)));
    a = a/max_a; %creates intensity values between 0 and 1 for colormap
    
    %creates background neglecting objects radius smaller than 2nd input
    background = imopen(a, strel('disk', 4)); %4
    
    a2 = a - background; % Remove background
    
    level = graythresh(a2); %threshold map, creates binary image
    bw = im2bw(a2, level);
    
    %remove background noise. removes from a binary image all connected
    %components that have fewer than P pixels (2nd input)
    bw = bwareaopen(bw, 40); %40
    
    cc = bwconncomp(bw, 8); %identify objects in image
    
    %**************************************************************
    %   Remove unwanted objects that may not be bacteria
    %**************************************************************
    stats = regionprops(cc, 'MajorAxisLength', 'MinorAxisLength');
    eccen = [stats.MajorAxisLength]./[stats.MinorAxisLength];
    idx = find(eccen > 2.0 & eccen < 7.5); %2.5
    bw2 = ismember(labelmatrix(cc), idx);
    
    cc = bwconncomp(bw2, 8); %identify objects in new image
    
    stats = regionprops(cc, 'Area','Perimeter');
    ratio = [stats.Area]./[stats.Perimeter];
    idx = find(ratio > 1.0);
    bw2 = ismember(labelmatrix(cc), idx);
%     
%     %hist(ratio)
%     
    cc = bwconncomp(bw2, 8); %obtain objects from edited image
%     
    %gives each object area, centroid, boundingbox
    graindata = regionprops(cc, 'basic');
    
    %vector which holds the area measurement for each grain
    % grain_areas = [graindata.Area];
    
    % Can play around with statistics min area etc.
    
    %   [min_area, idx] = min(grain_areas);
    %   grain = false(size(bw));
    %   grain(cc.PixelIdxList{idx}) = true;
    %   imshow(grain);
    
    % get centroids and reshape them into a list with two columns
    
    cents = [graindata.Centroid];
    n_cents = length(cents)/2;
    grain_cents = reshape(cents,2,n_cents)';
    
    % make a vector of time values the same length as the number of
    % centroids
    times = ones(n_cents,1)*t*(j-1);
    
    % join the centroids and the times into a list with three columns, x
    % centroid, y centroid, and time
    grain_cents = cat(2,grain_cents,times);
    
    % the first time through, estimate the size of pk. Assume that the
    % total number of centroids will be the number in this frame multiplied
    % by the number of frames
    
    if j==1
        pk = zeros(n_cents*n_frames,3);
    end
    
    % Put the current list of centroids and times into pk. Count keeps
    % track of how many centroids are in pk so that none are overwritten.
    % Even if there end up being more centroids than the above estimate,
    % they'll still be added to the list.
    
    pk(count+1:count+n_cents,:) = grain_cents;
    count = count+n_cents;
    
end

close(wb);

% since allocation was an estimate, there may be rows of zeros left
% this line removes them

pk(all(pk==0,2),:)=[];

%Parameters for track function
param.mem = 5;
param.good = 10; %20
param.quiet = 0;
param.dim = 2;

trajectories = track(pk,10,param); %10 This line will affect the error difficult combinatrix encountered

num_tracked = max(trajectories(:,4));

vel = zeros(length(trajectories) - num_tracked,1);

traj = zeros(num_tracked,4);
count = 0;
hold all

for i = 1:num_tracked
    traj = trajectories(trajectories(:,4) == i,:);
    dx = diff(traj(:,1));
    dy = diff(traj(:,2));
    dr = handles.mag.*sqrt(dx.^2 + dy.^2);
    v = dr./handles.dt;
    vel(count+1:count+size(traj,1)-1) = v;
    count = count+size(traj,1)-1;
    plot(traj(:,1),traj(:,2));
end
hold off

vel = vel(isfinite(vel)); %returns only finite velocities
figure()
axis()
hold on;
hist(vel,10000)
hold off;
avg_v = mean(vel)
std(vel) %standard deviation
cov(vel) %variance
mad(vel) %mean absolute deviation