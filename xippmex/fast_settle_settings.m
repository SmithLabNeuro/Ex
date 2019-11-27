function varargout = fast_settle_settings(varargin)
% FAST_SETTLE_SETTINGS MATLAB code for fast_settle_settings.fig
%      FAST_SETTLE_SETTINGS, by itself, creates a new FAST_SETTLE_SETTINGS or raises the existing
%      singleton*.
%
%      H = FAST_SETTLE_SETTINGS returns the handle to a new FAST_SETTLE_SETTINGS or the handle to
%      the existing singleton*.
%
%      FAST_SETTLE_SETTINGS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FAST_SETTLE_SETTINGS.M with the given input arguments.
%
%      FAST_SETTLE_SETTINGS('Property','Value',...) creates a new FAST_SETTLE_SETTINGS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fast_settle_settings_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fast_settle_settings_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help fast_settle_settings

% Last Modified by GUIDE v2.5 14-Feb-2013 17:40:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fast_settle_settings_OpeningFcn, ...
                   'gui_OutputFcn',  @fast_settle_settings_OutputFcn, ...
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

% --- Executes just before fast_settle_settings is made visible.
function fast_settle_settings_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to fast_settle_settings (see VARARGIN)

% Choose default command line output for fast_settle_settings
handles.output = hObject;

stim_main_input = find(strcmp(varargin, 'stimui'));

if (isempty(stim_main_input) ...
    || (length(varargin) <= stim_main_input) ...
    || (~ishandle(varargin{stim_main_input+1})))
    dontOpen = true;
else
    % Remember the handle, and adjust our position
    handles.stim_main = varargin{stim_main_input+1};
    % Obtain handles using GUIDATA with the caller's handle 
    main_handles = guidata(handles.stim_main);
    set(handles.found_list, 'String', ...
        get(main_handles.found_list, 'String'));
    labels = [];
    for iLabel=1:length(main_handles.fs_elec)
        labels = [labels; sprintf('%03d', main_handles.fs_elec(iLabel))];
    end
    set(handles.fs_list, 'String', labels);
end
% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = fast_settle_settings_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on selection change in found_list.
function found_list_Callback(hObject, eventdata, handles)
% hObject    handle to found_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes during object creation, after setting all properties.
function found_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to found_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), ...
%         get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end

function fs_list_CreateFcn(hObject, eventdata, handles)


% --- Executes on selection change in fs_list.
function fs_list_Callback(hObject, eventdata, handles)
% hObject    handle to fs_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in add_button.
function add_button_Callback(hObject, eventdata, handles)
% hObject    handle to add_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value_list = get(handles.found_list, 'Value');
found_labels = get(handles.found_list, 'String');
saved_labels = get(handles.fs_list, 'String');

for iValue=1:length(value_list)
    label = found_labels(value_list(iValue));
    if ~isempty(find(strcmp(label, saved_labels)))
        continue
    end
    saved_labels = [saved_labels; label];
end
set(handles.fs_list, 'String', saved_labels);
% Update handles structure
rt_electrodes(hObject, eventdata, handles);

guidata(hObject, handles);

function rt_electrodes(hObject, eventdata, handles)
saved_labels = get(handles.fs_list, 'String');
main_handles = guidata(handles.stim_main);
main_handles.fs_elec = zeros(1, length(saved_labels));

for iLabel=1:length(saved_labels)
    main_handles.fs_elec(iLabel) = str2double(saved_labels(iLabel));
end
guidata(handles.stim_main, main_handles);

% --- Executes on button press in remove_button.
function remove_button_Callback(hObject, eventdata, handles)
% hObject    handle to remove_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
value_list = get(handles.fs_list, 'Value');
found_labels = get(handles.fs_list, 'String');

for value_index=1:length(value_list)
    label = found_labels(value_index);
    found_labels = found_labels(~strcmp(label, found_labels));
end
set(handles.fs_list, 'String', found_labels);

set(handles.fs_list, 'Value', 1);

rt_electrodes(hObject, eventdata, handles);

% Update handles structure
guidata(hObject, handles);
