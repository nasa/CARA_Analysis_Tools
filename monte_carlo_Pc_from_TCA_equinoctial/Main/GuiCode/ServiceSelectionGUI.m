function [varargout] = ServiceSelectionGUI(varargin)
% SERVICESELECTIONGUI MATLAB code for ServiceSelectionGUI.fig
%      SERVICESELECTIONGUI, by itself, creates a new SERVICESELECTIONGUI or raises the existing
%      singleton*.
%
%      H = SERVICESELECTIONGUI returns the handle to a new SERVICESELECTIONGUI or the handle to
%      the existing singleton*.
%
%      SERVICESELECTIONGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SERVICESELECTIONGUI.M with the given input arguments.
%
%      SERVICESELECTIONGUI('Property','Value',...) creates a new SERVICESELECTIONGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ServiceSelectionGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ServiceSelectionGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ServiceSelectionGUI

% Last Modified by GUIDE v2.5 06-Jun-2018 08:06:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ServiceSelectionGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @ServiceSelectionGUI_OutputFcn, ...
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


% --- Executes just before ServiceSelectionGUI is made visible.
function ServiceSelectionGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ServiceSelectionGUI (see VARARGIN)

% Initialize Data Outpus
handles.t.Data = varargin{1};

% Choose default command line output for ServiceSelectionGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ServiceSelectionGUI wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function [varargout] = ServiceSelectionGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.t.Data;

% Delete Figure Window
delete(handles.figure1);




% --- Executes on button press in ProceedButton.
function [varargout] = ProceedButton_Callback(hObject, eventdata, handles)
% hObject    handle to ProceedButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Set Outputs
varargout{1} = handles.t.Data;

% Delete Figure Window
figure1_CloseRequestFcn(handles.figure1, eventdata, handles);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
if isequal(get(hObject,'waitstatus'),'waiting')
    % Wait for user input
    uiresume(hObject);
else
    % Close figure
    delete(hObject);
end
