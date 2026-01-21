function [fig, headerTable, signalTable] = header_gui(header, signal_header, varargin)
%HEADER_GUI  Display EDF header and signal header information in UI tables
%
%   Supports both UIFIGURE and UIPANEL parents.
%   Tables automatically resize if parent is resized.

%% ---------------- PARSE NAME-VALUE PAIRS ----------------
p = inputParser;
addParameter(p, 'Parent', [], @(x) isempty(x) || isgraphics(x));
parse(p, varargin{:});
parent = p.Results.Parent;

%% ---------------- VALIDATE INPUT STRUCTURES ----------------
if ~isstruct(header) || isempty(header)
    error('header_gui:InvalidInput', 'Header must be a non-empty structure.');
end
if ~isstruct(signal_header)
    error('header_gui:InvalidInput', 'Signal header must be a structure array.');
end

%% ---------------- HANDLE EMPTY SIGNAL HEADER ----------------
if isempty(signal_header)
    if isempty(parent)
        fig = uifigure('Name','EDF Header', 'Position',[100 100 800 200]);
        parent = fig;
    else
        fig = [];
    end
    uilabel(parent, ...
        'Text','No signal header data available.', ...
        'Position',[20 20 400 30]);
    headerTable = [];
    signalTable = [];
    return;
end

%% ---------------- CONVERT STRUCTS TO TABLES ----------------
H = struct2table(header, 'AsArray', true);
H.Properties.VariableNames = prettifyNames(H.Properties.VariableNames);

S = struct2table(signal_header);
S = addvars(S, (1:height(S))', 'Before', 1, 'NewVariableNames', 'Channel');
S.Properties.VariableNames = prettifyNames(S.Properties.VariableNames);

%% ---------------- CREATE UIFIGURE IF NEEDED ----------------
if isempty(parent)
    screenSize = get(0, 'ScreenSize');
    figWidth  = min(1500, screenSize(3));
    figHeight = min(700,  screenSize(4));
    x = (screenSize(3) - figWidth)/2;
    y = (screenSize(4) - figHeight)/2;
    fig = uifigure('Name','EDF Header Information', 'Position',[x y figWidth figHeight]);
    parent = fig;
else
    fig = [];
end

%% ---------------- POSITION TABLES IN PIXELS ----------------
originalUnits = parent.Units;
parent.Units = 'pixels';
parentPos = parent.Position;
parentWidth  = parentPos(3);
parentHeight = parentPos(4);

padding = 20;
rowHeight = 22;
headerHeight = rowHeight + 35;
signalHeight = parentHeight - headerHeight - 3*padding;

headerTable = uitable(parent, ...
    'Data', H, ...
    'Position', [padding, parentHeight - headerHeight - padding, ...
                 parentWidth - 2*padding, headerHeight], ...
    'ColumnWidth','auto');

signalTable = uitable(parent, ...
    'Data', S, ...
    'Position', [padding, padding, ...
                 parentWidth - 2*padding, signalHeight], ...
    'ColumnWidth','auto');

%% ---------------- CONVERT TO NORMALIZED UNITS ----------------
headerTable.Units = 'normalized';
signalTable.Units = 'normalized';

parent.Units = originalUnits; % restore parent units

end

%% =========================================================
function names = prettifyNames(names)
for k = 1:numel(names)
    name = strrep(names{k}, '_', ' ');
    parts = strsplit(lower(name), ' ');
    for p = 1:numel(parts)
        if ~isempty(parts{p})
            parts{p}(1) = upper(parts{p}(1));
        end
    end
    names{k} = strjoin(parts, ' ');
end
end
