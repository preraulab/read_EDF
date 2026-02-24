function varargout = read_EDF(edf_fname, varargin)
%READ_EDF  Load EDF or EDF+ file with full metadata, annotations, and MEX acceleration
%
%   READ_EDF reads European Data Format (EDF/EDF+) files using a compiled MEX
%   reader when available, and a pure MATLAB fallback otherwise. The function
%   provides full access to header metadata, per-signal_cellnal scaling (digital-to-
%   physical conversion), and EDF+ annotations.
%
%   Usage:
%       [header, signal_header, signal_cell, annotations] = read_EDF(filename, 'Channels', {'EEG Fpz-Cz'})
%
%   Inputs:
%       edf_fname      : string - EDF or EDF+ file path
%       'Channels'     : cell array - subset of channels to read (default: all)
%       'Epochs'       : 1x2 vector [start_epoch end_epoch] (0-indexed, default: all)
%       'Verbose'      : logical - print progress and status info (default: false)
%       'RepairHeader' : logical - correct invalid record counts and save with _fixed suffix (default: false)
%       'forceMATLAB'  : logical - disable MEX usage (default: false)
%       'debug'        : logical - debug mode for MEXt (default: false)
%
%   Outputs:
%       header   : structure containing EDF file-level metadata
%       signal_header    : structure array of per-signal_cellnal headers
%       signal_cell   : cell array containing each signal_cellnal vector (in physical units)
%       annotations   : structure array of EDF+ annotations with onset and text
%
%   Example:
%       [header, signal_header, signal_cell, annotations] = read_EDF('sleep.edf', 'Channels', {'EEG C3-A2'});
%
%   -------------------------------------------------------------------------
%   EDF File Specification Summary:
%       • Each EDF file begins with a fixed-length 256-byte main header
%       • Followed by per-signal_cellnal headers (16 fields × N signal_cellnals)
%       • Digital samples stored as int16 are scaled to physical units:
%
%             phys_val = phys_min + (dig_val - dig_min) * (phys_max - phys_min) / (dig_max - dig_min)
%
%       • EDF+ annotation channels (labelled 'EDF Annotations') contain onset
%         times and event texts in TAL (Time-Annotation List) format.

%% ---------------- INPUT PARSING ----------------
if nargin < 1
    error('read_EDF:InvalidInput', 'EDF filename required.');
end

edf_fname = char(edf_fname);
if ~isfile(edf_fname)
    error('read_EDF:FileNotFound', 'EDF file not found: %s', edf_fname);
end

p = inputParser;
addParameter(p, 'Channels', {}, @iscell);
addParameter(p, 'Epochs', [], @isnumeric);
addParameter(p, 'Verbose', false, @islogical);
addParameter(p, 'RepairHeader', false, @islogical);
addParameter(p, 'forceMATLAB', false, @islogical);
addParameter(p, 'debug', false, @islogical);
parse(p, varargin{:});

channels        = p.Results.Channels;
epochs          = p.Results.Epochs;
verbose         = p.Results.Verbose;
repair_header   = p.Results.RepairHeader;
force_matlab    = p.Results.forceMATLAB;
debug           = p.Results.debug;

%% ---------------- MEX HANDLING ----------------
script_dir = fileparts(mfilename('fullpath'));
mex_file = fullfile(script_dir, ['read_EDF_mex.' mexext]);
mex_exists = isfile(mex_file);

if ~force_matlab
    if mex_exists
        try
            [varargout{1:nargout}] = read_EDF_mex(edf_fname, channels, epochs, verbose, repair_header, debug);
            return
        catch ME
            if verbose
                fprintf('MEX failed (%s). Falling back to MATLAB reader.\n', ME.message);
            end
        end
    end
else
    if verbose
        fprintf('forceMATLAB = true, using MATLAB reader.\n');
    end
end

%% ---------------- MATLAB FALLBACK ----------------
[varargout{1:nargout}] = read_EDF_matlab(edf_fname, channels, epochs, verbose, repair_header);

end


%% =========================================================================
%  PURE MATLAB EDF READER
% =========================================================================
function varargout = read_EDF_matlab(edf_fname, channels, epochs, verbose, repair_header)

fid = fopen(edf_fname, 'r', 'ieee-le');
if fid < 0
    error('read_EDF:FileError', 'Cannot open file.');
end

%% ---------------- READ MAIN HEADER ----------------
A = fread(fid, 256, 'uint8=>char')';
field_sizes = [8,80,80,8,8,8,44,8,8,4];
fields = {'edf_ver','patient_id','local_rec_id',...
    'recording_startdate','recording_starttime',...
    'num_header_bytes','reserve_1','num_data_records',...
    'data_record_duration','num_signals'};

loc = [0; cumsum(field_sizes(:))];
header = struct();

for k = 1:numel(fields)
    str = strtrim(A(loc(k)+1:loc(k+1)));
    if any(strcmp(fields{k},{'num_header_bytes','num_data_records','data_record_duration','num_signals'}))
        header.(fields{k}) = str2double(str);
    else
        header.(fields{k}) = str;
    end
end

%% ---------------- READ SIGNAL HEADERS ----------------
fseek(fid, 256, 'bof');

num_signals = header.num_signals;
sig_sizes = [16,80,8,8,8,8,8,80,8,32];
sig_fields = {'signal_labels','transducer_type','physical_dimension',...
    'physical_min','physical_max','digital_min','digital_max',...
    'prefiltering','samples_in_record','reserve_2'};

A = fread(fid, header.num_header_bytes-256, 'uint8=>char')';
loc = [0; cumsum(sig_sizes(:)*num_signals)];

signal_header = struct();
for f = 1:numel(sig_fields)
    block = A(loc(f)+1:loc(f+1));
    for s = 1:num_signals
        idx1 = (s-1)*sig_sizes(f)+1;
        idx2 = s*sig_sizes(f);
        val = strtrim(block(idx1:idx2));
        if any(strcmp(sig_fields{f},{'physical_min','physical_max','digital_min','digital_max','samples_in_record'}))
            val = str2double(val);
        end
        signal_header(s).(sig_fields{f}) = val;
    end
end

%% ---------------- CHANNEL SELECTION ----------------
labels = strtrim({signal_header.signal_labels});

if isempty(channels)
    signal_indices = 1:num_signals;
else
    signal_indices = find(ismember(lower(labels), lower(strtrim(channels))));
    if isempty(signal_indices)
        error('No valid channels found.');
    end
end

%% ---------------- READ DATA ----------------
samples_per_record = [signal_header.samples_in_record];
record_size = sum(samples_per_record);

fseek(fid, header.num_header_bytes, 'bof');
raw = fread(fid, inf, 'int16=>double');
fclose(fid);

if isempty(raw)
    error('EDF contains no data samples.');
end

actual_records = floor(length(raw) / record_size);

if actual_records < 1
    error('File does not contain complete data records.');
end

if actual_records ~= header.num_data_records
    if verbose
        fprintf('Header record count mismatch. Using actual value: %d\n', actual_records);
    end

    if repair_header
        % Mirror MEX fix_num_records: write corrected count as left-justified
        % 8-char ASCII string at byte offset 236 (0-based), same as MEX does.
        fw = fopen(edf_fname, 'r+', 'ieee-le');
        if fw < 0
            warning('read_EDF:RepairFailed', ...
                'Could not open file for header repair: %s', edf_fname);
        else
            rec_str = sprintf('%-8d', actual_records);   % left-justified, space-padded to 8 chars
            fseek(fw, 236, 'bof');
            fwrite(fw, rec_str, 'char');
            fclose(fw);
            if verbose
                fprintf('Header repaired on disk: num_data_records written as %d.\n', actual_records);
            end
        end
    end

    header.num_data_records = actual_records;
end

%% ---------------- EPOCH VALIDATION ----------------
if isempty(epochs)
    start_epoch = 1;
    end_epoch = header.num_data_records;
else
    start_epoch = epochs(1) + 1;
    end_epoch   = epochs(2);
end

start_epoch = max(start_epoch,1);
end_epoch   = min(end_epoch, header.num_data_records);

if end_epoch < start_epoch
    error('Invalid epoch selection.');
end

num_epochs = end_epoch - start_epoch + 1;

%% ---------------- TOTAL DURATION ----------------
total_seconds = header.num_data_records * header.data_record_duration;
header.total_data_seconds = total_seconds;
header.total_data_hms = char(duration(0,0,total_seconds));

%% ---------------- SAMPLING FREQUENCY ----------------
for ii = 1:length(signal_header)
    signal_header(ii).sampling_frequency = ...
        signal_header(ii).samples_in_record / header.data_record_duration;
end

%% ---------------- DIGITAL TO PHYSICAL ----------------
signal_cells = cell(1,length(signal_indices));

for i = 1:length(signal_indices)

    sidx = signal_indices(i);
    offset = sum(samples_per_record(1:sidx-1));
    samples_per_epoch = samples_per_record(sidx);
    total_samples = samples_per_epoch * num_epochs;

    sig = zeros(1,total_samples);

    dig_min = signal_header(sidx).digital_min;
    dig_max = signal_header(sidx).digital_max;
    phy_min = signal_header(sidx).physical_min;
    phy_max = signal_header(sidx).physical_max;

    scale = (phy_max - phy_min) / (dig_max - dig_min);

    for r = 1:num_epochs
        rec = start_epoch + r - 1;
        idx1 = (rec-1)*record_size + offset + 1;
        idx2 = idx1 + samples_per_epoch - 1;

        raw_vals = raw(idx1:idx2);
        sig((r-1)*samples_per_epoch+1 : r*samples_per_epoch) = ...
            phy_min + (raw_vals - dig_min)*scale;
    end

    signal_cells{i} = sig;

end

signal_header = signal_header(signal_indices);

annotations = extractAnnotations(edf_fname, header, signal_header);

%% ---------------- OUTPUT ----------------
varargout{1} = header;
if nargout > 1, varargout{2} = signal_header; end
if nargout > 2, varargout{3} = signal_cells; end
if nargout > 3, varargout{4} = annotations; end

end

%% =========================================================================
%  EDF+ ANNOTATION PARSER — exclude blank annotation texts
% =========================================================================
function annotations = extractAnnotations(edf_fname, header, signal_header)
annotations = struct('onset', {}, 'text', {});
annotIdx = find(strcmp(strtrim({signal_header.signal_labels}), 'EDF Annotations'), 1);
if isempty(annotIdx), return; end

fid = fopen(edf_fname, 'rb');
if fid < 0, return; end

annotBytePos = sum([signal_header(1:annotIdx-1).samples_in_record]) * 2;
annotBytesPerRecord = signal_header(annotIdx).samples_in_record * 2;
totalBytesPerRecord = sum([signal_header.samples_in_record]) * 2;

fseek(fid, header.num_header_bytes, 'bof');
allAnnotations = [];
while ~feof(fid)
    data = fread(fid, totalBytesPerRecord, 'uint8=>uint8');
    if length(data) < totalBytesPerRecord, break; end
    annotData = data(annotBytePos + 1 : annotBytePos + annotBytesPerRecord);
    tals = parseTALs(annotData);
    allAnnotations = [allAnnotations; tals];
end
fclose(fid);

% Keep only annotations with non-empty text
if ~isempty(allAnnotations)
    mask = arrayfun(@(a) any(~cellfun(@isempty, strtrim(a.texts))), allAnnotations);
    valid = allAnnotations(mask);
    if ~isempty(valid)
        annotations = struct('onset', num2cell([valid.onset])', ...
            'text', {valid.texts}');
    end
end
end

%% =========================================================================
%  PARSE EDF+ TAL STRINGS
% =========================================================================
function tals = parseTALs(data)
tals = [];
idx = 1;
while idx <= length(data)
    while idx <= length(data) && data(idx) == 0
        idx = idx + 1;
    end
    if idx > length(data), break; end

    if data(idx) ~= 43 && data(idx) ~= 45  % '+' or '-'
        idx = idx + 1; continue;
    end
    onsetStart = idx;
    while idx <= length(data) && data(idx) ~= 20
        idx = idx + 1;
    end
    if idx > length(data), break; end
    onset = str2double(char(data(onsetStart:idx-1))');
    idx = idx + 1;

    texts = {};
    while idx <= length(data) && data(idx) ~= 0
        textStart = idx;
        while idx <= length(data) && data(idx) ~= 20 && data(idx) ~= 0
            idx = idx + 1;
        end
        texts{end+1} = strtrim(char(data(textStart:idx-1))');
        if idx <= length(data) && data(idx) == 20
            idx = idx + 1;
        end
    end
    if any(~cellfun(@isempty, strtrim(texts)))
        tals = [tals; struct('onset', onset, 'texts', {texts})];
    end
    idx = idx + 1;
end
end