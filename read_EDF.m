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
%       'RepairHeader' : logical - correct invalid record counts (default: false)
%       'forceMATLAB'  : logical - disable MEX usage (default: false)
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
parse(p, varargin{:});
channels = p.Results.Channels;
epochs = p.Results.Epochs;
verbose = p.Results.Verbose;
repair_header = p.Results.RepairHeader;
force_matlab = p.Results.forceMATLAB;

%% ---------------- MEX HANDLING ----------------
script_dir = fileparts(mfilename('fullpath'));
mex_file = fullfile(script_dir, ['read_EDF_mex.' mexext]);
mex_exists = isfile(mex_file);

if ~force_matlab
    if ~mex_exists
        try
            if verbose, fprintf('Compiling read_EDF_mex.c...\n'); end
            cd(script_dir);
            mex('read_EDF_mex.c');
            cd(script_dir);
            mex_exists = isfile(mex_file);
        catch
            if verbose, fprintf('MEX compilation failed. Using MATLAB fallback.\n'); end
            mex_exists = false;
        end
    end
    if mex_exists
        try
            [varargout{1:nargout}] = read_EDF_mex(edf_fname, channels, epochs, verbose, repair_header);
            % ---------------- Add total data length fields ----------------
            total_seconds = varargout{1}.num_data_records * varargout{1}.data_record_duration;
            varargout{1}.total_data_seconds = total_seconds;
            dur = seconds(total_seconds);
            dur.Format="hh:mm:ss.SSS";
            varargout{1}.total_data_hms = char(dur);


            if verbose
                fprintf('Total data length: %.2f sec (%s)\n', varargout{1}.total_data_seconds, varargout{1}.total_data_hms);
            end

            return;
        catch ME
            if verbose
                fprintf('MEX failed (%s). Falling back to MATLAB reader.\n', ME.message);
            end
        end
    end
else
    if verbose, fprintf('forceMATLAB = true, skipping MEX reader\n'); end
end

%% ---------------- MATLAB FALLBACK ----------------
[varargout{1:nargout}] = read_EDF_matlab(edf_fname, channels, epochs, verbose, repair_header);

% ---------------- Add total data length fields ----------------
total_seconds = varargout{1}.num_data_records * varargout{1}.data_record_duration;
varargout{1}.total_data_seconds = total_seconds;
dur = seconds(total_seconds);
dur.Format="hh:mm:ss.SSS";
varargout{1}.total_data_hms = char(dur);


if verbose
    fprintf('Total data length: %.2f sec (%s)\n', varargout{1}.total_data_seconds, varargout{1}.total_data_hms);
end

end


%% =========================================================================
%  PURE MATLAB EDF READER
% =========================================================================
function varargout = read_EDF_matlab(edf_fname, channels, epochs, verbose, repair_header)
fid = fopen(edf_fname, 'r', 'ieee-le');
if fid < 0
    error('read_EDF:FileError', 'Cannot open file: %s', edf_fname);
end

if verbose, fprintf('Reading EDF header...\n'); end

% ---------------- Read EDF Main Header ----------------
A = fread(fid, 256, 'uint8=>char')';
header_fields = {'edf_ver','patient_id','local_rec_id',...
    'recording_startdate','recording_starttime',...
    'num_header_bytes','reserve_1','num_data_records',...
    'data_record_duration','num_signals'};
field_sizes = [8,80,80,8,8,8,44,8,8,4];
loc = [0; cumsum(field_sizes(:))];
header = struct();
for k = 1:numel(header_fields)
    str = strtrim(A(loc(k)+1:loc(k+1)));
    switch header_fields{k}
        case {'num_header_bytes','num_data_records','data_record_duration','num_signals'}
            header.(header_fields{k}) = str2double(str);
        otherwise
            header.(header_fields{k}) = str;
    end
end

% Check for data record mismatches
if header.num_data_records == -1 || isnan(header.num_data_records)
    warning('read_EDF:InvalidDataRecords', 'num_data_records is -1 or invalid.');
end

% ---------------- Read Per-Signal Headers ----------------
fseek(fid, 256, 'bof');
num_signals = header.num_signals;
signal_cell_header_size = header.num_header_bytes - 256;
A = fread(fid, signal_cell_header_size, 'uint8=>char')';
signal_cell_fields = {'signal_labels','transducer_type','physical_dimension',...
    'physical_min','physical_max','digital_min','digital_max',...
    'prefiltering','samples_in_record','reserve_2'};
signal_cell_field_sizes = [16,80,8,8,8,8,8,80,8,32];
loc = [0; cumsum(signal_cell_field_sizes(:)*num_signals)];
signal_header = struct();
for f = 1:numel(signal_cell_fields)
    block = A(loc(f)+1:loc(f+1));
    for s = 1:num_signals
        start_idx = (s-1)*signal_cell_field_sizes(f)+1;
        end_idx = s*signal_cell_field_sizes(f);
        val = strtrim(block(start_idx:end_idx));
        if any(strcmp(signal_cell_fields{f},{'physical_min','physical_max','digital_min','digital_max','samples_in_record'}))
            val = str2double(val);
        end
        signal_header(s).(signal_cell_fields{f}) = val;
    end
end

% ---------------- Validate Channel List ----------------
labels = cellfun(@deblank,{signal_header.signal_labels},'UniformOutput',false);
if isempty(channels)
    signal_indices = 1:num_signals;
else
    channels_trimmed = cellfun(@strtrim, channels, 'UniformOutput', false);
    signal_indices = find(ismember(lower(labels), lower(channels_trimmed)));
    invalid_channels = setdiff(channels_trimmed, labels);
    if ~isempty(invalid_channels)
        invalid_str = sprintf('%s, ', invalid_channels{:});
        valid_str = sprintf('%s, ', labels{:});
        warning('read_EDF:InvalidChannel', ...
            'Some channels not found: %s\nValid channels: %s', ...
            invalid_str(1:end-2), valid_str(1:end-2));
    end
    if isempty(signal_indices)
        valid_str = sprintf('%s, ', labels{:});
        error('read_EDF:NoValidChannels', ...
            'No valid channels found.\nAvailable: %s', ...
            valid_str(1:end-2));
    end
end

% ---------------- Prepare Data Read ----------------
samples_per_record = [signal_header.samples_in_record];
record_size = sum(samples_per_record);
fseek(fid, header.num_header_bytes, 'bof');
raw = fread(fid, inf, 'int16=>double', 'ieee-le');
fclose(fid);

num_records = header.num_data_records;
if isempty(epochs)
    start_epoch = 1;
    end_epoch = num_records;
else
    start_epoch = epochs(1) + 1;
    end_epoch = epochs(2);
end
num_epochs = end_epoch - start_epoch + 1;

% Update header with actual record count and duration
actual_records = length(raw) / record_size;
if actual_records ~= num_records
    warning('read_EDF:HeaderUpdate', ...
        'Updating header: num_data_records changed from %d to %d', ...
        num_records, actual_records);
    header.num_data_records = actual_records;
end
total_seconds = header.num_data_records * header.data_record_duration;
header.total_data_seconds = total_seconds;
dur = seconds(total_seconds);
dur.Format="hh:mm:ss.SSS";
header.total_data_hms = char(dur);
if verbose
    fprintf('Header updated: Total data length: %.2f sec (%s)\n', header.total_data_seconds, header.total_data_hms);
end

% ---------------- Repair Header if Requested ----------------
record_size = sum(samples_per_record);  % samples per record
actual_records = length(raw) / record_size;

if repair_header
    if verbose
        fprintf('RepairHeader = true: correcting num_data_records from %d to %d\n', ...
            header.num_data_records, actual_records);
    end
    header.num_data_records = actual_records;

    % Optional: save repaired EDF copy
    [pathstr, name, ext] = fileparts(edf_fname);
    fixed_fname = fullfile(pathstr, [name '_fixed' ext]);
    fid_orig = fopen(edf_fname, 'r', 'ieee-le');
    fid_fixed = fopen(fixed_fname, 'w', 'ieee-le');
    if fid_orig > 0 && fid_fixed > 0
        hdr = fread(fid_orig, header.num_header_bytes, 'uint8=>uint8');
        % Overwrite num_data_records field (bytes 237-244 in header, zero-based index)
        num_records_str = sprintf('%-8d', actual_records);
        hdr(237:244) = uint8(num_records_str);
        fwrite(fid_fixed, hdr, 'uint8');
        % Copy remaining file data
        fseek(fid_orig, header.num_header_bytes, 'bof');
        data_bytes = fread(fid_orig, inf, 'uint8=>uint8');
        fwrite(fid_fixed, data_bytes, 'uint8');
        fclose(fid_orig); fclose(fid_fixed);
        if verbose
            fprintf('Repaired EDF saved as: %s\n', fixed_fname);
        end
    end
end


% ---------------- EDF Digital-to-Physical Conversion ----------------
signal_cells = cell(1, length(signal_indices));
for i = 1:length(signal_indices)
    sidx = signal_indices(i);
    signal_cell_offset = sum(samples_per_record(1:sidx-1));
    num_samples = samples_per_record(sidx) * num_epochs;
    signal_cell = zeros(1, num_samples, 'double');

    dig_min = double(signal_header(sidx).digital_min);
    dig_max = double(signal_header(sidx).digital_max);
    phy_min = double(signal_header(sidx).physical_min);
    phy_max = double(signal_header(sidx).physical_max);

    if dig_max < dig_min, [dig_min, dig_max] = deal(dig_max, dig_min); end
    if phy_max < phy_min, [phy_min, phy_max] = deal(phy_max, phy_min); end
    if dig_max == dig_min
        scale = 0;
    else
        scale = (phy_max - phy_min) / (dig_max - dig_min);
    end

    if verbose & ~isempty(signal_header), fprintf('     Reading %s...\n',signal_header(sidx).signal_labels); end

    for r = 1:num_epochs
        record_idx = start_epoch + r - 1;
        rec_start = (record_idx - 1) * record_size + signal_cell_offset + 1;
        rec_end = rec_start + samples_per_record(sidx) - 1;
        raw_vals = double(raw(rec_start:rec_end));
        signal_cell_start = (r - 1) * samples_per_record(sidx) + 1;
        signal_cell_end = r * samples_per_record(sidx);
        signal_cell(signal_cell_start:signal_cell_end) = phy_min + (raw_vals - dig_min) * scale;
    end
    signal_cells{i} = signal_cell;
end

signal_header = signal_header(signal_indices);
annotations = extractAnnotations(edf_fname, header, signal_header);

% ---------------- Output Assignal_cellnment ----------------
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