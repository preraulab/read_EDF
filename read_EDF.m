function varargout = read_EDF(edf_fname, varargin)
%READ_EDF  Load EDF / EDF+ / EDF.gz / EDF.zst file with full metadata, annotations, and MEX acceleration
%
%   READ_EDF reads European Data Format (EDF / EDF+) files using a compiled
%   MEX reader when available, and a pure MATLAB fallback otherwise. Plain
%   '.edf', gzip-compressed '.edf.gz', and zstd-compressed '.edf.zst' files
%   are all accepted; compressed files are streamed directly through zlib
%   or libzstd with no temp file (MEX path). The function provides full
%   access to header metadata, per-signal scaling (digital-to-physical
%   conversion), and EDF+ annotations.
%
%   Usage:
%       [header, signal_header, signal_cell, annotations] = ...
%           read_EDF(filename, 'Channels', {'EEG Fpz-Cz'})
%
%   Inputs:
%       edf_fname      : string - path to a .edf, .edf.gz, or .edf.zst file
%       'Channels'     : cell array of strings - which channels to emit
%                        (default {} = all). Each entry is parsed as
%                        either a passthrough label or an expression
%                        (see "Re-referencing on load" below).
%       'References'   : cell array of 'NAME = expr' strings (default {}).
%                        Defines named mean / linear-combination signals
%                        once, then makes the name usable inside any
%                        later 'References' or 'Channels' spec. Evaluated
%                        in declared order, so a later reference can
%                        build on an earlier one. References that are
%                        not also listed in 'Channels' are dropped from
%                        the returned outputs (they exist purely to
%                        support derivations).
%       'Epochs'       : 1x2 vector [start_epoch end_epoch] (0-indexed,
%                        default: all)
%       'Verbose'      : logical - print progress / status (default false)
%       'RepairHeader' : logical - correct an invalid num_data_records
%                        in the file header and save a '_fixed' copy
%                        (default false; ignored with a warning for
%                        compressed inputs).
%       'forceMATLAB'  : logical - disable the MEX backend and use the
%                        pure-MATLAB reader (default false). Compressed
%                        inputs are decompressed to a temp file in this
%                        path and auto-cleaned-up on return.
%       'debug'        : logical - debug mode for MEX (default false)
%       'deidentify'   : logical - blank PHI fields and save a
%                        '_deidentified' copy (default false; ignored
%                        with a warning for compressed inputs).
%
%   Outputs:
%       header        : struct of EDF file-level metadata
%       signal_header : struct array of per-signal headers
%       signal_cell   : cell array, one signal vector per channel,
%                       in physical units (scaled and offset)
%       annotations   : struct array of EDF+ annotations
%
%   -------------------------------------------------------------------------
%   Re-referencing on load
%
%   Each 'Channels' entry is preprocessed (greedy longest-match wraps
%   real labels in '$...$' tokens) and parsed as an expression in this
%   grammar:
%
%       spec   := [name '='] expr
%       expr   := [sign] term (sign term)*
%       term   := factor (('*'|'/') factor)*
%       factor := ['+'|'-'] (number | '$' label '$' | '(' expr ')'
%                            | mean '(' expr (',' expr)+ ')')
%       sign   := '+' | '-'
%
%   The output of parsing is always a flat linear combination
%   sum_j coef_j * signal[leaf_j] — number literals fold into the
%   leaf coefficients at parse time, so '(1/3)*C1 - 4*(C2-C3)/7'
%   becomes leaves={C1,C2,C3}, terms=[1/3, -4/7, 4/7].
%
%   Linearity is enforced: a single term may contain at most one
%   signal-valued factor. The following are rejected with
%   read_EDF:ParseError:
%       signal + scalar      (DC offsets not supported)
%       signal * signal      (not linear)
%       signal / signal      (not linear)
%       scalar / signal      (not linear)
%
%   The optional 'name =' prefix sets the output channel's label
%   (alias). Without it the output label is the raw spec string.
%   Aliases are mandatory for References. Aliased Channels entries
%   are visible by name to subsequent Channels entries (chaining);
%   unaliased entries are one-shot.
%
%   mean(arg1, ..., argN) requires N >= 2 and produces (1/N) * sum_i argi.
%   Args may be expressions, not just leaves: mean(C1, C2-C3) works.
%
%   Examples (each comment is the equivalent math):
%
%       % Inline linked-mastoid                  C3 - (A1+A2)/2
%       'C3 - mean(A1, A2)'
%
%       % Same, but bind a reusable name first
%       'References', {'LM = mean(A1, A2)'}, ...
%       'Channels',   {'C3-LM', 'C4-LM', 'Fpz_LM = -LM'}
%
%       % Sum of two channels                    A1 + A2
%       'A1 + A2'
%
%       % Reference building on a reference
%       'References', {'M = mean(A1,A2)', 'R = C3 - M'}, ...
%       'Channels',   {'R'}                            % emits C3 - M
%
%       % Arbitrary linear combination
%       '(1/3)*C1 - 4*(C2 - C3)/7'                     % weighted
%       '0.7*C3 + 0.3*C4'                              % weighted average
%       'mean(C1, 2*C2 + C3)'                          % expr inside mean
%
%   -------------------------------------------------------------------------
%   Resolution: how a spec becomes a signal
%
%   Every 'Channels' and 'References' spec is preprocessed before
%   parsing. A single greedy longest-match pass walks the augmented
%   label set (real EDF labels + earlier References + earlier aliased
%   Channels outputs) and wraps each occurrence in '$...$' tokens. The
%   parser then sees only '$LABEL$' tokens, operators ('+', '-', ',',
%   '(', ')', '='), and 'mean'.
%
%       EDF labels: {'EEG C3 - A2', 'EEG C3', 'A2', 'EMG'}
%
%       'EEG C3 - A2'         ->  '$EEG C3 - A2$'
%                                 (passthrough; longest match wraps the
%                                  whole spec)
%       'EEG C3 - A2 - A2'    ->  '$EEG C3 - A2$ - $A2$'
%                                 ((labeled C3-A2) minus A2; longest-
%                                  first then leftover)
%       'mean(EEG C3, A2)'    ->  'mean($EEG C3$, $A2$)'
%       'OUT = EEG C3 - A2'   ->  'OUT = $EEG C3 - A2$'
%                                 (alias name is left alone)
%
%   User-supplied '$...$' regions are preserved verbatim by the
%   preprocessor, so write '$LABEL$' explicitly to override the
%   longest-match interpretation:
%
%       % EDF has both 'EEG C3 - A2' (a recording-side reref) AND
%       % 'EEG C3' / 'A2' as separate channels.
%       read_EDF(f, 'Channels', {'EEG C3 - A2'})        % the labeled channel
%       read_EDF(f, 'Channels', {'$EEG C3$ - $A2$'})    % computed difference
%
%   Pass 'Verbose', true to see the wrapped form for each spec — useful
%   when you suspect a label-vs-expression collision.
%
%   '$LABEL$' is also the escape for labels containing characters the
%   parser would otherwise eat ('+', ',', '(', ')', '='), labels that
%   start with 'mean(', or labels that are substrings of others when you
%   want to force the shorter one:
%
%       read_EDF(f, 'Channels', {'$EEG A+B$ - mean($A1$, $A2$)'})
%       read_EDF(f, 'Channels', {'$mean(LF,RF)$'})           % literal label
%       read_EDF(f, 'Channels', {'$Label = X$'})             % '=' is literal
%       read_EDF(f, 'Channels', {'C3 - mean($A1$, $A2$)'})   % force short A1
%
%   -------------------------------------------------------------------------
%   Channels chaining
%
%   Aliased Channels entries become reusable names for any later
%   'Channels' entry, so multi-step derivations can be written inline:
%
%       'Channels', {'LM      = mean(A1, A2)', ...
%                    'C3_LM   = C3 - LM', ...
%                    'NewChan = C3_LM + 0.5'}
%
%   The difference vs 'References': aliased Channels entries are
%   *returned* as outputs; References outputs are hidden helpers (drop
%   out of the returned cell unless explicitly named in 'Channels').
%   Unaliased Channels entries (no '=') don't pollute the namespace.
%
%   -------------------------------------------------------------------------
%   Constraints. A per-entry resolution failure on a 'Channels' spec is
%   demoted to a warning (warn + skip that entry, keep loading the
%   rest); the same failure inside a 'References' entry is still a
%   hard error, since References are foundation pieces other channels
%   may depend on.
%
%   Demoted to warnings on a 'Channels' entry (errors elsewhere):
%       read_EDF:UnknownChannel  - a leaf in an expression doesn't
%                                  resolve to a label in the augmented
%                                  set. Also emitted by the legacy
%                                  backend when a plain label is not
%                                  found in the file.
%       read_EDF:ParseError      - malformed expression syntax (missing
%                                  '=' in a Reference, unterminated
%                                  '$' quote, dangling operator, ...).
%                                  Also fires when a plain label like
%                                  'C3-A2' is requested but the EDF
%                                  has only '[C3-A2 - B]' — nothing
%                                  matches, so the parser sees raw
%                                  'C3-A2' and trips at the first 'C'.
%       read_EDF:RefCollision    - an alias 'OUT = expr' tries to
%                                  redefine a name that already exists
%                                  (case-insensitive). Demoted so
%                                  fallback patterns like
%                                  {'A', 'A=B'} or {'A=X', 'A=Y'}
%                                  work the way you'd expect: first
%                                  match wins, later same-name entries
%                                  are quietly skipped.
%
%   Always errors:
%       read_EDF:RateMismatch    - leaves of one spec / reference span
%                                  different sampling rates. The
%                                  toolbox does not resample on the
%                                  fly; resample after read_EDF,
%                                  precompute compatible inputs, or
%                                  pass 'TargetFs' to resample raw
%                                  signals before the derived layer
%                                  runs.
%       read_EDF:BadMean         - mean(...) called with fewer than
%                                  2 arguments.
%
%   Validation runs whenever 'References' is non-empty or any
%   'Channels' entry contains 'mean(', '=', '+', or '$' — even if the
%   caller ignored read_EDF's outputs. Plain labels and legacy 'A-B'
%   strings stay on the existing fast backend (MEX or MATLAB) and are
%   bit-identical to prior versions; the MEX backend also emits
%   read_EDF_mex:UnknownChannel as a warning when a requested plain
%   label is not in the file (regardless of nargout).
%
%   -------------------------------------------------------------------------
%   Examples (full call sites):
%
%       % Plain EDF, all channels
%       [hdr, shdr, sc, ann] = read_EDF('sleep.edf');
%
%       % Subset by label, with a legacy A-B reref
%       [~, ~, sc] = read_EDF('sleep.edf', ...
%           'Channels', {'EEG C3-A2', 'EEG O2-A1'});
%
%       % Linked-mastoid via inline mean()
%       [~, ~, sc] = read_EDF('sleep.edf', ...
%           'Channels', {'C3 - mean(A1, A2)'});
%
%       % Linked-mastoid via a named reference reused across channels
%       [~, ~, sc] = read_EDF('sleep.edf', ...
%           'References', {'LM = mean(A1, A2)'}, ...
%           'Channels',   {'C3-LM', 'C4-LM', '-LM'});
%
%       % EDF where a label has a '+' in it (e.g. 'EEG A+B'). Bare
%       % 'EEG A+B' would be parsed as 'EEG A' + 'B' -- escape it:
%       [~, ~, sc] = read_EDF('weird.edf', ...
%           'Channels', {'$EEG A+B$ - mean($A1$, $A2$)'});
%
%       % EDF that has both 'C1', 'A2', AND a labeled 'C1-A2' channel.
%       % Bare 'C1-A2' returns the labeled channel; '$C1$-$A2$' forces
%       % the computed difference of the raw channels.
%       [~, ~, sc_labeled]  = read_EDF('overlap.edf', 'Channels', {'C1-A2'});
%       [~, ~, sc_computed] = read_EDF('overlap.edf', 'Channels', {'$C1$-$A2$'});
%
%       % Same-channel-different-name fallback: load whichever of A,
%       % X, Y, Z exists in this file under the canonical name 'A'.
%       % First match wins; misses warn-and-skip.
%       [~, ~, sc] = read_EDF('mixed.edf', ...
%           'Channels', {'A', 'A=X', 'A=Y', 'A=Z'});
%
%       % Compressed input (streamed via zlib, no temp file in MEX path)
%       [hdr, shdr, sc, ann] = read_EDF('sleep.edf.gz');
%
%   -------------------------------------------------------------------------
%   EDF File Specification Summary:
%       • Each EDF file begins with a fixed-length 256-byte main header
%       • Followed by per-signal headers (16 fields x N signals)
%       • Digital samples stored as int16 are scaled to physical units:
%
%             phys_val = phys_min + (dig_val - dig_min) * (phys_max - phys_min) / (dig_max - dig_min)
%
%       • EDF+ annotation channels (labelled 'EDF Annotations') contain onset
%         times and event texts in TAL (Time-Annotation List) format.
%
%   -------------------------------------------------------------------------
%   Compressed input notes:
%       • Files ending in '.gz' (case-insensitive) are read directly through
%         the bundled zlib (vendored in the 'zlib/' subdirectory).
%       • Peak memory is one record's worth of raw bytes plus the per-signal
%         output arrays — files do not need to fit decompressed on disk.
%       • RepairHeader and deidentify both modify the file in place, which is
%         not meaningful on a gzip archive; both are silently disabled (with a
%         warning) for .gz inputs. Decompress to .edf first if you need either.

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
addParameter(p, 'References', {}, @iscell);
addParameter(p, 'Epochs', [], @isnumeric);
addParameter(p, 'Verbose', false, @islogical);
addParameter(p, 'RepairHeader', false, @islogical);
addParameter(p, 'forceMATLAB', false, @islogical);
addParameter(p, 'debug', false, @islogical);
addParameter(p, 'deidentify', false, @islogical);
% TargetFs: resample every raw signal to this rate BEFORE running the
% derived pipeline. Required when References combine constituents at
% different native rates -- apply_channel_derivations rejects mixed
% rates by design, so the caller (typically load_data with resampling
% enabled) must hand us the target up front. [] = no resample.
addParameter(p, 'TargetFs', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x>0));
parse(p, varargin{:});

channels        = p.Results.Channels;
references      = p.Results.References;
epochs          = p.Results.Epochs;
verbose         = p.Results.Verbose;
repair_header   = p.Results.RepairHeader;
force_matlab    = p.Results.forceMATLAB;
debug           = p.Results.debug;
deidentify      = p.Results.deidentify;
target_fs       = p.Results.TargetFs;

% Decide whether to route through the derived-channels post-processing
% layer. The legacy backend (MEX in C, or read_EDF_matlab in this file)
% understands plain labels and 'A-B' rereference strings; nothing richer.
% So whenever the caller used new syntax (References, mean(), aliasing,
% '+' sums, or '$...$' escape-quoted labels), we ask the backend for
% EVERY channel via Channels = {} and filter / derive after it returns.
% The trade-off is decoding signals the user may not need -- acceptable
% for typical PSG layouts (8-30 channels) but worth knowing about for
% high-density EEG. See has_derived_syntax for the detector.
use_derived_pipeline = ~isempty(references) || has_derived_syntax(channels);
backend_channels     = channels;
if use_derived_pipeline
    backend_channels = {};
end

%% ---------------- COMPRESSED INPUT GATING ----------------
% RepairHeader and deidentify both require modifying the file in place,
% which is not meaningful on a compressed archive. Force them off and warn.
is_gz  = endsWith(edf_fname, '.gz',  'IgnoreCase', true);
is_zst = endsWith(edf_fname, '.zst', 'IgnoreCase', true);
is_compressed = is_gz || is_zst;
if is_compressed
    suffix = '.gz'; if is_zst, suffix = '.zst'; end
    if repair_header
        warning('read_EDF:RepairOnCompressed', ...
            'RepairHeader is not supported for %s inputs; in-memory header will still be corrected.', suffix);
        repair_header = false;
    end
    if deidentify
        warning('read_EDF:DeidentifyOnCompressed', ...
            'deidentify is not supported for %s inputs; decompress to .edf first if you need a deidentified copy.', suffix);
        deidentify = false;
    end
end

%% ---------------- MEX HANDLING ----------------
script_dir = fileparts(mfilename('fullpath'));
mex_file = fullfile(script_dir, ['read_EDF_mex.' mexext]);
mex_exists = isfile(mex_file);

% When the derived pipeline is active, validation and evaluation need
% (header, signal_header, signal_cell) regardless of how many outputs
% the caller asked for. Bumping backend_nargout up to 3 ensures that
% even a `read_EDF(f, 'Channels', {...})` call with no LHS captures
% will still raise UnknownChannel / RateMismatch / etc. on bad input.
% At the very end we trim varargout back down to the caller's nargout
% so the MATLAB return semantics are unchanged.
backend_nargout = nargout;
if use_derived_pipeline || ~isempty(target_fs)
    backend_nargout = max(nargout, 3);
end

mex_succeeded = false;
if ~force_matlab
    if ~mex_exists
        compile_edf_mex(script_dir, 'read_EDF_mex.c');
    end
    try
        [varargout{1:backend_nargout}] = read_EDF_mex(edf_fname, backend_channels, epochs, verbose, repair_header, debug);

        if backend_nargout>0
            %Add the total data in seconds
            total_seconds = varargout{1}.num_data_records * varargout{1}.data_record_duration;
            varargout{1}.total_data_seconds = total_seconds;
            varargout{1}.total_data_hms = char(duration(0,0,total_seconds));
        end

        %Remove the whitespace around the labels
        if backend_nargout>1
            new_signal_labels = cellfun(@strip,{varargout{2}.signal_labels},'UniformOutput', false);
            [varargout{2}.signal_labels] = deal(new_signal_labels{:});
        end

        mex_succeeded = true;
    catch ME
        if verbose
            fprintf('MEX failed (%s). Falling back to MATLAB reader.\n', ME.message);
        end
    end
else
    if verbose
        fprintf('forceMATLAB = true, using MATLAB reader.\n');
    end
end

%% ---------------- MATLAB FALLBACK ----------------
if ~mex_succeeded
    [varargout{1:backend_nargout}] = read_EDF_matlab(edf_fname, backend_channels, epochs, verbose, repair_header, deidentify);
end

%% ---------------- TARGET-FS RESAMPLE (pre-derive) ----------------
% Resample every raw signal to a common rate before References /
% derived-channel math runs. Without this, references that combine
% constituents at different native rates (e.g. EEG@128 + EOG@32) fail
% inside apply_channel_derivations. Done here -- on the raw,
% per-channel signal_cell -- so each row gets resampled with its own
% native Fs. After this loop every signal_header.sampling_frequency
% equals target_fs and the derived pipeline sees uniform rates.
if ~isempty(target_fs) && backend_nargout >= 3
    sh_raw = varargout{2};
    sc_raw = varargout{3};
    for ii = 1:numel(sc_raw)
        f_native = sh_raw(ii).sampling_frequency;
        if abs(f_native - target_fs) > 1e-9
            sc_raw{ii} = smartresample(sc_raw{ii}, f_native, target_fs);
            sh_raw(ii).sampling_frequency = target_fs;
        end
    end
    varargout{2} = sh_raw;
    varargout{3} = sc_raw;
end

%% ---------------- DERIVED PIPELINE (post-load) ----------------
% Both backends produced (header, signal_header, signal_cell) for the
% full channel set; this layer applies References and 'Channels'
% expressions identically regardless of which backend ran. It's
% deliberately OUTSIDE the MEX try/catch so a validation error
% (UnknownChannel, RateMismatch, RefCollision, BadMean, ParseError)
% propagates to the caller instead of triggering a silent fallback to
% the MATLAB reader, which would just hit the same error a moment
% later anyway.
if use_derived_pipeline
    [varargout{2}, varargout{3}] = ...
        apply_channel_derivations(varargout{2}, varargout{3}, channels, references, verbose);
end

%% ---------------- DEIDENTIFY (post-MEX only) ----------------
% deidentify modifies the file on disk; only meaningful when the MEX
% backend succeeded (the MATLAB fallback already handles it inline).
if mex_succeeded && deidentify
    deidentify_edf(edf_fname, repair_header, verbose);
end

% Restore the caller's nargout. backend_nargout may have been bumped
% above to make sh / sc available to the derived pipeline; trim off
% anything the caller didn't ask for so the function's apparent
% multiple-return signature matches MATLAB's expectations.
if numel(varargout) > nargout
    varargout = varargout(1:nargout);
end

end


%% =========================================================================
%  PURE MATLAB EDF READER
% =========================================================================
function varargout = read_EDF_matlab(edf_fname, channels, epochs, verbose, repair_header, deidentify)

% If input is compressed, decompress to a temp .edf and run the rest
% against it. Temp directory is removed when this function returns.
comp_cleanup = []; %#ok<NASGU>
if endsWith(edf_fname, '.gz', 'IgnoreCase', true)
    if verbose
        fprintf('Decompressing %s to temp file for MATLAB reader...\n', edf_fname);
    end
    tmpdir = tempname;
    mkdir(tmpdir);
    comp_cleanup = onCleanup(@() rmdir(tmpdir, 's'));
    gunzipped = gunzip(edf_fname, tmpdir);
    edf_fname = gunzipped{1};
elseif endsWith(edf_fname, '.zst', 'IgnoreCase', true)
    if verbose
        fprintf('Decompressing %s to temp file for MATLAB reader...\n', edf_fname);
    end
    tmpdir = tempname;
    mkdir(tmpdir);
    comp_cleanup = onCleanup(@() rmdir(tmpdir, 's'));
    [~, base, ~] = fileparts(edf_fname);   % strip .zst -> base.edf
    out_path = fullfile(tmpdir, base);
    cmd = sprintf('zstd -q -d -f -o %s %s', shell_quote(out_path), shell_quote(edf_fname));
    [rc, msg] = system(cmd);
    if rc ~= 0
        error('read_EDF:ZstdDecompress', ...
            'zstd command failed for %s: %s', edf_fname, strtrim(msg));
    end
    edf_fname = out_path;
end

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

new_signal_labels = cellfun(@strip,{signal_header.signal_labels},'UniformOutput', false);
[signal_header.signal_labels] = deal(new_signal_labels{:});

%% ---------------- CHANNEL PLAN ----------------
% Inspect all requested channels. Identify plain channel requests and A-B
% rereference pairs. Determine the superset of raw signals to decode.
all_labels = {signal_header.signal_labels};

if isempty(channels)
    signal_indices = 1:num_signals;
    plain_chs   = {};
    reref_pairs = {};
else
    [plain_chs, reref_pairs, load_labels] = parse_channel_plan(channels, all_labels);
    signal_indices = find(ismember(lower(all_labels), lower(load_labels)));
    if isempty(signal_indices)
        error('read_EDF:NoChannels', 'No valid channels found.');
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
            % uint8, not 'char': 'char' precision is encoding-dependent (can emit UTF-16).
            fwrite(fw, uint8(rec_str), 'uint8');
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
% Decode only the channels in signal_indices (plain + reref constituents).
signal_cells = cell(1, length(signal_indices));

for i = 1:length(signal_indices)

    sidx = signal_indices(i);
    offset = sum(samples_per_record(1:sidx-1));
    samples_per_epoch = samples_per_record(sidx);
    total_samples = samples_per_epoch * num_epochs;

    sig = zeros(1, total_samples);

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

% Trim signal_header to only the loaded subset (signal_cells{i} <-> signal_header(i)).
signal_header = signal_header(signal_indices);

%% ---------------- REREFERENCING ----------------
% Compute A-B difference signals and remove any constituent channels that
% were not independently requested.
if ~isempty(channels)
    [signal_header, signal_cells] = apply_channel_plan(signal_header, signal_cells, plain_chs, reref_pairs);
end

annotations = extractAnnotations(edf_fname, header, signal_header);

%% ---------------- DEIDENTIFY ----------------
if deidentify
    deidentify_edf(edf_fname, repair_header, verbose);

    % Also scrub the returned header struct so it matches the written file
    header.patient_id        = 'X X X X';
    header.local_rec_id      = 'Startdate X X X X';
    header.recording_startdate = '01.01.01';
end

%% ---------------- OUTPUT ----------------
varargout{1} = header;
if nargout > 1, varargout{2} = signal_header; end
if nargout > 2, varargout{3} = signal_cells; end
if nargout > 3, varargout{4} = annotations; end

end


%% =========================================================================
%  CHANNEL PLAN PARSER
% =========================================================================
function [plain_chs, reref_pairs, load_labels] = parse_channel_plan(channels, all_labels)
%PARSE_CHANNEL_PLAN  Classify requested channels as plain or A-B reref pairs.
%
%   For each entry in channels:
%     - If it directly matches a label (case-insensitive) → plain channel.
%     - Otherwise try each '-' as a split point. The first split where both
%       sides are valid label names is treated as a rereference request.
%     - Anything else is passed through as a plain name (will error on load).
%
%   Returns:
%     plain_chs   : cell array of directly-requested channel name strings
%     reref_pairs : cell array of {chA, chB, 'chA-chB'} rows
%     load_labels : unique set of labels to actually decode (plain + constituents)

all_labels_lower = lower(cellfun(@strtrim, all_labels, 'UniformOutput', false));
plain_chs   = {};
reref_pairs = {};

for k = 1:numel(channels)
    ch = strtrim(channels{k});

    % Direct label match → plain channel
    if ismember(lower(ch), all_labels_lower)
        plain_chs{end+1} = ch;
        continue;
    end

    % Try each '-' as a split point, leftmost first
    dashes = strfind(ch, '-');
    found  = false;
    for d = dashes
        chA = strtrim(ch(1:d-1));
        chB = strtrim(ch(d+1:end));
        if ~isempty(chA) && ~isempty(chB) && ...
                ismember(lower(chA), all_labels_lower) && ...
                ismember(lower(chB), all_labels_lower)
            reref_pairs{end+1} = {chA, chB, ch};
            found = true;
            break;
        end
    end

    if ~found
        % Unknown channel — warn and drop. Matches the MEX backend, which
        % emits read_EDF_mex:UnknownChannel and skips. Keeping the load
        % running on the channels that did resolve is intentional.
        warning('read_EDF:UnknownChannel', ...
            'Channel ''%s'' not found in file.', ch);
    end
end

% Superset of labels to decode = plain channels + both sides of every reref pair
load_labels = plain_chs;
for k = 1:numel(reref_pairs)
    load_labels{end+1} = reref_pairs{k}{1};
    load_labels{end+1} = reref_pairs{k}{2};
end
load_labels = unique(load_labels, 'stable');
end


%% =========================================================================
%  APPLY CHANNEL PLAN
% =========================================================================
function [sh_out, sc_out] = apply_channel_plan(signal_header, signal_cell, plain_chs, reref_pairs)
%APPLY_CHANNEL_PLAN  Filter to requested channels and compute A-B reref signals.
%
%   signal_header / signal_cell contain the loaded superset (plain + constituents).
%   Output keeps only the explicitly requested plain channels, then appends one
%   entry per reref pair. Constituent-only channels are dropped.

labels = {signal_header.signal_labels};

% Indices of plain channels in the loaded data (preserve request order).
% Anything that doesn't resolve gets dropped — parse_channel_plan should
% already have warned for unknown labels, but we guard here too in case
% a known channel didn't make it into the loaded subset for some reason.
plain_idx = zeros(1, numel(plain_chs));
for k = 1:numel(plain_chs)
    f = find(strcmpi(labels, strtrim(plain_chs{k})), 1);
    if ~isempty(f)
        plain_idx(k) = f;
    end
end
plain_idx = plain_idx(plain_idx > 0);

% Compute rereferenced signals
n = numel(reref_pairs);
reref_sc = cell(1, n);
if n > 0
    reref_sh(1:n) = signal_header(1);   % preallocate with matching struct fields
    for k = 1:n
        chA  = reref_pairs{k}{1};
        chB  = reref_pairs{k}{2};
        name = reref_pairs{k}{3};
        idxA = find(strcmpi(labels, chA), 1);
        idxB = find(strcmpi(labels, chB), 1);
        reref_sh(k)               = signal_header(idxA);
        reref_sh(k).signal_labels = name;
        reref_sc{k}               = signal_cell{idxA} - signal_cell{idxB};
    end
    sh_out = [signal_header(plain_idx), reref_sh];
else
    sh_out = signal_header(plain_idx);
end
sc_out = [signal_cell(plain_idx), reref_sc];
end


%% =========================================================================
%  DERIVED-CHANNELS PIPELINE
%  -----------------------------------------------------------------------
%  Post-load layer that adds named References (mean-of-N etc.) and a small
%  expression syntax inside 'Channels' to read_EDF without touching the
%  MEX. Activated by the public read_EDF dispatcher whenever 'References'
%  is non-empty or any 'Channels' entry contains a marker character that
%  the legacy backend wouldn't understand. When inactive, plain labels
%  and legacy 'A-B' strings flow straight through the MEX (or MATLAB
%  fallback) and the output is bit-identical to prior versions.
%
%  The pipeline operates on physical-units signals returned by the
%  backend, so re-referencing math runs in float and inherits the
%  scaling that read_EDF already applied.
% =========================================================================

function tf = has_derived_syntax(channels)
%HAS_DERIVED_SYNTAX  Detect whether any 'Channels' entry needs the new pipeline.
%
%   Returns true if any string contains a derived-syntax marker:
%     mean(  =  +  $  *  /  (
%
%   '+' covers sums; '*' and '/' cover scalar-weighted terms; '(' covers
%   grouping (and 'mean(' as a side effect); '=' covers aliased outputs;
%   '$' covers escape-quoted labels. Number literals always show up
%   next to one of these in a real arithmetic spec, so digits alone
%   don't need to trigger.
%
%   We trigger on '(' even though real EDF labels can contain parens —
%   the false-positive cost is just routing through the derived path,
%   which is still correct for plain labels.
%
%   These markers are absent from plain labels and legacy 'A-B' strings —
%   both of which the MEX / MATLAB backend handles directly — so this is
%   a cheap conservative gate. False positives are harmless (we just
%   route through the post-processing layer needlessly).
tf = false;
for k = 1:numel(channels)
    s = lower(channels{k});
    if contains(s, '=') || contains(s, '+') || contains(s, '$') ...
            || contains(s, '*') || contains(s, '/') || contains(s, '(')
        tf = true;
        return
    end
end
end




%% =========================================================================
%  DEIDENTIFY HELPER — copy file with PHI fields blanked in the header
% =========================================================================
function deidentify_edf(edf_fname, repair_header, verbose)
% Build output filename with appropriate suffix
[fdir, fname, fext] = fileparts(edf_fname);
if repair_header
    out_fname = fullfile(fdir, [fname '_fixed_deidentified' fext]);
else
    out_fname = fullfile(fdir, [fname '_deidentified' fext]);
end

% Copy the entire file first, then patch header fields in the copy
copyfile(edf_fname, out_fname);

fw = fopen(out_fname, 'r+', 'ieee-le');
if fw < 0
    warning('read_EDF:DeidentifyFailed', ...
        'Could not open output file for deidentification: %s', out_fname);
    return
end

% EDF main header field byte offsets and widths (0-based offsets):
%   patient_id          offset=8,   width=80
%   local_rec_id        offset=88,  width=80
%   recording_startdate offset=168, width=8
phi_fields = {
    'X X X X',          8,   80;   % patient_id
    'Startdate X X X X', 88,  80;   % local_rec_id
    '01.01.01',         168,  8;   % recording_startdate
    };

for k = 1:size(phi_fields, 1)
    val     = phi_fields{k,1};
    offset  = phi_fields{k,2};
    width   = phi_fields{k,3};

    % Left-justify value, space-pad to exact field width
    padded = sprintf('%-*s', width, val);
    padded = padded(1:width);   % truncate if somehow over width

    fseek(fw, offset, 'bof');
    % EDF header is 8-bit ASCII; 'char' precision is encoding-dependent and
    % can emit 2 bytes/char (UTF-16). uint8 forces one byte per character.
    fwrite(fw, uint8(padded), 'uint8');
end

fclose(fw);

if verbose
    fprintf('Deidentified file written: %s\n', out_fname);
end
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

%% =========================================================================
%  SHELL QUOTING (POSIX, for system() calls)
% =========================================================================
function s = shell_quote(s)
s = ['''' strrep(s, '''', '''\''''') ''''];
end
