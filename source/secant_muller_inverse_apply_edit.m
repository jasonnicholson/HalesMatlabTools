function [dataOut, eventOut] = secant_muller_inverse_apply_edit(dataIn, editedRow, editedCol, targetY, opts)
%SECANT_MULLER_INVERSE_APPLY_EDIT Apply row-n edit trigger and recompute row n+1 x.
% Inputs:
%   dataIn    - cell array with columns [Iteration, Method, Step Size, x, y]
%   editedRow - 1-based row index edited by user
%   editedCol - 1-based column index edited by user
%   targetY   - scalar desired y
%   opts      - optional struct:
%                 denominatorTol, duplicateXTol
%
% Outputs:
%   dataOut   - updated cell array
%   eventOut  - struct with trigger metadata and compute result

arguments
  dataIn cell
  editedRow (1,1) double {mustBeInteger, mustBePositive}
  editedCol (1,1) double {mustBeInteger, mustBePositive}
  targetY (1,1) double = 1.0
  opts struct = struct()
end

triggerColumns = [2, 3, 5]; % Method, Step Size, y

dataOut = dataIn;
eventOut = struct( ...
  "triggered", false, ...
  "triggerSourceRow", editedRow, ...
  "triggerTargetRow", editedRow + 1, ...
  "blockedReason", "", ...
  "ok", false, ...
  "xNext", NaN, ...
  "methodUsed", "");

if editedCol > size(dataIn, 2) || editedRow > size(dataIn, 1)
  eventOut.blockedReason = "out_of_range";
  return;
end

if ~ismember(editedCol, triggerColumns)
  eventOut.blockedReason = "non_trigger_field";
  return;
end

eventOut.triggered = true;
targetRow = editedRow + 1;
if targetRow > size(dataIn, 1)
  eventOut.blockedReason = "target_row_missing";
  return;
end

method = string(cell_to_default(dataIn{editedRow, 2}, ""));
stepSize = scalar_from_cell(dataIn{editedRow, 3}, 1.0);

points = collect_points_up_to_row(dataIn, editedRow);
out = secant_muller_inverse_next_x( ...
  points, targetY, stepSize, method, ...
  struct( ...
    "denominatorTol", get_opt(opts, "denominatorTol", 1e-12), ...
    "duplicateXTol", get_opt(opts, "duplicateXTol", 1e-12)));

eventOut.ok = out.ok;
eventOut.methodUsed = out.methodUsed;
eventOut.xNext = out.xNext;
eventOut.blockedReason = out.blockedReason;

if out.ok
  dataOut{targetRow, 4} = out.xNext;
end
end

function points = collect_points_up_to_row(data, rowIdx)
points = zeros(0, 2);
for k = 1:rowIdx
  xk = scalar_from_cell(data{k, 4}, NaN);
  yk = scalar_from_cell(data{k, 5}, NaN);
  if isfinite(xk) && isfinite(yk)
    points(end + 1, :) = [xk, yk]; %#ok<AGROW>
  end
end
end

function value = scalar_from_cell(raw, fallback)
if isempty(raw)
  value = fallback;
  return;
end

if isnumeric(raw) && isscalar(raw) && isfinite(raw)
  value = double(raw);
  return;
end

if isstring(raw) || ischar(raw)
  tmp = str2double(string(raw));
  if isfinite(tmp)
    value = tmp;
    return;
  end
end

value = fallback;
end

function value = cell_to_default(raw, fallback)
if isstring(raw) || ischar(raw)
  value = string(raw);
  return;
end
value = fallback;
end

function value = get_opt(opts, key, defaultValue)
if isfield(opts, key)
  value = opts.(key);
else
  value = defaultValue;
end
end
