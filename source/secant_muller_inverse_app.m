function secant_muller_inverse_app(opts)
  %SECANT_MULLER_INVERSE_APP Launch MATLAB GUI for Secant/Muller inverse workflow.
  % app = secant_muller_inverse_app() opens the UI and returns handles.

  arguments
    opts.NumRows (1,1) double {mustBeInteger, mustBePositive} = 20;
  end

  cfg = secant_muller_inverse_contract();
  data0 = build_default_data(opts.NumRows, cfg);

  fig = uifigure( ...
    'Name', 'Secant/Muller Inverse Solver', ...
    'Position', [100 100 900 500]);

  uilabel(fig, 'Position', [20 465 60 22], 'Text', 'Target y');
  targetField = uieditfield(fig, 'numeric', 'Position', [85 465 120 22], 'Value', cfg.targetYDefault);

  statusLabel = uilabel(fig, 'Position', [220 465 650 22], 'Text', 'Ready');

  tbl = uitable(fig, ...
    'Position', [20 20 860 430], ...
    'Data', data0, ...
    'ColumnName', cellstr(cfg.columns), ...
    'ColumnEditable', [false true true true true], ...
    'ColumnFormat', {'numeric', cellstr(cfg.methodOptions), 'numeric', 'numeric', 'numeric'}, ...
    'ColumnWidth', {'fit', 'fit', 'fit','auto','auto'}, ...
    'CellEditCallback', @(src, evt) on_table_edit(src, evt, targetField, statusLabel)); %#ok<NASGU>
end

function data = build_default_data(nRows, cfg)
  data = cell(nRows, 5);
  for i = 1:nRows
    data{i, 1} = i;
    if i >= 3
      data{i, 2} = char(cfg.defaultMethodByPointCount.threePoints);
    else
      data{i, 2} = char(cfg.defaultMethodByPointCount.twoPoints);
    end
    data{i, 3} = cfg.defaultStepSize;
    data{i, 4} = NaN;
    data{i, 5} = NaN;
  end
end

function on_table_edit(src, evt, targetField, statusLabel)
  idx = evt.Indices;
  row = idx(1);
  col = idx(2);

  % load cells that should not change in the columns that do change
  persistent STATIC_CELLS
  if isempty(STATIC_CELLS)
    STATIC_CELLS = [1 2; 1 3; 2 2];
  end

  % Restore the static cells
  if ismember(idx, STATIC_CELLS, "rows")
    src.Data{row, col} = evt.PreviousData;
    return;
  end

  [dataOut, eventOut] = secant_muller_inverse_apply_edit( src.Data, row, col, targetField.Value);

  src.Data = dataOut;

  if ~eventOut.triggered
    statusLabel.Text = 'Edit accepted (no recompute trigger).';
    return;
  end

  if eventOut.ok
    statusLabel.Text = sprintf('Computed row %d x = %.10g (%s)', eventOut.triggerTargetRow, eventOut.xNext, eventOut.methodUsed);
  else
    statusLabel.Text = sprintf('Row %d compute blocked: %s', eventOut.triggerTargetRow, eventOut.blockedReason);
  end
end
