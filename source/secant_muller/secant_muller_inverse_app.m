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

  statusLabel = uilabel(fig, 'Position', [20 438 860 22], 'Text', 'Ready');

  tbl = uitable(fig, ...
    'Position', [20 20 860 410], ...
    'Data', data0, ...
    'ColumnName', cellstr(cfg.columns), ...
    'ColumnEditable', [false true true true true], ...
    'ColumnFormat', {'numeric', cellstr(cfg.methodOptions), 'numeric', 'numeric', 'numeric'}, ...
    'ColumnWidth', {'fit', 'fit', 'fit','auto','auto'}, ...
    'CellEditCallback', @(src, evt) on_table_edit(src, evt, targetField, statusLabel));

  uibutton(fig, 'push', ...
    'Text', 'Load Example', ...
    'Position', [215 465 95 22], ...
    'ButtonPushedFcn', @(~, ~) on_load_example(tbl, targetField, cfg, statusLabel));

  uibutton(fig, 'push', ...
    'Text', 'Clear', ...
    'Position', [315 465 70 22], ...
    'ButtonPushedFcn', @(~, ~) on_clear_table(tbl, cfg, statusLabel));
end

function on_clear_table(tbl, cfg, statusLabel)
  nRows = size(tbl.Data, 1);
  tbl.Data = build_default_data(nRows, cfg);
  statusLabel.Text = 'Cleared. Enter x/y points. Secant needs 2 valid points; Muller needs 3.';
end

function on_load_example(tbl, targetField, cfg, statusLabel)
  nRows = size(tbl.Data, 1);
  tbl.Data = build_default_data(nRows, cfg);
  targetField.Value = 0;

  example = {
    1, 'secant', 1, 1.6, 1.0122;
    2, 'secant', 1, 1.5, 0.9828;
    3, 'muller', 1, -1.8417, -1.0734;
    4, 'muller', 1, -0.5038, -0.4667;
    5, 'muller', 1, 0.2463, 0.2415;
    6, 'muller', 1, 0.0200, 0.0200;
    7, 'muller', 1, -7.703e-04, NaN
  };

  nLoad = min(nRows, size(example, 1));
  for k = 1:nLoad
    r = example{k, 1};
    tbl.Data{r, 2} = example{k, 2};
    tbl.Data{r, 3} = example{k, 3};
    tbl.Data{r, 4} = example{k, 4};
    tbl.Data{r, 5} = example{k, 5};
  end

  statusLabel.Text = 'Example rows loaded (target y = 0). Edit y, step size, or method to compute next-row x.';
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
