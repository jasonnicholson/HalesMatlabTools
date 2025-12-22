function setupHalesToolsPath(isSetup)
  %setupHalesToolsPath sets up or breakdown the path
  
  arguments (Input)
    isSetup (1,1) logical = true;
  end

  thisFolder = fileparts(mfilename("fullpath"));
  if isSetup
    addpath(thisFolder);
  else
    rmpath(thisFolder);
  end
end