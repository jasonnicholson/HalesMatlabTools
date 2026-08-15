function setupHalesToolsPath(isSetup)
  %setupHalesToolsPath sets up or breakdown the path
  
  arguments (Input)
    isSetup (1,1) logical = true;
  end

  thisFolder = fileparts(mfilename("fullpath"));
  folderList = { 
    thisFolder;
    fullfile(thisFolder, "secant_muller")
    };
  if isSetup
    addpath(folderList{:});
  else
    rmpath(folderList{:});
  end
end