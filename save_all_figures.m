function save_all_figures(folder_name)

if isempty(folder_name) || false
    disp("Nothig was saved")
    return
end

if ~exist(folder_name, 'dir')
   mkdir(folder_name)
end


FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = sprintf('fig_%06d', get(FigHandle, 'Number'));
  
  disp(FigName)
  % set(0, 'CurrentFigure', FigHandle);
  filename = fullfile(folder_name, [FigName '.fig']);
  savefig(FigHandle, filename);
  filename = fullfile(folder_name, [FigName '.png']);
  saveas(FigHandle,filename);
end
end