function mkMovieGIF(Folder, gif_name)

filesInfo = dir(fullfile(Folder, '*final*.mat'));
fileNames = {filesInfo.name};
fileNumbers = cellfun(@(x) sscanf(x, 'Iteration_%d_final.mat'), fileNames);
[~, sortIdx] = sort(fileNumbers);
sortedFileNames = fileNames(sortIdx);
files = cellfun(@(f) fullfile(Folder, f), sortedFileNames, 'UniformOutput', false);

% Adjust figure size
% figureWidth = 720;
% figureHeight = 540;
figureWidth = 640;
figureHeight = 480;
% figureWidth = 480;
% figureHeight = 360;
figure('Position', [100, 100, figureWidth, figureHeight], 'Color', [1 1 1]);

gif_path = fullfile(Folder, [gif_name, '.gif']);
delayTime = 0.1; % seconds per frame (e.g. 5 fps)

for j = 1:2:length(files)*0.5
    loadedData = load(files{j});
    x = loadedData.x;
    y = loadedData.y;
    params = loadedData.params;

    clf;
    if j == 1
        cc = DrawCells(x, y, params);
    else
        DrawCells(x, y, params, cc);
    end
    axis off; % hide axes for clean GIF
    set(gca, 'Position', [0 0 1 1]);

    drawnow;
    frame = getframe(gcf);
    im = frame2im(frame);
    [A, map] = rgb2ind(im, 256);

    if j == 1
        imwrite(A, map, gif_path, 'gif', 'LoopCount', inf, 'DelayTime', delayTime);
    else
        imwrite(A, map, gif_path, 'gif', 'WriteMode', 'append', 'DelayTime', delayTime);
    end
end

close(gcf)
fprintf('Saved animated GIF to %s\n', gif_path);
