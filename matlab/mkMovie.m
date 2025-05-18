function mkMovie(Folder, movie_name)

%Folder = 'Data70';

filesInfo = dir(fullfile(Folder, '*final*.mat'));

% Extract file names
fileNames = {filesInfo.name};

% Extract numerical parts of file names and convert to numbers
% Assuming the format is always 'Iteration_X_final.mat'
fileNumbers = cellfun(@(x) sscanf(x, 'Iteration_%d_final.mat'), fileNames);

% Sort based on the numerical parts
[~, sortIdx] = sort(fileNumbers);
sortedFileNames = fileNames(sortIdx);

% Store in a cell array
files = cell(size(sortedFileNames));
for i = 1:length(sortedFileNames)
    files{i} = fullfile(Folder, sortedFileNames{i});
end

% files now contains sorted file paths
outputVideoName = [movie_name, '.mpeg'];

% Initialize the video writer object with MPEG-4 compression
video = VideoWriter([Folder,'\', outputVideoName], 'MPEG-4');
video.FrameRate = 5; % Set the frame rate of the video
video.Quality = 100; % Set the quality from 0 to 100 (higher means better quality)
open(video);

% Adjust figure size for desired resolution
figureWidth = 1280; % width in pixels
figureHeight = 720; % height in pixels
figure('Position', [100, 100, figureWidth, figureHeight], 'Color', [1 1 1]);


% Assuming files{} has been populated as per the previous example
for j = 1:length(files)
    loadedData = load(files{j});
    % Assuming x, y, and params are variables in your .mat file
    % Adjust these variable names as per your actual file contents
    x = loadedData.x;
    y = loadedData.y;
    params = loadedData.params;
    
    % Clear the current figure to avoid drawing over previous plots
    clf;
    
    % Your custom Draw function to plot data
    if j == 1
       cc = DrawCells(x, y, params);
    else
       DrawCells(x, y, params, cc);
    end

    
    % % Maximize figure to full screen - This is a trick to ensure the figure fully utilizes the screen
    set(gca, 'Position', [0 0 1 1]);
    % 
    % % Adjust axes to fill the figure window
    % ax = gca;
    % ax.Units = 'normalized';
    % ax.Position = [0 0 1 1]; % Make the axes fill the figure
    
    % After adjustments, capture the frame
    frame = getframe(gcf);

    % Write the frame to the video
    writeVideo(video, frame);
end

close(gcf)

% Close the video file
close(video);

% Your Draw function needs to be defined somewhere in your scripts
% It might look something like this:
% function Draw(x, y, params)
%     plot(x, y, params);
%     % Add any additional plotting commands here
% end


