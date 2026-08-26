% Auto-Format Active Figure for LaTeX (Full Page Width)
fig = gcf; % Get current figure

% 1. Set Figure Dimensions (17 cm width is standard A4 \textwidth)
fig.Units = 'centimeters';
currentPos = fig.Position;
aspectRatio = currentPos(4) / currentPos(3);
targetWidth = 17; % Centimeters
targetHeight = targetWidth * aspectRatio; 
fig.Position = [currentPos(1), currentPos(2), targetWidth, targetHeight];

% 2. Find all axes in the figure
allAxes = findall(fig, 'type', 'axes');

% 3. Apply formatting to each axis
for i = 1:length(allAxes)
    ax = allAxes(i);
    
    % Remove titles
    title(ax, '');
    
    % Set axes font size 10, bold, and thicken the axis lines
    set(ax, 'FontSize', 6, 'FontWeight', 'bold', 'LineWidth', 1.2);
    
    % Ensure X, Y, and Z labels are explicitly size 10 and bold
    set(ax.XLabel, 'FontSize', 6, 'FontWeight', 'bold');
    set(ax.YLabel, 'FontSize', 6, 'FontWeight', 'bold');
    if isprop(ax, 'ZLabel')
        set(ax.ZLabel, 'FontSize', 6, 'FontWeight', 'bold');
    end
end

% 4. Format all legends (Font size 10)
allLegends = findall(fig, 'type', 'legend');
if ~isempty(allLegends)
    set(allLegends, 'FontSize', 4, 'FontWeight', 'normal'); 
end

% 5. Format all colorbars (Font size 10, bold)
allColorbars = findall(fig, 'type', 'colorbar');
if ~isempty(allColorbars)
    set(allColorbars, 'FontSize', 6, 'FontWeight', 'bold', 'LineWidth', 1.2);
end

% 6. RETROACTIVE LINE AND MARKER FORMATTING
% Find all standard line plots and set thickness/marker size
allLines = findall(fig, 'Type', 'line');
if ~isempty(allLines)
    % Change 'MarkerSize' to your preferred value (e.g., 4 or 6)
    % Set 'Marker' to 'none' if you want them completely removed
    set(allLines, 'LineWidth', 1.5, 'MarkerSize', 4); 
end

% Find all scatter plots (your control nodes) and set dot size
allScatters = findall(fig, 'Type', 'scatter');
if ~isempty(allScatters)
    % 'SizeData' controls the area of the scatter circles
    set(allScatters, 'SizeData', 30); 
end

fprintf('Figure formatted: 17cm width, Size 10 font, Bold axes, Uniform markers.\n');