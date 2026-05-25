% Define the number of points
numpoints = 1000;

% Generate x values
x = linspace(0, 4*pi, numpoints);

% Define y values for two functions
y1 = sin(x);
y2 = cos(x);

% Create figure and set axes limits
figure;
axis([0, 4*pi, -1, 1]);

% Create animated line objects for each function
h1 = animatedline('Color', 'r', 'LineWidth', 2);
h2 = animatedline('Color', 'b', 'LineWidth', 2);

% Add points to the animated lines in a loop
for k = 1:numpoints
    addpoints(h1, x(k), y1(k));
    addpoints(h2, x(k), y2(k));
    drawnow; % Update the figure
    pause(0.01);  % Slow down the drawing (adjust time as needed)
end
