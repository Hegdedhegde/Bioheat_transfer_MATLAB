%% Convert to Cartesian coordinates for plotting
r = linspace(0, R, nr);
theta = linspace(0, pi, Ntheta);
[r_grid, theta_grid] = meshgrid(r, theta);
X = r_grid .* cos(theta_grid);
Y = r_grid .* sin(theta_grid);

%% Plot the results
figure;
contourf(X, Y, T', 20, 'LineColor', 'none');
colorbar;
title('Bioheat Conduction in a Semicircular Domain');
xlabel('x (cm)');
ylabel('y (cm)');
axis equal;
colormap('jet');

hold on;

%% plotting the tumor
% Circle parameters
r = rt; % radius of the circle
theta = linspace(0, 2*pi, 100); % angle values for the circle

% Specify the center of the circle (e.g., at (3, 4))
x_center = rc*cos(theta_c);
y_center = rc*sin(theta_c);

% Parametric equations for the circle with the specified center
x_circle = x_center + r * cos(theta); % x coordinates of the circle
y_circle = y_center + r * sin(theta); % y coordinates of the circle

% Plot the circle on top of the contour plot
plot(x_circle, y_circle, 'k-', 'LineWidth', 1); % 'k-' for black line

hold off;

disp(min(min(T)) + 273)
disp(max(max(T)) + 273)