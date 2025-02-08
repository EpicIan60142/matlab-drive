function makePlot(timeVector,StateMatrix)
% Plots the trajectory of the object
% Inputs: Time vector, state vector

figure % Open a new figure window
plot(StateMatrix(:,3),StateMatrix(:,4), 'g'); % Plot distance on the x axis and height on the y axis in green
hold on % Save the graph in the figure
grid on % Display the coordinate grid in the figure
title("Height vs. Distance"); % Give the figure an appropriate title
xlabel("Distance"); % Properly label the x axis
ylabel("Height"); % Properly label the y axis

figure % Open a new figure window
plot(timeVector,StateMatrix(:,2), 'b'); % Plot time on the x axis and vertical velocity on the y axis in blue
hold on % Save the graph in the figure
plot(timeVector,StateMatrix(:,4), 'm'); % Plot time on the x axis and height on the y axis in magenta
hold on % Save the graph in the figure
grid on % Display the coordinate grid in the figure
title("Vertical Velocity and Height vs. time"); % Give the figure an appropriate title
xlabel("Time"); % Properly label the x axis
ylabel("Vertical velocity, Height"); % Properly label the y axis
legend("Vertical velocity vs. time", "Height vs. time"); % Create a legend to distinguish between the two plots

figure % Create a new figure window
plot(timeVector,StateMatrix(:,1), 'c'); % Plot time on the x axis and horizontal velocity on the y axis in cyan
hold on % Save the graph in the figure
plot(timeVector,StateMatrix(:,3), 'r'); % Plot time on the x axis and height on the y axis in red
hold on % Save the graph in the figure
grid on % Display the coordinate grid in the figure
title("Horizontal Velocity and Distance vs. time"); % Give the figure an appropriate title
xlabel("Time"); % Properly label the x axis
ylabel("Horizontal velocity, Distance"); % Properly label the y axis
legend("Horizontal velocity vs. time", "Distance vs. time") % Create a legend to distinguish between the two plots

end
