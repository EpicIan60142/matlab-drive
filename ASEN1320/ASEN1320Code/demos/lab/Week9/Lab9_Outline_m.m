%% Lab 9 Outline
% 1. Print Formatting
% 2. Conditional Statements (syntax review)
% 3. Loops (syntax review)
% 4. Functions
% 5. Figures - plot and subplot

close all; clear; clc;


%% Print Formatting (fprintf)
% %s 	Format as a string.
% %d 	Format as an integer.
% %f 	Format as a floating point value.
% %e 	Format as a floating point value in scientific notation.
% %g 	Format in the most compact form: %f or %e.
% \n 	Insert a new line in the output string.
% \t 	Insert a tab in the output string.
clc

fprintf("Here is a string: %s\n", "string")
fprintf("Here is a int: %d\n", 1)
fprintf("Here is a double (floating point value): %f\n", 3.14)
fprintf("Here is scientific notation: %e\n", 3.14)
fprintf("Print a percent %% and a quote ""\n")


%% Conditional Statements
% go through this quickly - just remind them of the syntax
clc

% booleans 
true    % 1
false   % 0

% boolean logic
~true
~false
true && false
true || false
1 == 2
1 ~= 2
1 > 2
1 >= 2
1 < 2
1 <= 2

% if statement
if (true)
    % this will run
end

% if-else statement
if (false)
    % this won't run
else
    % this will run
end

% elseif statement
if (false)
    % this won't run
elseif (false)
    % this won't run
else
    % this will run
end


%% Loops
% go through this quickly as well - just to remind them of the syntax
clc

% note that 'i' and 'j' in matlab are the imaginary units, so consider a
% different variable for loop counters
i
j

% loop through specific numbers
for k = 1:5
    fprintf("%d ", k)
end
fprintf("\n")

% increement by 2 until the length of an array
array = [1 2 3 4 5 6 7 8 9 10];
for k = 1:2:length(array)
    fprintf("%d ", array(k))
end
fprintf("\n")

% loop through an array
numbers = ["one", "two", "three"];
for k = numbers
    fprintf("%s ", k)
end
fprintf("\n")

% while loop
k = 0;
while (k < 5)
    fprintf("%d ", k)
    k = k+1;  % can't do k++ or k+=1
end
fprintf("\n")


%% User-Defined Functions
% navigate to new -> function
clc

% using function outputs
quadratic(1, -4, 4)             % POLL - what prints? (only first output)
x1 = quadratic(2,4,-4)          % again, only the first output
[x1, x2] = quadratic(2,4,-4)    % gets both outputs - note these are separate
                                % variables, not an array


%% Plotting
close all; clc;

x = 0:0.5:100;
y1 = x;
y2 = 2*x;

% basics
figure          % creates a new figure
plot(x,y1)      % first arg is x, second is y
hold on         % plot two lines on same figure
plot(x,y2)
grid on         % shows grid lines
legend("y=x", "y=2x")
title("My Plots")
xlabel("value of x")
ylabel("value of y")

% subplots
figure

subplot(1,2,1) % first arg is rows, second arg is cols, third is plot num
plot(x,y1)
title("y=x")
xlabel("value of x")
ylabel("value of y")
xlim([0 100])  % show plots without limits first
ylim([0 200])

subplot(1,2,2) % first arg is rows, second arg is cols, third is plot num
plot(x,y2)
title("y=2x")
xlabel("value of x")
ylabel("value of y")
xlim([0 100])  % show plots without limits first
ylim([0 200])

