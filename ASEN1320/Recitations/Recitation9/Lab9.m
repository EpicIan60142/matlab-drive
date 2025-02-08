%% Array concatenation
% Putting multiple arrays together
clc
clear

array1 = 1:3:17
array2 = 0:-1:-5

length(array1)
length(array2)

% Concatenation
array = [array1, array2]
length(array)
newArray = [array1; array2]

%% Print formatting
% fprintf()
clc
% Format specifiers "%"
% %d - int
% %f - float/double
% %s - string

fprintf("here is an int: %d, here is a float: %f \n", 1, 3.14)

% Specify the number of decimals
fprintf("here is two decimal places: %.2f \n", 3.14)

%with variable
x = 9.9;
fprintf("x = %.1f \n", x)

%% Conditional statements
clc

% booleans
% C++: True, False
% Matlab: true, false
true        % = 1
false       % = 0


%Logic

% not
% C++: "!"
% Matlab: "~"
~true
~false

% "and" and "or"
true && false
true || false

% Comparisons
1 == 2
1 ~= 2
1 >= 2
1 > 2
1 <= 2
1 < 2

% If statements
% C++: {}
% Matlab: "end"

if(true)
   x = 4 
end

% else if
if(false)
    x = 6
elseif(false)
    x = 7
end % Only one end statement

% else
if(x > 7)
   x = 7
elseif(x < 7)
   x = 8
else
   x = 9 
end

%% Loops
clc
% For loop
% C++: 3 parts
% Matlab: different

for k = 1:5
    fprintf("%d ", k)
end
fprintf("\n")

array = [1, 3, 5, 7, 9, 11, 13];
for k = 1:2:length(array)
   fprintf("%d ", array(k)) 
end
fprintf("\n")

sum = 0;
array = ["one", "two", "three"];
for number = array
   fprintf("%s, ", number)
   x = 9;
   sum = sum + 9;
end
fprintf("\n")

% While loops
k = 0;
while(k < 3)
   fprintf("k = %d \n", k)
   k = k + 1; % No shorthand in Matlab :( 
    
end

%% User defined functions
% All functions are in separate files

% y = x^2 +3x + 2
[x1,x2] = quadratic(1, 3, 2)

%% Plotting
close all;
clc

x = 1:100;

% y = x
y1 = x;

% y = 2x
y2 = 2*x;

figure      % Creates a new figure window
plot(x, y1, "g--", "linewidth", 5)
hold on
plot(x, y2)
grid on
title("My plot")
xlabel("x value")
ylabel("y value")
legend("y = x", "y = 2x")

%% Test
clc

['cat' 'dog']
fprintf('The int is %d, the float is %f, the string is %s \n', 33 - 2, "xyz")
xor((1 < 5),~(1 == 1))
x = 2;
e = 2.713;
c = 5 + (x < 2 || (e < pi))

A = [0 1; 1 0];	
B=2;	
C = A + B

A = [2 3 0]; 
B = [4 0 8]; 
C = A.*B
A = [1:5];
A(3) = [];

ans + 5;
6 - ans
3+(5-3)*rand(4,6)

x = 2:3:7