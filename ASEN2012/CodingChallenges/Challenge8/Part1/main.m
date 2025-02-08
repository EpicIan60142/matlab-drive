%% Ian Faber
% SID: 108577813
% ASEN 2012
% Coding Challenge 8
% Last modified: 11/11/2021

%% Purpose
% The purpose of this code is to find the roots of a function using both mathematical approximation and MATLAB functions.

%% Housekeeping
clc
close all


%% Define error tolerance
s = 4; % number of decimals we want precision to
epsilon = 10^-(s+1); % tolerance (same for all methods)


%% Define constants
alpha = 0.008;
beta = 0.0002;
gamma = -0.26;

constants = [alpha beta gamma];

%% Define the function we want to set equal to zero and take the roots of
f = @(rho)((gamma*rho)/(1-(beta*rho)) - alpha*(rho).^2); % original function of energy density rho
f_prime = @(rho) (gamma*(1-(beta*rho)) - (-beta*gamma*rho))/(1-beta*rho)^2 - 2*alpha*rho; % Derivative of function with respect to rho


%% Determine a reasonable starting point for your recursive formulas
guess = [-32, 0]; % Initial guess(es) [hint: you should have as many distinct initial guesses you have roots - how can you approximate the values of these roots?]


%% Call function to calculate roots once per guess


% Reassign the outputs of your function  to the following variables
%rho_b: Roots for bisection method - should be a 2x1 column vector with both of your final roots, reported in ascending order
%rho_ss: Roots for successive substitution method - should be a 2x1 column vector with both of your final roots, reported in ascending order
%rho_n: Roots for Newton's method - should be a 2x1 column vector with both of your final roots, reported in ascending order
%rho_s: Roots for secant method - should be a 2x1 column vector with both of your final roots, reported in ascending order

[rho_b, rho_ss, rho_n, rho_s, rho_check, runtimes, iterations] = rootsfunc(epsilon, f, f_prime, guess, constants)

least = find(iterations == min(iterations));
fastest = find(runtimes == min(runtimes));

switch least
    case 1
        fprintf("Bisection Method ran the least iterations, which was %d!\n",iterations(least));
    case 2
        fprintf("Successive Substitution Method ran the least iterations, which was %d!\n",iterations(least));
    case 3
        fprintf("Newton's Method ran the least iterations, which was %d!\n",iterations(least));    
    case 4
        fprintf("Secant Method ran the least iterations, which was %d!\n",iterations(least));
end

answer1 = ['Newton''s method always converges the quickest, and does not change based on my initial guess. \n' ...
           'The number of iterations will only change depending on the function we are trying to find the root for.\n\n'];

fprintf(answer1);

switch fastest
    case 1
        fprintf("Bisection Method ran the fastest with a runtime of %f seconds!\n",runtimes(fastest));
    case 2
        fprintf("Successive Substitution Method ran the fastest with a runtime of %f seconds!\n",runtimes(fastest));
    case 3
        fprintf("Newton's method ran the fastest with a runtime of %f seconds!\n",runtimes(fastest));
    case 4
        fprintf("Secant Method ran the fastest with a runtime of %f seconds!\n",runtimes(fastest));
    case 5
        fprintf("fzero ran the fastest with a runtime of %f seconds!\n",runtimes(fastest));
    otherwise
        fprintf("All methods were equally fast!\n"); 
end

answer2 = ['Since the time each method takes to run appears to change every time I run the code, I can''t \n' ...
        'be confident about which method is the best without running a large number of simulations and \n' ...
        'seeing which methods tend to be the fastest more often. However, in general the closer the \n' ...
        'initial guess is to the actual root, the faster the method seems to run. I was expecting the number of\n' ...
        'iterations and runtime to match up, but this doesn''t seem to be the case. If I had to guess, the reason\n' ...
        'for this would be due to machine error and uncertainty, since computers aren''t perfect and are prone\n' ...
        'to random fluctuations that affect computation time\n\n'];

fprintf(answer2)

%% End of main script; start of sub-function
function [rho_b, rho_ss, rho_n, rho_s, rho_check, runtimes, iterations] = rootsfunc(epsilon, f, f_prime, guess, constants) % Define function to calculate roots for all methods
%% Note: since you want to check convergence iteratively, what type of loop might be useful to you for these methods?

alpha = constants(1);
beta = constants(2);
gamma = constants(3);

% Calculate how long each method takes to complete 
    % To do this, it is useful to note that "tic" starts a timer in MATLAB at the line where it is used.
    % Similarly, "toc" stops the timer wherever it appears. If you would like to report the total time elapsed 
        % between the two, you can assign "toc" to a variable.
        

        
%% Bisection method - the general process for this can be seen in slide 15 of lecture 20
% There is a wide range of initial bounds you could use, but if you would like some guidance on a starting point, 
% -50 to -15 might be a good range for your lowest (negative) root to fall in, and -10 to 10 might be good to check for the larger root.    
tic;

rho_b = zeros(2,1);
startBounds = [-50, -15; -10, 10];
iteration = 0;
for i = 1:2
   startInterval = startBounds(i,:);
   a = startInterval(1);
   b = startInterval(2);
   difference = b-a;
   while(abs(difference) > epsilon)
      c = (a + b)/2;
      if(sign(f(c)) ~= sign(f(a)))
          a = a;
          b = c;
          difference = b-a;
      elseif(sign(f(c)) ~= sign(f(b)))
          a = c;
          b = b;
          difference = b-a;
      end
      iteration = iteration + 1;
   end
   rho_b(i) = c;
end

iterations_b = iteration;
end_b = toc;

%% Successive substitution - the general process for this can be seen in slide 6 of lecture 21
% Manipulate original f(x) function to get x = g(x) [x, in this case, is rho (independent variable)]
tic;

rho_ss = zeros(2,1);
iteration = 0;
for i = 1:2
    x0 = guess(i);
    if i == 2
    % if you want the second-lowest root, solve for the rho in the numerator of the first term
        g = @(rho) 0;
    else
    % if you want the lowest (negative) root, use this equation: 
        g = @(rho) gamma / (alpha * (1 - (beta * rho)));
    end
    
    x1 = g(x0);
    while(abs(x1 - x0) >= epsilon)
        x0 = x1;
        x1 = g(x0);
        iteration = iteration + 1;
    end
    
    rho_ss(i) = x1;
end

iterations_ss = iteration;
end_ss = toc;

%% Newton's method - the general process for this can be seen on slide 9 of lecture 21
tic;

rho_n = zeros(2,1);
iteration = 0;
for i = 1:2
    x0 = guess(i);
    x1 = x0 - (f(x0)/f_prime(x0));
    while(abs(x1 - x0) > epsilon)
       x0 = x1; 
       x1 = x0 - (f(x0)/f_prime(x0));
       iteration = iteration + 1;
    end
    
    rho_n(i) = x1;
end

iterations_n = iteration;
end_n = toc;

%% Secant method - the general process for this can be seen on slide 14 of lecture 21
tic;

rho_s = zeros(2,1);
iteration = 0;
for i = 1:2
   x0 = guess(i);
   x1 = guess(i) - 1;
   x2 = x1 - ((f(x1)*(x1-x0))/(f(x1)-f(x0)));
   while(abs(x2 - x1) > epsilon)
       x0 = x1;
       x1 = x2;
       x2 = x1 - ((f(x1)*(x1-x0))/(f(x1)-f(x0)));
       iteration = iteration + 1;
   end
    
   rho_s(i) = x2; 
end

iterations_s = iteration;
end_s = toc;

%% Check using MATLAB built-in functions: fzero(), roots(), or fsolve()
% Don't forget "doc" and "help" are useful tools to find syntax of built-in functions
tic;

rho_check = zeros(2,1);
for i = 1:2
    rho_check(i) = fzero(f, [guess(i)-1, guess(i)+1]);
end
rho_check

end_check = toc;

runtimes = [end_b; end_ss; end_n; end_s; end_check];
iterations = [iterations_b; iterations_ss; iterations_n; iterations_s];

end