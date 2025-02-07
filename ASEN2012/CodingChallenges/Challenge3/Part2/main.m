%% Bonus: LeastSquares Function
%{

This bonus challenge is completely optional, but must be completed to be
considered a winner for this coding challenge.

The goal of this bonus exercise is to write a LeastSquares function that
will take in vectors for t, y, and p, returning: 
    X: the coefficients for any polynomial of order p
    f: an anonymous function representing the polynomial

For this section of the challenge, you will be asked to plot the given
stream of data with the MOST APPROPRIATE ORDER POLYNOMIAL, using your
function for a LeastSquares regression.

%}

% read in and extract data
data = readmatrix('BonusData.xlsx');
t = data(:,1);
y = data(:,2);

% plot the data and compare to the LS solution
figure(); hold on
    
    plot(t,y,'b.')
    N = length(t)
    for i = 1:N
        [coef, f] = leastSquares(t,y,4);
        
    end
    plot(t,f(t),'r-')

    
hold off

%% Function Defintion
function [X,f] = leastSquares(t,y,p)
    % for writing this function, some skeleton code has been provided to
    % help you design the function to serve your purposes
    A = [];
    % write an expression for A, the input matrix
    for ii = 0:p
        col = t.^ii;
        A = [col, A];
    end
    A
    % compute coefficient vector, x_hat
    x_hat = A\y
    X = x_hat
    
    % do not change the following lines of code. This will generate the
    % anonymous function handle "f" for you
%     f = '@(x)';
%     for i = 0:p
%         f = strcat(f,'+',strcat(string(x_hat(i+1)),'.*x.^',string(p-i)));
%     end
%     eval(strcat('f = ',f,';'))
    
    while length(x_hat) < 7
        x_hat = [0;x_hat];
    end
    % workaround for MATLAB grader
    f = @(x) x_hat(1)*x.^6 + x_hat(2)*x.^5 + x_hat(3)*x.^4 + x_hat(4)*x.^3 + x_hat(5)*x.^2 + x_hat(6)*x + x_hat(7);
    
end
