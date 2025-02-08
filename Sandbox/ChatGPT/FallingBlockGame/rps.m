% Rock-Paper-Scissors game in MATLAB

% Get user input
user_choice = input('Enter your choice (1 = Rock, 2 = Paper, 3 = Scissors): ');

% Generate computer's choice
comp_choice = randi([1 3]);

% Print computer's choice
if comp_choice == 1
    disp('Computer chooses Rock.');
elseif comp_choice == 2
    disp('Computer chooses Paper.');
else
    disp('Computer chooses Scissors.');
end

% Determine the winner
if user_choice == comp_choice
    disp('It is a tie!');
elseif user_choice == 1 && comp_choice == 3
    disp('You win! Rock beats Scissors.');
elseif user_choice == 2 && comp_choice == 1
    disp('You win! Paper beats Rock.');
elseif user_choice == 3 && comp_choice == 2
    disp('You win! Scissors beats Paper.');
else
    disp('Computer wins!');
end