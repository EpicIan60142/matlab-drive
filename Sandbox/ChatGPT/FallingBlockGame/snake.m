% Snake Game in MATLAB

% set up game parameters
board_size = [20, 20];  % board size
snake_init = [10, 10];  % initial position of snake
food = [0, 0];  % food location

% initialize the board with the snake and food
board = zeros(board_size);
board(snake_init(1), snake_init(2)) = 1;
place_food();

% set up game loop
game_over = false;
direction = 'right';
while ~game_over
    % display board
    clf;
    imagesc(board);
    colormap(gray);
    axis square;
    set(gca, 'XTick', [], 'YTick', []);
    drawnow;
    
    % move snake
    switch direction
        case 'up'
            snake_init(1) = snake_init(1) - 1;
        case 'down'
            snake_init(1) = snake_init(1) + 1;
        case 'left'
            snake_init(2) = snake_init(2) - 1;
        case 'right'
            snake_init(2) = snake_init(2) + 1;
    end
    
    % check for game over
    if snake_init(1) < 1 || snake_init(1) > board_size(1) || ...
            snake_init(2) < 1 || snake_init(2) > board_size(2) || ...
            board(snake_init(1), snake_init(2)) == 1
        game_over = true;
        continue;
    end
    
    % update board
    if isequal(snake_init, food)
        place_food();
    else
        tail = find(board == max(board(:)));
        board(tail(1)) = 0;
    end
    board(snake_init(1), snake_init(2)) = max(board(:)) + 1;
    
    % get input for direction change
    [keyIsDown, ~, keyCode] = KbCheck;
    if keyIsDown
        switch KbName(keyCode)
            case 'uparrow'
                if ~strcmp(direction, 'down')
                    direction = 'up';
                end
            case 'downarrow'
                if ~strcmp(direction, 'up')
                    direction = 'down';
                end
            case 'leftarrow'
                if ~strcmp(direction, 'right')
                    direction = 'left';
                end
            case 'rightarrow'
                if ~strcmp(direction, 'left')
                    direction = 'right';
                end
        end
    end
    
    % pause for animation
    pause(0.1);
end

% display game over message
msgbox('Game over!');

% function to place food randomly on the board
function place_food()
    global board food;
    food = randi(size(board));
    while board(food(1), food(2)) ~= 0
        food = randi(size(board));
    end
    board(food(1), food(2)) = -1;
end