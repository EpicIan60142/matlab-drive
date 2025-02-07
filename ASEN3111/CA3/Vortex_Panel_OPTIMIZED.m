function [CL] = Vortex_Panel_OPTIMIZED(XB,YB,VINF,ALPHA)
    % Overview: This function has been REWRITTEN by myself to heavily
    % optimize it by vectorizing all the operations. I also try my best to
    % reduce RAM usage.
    % Last Modified: October 17 2022


    % (R)am (L)imit (S)ize
    RLS = 900;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Input:                           %
    %                                  %
    % XB  = Boundary Points x-location %
    % YB  = Boundary Points y-location %
    % VINF  = Free-stream Flow Speed   %
    % ALPHA = AOA                      %
    %                                  %
    % Output:                          %
    %                                  %
    % CL = Sectional Lift Coefficient  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%
    % Convert to Radians %
    %%%%%%%%%%%%%%%%%%%%%%
    
    ALPHA = ALPHA*pi/180;
    
    %%%%%%%%%%%%%%%%%%%%%
    % Compute the Chord %
    %%%%%%%%%%%%%%%%%%%%%
    
    CHORD = max(XB)-min(XB);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Determine the Number of Panels %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % obtain length of various arrays
    M = length(XB) - 1;
    MP1 = M+1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Intra-Panel Relationships:                                  %
    %                                                             %
    % Determine the Control Points, Panel Sizes, and Panel Angles %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % add the X and Y together, vectorized
    X = 1/2 * (XB(1:end-1) + XB(2:end));
    Y = 1/2 * (YB(1:end-1) + YB(2:end));
    S = sqrt((XB(2:end) - XB(1:end-1)).^2 + (YB(2:end) - YB(1:end-1)).^2);
    
    % get angle, precomputed, with some funny business indexing
    THETA = atan2(YB(2:end) - YB(1:end-1), XB(2:end) - XB(1:end-1));

    % precompute sine and cosine, as well as RHS
    SINE = sin(THETA);
    COSINE = cos(THETA);
    RHS = sin(THETA - ALPHA);
    
    % get the theta array arranged vertically so the math does good things
    THETA_I = THETA';
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Inter-Panel Relationships:             %
    %                                        %
    % Determine the Integrals between Panels %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % chop off the end
    XB = XB(1:M);
    YB = YB(1:M);
    
    % math -- behave yourself!
    X = X';
    Y = Y';
    
    % RULE: I is horizontal rows are SAME
    % J is vertical colums are SAME
    % THIS RULE IS BROKEN WHEN NECESSARY
    % Basically it just "works" and the more you try to understand, the
    % more you become confused. I certainly am
    
    % calculate re-used calculations, at the expense of more RAM
    X_INTERMEDIATE = X - XB;
    Y_INTERMEDIATE = Y - YB;
    if M > RLS
        clear X Y XB YB;
    end
    SINE_INTERMEDIATE = sin(THETA_I - 2*THETA);
    COSINE_INTERMEDIATE = cos(THETA_I - 2*THETA);

    
    % Get the entire alphabet, et al
    A = -(X_INTERMEDIATE) .* COSINE - (Y_INTERMEDIATE) .* SINE;
    E = (X_INTERMEDIATE) .* SINE - (Y_INTERMEDIATE) .* COSINE;
    if M > RLS
        clear SINE COSINE;
    end
    P = (X_INTERMEDIATE) .* SINE_INTERMEDIATE + (Y_INTERMEDIATE) .* COSINE_INTERMEDIATE;
    Q = (X_INTERMEDIATE) .* COSINE_INTERMEDIATE - (Y_INTERMEDIATE) .* SINE_INTERMEDIATE;
    if M > RLS
        clear SINE_INTERMEDIATE COSINE_INTERMEDIATE;
    end
    B = (X_INTERMEDIATE).^2 + (Y_INTERMEDIATE).^2;
    if M > RLS
        clear X_INTERMEDIATE Y_INTERMEDIATE;
    end
    C = sin(THETA_I - THETA);
    D = cos(THETA_I - THETA);
    F = log(1 + S .* (S + 2 * A) ./ B);
    G = atan2(E .* S, B + A .* S);
    if M > RLS
        clear B THETA_I;
    end
    
    % diagonal logical index. This is fast!
    DI = eye(M, 'logical');
    
    CN2 = D + 1/2 * Q .* F ./ S - (A .* C + D .* E) .* G ./ S;
    CN2(DI) = 1;
    if M > RLS
        clear Q;
    end
    CN1 = 1/2 * D .* F + C .* G - CN2;
    CN1(DI) = -1;

    CT2 = C + 1/2 * P .* F ./ S + (A .* D - C .* E) .* G ./ S;
    CT2(DI) = pi/2;
    if M > RLS
        clear P A E;
    end
    CT1 = 1/2 * C .* F - D .* G - CT2;
    CT1(DI) = pi/2;
    if M > RLS
        clear C D F G DI;
    end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Inter-Panel Relationships:           %
    %                                      %
    % Determine the Influence Coefficients %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % initialize
    AN = zeros(MP1);
    AT = zeros(M, MP1);

    % Do NOT vectorize. This code is slow. Vectorized is slower. I tried, and I failed
    for J = 2:M
        AN(1:M, J) = CN1(:, J) + CN2(:, J-1);
        AT(1:M, J) = CT1(:, J) + CT2(:, J-1);
    end
    
    % extreme funny business
    AN(1:end-1, 1) = CN1(:, 1);
    AN(1:end-1, MP1) = CN2(:, M);
    AN(MP1,1) = 1;
    AN(MP1,MP1) = 1;
    AN(MP1, 2:M) = 0;
    
    AT(:, 1) = CT1(:, 1);
    AT(:, MP1) = CT2(:, M);
    
    % put in a zero for good measure
    RHS(MP1) = 0.0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % Solve for the gammas %
    %%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % thank you, matlab
    GAMA = linsolve(AN, RHS');
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solve for Tangential Veloity and Coefficient of Pressure %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % basic velocity stuff
    V = cos(THETA - ALPHA);
    V = V + sum(AT' .* GAMA, 1);
    
    %CP = 1 - V.^2; % not sure what this is used for??
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Solve for Sectional Coefficient of Lift %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % obtain circulation, and done!
    CIRCULATION = sum(S.*V);
    CL = 2*CIRCULATION/CHORD;
end