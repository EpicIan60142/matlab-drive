function BN = QUEST(Vb, Vn, w)
    % Vb is a list of body vectors, Vn is a list of inertial vectors, and w is
    % their weights

    % Normalize unit vectors
    for k = 1:size(Vb,2)
        Vb(:,k) = Vb(:,k)/norm(Vb(:,k));
    end
    
    for k = 1:size(Vn,2)
        Vn(:,k) = Vn(:,k)/norm(Vn(:,k));
    end
    
    % Define B, S, sigma, and Z
    B = zeros(3,3);
    
    for k = 1:length(w)
        mat = w(k)*Vb(:,k)*Vn(:,k)';
        B = B + mat;
    end
    S = B + B';
    sig = trace(B);
    Z = [
            B(2,3) - B(3,2);
            B(3,1) - B(1,3);
            B(1,2) - B(2,1);
        ];
    
    % Define K
    K = [
            sig, Z';
            Z, S - sig*eye(3)
        ];
    
    % Find optimal eigenvalue
    eig0 = sum(w);
    
    epsilon = 0.000001;
    steps = 10;
    f = @(s) det(K-s*eye(4));
    fPrime = @(s) (f(s+epsilon) - f(s))/epsilon; % Finite difference quotient
    
    eig = eig0;
    for k = 1:steps
        eigOpt = eig - f(eig)/fPrime(eig);
    end
    
    % Find CRP
    q = (((eigOpt + sig)*eye(3) - S)^-1)*Z;
    
    % Convert to EP's and DCM
    beta = (1/sqrt(1+q'*q))*[1;q];
    BN = EP2DCM(beta);

end