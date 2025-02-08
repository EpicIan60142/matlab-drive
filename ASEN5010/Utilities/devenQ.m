function BN = devenQ(Vb, Vn, w)
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
    
    % Solve eigenvalue problem
    eigVals = eig(K);
    [betas, ~] = eig(K);
    
    % Find max eigenvalue
    [~, orientIdx] = max(eigVals, [], 'all');
    
    % Choose orientation with max eigenvalue and make DCM
    optOrient = betas(:, orientIdx);
    
    BN = EP2DCM(optOrient);

end