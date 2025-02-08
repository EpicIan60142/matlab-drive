function [Q,R] = QRFac(key)
%Finds the QR factorization of a given key (from pg. 206 in the textbook)
%   Inputs: the key
%   Outputs: Q, R

[m,n] = size(key); %Find the key's size
R0 = zeros(n); %Reserve an n x n matrix for R

for j=1:m %Loop through each row of the key
    R0(j,j)=norm(key(:,j)); %Find the norm of the jth column vector of the key
    if (R0(j,j) == 0) %If the norm is 0
        fprintf('Key has linearly dependent columns, Q not possible');
        break; %Break out of the algorithm
    else
        for i = 1:n %Loop through each column of the key
            key(i,j)=key(i,j)/R0(j,j); %Calculate the unit vector of each key element in the row
        end
    end
    
    for k = j+1:n %Loop through each column vector, skipping the previous column (ex. if j = 1, start at column 2)
        R0(j,k)=key(:,j)'*key(:,k); %Calculate each R matrix element
        for i=1:n %Loop through each column of the key
            key(i,k) = key(i,k)-key(i,j)*R0(j,k); %Apply the Gram-Schmidt algorithm to calculate orthonormal basis vectors
        end
    end
        
end

Q = key; %Assign the Q matrix
R = R0; %Assign the R matrix

end

