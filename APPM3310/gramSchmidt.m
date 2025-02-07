% author: mostly stolen from online, so it doesn't really matter

% replace A with whatever vectors form the given basis. (in transpose)
A = [[-1 1 2]',[-1 -1 1]',[0 1 3]'];
n = rank(A);
Q = zeros(n);
R = Q;

for j=1:n
    v = A(:,j);
    for i=1:j-1
        R(i,j)=Q(:,j)'*A(:,j);
        v = v-R(i,j)*Q(:,i);
    end
    R(j,j)=norm(v);
    Q(:,j)=v/R(j,j);
end

Q