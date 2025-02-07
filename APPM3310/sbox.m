% S - Box 
% Kellen Martin 

close all; clear; clc;

modulo = 29;

% Creates a random S box 
S = randperm(modulo)-1

% Random Message 
key1 = [2 3 ; 1 4];
message1 = 'LINEARALGEBRA]';

message = message1;
key = key1;

inputASCII = double(message); %Convert string into ASCII
messageLength = length(inputASCII);



fprintf("Input message: %s \n",char(inputASCII)) 

input = inputASCII - 65;
% Encodes the Plaintext using the s box 
Scypher = zeros(1,messageLength);

for n = 1:messageLength

    Scypher(n) = S(input(n)+1);

end

sASCII = Scypher+64;

sMessage = char(sASCII);



encodedMessage = encode(sMessage, key, modulo);
%fprintf("Encoded message: %s \n",encodedMessage);

% Decodes the S Cypher message 

boxed = decode(encodedMessage, key, modulo);

unboxASCII = double(boxed); %Convert string into ASCII
box = unboxASCII - 65;


unboxed = zeros(1,messageLength);
for n = 1:messageLength
   check = box(n)+1;
   unboxed(n) = find(S == check);
end

decodedASCII = unboxed+64;

decoded = char(decodedASCII)








