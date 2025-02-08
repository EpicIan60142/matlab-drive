function [encodedMessage] = encode(inputMessage,key, modulo)
%Encodes an input message using a given key
%   Inputs: Message (string), key (square matrix), modulus (integer)
%   Outputs: Encoded message (string)

key %Display the current key
[m,n] = size(key); %Find the size of the key
if (m ~= n) %If the key isn't square
    fprintf("Key isn't square! Need to have a square key. \n");
    return
else
    fprintf("Key is square! Yay! \n");
end

inputASCII = double(inputMessage); %Convert string into ASCII codes
messageLength = length(inputASCII); %Find how long the original message is

while (mod(messageLength,n) ~= 0) %If the message length is not a multiple of the columns of the key
    inputASCII = append(char(inputASCII), char(modulo+65-1)); %Fill in the message with the last character of the alphabet
    messageLength = length(inputASCII); %Update the message length to avoid an infinite loop
end

fprintf("Input message: %s \n",char(inputASCII)) %Display the input message after filling

input = inputASCII - 65; %Convert ASCII into numbers of alphabet according 
                         %to the Hill Cipher. In our case, that's
                         %A B C D E F G H I J  K  L  M  N  O  P  Q  R  S 
                         % T  U  V  W  X  Y  Z  [  /  ]  (^  _) <- ^ and _ part of the Beylkinaut-Curryist alphabet
                         %into
                         %0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 
                         %19 20 21 22 23 24 25 26 27 28 (29 30)

A = reshape(input, n, []); %Convert the input message into a column vector of size n x 1

encoded = key*A; %Apply the key to the input message
encoded = mod(encoded, modulo); %Convert the values to their modular counterparts
encoded = reshape(encoded, 1, messageLength); %Convert the encoded column vector to a row vector

encodedASCII = encoded + 65; %Restore letters of the alphabet to their ASCII codes

encodedMessage = char(encodedASCII); %Convert ASCII codes back to characters

fprintf("Encoded message: %s \n",encodedMessage); %Display the encoded message

end

