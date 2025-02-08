function [decodedMessage] = decode(encodedMessage,key,modulo)
%Decodes an encoded message, as long as the key is known
%   Inputs: Encoded message (string), key (square matrix), modulus (integer)
%   Outputs: Decoded message (string)

[m,n] = size(key); %Find the size of the key
if (m ~= n) %If the key isn't square
    fprintf("Key isn't square! Need to have a square key. \n");
    return
else
    fprintf("Key is square! Yay! \n");
end

inverseKey = inv(key); %Take the normal inverse of the key
[numerator, denomenator] = numden(sym(inverseKey)); %Find the numerator and denomenator of each element in the key, as they are likely fractions

%Find the modular inverse for each value in the inverse key matrix by converting fractional elements to integers in the modulus
for k = 1:m %Loop through each row
    for l = 1:n %Loop through each column
        while true %Loop forever
            [currentNumerator, currentDenomenator] = numden(sym(inverseKey(k,l))); %Find the current numerator and denominator of an inverse key element
            if(currentDenomenator == 1) %If the element is no longer a fraction
                break; %Break out of the while loop
            elseif(inverseKey(k,1) >= modulo) %If the frational element never becomes an integer
                fprintf("Modular inverse not possible at %f %f, matrix element %f and modulo %f not coprime \n",k,l,inverseKey(k,l),modulo);
                break; %Break out of the while loop
            elseif(currentDenomenator >= 999999999999999999999) %If the denomenator is incredibly large (compensate for numerical instability and round-off error)
                inverseKey(k,l) = 0; %Set the element to 0
                break; %Break out of the while loop
            end
            inverseKey(k,l) = inverseKey(k,l) + (modulo/denomenator(k,l)); %Increment the element by a constant amount beased on the modulus and original denominator
        end
    end
end

inverseKey %Display the modular inverse

encodedASCII = double(encodedMessage); %Convert characters into their ASCII codes
messageLength = length(encodedASCII); %Find the length of the encoded message

encoded = encodedASCII - 65; %COnvert ASCII codes into their corresponding alphabet numbers

A = reshape(encoded, n, []); %Change the input message into a column vector of size n x 1

decoded = inverseKey*A; %Apply the modular inverse key to the encoded message
decoded = mod(decoded, modulo); %Convert the calculated values to their modular counterparts
decoded = reshape(decoded, 1, messageLength); %Convert the column vector to a row vector

decodedASCII = decoded + 65; %Convert letters of the alphabet to their ASCII codes

decodedMessage = char(decodedASCII); %Convert ASCII codes back to characters

fprintf("Decoded message: %s \n",decodedMessage); %Display the decoded message

end



