%% APPM 3310 Project: Improving the Hill Cipher
    % Code by Ian Faber
    % Teammates: Emanuel Quintana, Kellen Martin
    
clc; clear; close all; format;

charOffset = 65;

key1 = [2 3 ; 1 4];
message1 = 'LINEARALGEBRA';

key2 = ['IANEM';'ANUEL';'KELLE';'NIANE';'MANUE'] - charOffset;
message2 = 'ATTACKATDAWNNOSURVIVORS'; %ALL CAPS or all lowercase

key3 = [7 8 ; 11 11]; %"HILL"
message3 = 'SHORTEXAMPLE';

key4 = ['MA';'TH'] - charOffset;
message4 = 'GROOMSTER_BEING_DISRESPECTFUL^RECOMMEND_ACTION';
message5 = 'UNDERSTOOD^ATTACK_AT_DAWN';

message = message4; %Choose a message from above to encode
key = key2; %Choose a key from above to encode with
modulo = 31; %Choose a modulus to work in

encodedMessage = encode(message, key, modulo); %Encode the message with the chosen key in a specific modulus
decode(encodedMessage, key, modulo); %Decode the encoded message using the key in a specific modulus

[keyQ, keyR] = QRFac(key) %Find the QR factorization of the key

encodedMessage = encode(message, keyQ*keyR, modulo); %Encode the message with the QR factorization of the key, only works if both matrices are known
decode(encodedMessage, keyQ*keyR, modulo); %Decode the encoded message with the QR key, only works when both matrices are known



