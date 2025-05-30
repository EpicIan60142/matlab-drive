%% Given Constants
L = 12; %[in]
L_e = 4.5; %[in]
L_r = 5; %[in]
w = 1; %[in]
h = 1/8; %[in]
h_e = 1/4; %[in]
h_r = .040 ; %[in]
E = 10175000; %[psi]
rho = 0.0002505; %[lb-sec2/in4]
M_t = 1.131*rho; %[mass of assembly]
S_t = 0.5655*rho; %[first mass-moment of assembly]
I_t = 23.124*rho; %[seccond mass-moment of assembly]

%% Two Element Model
A = w*h;
I_zz = w*h^3/12;
C_M2 = rho*A*L/100800;
C_K2 = 4*E*I_zz/L^3;

M_2  = C_M2.*[19272 1458*L 5928 -642*L 0 0 ; 1458*L 172*L^2 642*L -73*L^2 0 0 ; 5928 642*L 38544 0 5928 -642*L ; -642*L -73*L^2 0 344*L^2 642*L -73*L^2 ; 0 0 5928 642*L 119272 -1458*L ; 0 0 -642*L -73*L^2 -1458*L 172*L^2] + [ 0 0 0 0 0 0 ;  0 0 0 0 0 0 ;  0 0 0 0 0 0 ; 0 0 0 0 0 0 ;  0 0 0 0 M_t S_t ;  0 0 0 0 S_t I_t ;];
M_2_hat = M_2(3:6,3:6);
K_2 = C_K2.*[24 6*L -24 6*L 0 0 ; 6*L 2*L^2 -6*L L^2 0 0; -24 -6*L 48 0 -24 6*L ; 6*L L^2 0 4*L^2 -6*L L^2 ; 0 0 -24 -6*L 24 -6*L; 0 0 6*L L^2 -6*L 2*L^2];
K_2_hat = K_2(3:6,3:6);
%% Four Element Model
C_M4 = rho*A*L/806400 ;
C_K4 = 8*E*I_zz/(L^3);

M_4 = C_M4.*[ 77088 2916*L 23712 -1284*L 0 0 0 0 0 0 ; 2916*L 172*L^2 1284*L -73*L^2 0 0 0 0 0 0; 23712 1284*L 154176 0 23712 -1284*L 0 0 0 0; -1284*L -73*L^2 0 344*L^2 1284*L -73*L^2 0 0 0 0 ;  0 0 23712 1284*L 154176 0 23712 -1284*L 0 0 ; 0 0 -1284*L -73*L^2 0 344*L^2 1284*L -73*L^2 0 0 ; 0 0 0 0 23712 1284*L 154176 0 23712 -1284*L ; 0 0 0 0 -1284*L -73*L^2 0 344*L^2 1284*L -73*L^2 ; 0 0 0 0 0 0 23712 1284*L 77088 -2916*L ; 0 0 0 0 0 0 -1284*L -73*L^2 -2916*L 172*L^2] + [ 0 0 0 0 0 0 0 0 0 0 ;  0 0 0 0 0 0 0 0 0 0 ;  0 0 0 0 0 0 0 0 0 0 ;  0 0 0 0 0 0 0 0 0 0 ;  0 0 0 0 0 0 0 0 0 0 ;  0 0 0 0 0 0 0 0 0 0 ;  0 0 0 0 0 0 0 0 0 0 ;  0 0 0 0 0 0 0 0 0 0 ;  0 0 0 0 0 0 0 0 M_t S_t ;  0 0 0 0 0 0 0 0 S_t I_t ];
M_4_hat = M_4(3:10 , 3:10);
K_4 = C_K4.*[96 12*L -96 12*L 0 0 0 0 0 0 ; 12*L 2*L^2 -12*L L^2 0 0 0 0 0 0 ; -96 -12*L 192 0 -96 12*L 0 0 0 0 ; 12*L L^2 0 4*L^2 -12*L L^2 0 0 0 0; 0 0 -96 -12*L 192 0 -96 12*L 0 0 ; 0 0 12*L L^2 0 4*L^2 -12*L L^2 0 0 ; 0 0 0 0 -96 -12*L 192 0 -96 12*L ; 0 0 0 0 12*L L^2 0 4*L^2 -12*L L^2 ; 0 0 0 0 0 0 -96 -12*L 96 -12*L ; 0 0 0 0 0 0 12*L L^2 -12*L 2*L^2];
K_4_hat =  K_4(3:6,3:6);

K_4_hat*U = w^2*M_4_hat*U


