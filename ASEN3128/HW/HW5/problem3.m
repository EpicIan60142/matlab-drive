%% Housekeeping
clc; clear; close all;

%% Original matrix, eigenvalues, and eigenvectors
A = [
        -0.02, 0.016, -0.65, -32.17;
        -0.13, -1.019, 454.21, 0;
        0, -0.005, -1.38, 0;
        0, 0, 1, 0
    ];

[V, D] = eig(A);

%% Phasor Setup
size = 25;
col = [[1 0 0];[0 1 0];[0 0 1];[1 0 1]];
variable = ["u", "\alpha", "q", "\theta"];

%% Short Period
VShort = V(:,1);
VShort = VShort/VShort(end)
realShort = real(VShort);
imagShort = imag(VShort);

figure
hold on
grid on
title("Short Period Mode Phasor Plot")
xline(0, 'k--')
yline(0, 'k--')
scatter(realShort, imagShort, size, col, 'filled')
for k = 1:length(VShort)
    shortPlot(k) = plot([0,realShort(k)], [0,imagShort(k)], 'Color', col(k,:));
    if k == 4
        label(k) = sprintf("\\Delta %s = 1", variable(k));
    else
        label(k) = sprintf("\\Delta %s", variable(k));
    end
end

subset = shortPlot;
legend(subset, label, 'Location', 'best')

%% Phugoid
VPhugoid = V(:,3);
VPhugoid = VPhugoid/VPhugoid(end)
realPhugoid = real(VPhugoid);
imagPhugoid = imag(VPhugoid);

figure
hold on
grid on
title("Phugoid Mode Phasor Plot")
xline(0, 'k--')
yline(0, 'k--')
scatter(realPhugoid, imagPhugoid, size, col, 'filled')
for k = 1:length(VPhugoid)
    phugoidPlot(k) = plot([0,realPhugoid(k)], [0,imagPhugoid(k)], 'Color', col(k,:));
    if k == 4
        label(k) = sprintf("\\Delta %s = 1", variable(k));
    else
        label(k) = sprintf("\\Delta %s", variable(k));
    end
end

subset = phugoidPlot;
legend(subset, label, 'Location', 'best')

%% Augmented matrix, eigenvalues, and eigenvectors
A_fp = [
            -0.02, 0.016, -0.65, -32.17, 0, 0;
            -0.13, -1.019, 454.21, 0, 0, 0;
            0, -0.005, -1.38, 0, 0, 0;
            0, 0, 1, 0, 0, 0;
            1, 0, 0, 0, 0, 0;
            0, 1, 0, -502, 0, 0
        ];

[V, D] = eig(A_fp)
