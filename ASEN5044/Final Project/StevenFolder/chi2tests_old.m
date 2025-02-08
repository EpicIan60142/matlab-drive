function [epsilon_x_bar, epsilon_y_bar, NEESbounds, NISbounds, filterData] = chi2tests(filterFunc, N, alpha, filterParams, Qtrue, Rtrue, filterInit, timesteps, x0)
    
    if timesteps > 533
        figure
        hold on; grid on; axis equal;
        title("All Monte Carlo runs")
        points = [55, 102, 149, 196, 243, 290, 346, 391, 437, 485, 533];
        % diff(points)
        colors = colororder('glow12');
    end

    epsilon_mat_x = zeros(N, timesteps);
    epsilon_mat_y = zeros(N, timesteps);
    NISbounds = zeros(2,timesteps);
    for i =1:N
        [xtrue, ytrue] = orbitMonteCarlo(Qtrue, Rtrue, filterInit.dt, timesteps, x0);
    
        filterData = filterFunc(ytrue, filterParams, filterInit);
        state_error = xtrue - filterData.x_KF;
        
        if timesteps > 533
            plot(filterData.x_KF(1,:), filterData.x_KF(3,:));
            % plot(xtrue(1,:), xtrue(3,:), '-.')
            if i == 1
                start = plot(filterData.x_KF(1,1), filterData.x_KF(3,1), 'g.', 'MarkerSize', 15);
                spikes = plot(filterData.x_KF(1,points), filterData.x_KF(3,points), 'k.', 'MarkerSize', 15);
            else
                plot(filterData.x_KF(1,1), filterData.x_KF(3,1), 'g.', 'MarkerSize', 15);
                plot(filterData.x_KF(1,points), filterData.x_KF(3,points), 'k.', 'MarkerSize', 15);
                for kk = points
                    plot(filterData.stations(filterData.IDs_KF{kk}).x(kk), filterData.stations(filterData.IDs_KF{kk}).y(kk), 'Color', colors(filterData.IDs_KF{kk},:), 'Marker', '.', 'MarkerSize', 15)
                end
            end
        end

        for j = 1:timesteps
            if(isempty(filterData.innovation_KF{j}))
                epsilon_mat_x(i,j) = nan;
                epsilon_mat_y(i,j) = nan;
            else
                epsilon_mat_x(i,j) = state_error(:,j)' * inv(filterData.P_KF(:,:,j)) * state_error(:,j);
                epsilon_mat_y(i,j) = filterData.innovation_KF{j}' * inv(filterData.Sk_KF{j}) * filterData.innovation_KF{j};
            end
        end
        
    end
    epsilon_x_bar = mean(epsilon_mat_x);
    epsilon_y_bar = mean(epsilon_mat_y);

    NEESbounds = [chi2inv(alpha/2, N*length(filterInit.dx0))./N, chi2inv(1-alpha/2, N*length(filterInit.dx0))./N];

    for i = 1:length(NISbounds)
        NISbounds(:,i) =  [chi2inv(alpha/2, N*length(filterData.innovation_KF{i}))./N; chi2inv(1-alpha/2, N*length(filterData.innovation_KF{i}))./N];
    end

    if timesteps > 533
        nom = plot(filterInit.xNom(1,:), filterInit.xNom(3,:), 'k--');
        legend([start, spikes, nom], ["Start point", "NEES Spike locations", "Nominal solution"])
    end

end