function plotOrbitStatesAndMeasurements(thist, state_hist, yhist)
    figure()
    subplot(411)
    plot(thist, state_hist(1,:))
    
    subplot(412)
    plot(thist, state_hist(2,:))
    
    subplot(413)
    plot(thist, state_hist(3,:))
    
    subplot(414)
    plot(thist, state_hist(4,:))
    
    figure()
    for i=1:12
        subplot(311)
        hold on
        plot(thist, yhist(3*i-2,:),"X")
        subplot(312)
        hold on
        plot(thist, yhist(3*i-1,:),"o")
        subplot(313)
        hold on
        plot(thist, yhist(3*i,:),"square")
    end
end