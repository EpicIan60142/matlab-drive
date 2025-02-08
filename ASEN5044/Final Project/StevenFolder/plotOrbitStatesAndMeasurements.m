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
        subplot(411)
        hold on
        plot(thist, yhist(3*i-2,:),"x")
        subplot(412)
        hold on
        plot(thist, yhist(3*i-1,:),"x")
        subplot(413)
        hold on
        plot(thist, yhist(3*i,:),"x")
        subplot(414)
        hold on
        timeSeen = thist(~isnan(yhist(3*i,:)));
        plot(timeSeen, i*ones(length(timeSeen), 1), "x")
    end
end