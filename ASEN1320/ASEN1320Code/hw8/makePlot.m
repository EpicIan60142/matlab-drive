function makePlot(timeVector,StateMatrix)

vxVector = StateMatrix(:,1);    
vyVector = StateMatrix(:,2);   
xVector  = StateMatrix(:,3); 
yVector  = StateMatrix(:,4);

figure(1)
plot(xVector,yVector) 
xlabel('Horizontal Position (m)')
ylabel('Vertical Position (m)') 

figure(2)
plot(timeVector,yVector,timeVector,vyVector) 
xlabel('Time (s)')
ylabel('Position (m)') 
legend('Vertical Position (m)','Vertical Velocity (m)');

figure(3)
plot(timeVector,xVector,timeVector,vxVector) 
xlabel('Time (s)')
ylabel('Position (m)') 
legend('Horizontal Position (m)','Horizontal Velocity (m)');

end