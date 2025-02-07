 f = @(t,y) (0.65*t-(0.65/8.1)*t^2-((1.2*t^2)/(1+t^2)))
 tval = t1:dt:t2; %(find t1,t2 and dt)
 yval = y1:dy:y2; %(find y1,y2,dy) 
 dirfield(f, t1:dt:t2, y1:dy:y2, 'Dirfield graph')