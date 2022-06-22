function h = MomentArrow(r,t1,t2,location, style, thickness)
t      = linspace(t1,t2,50)';
loc2   = [r*cos(t2),r*sin(t2)] + location;
r2     = r/2;
ts1    = t2+sign(t1-t2)*3*pi/4;
ts2    = t2+sign(t1-t2)*5*pi/16;
p1     = [r2*cos(ts1),r2*sin(ts1)] + loc2;
p2     = [r2*cos(ts2),r2*sin(ts2)] + loc2;
L      = [r*cos(t),r*sin(t)] + ones(50,1)*location;
         plot(L(:,1),L(:,2),style, 'linewidth',thickness);
         plot([L(50,1),p1(1)],[L(50,2),p1(2)],style, 'linewidth',thickness);
h      = plot([L(50,1),p2(1)],[L(50,2),p2(2)],style, 'linewidth',thickness);




