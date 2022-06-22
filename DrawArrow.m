function h = DrawArrow(x,y,arrowhead, varargin)
ax = gca;
holding = ishold(ax); 
hold on;
plot([x,x],y,varargin{:});
Xhead = [x - arrowhead(1),x,x + arrowhead(1)];
Yhead = [y(2) + sign(y(1))*arrowhead(2), y(2), y(2) + sign(y(1))*arrowhead(2)];
h = plot(Xhead,Yhead,varargin{:});
if(~holding) 
    hold off
end


