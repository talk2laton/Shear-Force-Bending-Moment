function y = Interpolate(X,Y,x)
y = [];
for n = 1:numel(x)-1
    xn = x(n);
    for i = 1:numel(X)
        xlow = X(i); xhigh = X(i + 1);
        if(xn >= xlow && xn <= xhigh) || (xn <= xlow && xn >= xhigh)
            yn = Y(i) + (xn - xlow)*(Y(i+1) - Y(i))/(xhigh - xlow);
            break;
        end
    end
    if (xlow ~= xhigh) y = [y;yn]; end
end