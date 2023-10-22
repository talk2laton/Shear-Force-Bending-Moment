%class definition
classdef SFBMProb < handle
    properties 
        Name              string
        Source            string
        Length            (1, 1) double
        Supports          
        MomentLoad        (1, :) Moment
        ConcentratedLoad  (1, :) PointLoad
        DistributedLoad   (1, :) DistLoad
        AddMoment
        AddPointLoad
        AddDistLoad
        RemoveMoment
        RemovePointLoad
        RemoveDistLoad
        ForceUnit        string
        LengthUnit       string
        Solve
    end

    properties (SetAccess = private)
        axhandle
    end

    methods
        function prob = SFBMProb(name, length, supports)
            prob.Name = name;
            prob.Length = length;
            prob.Supports = supports;
            prob.ForceUnit = "KN";
            prob.LengthUnit = "m";
            prob.AddMoment = @(val, loc) AddMomentLoad(prob, val, loc);
            prob.AddPointLoad = @(val, loc) AddConcentratedLoad(prob, val, loc);
            prob.AddDistLoad = @(val, loc) AddDistributedLoad(prob, val, loc);
            prob.RemoveMoment = @(indx) RemoveMomentLoad(prob, indx);
            prob.RemovePointLoad = @(indx) RemoveConcentratedLoad(prob, indx);
            prob.RemoveDistLoad = @(indx) RemoveDistributedLoad(prob, indx);
            prob.axhandle = FreeBody(prob.Length, prob.Supports);
            prob.Solve = @() Solver(prob);
        end
    
        function prob =  AddMomentLoad(prob, val, loc)
            axes(prob.axhandle);
            M = Moment(val, loc);
            t = sign(val)*[-3*pi/4,3*pi/4];
            gs = numel(prob.axhandle.Children);
            MomentArrow(0.2,t(1),t(end),[loc*5/prob.Length,0],'r', 2);
            text(loc*5/prob.Length,0.3*sign(val),num2str(abs(val),'%.2f') +...
                prob.ForceUnit + "/" + prob.LengthUnit, 'FontWeight','bold',...
                'HorizontalAlignment','center', 'interpreter','latex')
            ge = numel(prob.axhandle.Children);
            M.num = ge-gs;
            M.sh = prob.axhandle.Children(1);
            prob.MomentLoad = [prob.MomentLoad, M];
        end

        function prob =  AddConcentratedLoad(prob, val, loc)
            Arrowhead = [0.05, 0.1];
            axes(prob.axhandle);
            P = PointLoad(val, loc);
            xl  = loc*5/prob.Length;
            yl = sign(-val)*[0.6, 0.05];
            gs = numel(prob.axhandle.Children);
            DrawArrow(xl,yl, Arrowhead,'linewidth', 2,'color','k');
            text(xl,1.1*yl(1),num2str(abs(val),'%.2f') + prob.ForceUnit, ...
                'FontWeight','bold', 'HorizontalAlignment','center', ...
                'interpreter','latex')
            ge = numel(prob.axhandle.Children);
            P.num = ge-gs;
            P.sh = prob.axhandle.Children(1);
            prob.ConcentratedLoad = [prob.ConcentratedLoad, P];
        end

        function prob =  AddDistributedLoad(prob, val, loc)
            Arrowhead = [0.05, 0.1];
            axes(prob.axhandle);
            D = DistLoad(val, loc);
            xl     = loc*5/prob.Length;
            num    = ceil((xl(end) - xl(1))/0.13);
            xdist  = linspace(xl(1),xl(end),num);
            ps     = polyfit(xl,val*0.3/max(abs(val)),numel(xl)-1);
            ydist1 = polyval(ps,xdist); 
            ydist2 = 0.05*sign(-ydist1);
            ydist1 = -ydist1 + ydist2;
            ydist = [ydist1; ydist2];
            gs = numel(prob.axhandle.Children);
            for m = 1:num
                DrawArrow(xdist(m), ydist(:,m),Arrowhead,'linewidth',1.5,'color','m');
            end
            plot(xdist,ydist1,'m','linewidth',1.5);
            text(xdist(1),1.1*ydist1(1), num2str(abs(val(1)),'%.2f') + ....
                prob.ForceUnit + "/" + prob.LengthUnit, 'FontWeight','bold', ....
                'HorizontalAlignment','center', 'interpreter','latex')
            text(xdist(end),1.1*ydist1(end), num2str(abs(val(end)),'%.2f') +...
                prob.ForceUnit + "/" + prob.LengthUnit, 'FontWeight','bold',...
                'HorizontalAlignment','center', 'interpreter','latex')
            ge = numel(prob.axhandle.Children);
            D.num = ge-gs;
            D.sh = prob.axhandle.Children(1);
            prob.DistributedLoad = [prob.DistributedLoad, D];
        end

        function prob =  RemoveMomentLoad(prob, index)
            i = find(prob.axhandle.Children == prob.MomentLoad(index).sh);
            for j = 1:prob.MomentLoad(index).num
                delete(prob.axhandle.Children(i));
            end
            prob.MomentLoad(:,index) = [];
        end

        function prob =  RemoveConcentratedLoad(prob, index)
            i = find(prob.axhandle.Children == prob.ConcentratedLoad(index).sh);
            for j = 1:prob.ConcentratedLoad(index).num
                delete(prob.axhandle.Children(i));
            end
            prob.ConcentratedLoad(:,index) = [];
        end

        function prob =  RemoveDistributedLoad(prob, index)
            i = find(prob.axhandle.Children == prob.DistributedLoad(index).sh);
            for j = 1:prob.DistributedLoad(index).num
                delete(prob.axhandle.Children(i));
            end
            prob.DistributedLoad(:,index) = [];
        end
        
        function Solver(prob)
            global ForceUnit LengthUnit
            ForceUnit = prob.ForceUnit;
            LengthUnit = prob.LengthUnit;
            input = {prob.Name, [prob.Length, prob.Supports]};
            for M = prob.MomentLoad
                input = [input, {{'M', M.Val, M.Loc}}];
            end
            for P = prob.ConcentratedLoad
                input = [input, {{'CF', P.Val, P.Loc}}];
            end
            for D = prob.DistributedLoad
                input = [input, {{'DF', D.Val, D.Loc}}];
            end
            f1 = getframe(gcf);
            [X, ShearF, BendM, f2, f3] = SFBM(input{:});
            sz1 = size(f1.cdata,1);
            pad1 = floor((size(f2.cdata,2) - size(f3.cdata,2))/2);
            pad2 = size(f2.cdata,2) - size(f3.cdata,2)-pad1;
            v1 = repmat(255, sz1, pad1, 3);
            v2 = repmat(255, sz1, pad2, 3);
            F = [f1.cdata; f2.cdata; 
                [uint8(v1), f3.cdata, uint8(v2)]]; 
            if(~isempty(prob.Source))
                imtext=text2im(char("   source = "+prob.Source));
                [m, n] = size(imtext);
                [~, p, ~] = size(F);
                pad = 255*uint8(ones(3*m, p, 3));
                pad(m+1:2*m,1:n,:) = 255*uint8(cat(3, imtext,imtext,imtext));
                F = [F;pad(:,1:p,:)];
            end
            imwrite(F, prob.Name+".tiff");
        end
    end
end

function axhandle = FreeBody(Length, Supports)
    figure(Color = 'w', Units = 'normalized', Outerposition = [0.01 0.05 0.647 0.95])
    axhandle = axes('Visible','off'); axes(axhandle);
    axis([-0.5, 5.5, -1, 1]); axis equal; axis off;  hold on;
    xs = Supports*5/Length;
    if(numel(xs) == 2)
        [Xt,Yt, Xr, Yr] = makesupport1([xs(1), 0], 0.1, 0);
        fill(Xt,Yt,'k','LineWidth',1.5, 'FaceAlpha', 0.5);
        fill(Xr,Yr,'k','LineWidth',1.5, 'FaceAlpha', 0.3);
        [Xt,Yt, Xr, Yr, Xl, Yl, Xc, Yc] = makesupport2([xs(2), 0], 0.1, 0);
        fill(Xt,Yt,'k','LineWidth',1.5, 'FaceAlpha', 0.5);
        fill(Xr,Yr,'k','LineWidth',1.5, 'FaceAlpha', 0.3);
        fill(Xl,Yl,'k','LineWidth',1.5, 'FaceAlpha', 0.3);
        arrayfun(@(n) fill(Xc(n,:),Yc(n,:),'k','LineWidth',1.5, ...
            'FaceAlpha', 0.3), 1:size(Xc,1));
        plot(xs, 0, 'or', MarkerFaceColor = 'r', MarkerSize = 2); 
    elseif(numel(xs) == 1)
        tt = pi*linspace(0.5,1.5,100); R = 0.5+0.05*rand(1,100);
        xx = R.*cos(tt)+0.1; yy = R.*sin(tt);
        if(xs == 0)
            fill(xx,yy,'k','LineWidth',0.01, 'FaceAlpha', 0.2);
            plot(0.1+[xs, xs], [0.05, 0.5],'k', 'LineWidth',1.5);
            plot(0.1+[xs, xs], -[0.05, 0.5],'k', 'LineWidth',1.5);
        else
            fill(xs-xx,yy,'k','LineWidth',0.01, 'FaceAlpha', 0.2);
            plot(-0.1+[xs, xs], [0.05, 0.5],'k', 'LineWidth',1.5);
            plot(-0.1+[xs, xs], -[0.05, 0.5],'k', 'LineWidth',1.5);
        end
    end
    fill(5*[0,0,1,1],0.05*[-1,1,1,-1],'b','LineWidth',1.5, 'FaceAlpha', 0.5);
end

function [Xt, Yt, Xr, Yr] = makesupport1(Node, t, n)
    t    = 0.7*t; y = 4*t; x = t*cos(pi/6)*5.5/1.5;
    t1   = linspace(5*pi/6, pi/6, 50); angle = n*pi/2;
    XYpt = [t*cos(t1), x, -x; t*sin(t1),-y,-y];
    XYpr = [1.2*x, 1.2*x, -1.2*x, -1.2*x; -y, -1.2*y, -1.2*y, -y];
    XYt  = [cos(angle), -sin(angle); sin(angle), cos(angle)]*XYpt;
    XYr  = [cos(angle), -sin(angle); sin(angle), cos(angle)]*XYpr;
    Xt   = Node(1) + XYt(1,:);    Yt     = Node(2) + XYt(2,:);
    Xr   = Node(1) + XYr(1,:);    Yr     = Node(2) + XYr(2,:);
end

function [Xt, Yt, Xr, Yr, Xl, Yl, Xc, Yc] = makesupport2(Node, t, n)
    t    = 0.7*t; y = 3*t; x = t*cos(pi/6)*4.5/1.5;
    t1   = linspace(5*pi/6, pi/6, 50); angle = n*pi/2;
    XYpt = [t*cos(t1), x, -x; t*sin(t1),-y,-y];
    y2   = 4*t; x2 = t*cos(pi/6)*5.5/1.5;
    XYpl = [x2, -x2; -y, -y];
    XYpr = [1.2*x2, 1.2*x2, -1.2*x2, -1.2*x2; -y2, -1.2*y2, -1.2*y2, -y2];
    XYt  = [cos(angle), -sin(angle); sin(angle), cos(angle)]*XYpt;
    XYr  = [cos(angle), -sin(angle); sin(angle), cos(angle)]*XYpr;
    XYl  = [cos(angle), -sin(angle); sin(angle), cos(angle)]*XYpl;
    Xt   = Node(1) + XYt(1,:);    Yt     = Node(2) + XYt(2,:);
    Xr   = Node(1) + XYr(1,:);    Yr     = Node(2) + XYr(2,:);
    Xl   = Node(1) + XYl(1,:);    Yl     = Node(2) + XYl(2,:);
    Xc   = []; Yc = [];
    for n = 1:6
        xc = (-3.5 + n)*t + 0.5*t*cos(linspace(0,2*pi,100));
        yc = -3.5*t + 0.5*t*sin(linspace(0,2*pi,100));
        XYpc = [xc;yc];
        XYc  = [cos(angle), -sin(angle); sin(angle), cos(angle)]*XYpc;
        Xc   = [Xc; Node(1) + XYc(1,:)] ;
        Yc   = [Yc; Node(2) + XYc(2,:)] ;
    end
end
