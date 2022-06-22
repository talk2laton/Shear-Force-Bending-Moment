function Diagrams(x,sf,bm)
close all
clc
global Name nc nd nm Cload Dload Mload Xload Cloc Dloc Mloc Xtick YtickSF ....
       YtickBM PSF PBM DeciPlace TypeF TypeM TypeD EqRange



%% Display Information Processing
Loads = [];
for m = 1:numel(Cload)
    Loads = [Loads,Cload{m}];
end
for m = 1:numel(Dload)
    Loads = [Loads,Dload{m}];
end

% Factors needed for proper units so Matlab axis is not adjusted after creation
Xfactor     = 10^floor(log10(max(abs(x))));
SFfactor    = 10^floor(log10(max(abs(sf))));
BMfactor    = 10^floor(log10(max(abs(bm))));
maxforce    = max(Loads); minforce = min(Loads);

% To adjust length of arrows and number of arrows
lenperforce = 2/max([maxforce, -minforce]);
numperlen   = 30/(max(x) - min(x));


% Tick Points and Values
Xtick   = unique(Round(Xtick,DeciPlace));
YtickSF = unique(Round(YtickSF,DeciPlace));
YtickBM = unique(Round(YtickBM,DeciPlace));
XTICK   = {};
YTICKSF = {};
YTICKBM = {};
XLOAD   = {};
format  = ['%.',num2str(DeciPlace),'f'];
for n = 1:numel(Xtick)
    XTICK{n} = num2str(Xtick(n)/Xfactor,format);
end

for n = 1:numel(Xload)
    XLOAD{n} = num2str(Xload(n)/Xfactor,format);
end

for n = 1:numel(YtickSF)
    YTICKSF{n} = num2str(YtickSF(n)/SFfactor,format);
end

for n = 1:numel(YtickBM)
    YTICKBM{n} = num2str(YtickBM(n)/BMfactor,format);
end


%% Creating Freebody Diagram and Equation Figure
h1 = figure(1);
set(h1,'units','normalized','outerposition',[0 0.05 0.6 0.95],'Color', [1.00 1.00 1.00])
ax1 = axes('Position',[0 0 1 1],'Visible','off');
ax2 = axes('Position',[0.05,0.50,0.90,0.45],'Visible','off');
axes(ax2)
t = 0.2; 
fill([0,0,x(end),x(end),0],[-t/2,t/2,t/2,-t/2,-t/2],'b','linewidth',1.5)
minx = min(x); maxx = max(x);
Xmin = minx - 0.05*maxx;
Xmax = maxx + 0.05*maxx;
ymin = 0; ymax = 0;
box off
hold on
PlotHandles = [];
PlotNames   = {};

Arrowhead = [(Xmax - Xmin)/150, (maxforce-minforce)/(75*SFfactor)];
% Concentrated Forces
nca = 0; ncr = 0;
for n   = 1:nc
    xl  = Cloc{n};
    yl2 = sign(-Cload{n})*t/2;
    yl1 = -Cload{n}*lenperforce + yl2;
    ymin = min([ymin,yl1]); ymax = max([ymax,yl1]);
    text(xl,yl1+0.2*sign(yl1),[num2str(abs(Cload{n}),'%.2f'),'KN'], 'FontSize',12 ,'FontWeight','bold', 'HorizontalAlignment','center', 'interpreter','latex')
    if yl1 ~= yl2
        if TypeF(n) == 'a';
            if (nca == 0)
                nca = 1;
                hca = DrawArrow(xl,[yl1,yl2],Arrowhead,'linewidth',3,'color','k');
                PlotHandles = [PlotHandles; hca]; PlotNames = [PlotNames;{'Applied Concentrated Force'}];
            else
                DrawArrow(xl,[yl1,yl2],Arrowhead,'linewidth',1.5,'color','k');
            end
        else
            if (ncr == 0)
                ncr = 1;
                hcr = DrawArrow(xl,[yl1,yl2],Arrowhead,'linewidth',3,'color','k','linestyle','-.');
                PlotHandles = [PlotHandles; hcr];  PlotNames = [PlotNames;{'Reacting Concentrated Force'}];
            else
                DrawArrow(xl,[yl1,yl2],Arrowhead,'linewidth',1.5,'color','k','linestyle','-.');
            end
        end
    end
end

% Distributed Forces
nda = 0; ndr = 0;
for n  = 1:nd
    xl     = Dloc{n};
    num    = ceil(numperlen * (xl(end) - xl(1)));
    xdist  = linspace(xl(1),xl(end),num);
    ps     = polyfit(xl,-Dload{n}*lenperforce,numel(xl)-1);
    ydist1 = polyval(ps,xdist);
    IDzero = find(ydist1==0);    
    ydist2 = sign(ydist1)*t/2;
    if IDzero == 1
        ydist2(IDzero) = sign(ydist1(end))*t/2; 
    elseif IDzero == numel(ydist1)
        ydist2(IDzero) = sign(ydist1(1))*t/2;
    end
    ydist1 = ydist1 + ydist2;
    ymin   = min([ymin,min(ydist1)]); ymax = max([ymax,max(ydist1)]);
    
    text(xdist(1),ydist1(1) + 0.2*sign(ydist1(1)),[num2str(abs(Dload{n}(1)),'%.2f'),'KN/m'], 'FontSize',12 ,'FontWeight','bold', 'HorizontalAlignment','center', 'interpreter','latex')
    text(xdist(end),ydist1(end) + 0.2*sign(ydist1(end)),[num2str(abs(Dload{n}(end)),'%.2f'),'KN/m'], 'FontSize',12 ,'FontWeight','bold', 'HorizontalAlignment','center', 'interpreter','latex')
        
    for nn = 1:num
        if ydist2(nn) ~= ydist1(nn)
            if TypeD(n) == 'a';
                if (nda == 0)
                    nda = 1;
                    hda = DrawArrow(xdist(nn),[ydist1(nn), ydist2(nn)],Arrowhead,'linewidth',1.5,'color','m');
                    PlotHandles = [PlotHandles; hda];  PlotNames = [PlotNames;{'Applied Distributed Force'}];
                else
                    DrawArrow(xdist(nn),[ydist1(nn), ydist2(nn)],Arrowhead,'linewidth',1.5,'color','m');
                end
            else
                if (ndr == 0)
                    ndr = 1;
                    hdr = DrawArrow(xdist(nn),[ydist1(nn), ydist2(nn)],Arrowhead,'linewidth',1.5,'color','m','LineStyle','-.');
                    PlotHandles = [PlotHandles; hdr];  PlotNames = [PlotNames;{'Reacting Distributed Force'}];
                else
                    DrawArrow(xdist(nn),[ydist1(nn), ydist2(nn)],Arrowhead, 'linewidth',1.5,'color','m','LineStyle','-.');
                end
            end
        end
    end
    plot(xdist,ydist1,'m','linewidth',1.5)
end

% Moments
nma = 0; nmr = 0;
for n = 1:nm
    radius = 0.05*x(end);
    t      = sign(Mload{n})*[-3*pi/4,3*pi/4];
    t1     = t(1);
    t2     = t(2);
    if TypeM(n) == 'a';
        if(nma == 0)
            nma = 1;
            hma = MomentArrow(radius,t1,t2,[Mloc{n},0],'r', 1.5);
            PlotHandles = [PlotHandles; hma];  PlotNames = [PlotNames;{'Applied Torque'}];
        else
            MomentArrow(radius,t1,t2,[Mloc{n},0],'r', 1.5);
            
        end
        plot(Mloc{n},0,'or','markersize',5,'markerfacecolor','r')
    else
        if(nmr == 0)
            nmr = 1;
            hmr = MomentArrow(radius,t1,t2,[Mloc{n},0],'-.r', 1.5);
            PlotHandles = [PlotHandles; hmr];  PlotNames = [PlotNames;{'Reacting Torque'}];
        else
            MomentArrow(radius,t1,t2,[Mloc{n},0],'-.r', 1.5);
        end
        plot(Mloc{n},0,'ok','markersize',5,'markerfacecolor','r')
    end
    text(Mloc{n},0.3*sign(Mload{n}),[num2str(abs(Mload{n}),'%.2f'),'KN-m'], 'FontSize',12, 'FontWeight','bold', 'HorizontalAlignment','center', 'interpreter','latex')
end

Ymin = ymin - 0.5;
Ymax = ymax + 2.5;
axis([Xmin, Xmax,Ymin,Ymax])
title('$Free~Body~Diagram$','FontSize',12, 'interpreter','latex')
hold off
box off
ax2.YTick = [];
ax2.XTick = Xload;
ax2.TickLabelInterpreter = 'Latex';
ax2.FontSize = 15;
xlabel('$m$', 'Interpreter','latex');

[xa0, ya0] = ConvertCoordinates(ax2, Xload,zeros(size(Xload)));
[xb0, yb0] = ConvertCoordinates(ax2, Xload,repmat(Ymin,size(Xload)));

for n = 1:numel(Xload)
    h4 = annotation('line',[xa0(n) xb0(n)],[ya0(n) yb0(n)],'Tag' , 'connect1');
    set(h4,'LineStyle','--'); set(h4,'Color','b'); 
end 
legend(PlotHandles,PlotNames, 'FontSize',10,'interpreter', 'latex')

% Equations
Range = {'$Range$'};
Equa1 = {'$ Equations~of~Shear~Force: $'};
Equa2 = {'$ Equations~of~Bending~Moment: $'};

if numel(Xtick) > (numel(PSF)  + 1)
    Rangex = max(Xtick) - min(Xtick);
    dx = Rangex/50;
    ind = 0;
    for i = 1:numel(Xtick)
        if(any(abs(Xtick(i)-Xload)<dx))
            ind = i; TobeRemoved = 1;
            break;
        end
    end 
    while(TobeRemoved)
        Xtick(ind) = []; 
        TobeRemoved = 0;
        for i = 1:numel(Xtick)
            if(any(abs(Xtick(i)-Xload)<dx))
                ind = i; TobeRemoved = 1;
                break; 
            end
        end 
    end
end
Xtick = union(Xtick, Xload);

for n = 2:numel(EqRange)
    range = ['$',num2str(EqRange(n - 1)), '~to~',num2str(EqRange(n)),'$'];
    equa1 = ['$',makeequation(PSF{n - 1}),'$'];
    equa2 = ['$',makeequation(PBM{n - 1}),'$'];
    Range = [Range;{range}]; Equa1 = [Equa1;{equa1}]; Equa2 = [Equa2;{equa2}]; 
end
axes(ax1)
text(0.05, 0.40, Range, 'VerticalAlignment', 'cap', 'FontSize', 12, 'interpreter', 'latex')
text(0.20, 0.40, Equa1, 'VerticalAlignment', 'cap', 'FontSize', 12, 'interpreter', 'latex')
text(0.50, 0.40, Equa2, 'VerticalAlignment', 'cap', 'FontSize', 12, 'interpreter', 'latex')

f1 = getframe(gcf);

%%
h2 = figure(2);
set(h2,'units','normalized','outerposition',[0.6 0.05 0.4 0.95],'Color', [0.98 0.98 0.98]);
y = zeros(size(x)); % To be used for fill
% Shear Force
subplot(2,1,1);
fill([x;flipud(x)],[y;flipud(sf)],'b','facealpha',0.15,'edgecolor','none');hold on
plot(x,sf,'b','Linewidth',2);
offset = 0.01 + 0.05*(max(sf) - min(sf));
Vmin = min(sf) - offset;
Vmax = max(sf) + offset;
axis([Xmin, Xmax, Vmin, Vmax]);
plot([Xmin, Xmax], [0,0], 'k'); hold off
% Latex Formatting
ax3 = gca;
ax3.YTick = unique(round(YtickSF, 3, 'significant'));
ax3.XTick= unique(round(Xtick, 3, 'significant'));
ax3.TickLabelInterpreter = 'latex';

ylabel('$V(KN)$', 'Interpreter','latex');
xlabel('$m$', 'Interpreter','latex');

title('$Shear~Force~Diagram$','FontSize',12, 'interpreter','latex')


%%
subplot(2,1,2);
fill([x;flipud(x)],[y;flipud(bm)],'b','facealpha',0.15,'edgecolor','none');hold on
plot(x,bm,'b','Linewidth',2);
offset = 0.01 + 0.05*(max(bm) - min(bm));
Mmin = min(bm) - offset;
Mmax = max(bm) + offset;
axis([Xmin, Xmax, Mmin, Mmax]);
plot([Xmin, Xmax], [0,0], 'k'); hold off
% Latex Formatting
ax4 = gca;
ax4.YTick = unique(round(YtickBM, 3, 'significant'));
ax4.XTick = unique(round(Xtick, 3, 'significant'));
ax4.TickLabelInterpreter = 'latex';
ylabel('$M(KN-m)$', 'Interpreter','latex');
xlabel('$m$', 'Interpreter','latex');
title('$Bending~Monent~Diagram$','FontSize',12, 'interpreter','latex')

%% annotation
[xa1, ya1] = ConvertCoordinates(ax3, Xtick,repmat(Vmax,size(Xtick)));
[xb1, yb1] = ConvertCoordinates(ax3, repmat(Xmin,size(YtickSF)),YtickSF);
[xc1, yc1] = ConvertCoordinates(ax3, repmat(Xmax,size(YtickSF)),YtickSF);

[xa2, ya2] = ConvertCoordinates(ax4, Xtick,repmat(Mmin,size(Xtick)));
[xb2, yb2] = ConvertCoordinates(ax4, repmat(Xmin,size(YtickBM)),YtickBM);
[xc2, yc2] = ConvertCoordinates(ax4, repmat(Xmax,size(YtickBM)),YtickBM);

for n = 1:numel(Xtick)
    h4 = annotation('line',[xa1(n) xa2(n)],[ya1(n) ya2(n)],'Tag' , 'connect1');
    set(h4,'LineStyle','--'); set(h4,'Color','b'); 
end

for n = 1:numel(YtickSF)
    h4 = annotation('line',[xb1(n) xc1(n)],[yb1(n) yc1(n)],'Tag' , 'connect1');
    set(h4,'LineStyle','--'); set(h4,'Color','b'); 
end

for n = 1:numel(YtickBM)
    h4 = annotation('line',[xb2(n) xc2(n)],[yb2(n) yc2(n)],'Tag' , 'connect1');
    set(h4,'LineStyle','--'); set(h4,'Color','b'); 
end

f2 = getframe(gcf);
F = [f1.cdata,f2.cdata];
imwrite(F,[Name,'.png'])

function equation = makeequation(p)
    equation = [];
    for m = 1:numel(p)
        vv = p(m); power = numel(p) - m;
        avv = abs(vv);
        vvstr = [];
        if (avv ~= 0)
            if (avv == 1) 
                if(isempty(equation))
                    if(vv < 0)
                        if power == 0
                            vvstr = ' - 1';
                        else
                            vvstr = ' - ';
                        end
                    else
                        if power == 0
                            vvstr = ' 1';
                        end
                    end
                else
                    if(vv < 0)
                        vvstr = ' - ';
                    else
                        vvstr = ' + ';
                    end
                end
            else
                if(isempty(equation))
                    vvstr = num2str(vv);
                else
                    if(vv < 0)
                        vvstr = num2str(vv);
                    else
                        vvstr = [' + ',num2str(vv)];
                    end
                end
            end
            if (power > 0)
                if (power > 1) vvstr = [vvstr,'x^',num2str(power)];
                else vvstr = [vvstr,'x'];
                end
            end
        end
        equation = [equation, vvstr];
    end
    if(isempty(equation))
        equation = '0';
    end






