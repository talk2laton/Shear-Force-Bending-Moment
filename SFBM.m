function  [varargout] = SFBM( varargin)
%SFBM calculates the shear force and bending moment, plot the diagram and
%computes the equations of the lines.
global Name nc nd nm Cload Dload Mload Cloc Dloc Mloc Xload DistFfun ....
       DistPoly Xtick YtickSF YtickBM DeciPlace TypeF TypeM TypeD EqRange

N         = nargin;
NP        = 10000;
DeciPlace = 2;
Name      = varargin{1};
Length    = varargin{2}(1);
if numel(varargin{2}) == 1
    Supports = 0;
else
    Supports  = varargin{2}(2:end);
end
Support1  = Supports(1);
M         = 0; % sum moment about support 1.
F         = 0; % sum vertical forces .
nc        = numel(varargin{2}) - 1; nd    = 0; nm  = 0; 
DistFfun  = {};  DistPoly = {};  Xload    = [0,varargin{2}];
Cload     = {};  Dload    = {};  Mload    = {};
Cloc      = {};  Dloc     = {};  Mloc     = {};
Xtick     = [];  YtickSF  = [];  YtickBM  = [];
TypeF     = [];  TypeM    = [];  TypeD    = [];
for n = 3:N
    LoadCell = varargin{n};
    type     = LoadCell{1};
    load     = LoadCell{2};
    loc      = LoadCell{3};
    Xload    = union(Xload, loc);
    if strcmp(type,'CF')
        Mcl      = load*(loc - Support1);   
        M        = M + Mcl;
        Fcl      = load;   
        F        = F + Fcl;
        Cload    = [Cload;{load}];
        Cloc     = [Cloc;{loc}];
        TypeF    = [TypeF,'a'];
        nc       = nc + 1;
    elseif strcmp(type,'DF')
        if (numel(load) == 1) load = repmat(load,size(loc)); end
        p        = polyfit(loc - Support1, load, numel(loc) - 1); 
        DistPoly = [DistPoly, {p}];
        power    = numel(p) - 1: -1: 0;
        Ffun     = @(x) sum(p.*(x - Support1).^(power + 1)./(power + 1)); 
        DistFfun = [DistFfun, {Ffun}];
        Fbl      = Ffun(loc(end)) - Ffun(loc(1));
        F        = F + Fbl;
        Mfun     = @(x) sum(p.*(x - Support1).^(power + 2)./(power + 2)); 
        Mbl      = Mfun(loc(end)) - Mfun(loc(1));
        M        = M + Mbl;
        Dload    = [Dload;{load}];
        Dloc     = [Dloc;{loc}];
        TypeD    = [TypeD,'a'];
        nd       = nd + 1;
    elseif strcmp(type,'M')
        M        = M + load;
        Mload    = [Mload;{load}];
        Mloc     = [Mloc;{loc}];
        TypeM    = [TypeM,'a'];
        nm       = nm + 1;
    end
end
Xload = Round(Xload,2);
if numel(varargin{2}) == 1
    TypeD    = [TypeD,'r'];
    p        = -F/Length;
    power    = numel(p) - 1: -1: 0;
    Dload    = [Dload;{[p, p]}];
    Dloc     = [Dloc;{[0,Length]}]; 
    DistPoly = [DistPoly, {p}];
    Ffun     = @(x) sum(p.*(x - Support1).^(power + 1)./(power + 1)); 
    DistFfun = [DistFfun, {Ffun}];
    Xload    = [Xload,Length];
    nd       = nd + 1;
elseif (numel(Supports) > 1)
    B        = -M/(Supports(2) - Support1); % Reaction B
    A        = -(F + B); % Reaction A
    Cload    = [{A};Cload;{B}];
    TypeF    = ['r',TypeF,'r'];
    Cloc     = [{Supports(1)};Cloc;{Supports(2)}];
else
    if Support1 > 0
        A     = -F; % Reaction A
        Cload = [Cload;{A}];
        TypeF = [TypeF,'r'];
        Cloc  = [Cloc;{Support1}];
        Mload = [Mload;{-M}];
        TypeM = [TypeM,'r'];
        Mloc  = [Mloc;{Support1}];
    else
        A     = -F; % Reaction A
        Cload = [{A};Cload];
        TypeF = ['r',TypeF];
        Cloc  = [{Support1};Cloc];
        Mload = [{-M};Mload];
        TypeM = ['r',TypeM];
        Mloc  = [{Support1};Mloc];
    end
    nm    = nm + 1;
end
NewCcrossed = 0; NewMcrossed = 0;
    function [sf, bm] = SF(x)
        sf = 0; bm = 0;
        NewCcrossed = 0; NewMcrossed = 0;
        for i = 1:nc
            if x > Cloc{i}
                if i > Ccrossed
                    Ccrossed = i;  NewCcrossed = 1;
                end
                sf = sf + Cload{i};
            end
        end
        
        for i = 1:nd
            loc = Dloc{i};
            if x > loc(1)
                Ffun = DistFfun{i};
                sf = sf + Ffun(min([x,loc(end)])) - Ffun(loc(1));
            end
        end
        
        for i = 1:nm
            if x > Mloc{i}
                if i > Mcrossed
                    Mcrossed = i;  NewMcrossed = 1;
                end
                bm = bm - Mload{i};
            end
        end
    end
EqRange  = [0;cell2mat([Cloc(:); Mloc(:)]);Length]';
for n = 1:numel(Dloc)
    EqRange = [EqRange,Dloc{n}];
end
EqRange  = unique(EqRange);
dx       = Length/NP;
X        = 0; base = 0;
Ccrossed = 0; Mcrossed = 0;
ShearF   = 0; 
BendM    = 0;
for n = 1:NP + 2
    xx = (n - 1)*dx;
    [sf,bm] = SF(xx);
    % handling discontinuities
    if(NewCcrossed || NewMcrossed)
        if (NewCcrossed && NewMcrossed)
            [sfc,bmc] = SF(Cloc{Ccrossed}); xc = Cloc{Ccrossed};
            sfc2 = sfc + Cload{Ccrossed}; bmc2 = bmc - Mload{Mcrossed};
        elseif (NewCcrossed && ~NewMcrossed)
            [sfc,bmc] = SF(Cloc{Ccrossed}); xc = Cloc{Ccrossed};
            sfc2 = sfc + Cload{Ccrossed}; bmc2 = bmc;
        elseif (~NewCcrossed && NewMcrossed)
            [sfc,bmc] = SF(Mloc{Mcrossed}); xc = Mloc{Mcrossed};
            sfc2 = sfc; bmc2 = bmc - Mload{Mcrossed};
        end
        X = [X; xc; xc]; Xtick = [Xtick; xc; xc]; 
        ShearF = [ShearF; sfc; sfc2]; YtickSF = [YtickSF; sfc; sfc2];
        BendM = [BendM; base + bmc; base + bmc2]; YtickBM = [YtickBM; base + bmc; base + bmc2];
    end
    if n > 1 && n < NP + 2 
        X = [X;xx];
        ShearF = [ShearF; sf];
        base = base + 0.5*dx*(sf + ShearF(end - 1));
        BendM = [BendM; bm + base];
    end
end
SpecialPoints(X, ShearF, BendM);
Diagrams(X, ShearF, BendM)
varargout = {X, ShearF, BendM};
end