classdef PointLoad < handle
    properties 
        Val
        Loc
        sh
        num
    end

    methods
        function p = PointLoad(val, loc)
            p.Val = val;
            p.Loc = loc;
        end
    end
end