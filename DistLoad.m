classdef DistLoad < handle
    properties 
        Val
        Loc
        sh
        num
    end

    methods
        function d = DistLoad(val, loc)
            d.Val = val;
            d.Loc = loc;
        end
    end
end