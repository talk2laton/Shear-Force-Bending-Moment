classdef Moment < handle
    properties 
        Val
        Loc
        sh
        num
    end

    methods
        function m = Moment(val, loc)
            m.Val = val;
            m.Loc = loc;
        end
    end
end