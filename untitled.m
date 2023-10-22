        %Problem Name
        Name = 'Example221';
        
        %Length and Supports
        Length = 10; Supports = [2, 10]; % length  = 10, supports at 2 and 10;
        prob = SFBMProb(Name, Length, Supports);
        
        %Set Unit
        prob.ForceUnit = "lb";
        prob.LengthUnit = "inch";
        
        %Concetrated Loads
        prob.AddPointLoad(-5, 0); % 5N downward at point 0
        prob.AddPointLoad(-10, 8); % 10N downward at point 8
        
        %Torques
        prob.AddMoment(10, 3);  % ACW 10Nm at point 3
        prob.AddMoment(-10, 7); % CW 10Nm at point 7
        
        %Solve the problem
        prob.Solve()
        