classdef Local
    %LOCAL 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        beta;                       % the temperature
        M;                          % the Troter number
        StatesKeeped;               % the keeped state number
    end
    
    properties (Constant)
        MyRank = 1;
        TotNodeNo = 1;
        EvaluateInternalEnergy = 1;
        StateNoPerSite = 2;
        vSpin = [-1,1];
        vDim = [1,1];
        TargetSpin = 0;
    end
    
    properties (Dependent)      
        lambda;    
        LocalTensor;
    end
    
    methods
        function obj = Local(varargin)
             if length(varargin) == 3
                obj.beta = varargin{1};
                obj.M = varargin{2};
                obj.StatesKeeped = varargin{3};
             end         
        end
        
        function lambda = get.lambda(theLocal)
           lambda = theLocal.beta/theLocal.M*2; 
        end
        
        function LocalTensor = get.LocalTensor(theLocal)
            Sz = [1/2 0;0 -1/2];    %sigma_z
            Sp = [0 1;0 0];         %s^+， Sp=Sx+Im*Sy;
            Sm = [0 0;1 0];         %s^-， Sm=Sx-Im*Sy;
            SI = eye(2);
            s1z = kron(Sz,SI); s2z = kron(SI,Sz);
            s1p = kron(Sp,SI); s2p = kron(SI,Sp);
            s1m = kron(Sm,SI); s2m = kron(SI,Sm);

            local_horizontal = ...
                2*(s1z*s2z + 1/2*(s1p*s2m+s1m*s2p))-0.5*eye(4);
            tm = expm(-theLocal.lambda*local_horizontal);
            a = 1/2*(tm(1,1)+tm(4,4));
            b = 1/2*(tm(2,2)+tm(3,3));
            c = 1/2*(tm(2,3)+tm(3,2));
            cS1 = cell(1,4);
            cS1{1} = kron(Sz,SI); cS1{2} = kron(Sp,SI);
            cS1{3} = kron(Sm,SI); cS1{4} = kron(SI,SI);
            cS2 = cell(1,4);
            cS2{1} = kron(SI,Sz); cS2{2} = kron(SI,Sp);
            cS2{3} = kron(SI,Sm); cS2{4} = kron(SI,SI);
            
            g = @(x,y) (c+b)/2*x{4}*y{4}+2*(c-b)*x{1}*y{1}...
                +a*(x{2}*y{3}+x{3}*y{2});
            LocalTensor = permute(reshape(g(cS1,cS2),[2,2,2,2]),[1,3,2,4]);
            
        end
    end
    
end

