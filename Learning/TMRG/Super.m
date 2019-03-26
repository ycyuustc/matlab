classdef Super
    %SUPER 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        Sys;
        Env;
        SysDown;
        SysUp;
        EnvDown;
        EnvUp;
        SysTm;
        EnvTm;
        FindQuanNos2e;
        FindQuanNoe2s;
        mItoF;
        cFtoI;
        SysNew;
        EnvNew;
        DmSys;
        DmEnv;
        PsiL;
        PsiR;
        Matrix;
        EnvNewRight;
        EnvNewLeft;
        EnvNewEigValue;
        SysNewRight;
        SysNewLeft;
        SysNewEigValue;
    end
    
    methods
        %% overloaded construction method
        function obj = Super(varargin)
            if length(varargin) == 7
                % Super(local,sys,env,sysdown,sysup,envdown,envup)
                obj = obj.Super_construct_ABCD(obj,varargin);
            elseif length(varargin) == 5
                % Super(local,sys,env,systm,envtm)
                obj = obj.Super_construct_AB(obj,varargin);
            else
                error('error of initializing Sub object');
            end
        end
        
        function obj = Super_construct_ABCD(obj,varargin)
            vc = varargin{2};
            local = vc{1}; sys = vc{2}; env = vc{3};
            sysdown = vc{4}; sysup = vc{5};
            envdown = vc{6}; envup = vc{7};
            obj.Sys = sys;
            obj.Env = env;
            obj.SysDown = sysdown;
            obj.SysUp = sysup;
            obj.EnvDown = envdown;
            obj.EnvUp = envup;
            submid1 = sysdown.Sub_up;
            submid2 = envdown.Sub_up;
            
            % calculate the findquanno tensor
            vS1 = sys.vSpin; vS2 = env.vSpin;
            vD1 = sys.vDim; vD2 = env.vDim;
            vSa = local.vSpin; vDa = local.vDim;
            vS1mid = submid1.vSpin;
            vD1mid = submid1.vDim;
            vS2mid = submid2.vSpin;
            vD2mid = submid2.vDim;
            dim0 = length(vSa);
            dim1 = length(vS1);
            dim2 = length(vS2);
            fqn = zeros(dim1,dim0,dim0);
            fqn2 = zeros(dim1,dim0,dim0);
            for a1 = 1:dim0
                for a2 = 1:dim0
                    for m1 = 1:dim1
                        for m2 = 1:dim2
                            ts = vSa(a1)+vSa(a2)+vS1(m1)+vS2(m2);
                            if ts==local.TargetSpin
                                fqn(m1,a1,a2) = m2;
                                fqn2(m2,a1,a2) = m1;
                            end
                        end
                    end
                end
            end
            obj.FindQuanNos2e = fqn;
            obj.FindQuanNoe2s = fqn2;
            
            % find the mapping between I and F
            tmItoF = zeros(dim0,dim0,dim1,dim1,max(vS1),max(vS2));
            tcFtoI = cell(1,dim0*dim0*dim1*dim1*max(vS1)*max(vS2));
            point = 0;
            for a1 = 1:dim0
                for a2 = 1:dim0
                    for m1 = 1:dim1
                        m2 = fqn(m1,a1,a2);
                        if m2~=0
                           for mi1 = 1:vD1(m1)
                               for mi2 = 1:vD2(m2)
                                   point = point + 1;
                                   tmItoF(a1,a2,m1,m2,mi1,mi2) = point;
                                   tcFtoI{point} = [a1,a2,m1,m2,mi1,mi2];
                               end
                           end
                        end
                    end
                end
            end
            superdim = point;
            tcFtoI = tcFtoI(1:superdim);
            obj.mItoF = tmItoF;
            obj.cFtoI = tcFtoI;
                  
            % find the super block matrix
            matrix = zeros(superdim,superdim);
            for p = 1:superdim
                for q = 1:superdim
                    cp = tcFtoI{p};
                    cq = tcFtoI{q};
                    a1 = cp(1); a2 = cp(2); m1 = cp(3); m2 = cp(4);
                    mi1 = cp(5); mi2 = cp(6);
                    b1 = cq(1); b2 = cq(2); n1 = cq(3); n2 = cq(4);
                    ni1 = cq(5); ni2 = cq(6);
                    pointsum = 0;
                    for c1 = 1:dim0
                        for c2 = 1:dim0
                            x = sysdown.FindQuanNo_up(n1,c1,b1);
                            y = envdown.FindQuanNo_up(n2,c2,b2);
                            if x~=0 && y~=0
                            if x == sysup.FindQuanNo_down(m1,a2,c2) && ...
                                    y == envup.FindQuanNo_down(m2,a1,c1)
                                for xi1 = 1:vD1mid(x)
                                    for yi2 = 1:vD2mid(y)
                                        
                                        tc = sysdown.Tm;
                                        tm = tc{x,c1,b1};
                                        if isempty(tm)
                                            p1 = 0;
                                        else
                                            p1 = tm(xi1,ni1);
                                        end
                                        
                                        tc = sysup.Tm;
                                        tm = tc{m1,a2,c2};
                                        if isempty(tm)
                                            p2 = 0;
                                        else
                                            p2 = tm(mi1,xi1);
                                        end
                                        
                                        tc = envdown.Tm;
                                        tm = tc{y,c2,b2};
                                        if isempty(tm)
                                            p3 = 0;
                                        else
                                            p3 = tm(yi2,ni2);
                                        end
                                        
                                        tc = envup.Tm;
                                        tm = tc{m2,a1,c1};
                                        if isempty(tm)
                                            p4 = 0;
                                        else
                                            p4 = tm(mi2,yi2);
                                        end
                                        
                                        pointsum = pointsum + p1*p2*p3*p4;
                                        
                                    end
                                end              
                            end
                            end
                        end
                    end
                    matrix(p,q) = pointsum;
                end
            end
            
            % find the density matrix for system and enviroment
            % we should first generate a expanded sub site.
            [psi_L,~] = eigs(matrix,1);
            [psi_R,~] = eigs(transpose(matrix),1);
            obj.PsiL = psi_L;
            obj.PsiR = psi_R;
            obj.Matrix = matrix;
            % the expanded site: 
            sysnew = Sub(sys,local,'add');
            dms = cell(1,sysnew.totSubSetNo);
            for M = 1:sysnew.totSubSetNo
                tmatrix = zeros(sysnew.vDim(M),sysnew.vDim(M));
               for Mi = 1:sysnew.vDim(M)
                   for Ni = 1:sysnew.vDim(M)             
                       tv = sysnew.cItoOldI{M,Mi};
                       m1 = tv(1); mi1 = tv(2); a2 = tv(3);
                       n1 = tv(1); ni1 = tv(2); b2 = tv(3);
                       sumup = 0;
                       for a1 = 1:dim0
                           m2 = fqn(m1,a1,a2);
                           n2 = fqn(n1,a1,b2);
                           if m2 == n2 && m2~=0
                               for mi2 = 1:vD2(m2)
                                   p = tmItoF(a1,a2,m1,m2,mi1,mi2);
                                   q = tmItoF(a1,b2,n1,m2,ni1,mi2);
                                   sumup = sumup + psi_L(p)*psi_R(q);
                               end
                           end
                       end
                       tmatrix(Mi,Ni) = sumup;
                   end
               end
               dms{M} = tmatrix;
            end
            obj.SysNew = sysnew;
            obj.DmSys = dms;
            
            % the expanded site:
            envnew = Sub(env,local,'add');
            dme = cell(1,envnew.totSubSetNo);
            for M = 1:envnew.totSubSetNo
                tmatrix = zeros(envnew.vDim(M),envnew.vDim(M));
                for Mi = 1:envnew.vDim(M)
                    for Ni = 1:envnew.vDim(M)
                        tv = envnew.cItoOldI{M,Mi};
                        m2 = tv(1); mi2 = tv(2); a2 = tv(3);
                        n2 = tv(1); ni2 = tv(2); b2 = tv(3);
                        sumup = 0;
                        for a1 = 1:dim0
                            m1 = fqn2(m2,a1,a2);
                            n1 = fqn2(n2,a1,b2);
                            if m1 == n1 && m1~=0
                                for mi1 = 1:vD1(m1)
                                    p = tmItoF(a1,a2,m1,m2,mi1,mi2);
                                    q = tmItoF(a1,b2,m1,n2,mi1,ni2);
                                    sumup = sumup + psi_L(p)*psi_R(q);
                                end
                            end
                        end
                        tmatrix(Mi,Ni) = sumup;
                    end
                end
                dme{M} = tmatrix;         
            end
            obj.EnvNew = envnew;
            obj.DmEnv = dme;
            
            % calculate the right vector and left vector for Dm
            % the system
            sysnewright = cell(1,sysnew.totSubSetNo);
            sysnewleft = cell(1,sysnew.totSubSetNo);
            sysneweigvalue = cell(1,sysnew.totSubSetNo);
            for M = 1:sysnew.totSubSetNo
               matrix = dms{M};
               [mR,mL,vd] = fn_eig(matrix,size(matrix,1));
               sysnewright{M} = mR;
               sysnewleft{M} = mL;
               sysneweigvalue{M} = vd;
            end
            
            % the enviroment
            envnewright = cell(1,envnew.totSubSetNo);
            envnewleft = cell(1,envnew.totSubSetNo);
            envneweigvalue = cell(1,envnew.totSubSetNo);
            for M = 1:envnew.totSubSetNo
               matrix = dme{M};
               [mR,mL,vd] = fn_eig(matrix,size(matrix,1));
               envnewright{M} = mR;
               envnewleft{M} = mL;
               envneweigvalue{M} = vd;
            end
            
            obj.SysNewRight = sysnewright;
            obj.SysNewLeft = sysnewleft;
            obj.SysNewEigValue = sysneweigvalue;
            obj.EnvNewRight = envnewright;
            obj.EnvNewLeft = envnewleft;
            obj.EnvNewEigValue = envneweigvalue;
            
            
                   
        end
        
         function obj = Super_construct_AB(obj,varargin)
            vc = varargin{2};
            local = vc{1}; sys = vc{2}; env = vc{3};
            systm = vc{4}; envtm = vc{5};
            obj.Sys = sys;
            obj.Env = env;
            obj.SysTm = systm;
            obj.EnvTm = envtm;
            
%             submid1 = sysdown.Sub_up;
%             submid2 = envdown.Sub_up;
            
            % calculate the findquanno tensor
            vS1 = sys.vSpin; vS2 = env.vSpin;
            vD1 = sys.vDim; vD2 = env.vDim;
            vSa = local.vSpin; vDa = local.vDim;
           
            dim0 = length(vSa);
            dim1 = length(vS1);
            dim2 = length(vS2);
            fqn = zeros(dim1,dim0,dim0);
            fqn2 = zeros(dim1,dim0,dim0);
            for a1 = 1:dim0
                for a2 = 1:dim0
                    for m1 = 1:dim1
                        for m2 = 1:dim2
                            ts = vSa(a1)+vSa(a2)+vS1(m1)+vS2(m2);
                            if ts==local.TargetSpin
                                fqn(m1,a1,a2) = m2;
                                fqn2(m2,a1,a2) = m1;
                            end
                        end
                    end
                end
            end
            obj.FindQuanNos2e = fqn;
            obj.FindQuanNoe2s = fqn2;
            
            % find the mapping between I and F
            tmItoF = zeros(dim0,dim0,dim1,dim1,max(vS1),max(vS2));
            tcFtoI = cell(1,dim0*dim0*dim1*dim1*max(vS1)*max(vS2));
            point = 0;
            for a1 = 1:dim0
                for a2 = 1:dim0
                    for m1 = 1:dim1
                        m2 = fqn(m1,a1,a2);
                        if m2~=0
                           for mi1 = 1:vD1(m1)
                               for mi2 = 1:vD2(m2)
                                   point = point + 1;
                                   tmItoF(a1,a2,m1,m2,mi1,mi2) = point;
                                   tcFtoI{point} = [a1,a2,m1,m2,mi1,mi2];
                               end
                           end
                        end
                    end
                end
            end
            superdim = point;
            tcFtoI = tcFtoI(1:superdim);
            obj.mItoF = tmItoF;
            obj.cFtoI = tcFtoI;
                  
            % find the super block matrix
            matrix = zeros(superdim,superdim);
            for p = 1:superdim
                for q = 1:superdim
                    cp = tcFtoI{p};
                    cq = tcFtoI{q};
                    a1 = cp(1); a2 = cp(2); m1 = cp(3); m2 = cp(4);
                    mi1 = cp(5); mi2 = cp(6);
                    b1 = cq(1); b2 = cq(2); n1 = cq(3); n2 = cq(4);
                    ni1 = cq(5); ni2 = cq(6);
                    pointsum = 0;
                    for c1 = 1:dim0
                        for c2 = 1:dim0
                            if systm.FindQuanNo_down(m1,c1,b1,c2,b2) == n1
                            if envtm.FindQuanNo_down(m2,a2,c2,a1,c1) == n2
                               ttm = systm.Tm{m1,c1,b1,c2,b2};
                               p1 = ttm(mi1,ni1);
                               ttm = envtm.Tm{m2,a2,c2,a1,c1};
                               p2 = ttm(mi2,ni2);
                               pointsum = pointsum + p1*p2;
                            end
                            end
                        end
                    end
                    matrix(p,q) = pointsum;
                end
            end
            
            % find the density matrix for system and enviroment
            % we should first generate a expanded sub site.
            [psi_L,~] = eigs(matrix,1);
            [psi_R,~] = eigs(transpose(matrix),1);
            % normalize
            fac = transpose(psi_L)*psi_R;
            if fac>0
                psi_L = psi_L/sqrt(fac);
                psi_R = psi_R/sqrt(fac);
            else
               psi_L = -psi_L/sqrt(-fac);
               psi_R = psi_R/sqrt(-fac);
            end
                  
            obj.PsiL = psi_L;
            obj.PsiR = psi_R;
            obj.Matrix = matrix;
            % the expanded site: 
            sysnew = Sub(sys,local,'add');
            dms = cell(1,sysnew.totSubSetNo);
            for M = 1:sysnew.totSubSetNo
                tmatrix = zeros(sysnew.vDim(M),sysnew.vDim(M));
               for Mi = 1:sysnew.vDim(M)
                   for Ni = 1:sysnew.vDim(M)             
                       tv = sysnew.cItoOldI{M,Mi};
                       m1 = tv(1); mi1 = tv(2); a2 = tv(3);
                       n1 = tv(1); ni1 = tv(2); b2 = tv(3);
                       sumup = 0;
                       for a1 = 1:dim0
                           m2 = fqn(m1,a1,a2);
                           n2 = fqn(n1,a1,b2);
                           if m2 == n2 && m2~=0
                               for mi2 = 1:vD2(m2)
                                   p = tmItoF(a1,a2,m1,m2,mi1,mi2);
                                   q = tmItoF(a1,b2,n1,m2,ni1,mi2);
                                   sumup = sumup + psi_L(p)*psi_R(q);
                               end
                           end
                       end
                       tmatrix(Mi,Ni) = sumup;
                   end
               end
               dms{M} = tmatrix;
            end
            obj.SysNew = sysnew;
            obj.DmSys = dms;
            
            % the expanded site:
            envnew = Sub(env,local,'add');
            dme = cell(1,envnew.totSubSetNo);
            for M = 1:envnew.totSubSetNo
                tmatrix = zeros(envnew.vDim(M),envnew.vDim(M));
                for Mi = 1:envnew.vDim(M)
                    for Ni = 1:envnew.vDim(M)
                        tv = envnew.cItoOldI{M,Mi};
                        m2 = tv(1); mi2 = tv(2); a2 = tv(3);
                        n2 = tv(1); ni2 = tv(2); b2 = tv(3);
                        sumup = 0;
                        for a1 = 1:dim0
                            m1 = fqn2(m2,a1,a2);
                            n1 = fqn2(n2,a1,b2);
                            if m1 == n1 && m1~=0
                                for mi1 = 1:vD1(m1)
                                    p = tmItoF(a1,a2,m1,m2,mi1,mi2);
                                    q = tmItoF(a1,b2,m1,n2,mi1,ni2);
                                    sumup = sumup + psi_L(p)*psi_R(q);
                                end
                            end
                        end
                        tmatrix(Mi,Ni) = sumup;
                    end
                end
                dme{M} = tmatrix;         
            end
            obj.EnvNew = envnew;
            obj.DmEnv = dme;
            
            % calculate the right vector and left vector for Dm
            % the system
            sysnewright = cell(1,sysnew.totSubSetNo);
            sysnewleft = cell(1,sysnew.totSubSetNo);
            sysneweigvalue = cell(1,sysnew.totSubSetNo);
            for M = 1:sysnew.totSubSetNo
               matrix = dms{M};
               [mR,mL,vd] = fn_eig(matrix,size(matrix,1));
               sysnewright{M} = mR;
               sysnewleft{M} = mL;
               sysneweigvalue{M} = vd;
            end
            
            % the enviroment
            envnewright = cell(1,envnew.totSubSetNo);
            envnewleft = cell(1,envnew.totSubSetNo);
            envneweigvalue = cell(1,envnew.totSubSetNo);
            for M = 1:envnew.totSubSetNo
               matrix = dme{M};
               [mR,mL,vd] = fn_eig(matrix,size(matrix,1));
               envnewright{M} = mR;
               envnewleft{M} = mL;
               envneweigvalue{M} = vd;
            end
            
            obj.SysNewRight = sysnewright;
            obj.SysNewLeft = sysnewleft;
            obj.SysNewEigValue = sysneweigvalue;
            obj.EnvNewRight = envnewright;
            obj.EnvNewLeft = envnewleft;
            obj.EnvNewEigValue = envneweigvalue;
                               
        end
        
    end
    
end

