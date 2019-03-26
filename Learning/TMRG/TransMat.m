classdef TransMat
    %TRANSMAT 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        Local;
        Sub_up;
        Sub_down;
        Tm;
        FindQuanNo_down;
        FindQuanNo_up;
        Type;
    end
    
    methods
        function obj = TransMat(varargin)
            if length(varargin) == 2
                % TransMat(local,sub);
                obj = obj.TransMat_initial(obj,varargin);
            elseif length(varargin) == 5
                % TransMat(local,subin,tm_down,tm_up,'sys')
                if strcmp(varargin{5},'sys')
                    obj = obj.TransMat_systm(obj,varargin);
                end
                if strcmp(varargin{5},'env')
                    obj = obj.TransMat_envtm(obj,varargin);
                end
            elseif length(varargin) == 4
               % TransMat(local,subin,systm,'sysup')
               if strcmp(varargin{4},'sysup')
                   obj = obj.TransMat_sysup(obj,varargin);
               end
               if strcmp(varargin{4},'sysdown')
                   obj = obj.TransMat_sysdown(obj,varargin);
               end
               if strcmp(varargin{4},'envup')
                   obj = obj.TransMat_envup(obj,varargin);
               end
               if strcmp(varargin{4},'envdown')
                   obj = obj.TransMat_envdown(obj,varargin);
               end
               
               % trancation 
               % TransMat(local,tm0,subnew,up_or_down)
               if length(varargin) ==4
               if strcmp(varargin{4},'up') || strcmp(varargin{4},'down')
                   obj = obj.TransMat_truncate_type1(obj,varargin);
               end
               end
               
               if length(varargin) == 4 && isa(varargin{2},'TransMat')...
                   && isa(varargin{4},'Sub')
               % TransMat(local,tm0,subup,subdown)
                   obj = obj.TransMat_truncate_type2(obj,varargin);
               end
               
            end
        end
        
        
        function obj = TransMat_initial(obj,varargin)
            vc = varargin{2}; local = vc{1}; subin = vc{2};
            if local.StateNoPerSite ~= subin.StateNoPerSite
                error('wrong initialization');
            end
            obj.Type = 1;
            tm_ori = local.LocalTensor;
            obj.Sub_up = subin;
            obj.Sub_down = subin;
            tmItoF = subin.mItoF;
            
            vS = subin.vSpin; vD = subin.vDim;
            vSa = local.vSpin; vDa = local.vDim;
            DimSmall = local.StateNoPerSite;
            
            % calculate the FindQuanNo
            findquanno_down = zeros(subin.totSubSetNo,DimSmall,DimSmall);
            findquanno_up = zeros(subin.totSubSetNo,DimSmall,DimSmall);
            for m = 1:subin.totSubSetNo
                for a1 = 1:DimSmall
                    for b1 = 1:DimSmall
                        for n = 1:subin.totSubSetNo
                            if vSa(a1)+vS(m)==vSa(b1)+vS(n)
                                findquanno_down(m,a1,b1) = n;
                                findquanno_up(n,a1,b1) = m;
                            end
                        end
                    end
                end
            end
            obj.FindQuanNo_down = findquanno_down;
            obj.FindQuanNo_up = findquanno_up;
            
            % calculate TM, stored by cell Tm[n][a1][b1], 
            tm = cell(subin.totSubSetNo,DimSmall,DimSmall);
            for a1 = 1:DimSmall
                for b1 = 1:DimSmall
                    for m = 1:subin.totSubSetNo
                        n = findquanno_down(m,a1,b1);
                        if n~=0
                            tmatrix = zeros(vD(m),vD(n));
                            for ni = 1:vD(n)
                                for mi = 1:vD(m)
                                    q = tmItoF(n,ni);
                                    p = tmItoF(m,mi);
                                    tmatrix(mi,ni) = tm_ori(a1,b1,p,q);
                                end
                            end
                            tm{m,a1,b1} = tmatrix;
                        end
                    end
                end
            end
            obj.Tm = tm;
            
        end
        
        
        function obj = TransMat_systm(obj,varargin)
            vc = varargin{2};local = vc{1};
            subin = vc{2}; tm_down = vc{3}; tm_up = vc{4};  
            obj.Type = 2;
            obj.Sub_up = subin;
            obj.Sub_down = subin;
            submid = tm_down.Sub_up;
            
            vS = subin.vSpin; vD = subin.vDim;
            vSa = local.vSpin; vDa = local.vDim;
            DimSmall = local.StateNoPerSite;
            
            % calculate the FindQuanNo
            findquanno_down = zeros(subin.totSubSetNo,DimSmall,DimSmall,...
                DimSmall,DimSmall);
            findquanno_up = zeros(subin.totSubSetNo,DimSmall,DimSmall,...
                DimSmall,DimSmall);
            for m = 1:subin.totSubSetNo
                for a1 = 1:DimSmall
                    for b1 = 1:DimSmall
                        for a2 = 1:DimSmall
                            for b2 = 1:DimSmall
                                for n = 1:subin.totSubSetNo
                                    if vSa(a1)+vSa(a2)+vS(m)...
                                            == vSa(b1)+vSa(b2)+vS(n)
                                        findquanno_down(m,a1,b1,a2,b2) = n;
                                        findquanno_up(n,a1,b1,a2,b2) = m;
                                    end
                                end
                            end
                        end
                    end
                end
            end
            obj.FindQuanNo_down = findquanno_down;
            obj.FindQuanNo_up = findquanno_up;
            
            % calculate TM, stored by cell Tm[m][a1][b1][a2][b2], 
            tm = cell(subin.totSubSetNo,DimSmall,DimSmall,...
                DimSmall,DimSmall);
            for a1 = 1:DimSmall
            for b1 = 1:DimSmall
            for a2 = 1:DimSmall
            for b2 = 1:DimSmall
            for m = 1:subin.totSubSetNo
                n = findquanno_down(m,a1,b1,a2,b2);
                if n~=0
                    tmatrix = zeros(vD(m),vD(n));
                    for ni = 1:vD(n)
                        for mi = 1:vD(m)
                           tv = subin.cItoOldI{m,mi};
                           m0 = tv(1); mi0 = tv(2); a0 = tv(3);
                           tv = subin.cItoOldI{n,ni};
                           n0 = tv(1); ni0 = tv(2); b0 = tv(3);
                           sumup = 0;
                           for amid = 1:DimSmall
                              mmid = tm_up.FindQuanNo_down(m0,a0,amid);
                              mmid2 = tm_down.FindQuanNo_up(n0,a1,b1);
                              if mmid~=0 && mmid==mmid2
                                  for mmidi = 1:submid.vDim(mmid)
                                     tm1 = tm_down.Tm{mmid,a1,b1};
                                     p1 = tm1(mmidi,ni0);
                                     
                                     tm2 = tm_up.Tm{m0,a0,amid};
                                     p2 = tm2(mi0,mmidi);
                                     
                                     tm3 = local.LocalTensor;
                                     p3 = tm3(amid,b0,a2,b2);
                                     
                                     sumup = sumup + p1*p2*p3;
                                     
                                  end
                              end
                           end
                           tmatrix(mi,ni) = sumup;
                        end
                    end
                    tm{m,a1,b1,a2,b2} = tmatrix;
                end
            end
            end
            end
            end
            end
            obj.Tm = tm;
            
        end
        
        
        function obj = TransMat_envtm(obj,varargin)
            vc = varargin{2};local = vc{1};
            subin = vc{2}; tm_down = vc{3}; tm_up = vc{4};  
            obj.Type = 2;
            obj.Sub_up = subin;
            obj.Sub_down = subin;
            submid = tm_down.Sub_up;
            
            vS = subin.vSpin; vD = subin.vDim;
            vSa = local.vSpin; vDa = local.vDim;
            DimSmall = local.StateNoPerSite;
            
            % calculate the FindQuanNo
            findquanno_down = zeros(subin.totSubSetNo,DimSmall,DimSmall,...
                DimSmall,DimSmall);
            findquanno_up = zeros(subin.totSubSetNo,DimSmall,DimSmall,...
                DimSmall,DimSmall);
            for m = 1:subin.totSubSetNo
                for a1 = 1:DimSmall
                    for b1 = 1:DimSmall
                        for a2 = 1:DimSmall
                            for b2 = 1:DimSmall
                                for n = 1:subin.totSubSetNo
                                    if vSa(a1)+vSa(a2)+vS(m)...
                                            == vSa(b1)+vSa(b2)+vS(n)
                                        findquanno_down(m,a1,b1,a2,b2) = n;
                                        findquanno_up(n,a1,b1,a2,b2) = m;
                                    end
                                end
                            end
                        end
                    end
                end
            end
            obj.FindQuanNo_down = findquanno_down;
            obj.FindQuanNo_up = findquanno_up;
            
            % calculate TM, stored by cell Tm[m][a1][b1][a2][b2], 
            tm = cell(subin.totSubSetNo,DimSmall,DimSmall,...
                DimSmall,DimSmall);
            for a1 = 1:DimSmall
            for b1 = 1:DimSmall
            for a2 = 1:DimSmall
            for b2 = 1:DimSmall
            for m = 1:subin.totSubSetNo
                n = findquanno_down(m,a1,b1,a2,b2);
                if n~=0
                    tmatrix = zeros(vD(m),vD(n));
                    for ni = 1:vD(n)
                        for mi = 1:vD(m)
                           tv = subin.cItoOldI{m,mi};
                           m0 = tv(1); mi0 = tv(2); a0 = tv(3);
                           tv = subin.cItoOldI{n,ni};
                           n0 = tv(1); ni0 = tv(2); b0 = tv(3);
                           sumup = 0;
                           for amid = 1:DimSmall
                              mmid = tm_up.FindQuanNo_down(m0,a2,b2);
                              mmid2 = tm_down.FindQuanNo_up(n0,amid,b0);
                              if mmid~=0 && mmid==mmid2
                                  for mmidi = 1:submid.vDim(mmid)
                                     tm1 = local.LocalTensor;
                                     p1 = tm1(a1,b1,a0,amid);
                                     
                                     tm2 = tm_down.Tm{mmid,amid,b0};
                                     p2 = tm2(mmidi,ni0);
                                     
                                     tm3 = tm_up.Tm{m0,a2,b2};
                                     p3 = tm3(mi0,mmidi);
                                     
                                     sumup = sumup + p1*p2*p3;
                                     
                                  end
                              end
                           end
                           tmatrix(mi,ni) = sumup;
                        end
                    end
                    tm{m,a1,b1,a2,b2} = tmatrix;
                end
            end
            end
            end
            end
            end
            obj.Tm = tm;
            
        end
            
        
        function obj = TransMat_sysdown(obj,varargin)
            vc = varargin{2}; local = vc{1};
            subin = vc{2}; tm0 = vc{3};
            obj.Type = 1;
            obj.Sub_up = subin;
            obj.Sub_down = subin;
            
            vS = subin.vSpin; vD = subin.vDim;
            vSa = local.vSpin; vDa = local.vDim;
            DimSmall = local.StateNoPerSite;
            
            % calculate the FindQuanNo
            findquanno_down = zeros(subin.totSubSetNo,DimSmall,DimSmall);
            findquanno_up = zeros(subin.totSubSetNo,DimSmall,DimSmall);
            for m = 1:subin.totSubSetNo
                for a1 = 1:DimSmall
                    for b1 = 1:DimSmall
                        for n = 1:subin.totSubSetNo
                            if vSa(a1)+vS(m) == vSa(b1)+vS(n)
                               findquanno_down(m,a1,b1) = n;
                               findquanno_up(n,a1,b1) = m;
                            end
                        end
                    end
                end
            end
            obj.FindQuanNo_down = findquanno_down;
            obj.FindQuanNo_up = findquanno_up;
            
            tm = cell(subin.totSubSetNo,DimSmall,DimSmall);
            for a1 = 1:DimSmall
            for b1 = 1:DimSmall
            for m = 1:subin.totSubSetNo
                n = findquanno_down(m,a1,b1);
                if n~=0
                    tmatrix = zeros(vD(m),vD(n));
                    for ni = 1:vD(n)
                        for mi = 1:vD(m)
                            tv = subin.cItoOldI{m,mi};
                            m0 = tv(1); mi0 = tv(2); a2 = tv(3);
                            tv = subin.cItoOldI{m,mi};
                            n0 = tv(1); ni0 = tv(2); b2 = tv(3);
                            ttm = tm0.Tm{m0,a1,b1,a2,b2};
                            if n0 == tm0.FindQuanNo_down(m0,a1,b1,a2,b2)
                                tmatrix(mi,ni) = ttm(mi0,ni0);
                            end
                        end
                    end
                    tm{m,a1,b1} = tmatrix;
                end         
            end
            end
            end
            
            obj.Tm = tm;
            
        end
            
        
        function obj = TransMat_sysup(obj,varargin)
            vc = varargin{2}; local = vc{1};
            subin = vc{2}; tm0 = vc{3};
            obj.Type = 1;
            obj.Sub_up = subin;
            obj.Sub_down = subin;
            
            vS = subin.vSpin; vD = subin.vDim;
            vSa = local.vSpin; vDa = local.vDim;
            DimSmall = local.StateNoPerSite;
            
            % calculate the FindQuanNo
            findquanno_down = zeros(subin.totSubSetNo,DimSmall,DimSmall);
            findquanno_up = zeros(subin.totSubSetNo,DimSmall,DimSmall);
            for m = 1:subin.totSubSetNo
                for a1 = 1:DimSmall
                    for b1 = 1:DimSmall
                        for n = 1:subin.totSubSetNo
                            if vSa(a1)+vS(m) == vSa(b1)+vS(n)
                               findquanno_down(m,a1,b1) = n;
                               findquanno_up(n,a1,b1) = m;
                            end
                        end
                    end
                end
            end
            obj.FindQuanNo_down = findquanno_down;
            obj.FindQuanNo_up = findquanno_up;
            
            tm = cell(subin.totSubSetNo,DimSmall,DimSmall);
            for a1 = 1:DimSmall
            for b1 = 1:DimSmall
            for m = 1:subin.totSubSetNo
                n = findquanno_down(m,a1,b1);
                if n~=0
                    tmatrix = zeros(vD(m),vD(n));
                    for ni = 1:vD(n)
                        for mi = 1:vD(m)
                            tv = subin.cItoOldI{m,mi};
                            m0 = tv(1); mi0 = tv(2); a2 = tv(3);
                            tv = subin.cItoOldI{m,mi};
                            n0 = tv(1); ni0 = tv(2); b2 = tv(3);
                            ttm = local.LocalTensor;
                            if m0 == n0 && mi0 == ni0
                                tmatrix(mi,ni) = ttm(a2,b2,a1,b1);
                            end
                        end
                    end
                    tm{m,a1,b1} = tmatrix;
                end
                
            end
            end
            end
            
            obj.Tm = tm;
            
        end
            
        
        function obj = TransMat_envdown(obj,varargin)
            vc = varargin{2}; local = vc{1};
            subin = vc{2}; tm0 = vc{3};
            obj.Type = 1;
            obj.Sub_up = subin;
            obj.Sub_down = subin;
            
            vS = subin.vSpin; vD = subin.vDim;
            vSa = local.vSpin; vDa = local.vDim;
            DimSmall = local.StateNoPerSite;
            
            % calculate the FindQuanNo
            findquanno_down = zeros(subin.totSubSetNo,DimSmall,DimSmall);
            findquanno_up = zeros(subin.totSubSetNo,DimSmall,DimSmall);
            for m = 1:subin.totSubSetNo
                for a1 = 1:DimSmall
                    for b1 = 1:DimSmall
                        for n = 1:subin.totSubSetNo
                            if vSa(a1)+vS(m) == vSa(b1)+vS(n)
                               findquanno_down(m,a1,b1) = n;
                               findquanno_up(n,a1,b1) = m;
                            end
                        end
                    end
                end
            end
            obj.FindQuanNo_down = findquanno_down;
            obj.FindQuanNo_up = findquanno_up;
            
            tm = cell(subin.totSubSetNo,DimSmall,DimSmall);
            for a1 = 1:DimSmall
            for b1 = 1:DimSmall
            for m = 1:subin.totSubSetNo
                n = findquanno_down(m,a1,b1);
                if n~=0
                    tmatrix = zeros(vD(m),vD(n));
                    for ni = 1:vD(n)
                        for mi = 1:vD(m)
                            tv = subin.cItoOldI{m,mi};
                            m0 = tv(1); mi0 = tv(2); a2 = tv(3);
                            tv = subin.cItoOldI{m,mi};
                            n0 = tv(1); ni0 = tv(2); b2 = tv(3);
                            ttm = local.LocalTensor;
                            if m0 == n0 && mi0 == ni0
                                tmatrix(mi,ni) = ttm(a1,b1,a2,b2);
                            end
                        end
                    end
                    tm{m,a1,b1} = tmatrix;
                end
                
            end
            end
            end
            
            obj.Tm = tm;
            
        end
        
           
        function obj = TransMat_envup(obj,varargin)
            vc = varargin{2}; local = vc{1};
            subin = vc{2}; tm0 = vc{3};
            obj.Type = 1;
            obj.Sub_up = subin;
            obj.Sub_down = subin;
            
            vS = subin.vSpin; vD = subin.vDim;
            vSa = local.vSpin; vDa = local.vDim;
            DimSmall = local.StateNoPerSite;
            
            % calculate the FindQuanNo
            findquanno_down = zeros(subin.totSubSetNo,DimSmall,DimSmall);
            findquanno_up = zeros(subin.totSubSetNo,DimSmall,DimSmall);
            for m = 1:subin.totSubSetNo
                for a1 = 1:DimSmall
                    for b1 = 1:DimSmall
                        for n = 1:subin.totSubSetNo
                            if vSa(a1)+vS(m) == vSa(b1)+vS(n)
                               findquanno_down(m,a1,b1) = n;
                               findquanno_up(n,a1,b1) = m;
                            end
                        end
                    end
                end
            end
            obj.FindQuanNo_down = findquanno_down;
            obj.FindQuanNo_up = findquanno_up;
            
            tm = cell(subin.totSubSetNo,DimSmall,DimSmall);
            for a1 = 1:DimSmall
            for b1 = 1:DimSmall
            for m = 1:subin.totSubSetNo
                n = findquanno_down(m,a1,b1);
                if n~=0
                    tmatrix = zeros(vD(m),vD(n));
                    for ni = 1:vD(n)
                        for mi = 1:vD(m)
                            tv = subin.cItoOldI{m,mi};
                            m0 = tv(1); mi0 = tv(2); a2 = tv(3);
                            tv = subin.cItoOldI{m,mi};
                            n0 = tv(1); ni0 = tv(2); b2 = tv(3);
                            ttm = tm0.Tm{m0,a2,b2,a1,b1};
                            if n0 == tm0.FindQuanNo_down(m0,a2,b2,a1,b1)
                                tmatrix(mi,ni) = ttm(mi0,ni0);
                            end
                        end
                    end
                    tm{m,a1,b1} = tmatrix;
                end
                
            end
            end
            end
            
            obj.Tm = tm;
            
        end
        
        
        function obj = TransMat_truncate_type1(obj,varargin)
            % TransMat(local,tm0,subnew,up_or_down)
            vc = varargin{2}; local = vc{1};
            tm0 = vc{2}; subnew = vc{3}; up_or_down = vc{4};
            
            obj.Local = local;
            
            vSa = local.vSpin; vDa = local.vDim;
            vS1 = tm0.Sub_up.vSpin;
            vS2 = tm0.Sub_down.vSpin;
            vD1 = tm0.Sub_up.vDim;
            vD2 = tm0.Sub_down.vDim;
            vS = subnew.vSpin;
            vD = subnew.vDim;
            DimSmall = local.StateNoPerSite;
            
            if strcmp(up_or_down,'up')
                obj.Sub_up = subnew;
                obj.Sub_down = tm0.Sub_down;
            % calculate the FindQuanNo
            findquanno_down = zeros(subnew.totSubSetNo,DimSmall,DimSmall);
            findquanno_up = ...
                zeros(tm0.Sub_down.totSubSetNo,DimSmall,DimSmall);
            for m = 1:subnew.totSubSetNo
                for a1 = 1:DimSmall
                    for b1 = 1:DimSmall
                        for n = 1:tm0.Sub_down.totSubSetNo
                            if vSa(a1)+vS(m) == vSa(b1)+vS2(n)
                               findquanno_down(m,a1,b1) = n;
                               findquanno_up(n,a1,b1) = m;
                            end
                        end
                    end
                end
            end
            obj.FindQuanNo_down = findquanno_down;
            obj.FindQuanNo_up = findquanno_up;
            tm = cell(subnew.totSubSetNo,DimSmall,DimSmall);
            for a1 = 1:DimSmall
            for b1 = 1:DimSmall
            for m = 1:subnew.totSubSetNo
                n = findquanno_down(m,a1,b1);
                if n~=0
                   ttm = tm0.Tm{subnew.vNewToOld(m),a1,b1};
                   convm = subnew.cTruncationMatrix{m};
                   matrix = transpose(convm)*ttm;
                   tm{m,a1,b1} = matrix;
                end
            end
            end
            end
            
            
            elseif strcmp(up_or_down,'down')
                obj.Sub_up = tm0.Sub_up;
                obj.Sub_down = subnew;
                
                % calculate the FindQuanNo
            findquanno_down = ...
                zeros(tm0.Sub_up.totSubSetNo,DimSmall,DimSmall);
            findquanno_up = zeros(subnew.totSubSetNo,DimSmall,DimSmall);
            for m = 1:tm0.Sub_up.totSubSetNo
                for a1 = 1:DimSmall
                    for b1 = 1:DimSmall
                        for n = 1:subnew.totSubSetNo
                            if vSa(a1)+vS1(m) == vSa(b1)+vS(n)
                               findquanno_down(m,a1,b1) = n;
                               findquanno_up(n,a1,b1) = m;
                            end
                        end
                    end
                end
            end
            obj.FindQuanNo_down = findquanno_down;
            obj.FindQuanNo_up = findquanno_up;
            tm = cell(tm0.Sub_up.totSubSetNo,DimSmall,DimSmall);
            for a1 = 1:DimSmall
            for b1 = 1:DimSmall
            for m = 1:tm0.Sub_up.totSubSetNo
                n = findquanno_down(m,a1,b1);
                if n~=0
                   ttm = tm0.Tm{m,a1,b1};
                   convm = subnew.cTrancationMatrix{m};
                   matrix = ttm * convm;
                   tm{m,a1,b1} = matrix;
                end
            end
            end
            end
                
            end
            
        end
        
        
        function obj = TransMat_truncate_type2(obj,varargin)
            % TransMat(local,tm0,subup,subdown)
            vc = varargin{2}; local = vc{1};
            tm0 = vc{2}; subup = vc{3}; subdown = vc{4};
            
            obj.Local = local;
            
            vSa = local.vSpin; vDa = local.vDim;
            vS1 = tm0.Sub_up.vSpin;
            vS2 = tm0.Sub_down.vSpin;
            vD1 = tm0.Sub_up.vDim;
            vD2 = tm0.Sub_down.vDim;
            vSnew1 = subup.vSpin;
            vSnew2 = subdown.vSpin;
            DimSmall = local.StateNoPerSite;
            
            
            % calculate the FindQuanNo
            findquanno_down = ...
                zeros(subup.totSubSetNo,DimSmall,DimSmall,DimSmall,DimSmall);
            findquanno_up = ...
                zeros(subdown.totSubSetNo,DimSmall,DimSmall);
            for m = 1:subup.totSubSetNo
            for a1 = 1:DimSmall
            for b1 = 1:DimSmall
            for a2 = 1:DimSmall
            for b2 = 1:DimSmall
            for n = 1:subdown.totSubSetNo
                if vSa(a1)+vSa(a2)+vSnew1(m) == vSa(b1)+vSa(b2)+vSnew2(n)
                    findquanno_down(m,a1,b1,a2,b2) = n;
                    findquanno_up(n,a1,b1,a2,b2) = m;
                end
            end
            end
            end
            end
            end
            end
            obj.FindQuanNo_down = findquanno_down;
            obj.FindQuanNo_up = findquanno_up;
            
            tm = cell(subup.totSubSetNo,DimSmall,DimSmall,DimSmall,DimSmall);
            for a1 = 1:DimSmall
            for b1 = 1:DimSmall
            for a2 = 1:DimSmall
            for b2 = 1:DimSmall
            for m = 1:subup.totSubSetNo
                n = findquanno_down(m,a1,b1,a2,b2);
                if n~=0
                   ttm = tm0.Tm{subup.vNewToOld(m),a1,b1,a2,b2};
                   convmup = subup.cTruncationMatrix{m};
                   convmdown = subdown.cTruncationMatrix{n};
                   matrix = transpose(convmup)*ttm*convmdown;
                   tm{m,a1,b1,a2,b2} = matrix;
                end
            end
            end
            end
            end
            end
            
            obj.Tm = tm;
            obj.Sub_up = subup;
            obj.Sub_down = subdown;
            obj.Type = 2; 
            
        end
           
        
    end  % End of methods
    
    
end  % End of class

