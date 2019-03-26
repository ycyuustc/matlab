classdef Sub
    % SUB: the class of quantum structure
    % the Sub class records the information of the quantum structure
    % of the spin (charge) space. 
    
    properties
        vDim;   % the vectore stored the dimension of each subspace
        vSpin;  % the spin value of each subspace
        vDimOld;
        vSpinOld;
        mOldQuanNo;     % the matrix used to 
        StateNoPerSite; % the state number per site 
        totSubSetNo;    % the total number of subsets 
        totStateNo;     % the total number of all states
        mItoF;
        cFtoI;
        mItoF_old;
        cFtoI_old;
        cItoOldI;
        cOldItoI;
        vNewToOld;
        cTruncationMatrix;
    end
    
    %%
    methods
        %% overloaded construction method
        function obj = Sub(varargin)
            if length(varargin) == 1
                obj = obj.Sub_initial(obj,varargin);
            elseif length(varargin) == 5
                % Sub(local,sub,super,'sys/env','up/down')
                obj = obj.Sub_truncate(obj,varargin);
            elseif strcmp(varargin{3}, 'add')
                obj = obj.Sub_add(obj,varargin);
            else
                error('error of initializing Sub object');
            end
        end
        
        %% initial the Sub object from local object, Sub(local);
        function obj = Sub_initial(obj,varargin)
            local = varargin{2}; local = local{1};
            vS = local.vSpin;
            vD = local.vDim;
            
            k = 1;
            while k<=length(vS)
               v_find = find(vS == vS(k)); 
               vS(v_find(2:end)) = [];
               vD(k) = sum(vD(v_find));
               vD(v_find(2:end)) = [];
               k = k + 1;
            end
            
            obj.vSpin = vS;
            obj.vDim = vD;
            obj.totSubSetNo = length(vS);
            obj.StateNoPerSite = local.StateNoPerSite;
            obj.mOldQuanNo = 1:length(vS);
            obj.mItoF = transpose(1:length(vS));
            obj.totStateNo = sum(vD);
            
        end
        
       %% add one site to the Sub object, Sub(subold,local)
        function obj = Sub_add(obj,varargin)
            vc = varargin{2}; oldsub = vc{1}; local = vc{2};
            vS = oldsub.vSpin; vD = oldsub.vDim;
            vS_add = local.vSpin; vD_add = local.vDim;
            
            %%%%%%%  try to combine two spin 
            mS = ones(length(vS_add),1)*vS;
            mS_add = transpose(vS_add)*ones(1,length(vS)); 
            mD = ones(length(vD_add),1)*vD;
            mD_add = transpose(vD_add)*ones(1,length(vD));
            mS_index = ones(length(vS_add),1)*(1:length(vS));
            mS_add_index = transpose(1:length(vS_add))*ones(1,length(vS));
            mOldQN = zeros(length(vS)*length(vS_add),length(vS_add));
            
            m_comb = mS + mS_add;
            m_comb_D = mD.*mD_add;
            m_comb = reshape(m_comb,[1,numel(m_comb)]);
            m_comb_D = reshape(m_comb_D,[1,numel(m_comb_D)]);
            mS_index = reshape(mS_index,[1,numel(mS_index)]);
            mS_add_index = reshape(mS_add_index,[1,numel(mS_add_index)]);
            
            k = 1;
            while k<=length(m_comb)
               v_find = find(m_comb == m_comb(k)); 
               m_comb_D(k) = sum(m_comb_D(v_find));
               mOldQN(mS_add_index(v_find),k) = mS_index(v_find);
               % delete the sites with the same spin, keep the first
               m_comb(v_find(2:end)) = [];
               m_comb_D(v_find(2:end)) = [];
               mS_index(v_find(2:end)) = [];
               mS_add_index(v_find(2:end)) = [];
               k = k + 1;
            end
            
            % find the relation between oldI and newI
            c_OldItoI = cell(length(vD),max(vD),length(vS_add));
            vS_comb = m_comb;
            vD_comb = zeros(1,length(vS_comb));
            c_ItoOldI = cell(length(vS_comb),max(vD)*max(vD_add));
            for n0 = 1:length(vD)
                for ni0 = 1:vD(n0)
                    for s = 1:length(vS_add)
                        ts = vS(n0) + vS_add(s);
                        vfind = find(vS_comb==ts);
                        if ~isempty(vfind)
                        vD_comb(vfind) = vD_comb(vfind) + 1;   
                        c_OldItoI{n0,ni0,s} = [vfind, vD_comb(vfind)];   
                        c_ItoOldI{vfind,vD_comb(vfind)} = [n0,ni0,s];
                        end
                    end
                end
            end
            vS = m_comb; vD = m_comb_D;
            obj.StateNoPerSite = local.StateNoPerSite;
            obj.vSpinOld = oldsub.vSpin;
            obj.vDimOld = oldsub.vDim;
            obj.vSpin = vS;
            obj.vDim = vD;
            obj.mOldQuanNo = mOldQN;
            obj.totSubSetNo = length(m_comb);
            obj.cOldItoI = c_OldItoI;
            obj.cItoOldI = c_ItoOldI;
            
            % update the mapping between I and F
            obj.mItoF_old = obj.mItoF;
            obj.cFtoI_old = obj.cFtoI;
            tmItoF = zeros(length(vS),max(vD));
            tcFtoI = cell(1,length(vS)*max(vD));
            tot_num = 0;
            for k = 1:length(vS)
                for ki = 1:vD(k)
                    tot_num = tot_num + 1;
                    tmItoF(k,ki) = tot_num;
                    tcFtoI{tot_num} = [k,ki];
                end
            end
            % shrinking tcFtoI
            tcFtoI = tcFtoI(1:tot_num);
            
            obj.mItoF = tmItoF;
            obj.cFtoI = tcFtoI;
            obj.totStateNo = sum(obj.vDim);
            
            
        end
        
        %% truncation of the Sub
        function obj = Sub_truncate(obj,varargin)
            vc = varargin{2}; local = vc{1};
            subold = vc{2}; super = vc{3}; 
            sys_or_env = vc{4}; up_or_down = vc{5};
            
            if strcmp(sys_or_env,'sys')
               right = super.SysNewRight;
               left = super.SysNewLeft;
               eig_value = super.SysNewEigValue;
            end
            
            if strcmp(sys_or_env,'env')
                right = super.EnvNewRight;
                left = super.EnvNewLeft;
                eig_value = super.EnvNewEigValue;
            end
            
            msort = zeros(3,subold.totStateNo);
            point = 0;
            for m = 1:subold.totSubSetNo
                for mi = 1:subold.vDim(m)
                    point = point + 1;
                    msort(1,point) = m;
                    msort(2,point) = mi;
                    tv = eig_value{m};
                    msort(3,point) = tv(mi);
                end
            end
            
            [~,vsort] = sort(abs(msort(3,:)),'descend');
            msort = msort(:,vsort);
            vdimnew = zeros(1,subold.totSubSetNo);
            vspinnew = subold.vSpin;
            for k = 1:local.StatesKeeped
                vdimnew(msort(1,k)) = vdimnew(msort(1,k)) + 1;
            end
            vfind = find(vdimnew == 0);
            vdimnew(vfind) = [];
            vspinnew(vfind) = [];
            vnewtoold = 1:subold.totSubSetNo;
            vnewtoold(vfind) = [];
            obj.vDim = vdimnew;
            obj.vSpin = vspinnew;
            obj.vNewToOld = vnewtoold;
            obj.totSubSetNo = length(vdimnew);
            truncationmatrix = cell(1,obj.totSubSetNo);
            for k = 1:obj.totSubSetNo
                if strcmp(up_or_down,'up')
                    tt = left{vnewtoold(k)};
                    tt = tt(:,1:vdimnew(k));
                    truncationmatrix{k} = tt;
                end
                if strcmp(up_or_down,'down')
                    tt = right{vnewtoold(k)};
                    tt = tt(:,1:vdimnew(k));
                    truncationmatrix{k} = tt;
                end          
            end
            obj.cTruncationMatrix = truncationmatrix;
            
            
        end
       
        
    end
    
end

