function [X,numindX]=fn_contract(X,numindX,v_indX,Y,numindY,v_indY)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明

v_Xsize=ones(1,numindX);
v_Xsize(1:length(size(X)))=size(X);
v_Ysize=ones(1,numindY);
v_Ysize(1:length(size(Y)))=size(Y);

v_indX_left=1:numindX;
v_indX_left(v_indX)=[];

v_indY_right=1:numindY;
v_indY_right(v_indY)=[];

v_sizeX_left=v_Xsize(v_indX_left);
v_sizeY_right=v_Ysize(v_indY_right);
v_sizeX=v_Xsize(v_indX);
v_sizeY=v_Ysize(v_indY);

if prod(v_sizeX)~=prod(v_sizeY)
    error('indX and indY are not of the same dimension');
end

if isempty(v_indY_right)
    if isempty(v_indX_left)
        
        X=permute(X,v_indX);
        X=reshape(X,[1,prod(v_sizeX)]);
        
        Y=permute(Y,v_indY);
        Y=reshape(Y,[prod(v_sizeY),1]);
        
        X=X*Y;
        
        v_Xsize=1;
        X=reshape(X,[v_Xsize,1]);
        
        return
        
    else
        
        X=permute(X,[v_indX_left,v_indX]);
        X=reshape(X,[prod(v_sizeX_left),prod(v_sizeX)]);
        Y=permute(Y,v_indY);
        Y=reshape(Y,[prod(v_sizeY),1]);
        
        X=X*Y;
        v_Xsize=v_sizeX_left;
        
        X=reshape(X,[v_Xsize,1]);
        
        return
        
    end
    
end

X=permute(X,[v_indX_left,v_indX]);
X=reshape(X,[prod(v_sizeX_left),prod(v_sizeX)]);

Y=permute(Y,[v_indY,v_indY_right]);
Y=reshape(Y,[prod(v_sizeY),prod(v_sizeY_right)]);

X=X*Y;

v_Xsize=[v_sizeX_left,v_sizeY_right];

numindX=length(v_Xsize);

X=reshape(X,[v_Xsize,1]);

end



        
        
        
        
        
        
        
      

