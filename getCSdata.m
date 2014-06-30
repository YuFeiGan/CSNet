function [img_cs_1d,Theta_1d,Phi,mat_dct_1d]=getCSdata(Src,scale)
% if model==0
%     load temp
    [height,width]=size(Src);
    mat_dct_1d=zeros(height,width);  
    for k=0:1:floor(height/scale)-1 
        dct_1d=cos([0:1:width-1]'*k*pi/width);
        if k>0
            dct_1d=dct_1d-mean(dct_1d); 
        end;
        mat_dct_1d(:,k+1)=dct_1d/norm(dct_1d);
    end

    Phi=randn(floor(width/scale),height);  % only keep one third of the original data  
    Phi = Phi./repmat(sqrt(sum(Phi.^2,1)),[floor(width/scale),1]); % normalize each column
    img_cs_1d=Phi*Src;
    Theta_1d=Phi*mat_dct_1d;
% else
% %     load temp1
%     [height,width]=size(Src);
%     mat_dct_1d=zeros(height,width);  
%     for k=0:1:width-1 
%         dct_1d=cos([0:1:height-1]'*k*pi/height);
%         if k>0
%             dct_1d=dct_1d-mean(dct_1d); 
%         end
%             mat_dct_1d(:,k+1)=dct_1d/norm(dct_1d);
%     end
% %     save temp1
%     Phi=randn(floor(width/scale),height);  % only keep one third of the original data  
%     Phi = Phi./repmat(sqrt(sum(Phi.^2,1)),[floor(width/scale),1]); % normalize each column
%     img_cs_1d=Phi*Src;
%     Theta_1d=Phi*mat_dct_1d;
% %     save temp1
% end
% save temp