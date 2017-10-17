function [AngMtx,SelMtx]=AngSelFcnFCN(yvsdFour)
N=size(yvsdFour,3);
COMX=reshape(yvsdFour(1,:,:),N,N)-reshape(yvsdFour(3,:,:),N,N);% H-V
COMY=reshape(yvsdFour(2,:,:),N,N)-reshape(yvsdFour(4,:,:),N,N);% A-D
AngMtx=atan2(COMY,COMX);
SelMtx=sqrt(COMX.^2+COMY.^2);
end
