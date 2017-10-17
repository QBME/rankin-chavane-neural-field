function [yvsdFour,yvsdMax]=utoVSD_FCN(L,N,S,mu,svStruct,VKHatLoc,VKHatLat,SKHat,JStructLat,JStructLoc)

ufAllH=   svStruct.ufAllH;
ufAllA=   svStruct.ufAllA;
ufAllV=   svStruct.ufAllV;
ufAllD=   svStruct.ufAllD;
yvsdFour=zeros(size(ufAllH));
yvsdMax=zeros(size(ufAllH,1),1);
for i=1:size(ufAllH,1)% run through by OR
    vsdH=(2*L/N)^2*ifftshift(real(ifft2(VKHatLat .* fft2(S(mu*reshape(ufAllH(i,:,:),N,N)))))).*JStructLat{1}+...
        (2*L/N)^2*ifftshift(real(ifft2(VKHatLoc .* fft2(S(mu*reshape(ufAllH(i,:,:),N,N)))))).*JStructLoc{1};
    vsdA=(2*L/N)^2*ifftshift(real(ifft2(VKHatLat .* fft2(S(mu*reshape(ufAllA(i,:,:),N,N)))))).*JStructLat{2}+...
        (2*L/N)^2*ifftshift(real(ifft2(VKHatLoc .* fft2(S(mu*reshape(ufAllA(i,:,:),N,N)))))).*JStructLoc{2};
    vsdV=(2*L/N)^2*ifftshift(real(ifft2(VKHatLat .* fft2(S(mu*reshape(ufAllV(i,:,:),N,N)))))).*JStructLat{3}+...
        (2*L/N)^2*ifftshift(real(ifft2(VKHatLoc .* fft2(S(mu*reshape(ufAllV(i,:,:),N,N)))))).*JStructLoc{3};
    vsdD=(2*L/N)^2*ifftshift(real(ifft2(VKHatLat .* fft2(S(mu*reshape(ufAllD(i,:,:),N,N)))))).*JStructLat{4}+...
        (2*L/N)^2*ifftshift(real(ifft2(VKHatLoc .* fft2(S(mu*reshape(ufAllD(i,:,:),N,N)))))).*JStructLoc{4};
    yvsdFour(i,:,:)=(2*L/N)^2*ifftshift(real(ifft2(SKHat .* fft2(vsdH+vsdA+vsdV+vsdD))));
    yvsdMax(i)=max(max((yvsdFour(i,:,:))));
end
%% 

return
% function [yvsdFour,yvsdMax]=SutoVSD_BR_LL_FCN(L,N,Su,VKHatLoc,VKHatLat,SKHat,JStruct)
% yvsdFour=zeros(size(Su));
% yvsdMax=zeros(size(Su,1),1);
% for i=1:size(Su,1)
%     SuHat=fft2(reshape(Su(i,:,:),N,N));
%     yvsdtmp=(2*L/N)^2*ifftshift(real(ifft2(VKHatLat .* SuHat))).*JStruct{i}+...
%         (2*L/N)^2*ifftshift(real(ifft2(VKHatLoc .* SuHat)));
%     yvsdFour(i,:,:)=(2*L/N)^2*ifftshift(real(ifft2(SKHat .* fft2(yvsdtmp))));
%     yvsdMax(i)=max(max((yvsdFour(i,:,:))));
% end
% end
