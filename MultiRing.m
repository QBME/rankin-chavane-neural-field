function [out,varargout]=MultiRing(r,RingDis,RingAmp,RingWidth)

if size(RingDis)~=size(RingAmp) 
    error('bad ring dis or amp')
end

out=zeros(size(r));
GaussRing=@(r,s,lam) exp(-((r-lam).^2)/2/s^2);
GaussRingZeroMode=@(s,lam)2*pi*s^2*exp(-lam^2/2/s^2)+pi*s*lam*sqrt(2*pi)*(1+erf(lam*sqrt(2)/2/s));
rescalefac=0;
for i=1:length(RingDis)
    rescalefac=rescalefac+RingAmp(i)*GaussRingZeroMode(RingWidth,RingDis(i));
    out=out+RingAmp(i)*GaussRing(r,RingWidth,RingDis(i));    
end
out=out/rescalefac;
if nargout==2
    varargout{1}=rescalefac;
end

