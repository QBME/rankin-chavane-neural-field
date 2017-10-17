function [w,wLoc,wLat,WEFun,WEFunLat,WEFunLoc,WIFun,varargout]=BuildMultiRingFcns(L,N,R,X,Y,lambda,C,EnvelopeDecayConstEx,RingWidthEx,RingWidthIn,ifplot)

EnvelopeFcn=@(x,s)exp(-x./s);
EnvelopeDecayConstIn=lambda;
RingDisEx=0:lambda:2*lambda;

RingDisExLoc=[0];
RingDisExLat=[lambda,2*lambda];
RingDisIn=0;


RingAmpEx=EnvelopeFcn(RingDisEx,EnvelopeDecayConstEx);
RingAmpExLoc=EnvelopeFcn(RingDisExLoc,EnvelopeDecayConstEx);
RingAmpExLat=EnvelopeFcn(RingDisExLat,EnvelopeDecayConstEx);
RingAmpIn=EnvelopeFcn(RingDisIn,EnvelopeDecayConstIn);
% RingWidthIn=lambda;
WEFun=@(r)MultiRing(r,RingDisEx,RingAmpEx,RingWidthEx);
[~,WEFunRS]=MultiRing(R,RingDisEx,RingAmpEx,RingWidthEx);
[~,WEFunLocRS]=MultiRing(R,RingDisExLoc,RingAmpExLoc,RingWidthEx);
[~,WEFunLatRS]=MultiRing(R,RingDisExLat,RingAmpExLat,RingWidthEx);
WEFunLoc=@(r)WEFunLocRS/WEFunRS*MultiRing(r,RingDisExLoc,RingAmpExLoc,RingWidthEx);
WEFunLat=@(r)WEFunLatRS/WEFunRS*MultiRing(r,RingDisExLat,RingAmpExLat,RingWidthEx);
WIFun=@(r)MultiRing(r,RingDisIn,RingAmpIn,RingWidthIn);
x=linspace(0,L,N);        kmax=5;        k=linspace(0,kmax,N);
WRadial=@(r)(WEFunLoc(r)-(1-C)*WIFun(r)+WEFunLat(r));
WRadHat=ComputeFTRadialHankel(WRadial,x,k);
max_mode=max(WRadHat);
P=  8.8647/max_mode;
w=@(x,y)P*(WEFun(sqrt(x.^2+y.^2))-(1-C)*WIFun(sqrt(x.^2+y.^2)));
wLoc=@(x,y)P*(WEFunLoc(sqrt(x.^2+y.^2))-(1-C)*WIFun(sqrt(x.^2+y.^2)));
wLat=@(x,y)P*(WEFunLat(sqrt(x.^2+y.^2)));


if nargout>=8
    varargout{1}=P;
end
    
    
if max(max(WEFun(R)-(WEFunLoc(R)+WEFunLat(R))))>1e-6 || max(max(w(X,Y)-(wLoc(X,Y)+wLat(X,Y))))>1e-6
    error('Loc-lat decomposition broken')
end

if ifplot
    figure(23);clf; hold on
    subplot(2,1,1); hold on
plot(x,WRadial(x),'k','linewidth',3)
title('Spatial domain')
plot(x,(WEFunLoc(x)+WEFunLat(x)),'b','linewidth',3)
plot(x,(WIFun(x)),'r','linewidth',3)
legend('W=Exc-Inh','Exc','Inh')
plot([0 60],[0 0],'k-','linewidth',2)
set(gca,'xlim',[0 15])
set(gca,'xtick',[0 2*pi 4*pi,6*pi],'xticklabel',[0 1 2 3])
set(gca,'fontname','helvetica','fontsize',16,'linewidth',2);
% set(gcf,'units','centimeters','position',[20 24,10,8]);
set(gcf,'color','w');box on
% set(gca,'position',[0.15 0.1,.80,.85]);
% set(gca,'ytick',[-0.02:0.01:0.04]);
% set(gca,'ylim',[-0.01 0.04]);
% plot([lambda/2 lambda/2],[-5 5],'k-','linewidth',2)
% plot([lambda lambda],[-5 5],'k--','linewidth',2)
% plot([lambda*2 lambda*2],[-5 5],'k--','linewidth',2)
% figure(24);clf; hold on
    subplot(2,1,2);hold on
    title('Fourier domain')
plot(k,WRadHat,'k','linewidth',3)
% plot(k,WEHat-(1-C)*WIHat,'r--','linewidth',3)
plot([0 5],[0 0],'k-','linewidth',2)
set(gca,'xlim',[0 3])
set(gca,'fontname','helvetica','fontsize',16,'linewidth',2);
set(gcf,'color','w');box on
% set(gcf,'units','centimeters','position',[20 10,10,8]);
% set(gca,'position',[0.15 0.1,.80,.85]);

    
end
