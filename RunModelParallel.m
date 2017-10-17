function RunModel(varargin)

switch nargin
    case 1
        FigString=varargin{1};
    case 0
        FigString='7E';
    otherwise
        error('Accepts 0 or 1 arguments, e.g. RunModel or RunModel(''6B'')')
end

clc;close all;

%% Values for parameters that change in the paper
% n sets which if the 5 map locations used in the paper; n=0 gives a
% random location

p0=zeros(1,24);
switch FigString 
     case '5B' % also panels C and D
        n=4;
        beta_rec=0.0; p0(24)=beta_rec;
        C=-.4;
        RWEx=0.1;
     case '5E' % also panels C and D
        n=4;
        beta_rec=0.0; p0(24)=beta_rec;
        C=-.4;
        RWEx=0.25;
    case '6B'
        n=2;
        beta_rec=0.0; p0(24)=beta_rec;
        C=-.4;
        RWEx=0.25;
    case '6C'
        n=2;
        beta_rec=0.3; p0(24)=beta_rec;
        C=-.4;
        RWEx=0.25;
    case '6D'
        n=2;
        beta_rec=0.5; p0(24)=beta_rec;
        C=-.4;
        RWEx=0.25;
    case '6E'
        n=2;
        beta_rec=0.9; p0(24)=beta_rec;
        C=-.4;
        RWEx=0.25;
    case '7D'
        n=5;
        beta_rec=0.9; p0(24)=beta_rec;
        C=-.4;
        RWEx=0.25;
    case '7E'
        n=5;
        beta_rec=0.6; p0(24)=beta_rec;
        C=-.4;
        RWEx=0.225;
    case '7F'
        n=5;
        beta_rec=0.4; p0(24)=beta_rec;
        C=-.4;
        RWEx=0.15;    
    case '9Btop' % zoomout: select figure and set axis([-30,30,-30,30])
        n=5;
        beta_rec=0.9; p0(24)=beta_rec;
        C=-.2;
        RWEx=0.25;     
    case '9Bbottom' % zoomout: select figure and set axis([-30,30,-30,30])
        n=5;
        beta_rec=0.4; p0(24)=beta_rec;
        C=-.2;
        RWEx=0.225;        
    otherwise % specify your own values 
        % random map location for n=0, set to a value between 1 to 5 for
        % locations used in study
        n=0;
        beta_rec=0.6; p0(24)=beta_rec;
        C=-.4;
        RWEx=0.225;
end


%% Define discretisation 
L = 30; N =  128; h = 2*L/N; x  = (-L+(0:N-1)*h)'; 
[X,Y] = meshgrid(x,x); 
R = sqrt(X.^2+Y.^2);
LinX=X(1,:);LinY=LinX;
% Lambda in the paper
lambda=2*pi;

%% Set model paramters
% Firing rate function; see ComputeFiringRate.m
th=5.6;                 p0(1)=th;% \theta ; sigmoid threshold
mu=2.3;                 p0(2)=mu;% \mu ; sigmoid slope
% Input parameters
IX0=0;                  p0(3)=IX0;% stim location for H (0 deg)
IY0=0;                  p0(4)=IY0;% stim location for H (0 deg)
RI=0.70*lambda;         p0(5)=RI; % no longer used
RIbase=0.025*lambda; 
beta=0.25;              p0(6)=beta;% beta_inp in paper
                        p0(7)=L; % passes domain size to RHS/ODE
I0=2.8;                 p0(8)=I0;% k_1 in the paper
IL=.3*lambda;           p0(9)=IL;% in Eq (4), Gaussian decay at stimulus border 

                        p0(10)=0;% no longer used
p0(13:18)=0; % stimulus locations for A,V,D (45,90,135 deg), always at origin
% these are used for plotting and would be non-zero if stimuli were not at
% origin
IX0H=0;IY0H=0;IX0A=0;IY0A=0;IX0V=0;IY0V=0;IX0D=0;IY0D=0;
Ilil=0.5;               p0(19)=Ilil; %k_2=Ilil*k_1=1.4


rho = 0.1;            p0(11)=rho; % \rho; cross inhibition strength
taup=10;              p0(20)=taup; % timescale tau=10ms
                      p0(21)=0;% no longer used
                      p0(22)=0;% no longer used
                      p0(23)=0;% no longer used

% connectivity parameters
EDCEx=0.625; %\zeta in paper; multiplied by lambda later
RWIn=0.55; % RW_in in paper

% set up ODE integration interval and which points to save for plotting
tfinal=550;tint=0:tfinal;
tsaveIdx=1:10:length(tint);

% post-processing parameters
SelThreshFrac=0.5; % eta_sel 
ActThreshFrac=0.2; % eta_act
ADThresh=30/180*pi; % +-30 degrees described for yellow contours in Fig 5 caption
brtn=1.3; % brighten colormap
gscale=0.075*lambda; % sigma_oi in the paper; gaussian blurring of VSD signal
vsdwi=0.15/0.85; % p_I in paper; in eq (20)
ActivationRadius=1.1*lambda;
taump=240; % tau_lat in paper, used in eq (20)
S = @(x) ComputeFiringRate(x,th); % 1./(1+exp(-x+th))-1/(1+exp(th));


%% Load orientation preferences maps 
Ji=load('OrientationMapsJi.mat');
JHdef=Ji.JHdef;
JAdef=Ji.JAdef;
JVdef=Ji.JVdef;
JDdef=Ji.JDdef;

%% Shift the location in the orientation preference map 
% 5 locations stored in MapShiftVals.mat were used in the paper
if n==0
    MpShMtx=zeros(2,1);
    rng('shuffle')
    MpShMtx(1,:)=randperm(N,1);
    MpShMtx(2,:)=randperm(N,1);
    MapShiftX=MpShMtx(1,1);
    MapShiftY=MpShMtx(2,1);
else
    load(['MapShiftVals.mat']);
    MpShMtx=[MapShiftXVals';MapShiftYVals'];
    MapShiftX=MpShMtx(1,n);
    MapShiftY=MpShMtx(2,n);
end
Rval=(RI+RIbase);
JH=circshift(JHdef,[MapShiftX,MapShiftY]);
JA=circshift(JAdef,[MapShiftX,MapShiftY]);
JV=circshift(JVdef,[MapShiftX,MapShiftY]);
JD=circshift(JDdef,[MapShiftX,MapShiftY]);
COMX=JH-JV;
COMY=JA-JD;
MapMtx=atan2(COMY,COMX);

%% Connectivity
plotConnFlag=0;
[w,wLoc,wLat,WEFun,WEFunLat,WEFunLoc,WIFun,P]=BuildMultiRingFcns(...
    L,N,R,X,Y,lambda,C,EDCEx*lambda,RWEx*lambda,RWIn*lambda,plotConnFlag);
% The three components of the connectivity in real domain
WELoc=P*WEFunLoc(R);
WELat=P*WEFunLat(R);
WI=P*(1-C)*WIFun(R);
% In Fourier domain (this are passed to ODE solver)
wELocHat=fft2(WELoc);
wELatHat=fft2(WELat);
wIHat=fft2(WI);

%% Pre-allocate a struct for plotting later
sumSu=zeros(4,N,N);
svStruct=struct([]);
for i=1:length(tsaveIdx)
    svStruct(i).ufAllH=sumSu;
    svStruct(i).ufAllA=sumSu;
    svStruct(i).ufAllV=sumSu;
    svStruct(i).ufAllD=sumSu;
end

%% Main for loop where simulations are run via ODE solver
yBig={};% initialize data structure for be filled by ODE solver
% Run the ode solver for each orientation, parallelised as these are
% independent (could be run in serial 
parfor OR=1:4 % ***for serial implementation replace parfor with for***
    % initialize from a small random initial condition in each layer
    u0amp=0.1;
uH0=u0amp*randn(size(X));
uA0=u0amp*randn(size(X));
uV0=u0amp*randn(size(X));
uD0=u0amp*randn(size(X));
u0=[uH0(:);uA0(:);uV0(:);uD0(:)];

% Sets the stimulus radius to zero for all except the present OR. Within
% ModelRHSForODE any OR with radius 0 gets zero input
RIFour=zeros(4,1);
RIFour(OR)=Rval;

opts=odeset;
odehandle = @(t,u) ModelRHSForODE(t,u,p0,wELocHat,wELatHat,wIHat,JH,JA,JV,JD,RIFour,X,Y);

% run the ODE solver
tic
[t,yBig{OR}]=ode113(odehandle,tint,u0,opts);
toc
end

% Indices for sub-populations
HIdx=1:N*N;
AIdx=N*N+1:2*N*N;
VIdx=2*N*N+1:3*N*N;
DIdx=3*N*N+1:4*N*N;
% Fill up struct for plotting later
for ORIdx=1:4
for i=1:length(tsaveIdx)
uf=yBig{ORIdx}(tsaveIdx(i),:);
ufH=reshape(uf(HIdx),N,N);
ufA=reshape(uf(AIdx),N,N);
ufV=reshape(uf(VIdx),N,N);
ufD=reshape(uf(DIdx),N,N);
svStruct(i).ufAllH(ORIdx,:,:)=ufH;
svStruct(i).ufAllA(ORIdx,:,:)=ufA;
svStruct(i).ufAllV(ORIdx,:,:)=ufV;
svStruct(i).ufAllD(ORIdx,:,:)=ufD;
end
end

%% Post-processing; see "Conversion of model output to VSD-like signal"
% call make_colors.m to define some nice colours
make_colors
xymax=12; % axis limits

FFRadius=ActivationRadius;
FFArea=numel(find(R<FFRadius));

% Gaussian smoothing with sig_oi
sig=gscale;w=@(r) 1/2/pi/sig.^2*exp(-r.^2./2/sig^2); 
SmoothKernel=w(R);SKHat=fft2(SmoothKernel);

% process final frame
Idx=length(tsaveIdx);
tval=tint(tsaveIdx(Idx));
svTmp=svStruct(Idx);

% Used as part of equation (20) 
VSDKernelLoc=WEFunLoc(R)-vsdwi*WIFun(R);% independent of time 
fmp=@(t) 1-exp(-t/taump);  % appears in eq (20)
VSDKernelLat=@(t)fmp(t)*WEFunLat(R);
VKHatFcnLoc=fft2(VSDKernelLoc);
VKHatFcnLat=@(t)fft2(VSDKernelLat(t));

% Lateral excitatory connections modulated by orientation map
JStructLat={1+(beta_rec)*JH,1+(beta_rec)*JA,1+(beta_rec)*JV,1+(beta_rec)*JD};
% Local excitatory connections are not 
JStructLoc={ones(size(JH)),ones(size(JH)),ones(size(JH)),ones(size(JH))};

% Does the computation of eq (20)
[yvsdFour,yvsdMax]=utoVSD_FCN(L,N,S,mu,svTmp,VKHatFcnLoc,VKHatFcnLat(tval),SKHat,JStructLat,JStructLoc);

% Before computing the preference Pref and selectivity Sel, the VSD signals
% are normalized by a scale factor that accounts for differences in the
% maximum value over ( x, y ) across the four simulations with different
% orientations.
maxAct=max(yvsdMax);
yvsdFour=yvsdFour/maxAct;
yvsdFourSel=zeros(size(yvsdFour));
meanysdFourSel=max(reshape(yvsdFour(:,:,:),4,N*N),[],2);
for i=1:4
    scalefac=(1+(mean(meanysdFourSel)-meanysdFourSel(i)));
    yvsdFourSel(i,:,:)=yvsdFour(i,:,:)*scalefac;
end

% Eq (24)-(25) applied in here to get Sel(x,y), Pref(x,y)
[AngMtx,SelMtx]=AngSelFcnFCN(yvsdFourSel);
% The intensity of the colormaps shown is scaled by the final point
MaxSel=max(SelMtx(:))*1.1;

% Gets a colormap scaled by secltivity Sel(x,y) at each point
AngMtxSel=FCNHSVAngMtx(AngMtx,SelMtx,MaxSel,brtn);
meanAct=reshape(sum(yvsdFour,1)/4,N,N);

% used to compute mean in eq (26)--(27)
MeanPlateauIdx=find(R<Rval);

% \bar{Act}_{FF} in eq (28)
MeanPlateauAct=mean(meanAct(MeanPlateauIdx));
% \bar{Sel}_{FF} in eq (28)
MeanPlateauSel=mean(SelMtx(MeanPlateauIdx));

% Eq (28)
SelThresh=MeanPlateauSel*SelThreshFrac;
ActThresh=MeanPlateauAct*ActThreshFrac;

%% Fits to NakaRushton function for radial profiles
% Transform into polar coordinates
RadialCoord=(.4:.025:3)*lambda; % Radial discretisation
AngularCoord=linspace(0,2*pi,100); % Angular discretisation
[rr,tt]=meshgrid(RadialCoord,AngularCoord); % Create grid 
[X2,Y2]=pol2cart(tt,rr); 
RadialMA=griddata(X,Y,meanAct,X2,Y2); %#ok<GRIDD>
RadialSM=griddata(X,Y,SelMtx,X2,Y2); %#ok<GRIDD>

% Take averages in angular coordinate
DecayOfAct=mean(RadialMA,1)/MeanPlateauAct;
DecayOfSel=mean(RadialSM,1)/MeanPlateauSel;
% Eq (29)
NRFun=@(r,n,rhh,M) 1-(1-M)*r.^n./(r.^n+rhh.^n);
% Parameters to make sure fitting to Naka-Rushton behaves
slopeWidth=3;
MinRadius=4.5;
[SteepestIdxAct,PolySlopeAct,linIdxAct,NRPAct]=RadFindSlope(RadialCoord,DecayOfAct,slopeWidth,MinRadius,NRFun);
[SteepestIdxSel,PolySlopeSel,linIdxSel,NRPSel]=RadFindSlope(RadialCoord,DecayOfSel,slopeWidth,MinRadius,NRFun);

%% -------------- Everything after this point is plotting -----------------

%% Plot radial decay if active and selective region e.g. Fig 5D
%  show values of n_act and n_sel on plot
figure(1);clf;hold on
SelPlotFac=0.5;
drawnow;pause(0.1);
plot([FFRadius,FFRadius],[0,1.5],'r-','linewidth',1)
plot(RadialCoord,DecayOfAct,'color',dark_green)
plot(RadialCoord,SelPlotFac*DecayOfSel,'color',blue)
set(gca,'ylim',[0,1.2],'ytick',[0,1])

text(9,1.0,['n_{Act}=',num2str(NRPAct(1),'%.2f')],'color',dark_green)
text(9,0.8,['n_{Sel}=',num2str(NRPSel(1),'%.2f')],'color',blue)
set(findall(gcf,'-property','Fontname'),'Fontname','helvetica')
set(findall(gcf,'-property','FontSize'),'FontSize',8)

xlabel('r (a.u.)')
ylabel('activity')
plot(RadialCoord,NRFun(RadialCoord,NRPAct(1),NRPAct(2),NRPAct(3)),'-.','color',dark_green)
plot(RadialCoord,SelPlotFac*NRFun(RadialCoord,NRPSel(1),NRPSel(2),NRPSel(3)),'-.','color',blue)

% set(gcf,'units','centimeters','position',[6,10,6.4,5.6]);
drawnow;pause(0.1);

    %%
cc=linspace(0,2*pi,100);

%% Go through each of 4 sub-populations and plot ACT
figure (2);clf;hold on
for i=1:4 %
subplot(2,2,i);hold on
vsdH=reshape(yvsdFourSel(i,:,:)/max(yvsdFourSel(:)),N,N);
colorbar
imagesc([-L L],[-L L],vsdH,[0,1]);axis equal;
plot(ActivationRadius*(cos(cc)),ActivationRadius*(sin(cc)),'color',red,'linewidth',1.5)
contour(LinX,LinY,vsdH,[0.07 0.45 0.8],'color','w','linewidth',1) ;
axis off
% axis([-axlim axlim -axlim axlim]);
end

colormap bone

SelAreaTH=zeros(size(svStruct,2),1);
SelTH=zeros(size(svStruct,2),1);
ActAreaTH=zeros(size(svStruct,2),1);
SelOutTH=zeros(size(svStruct,2),1);
SelPrefTH=zeros(size(svStruct,2),1);
SelOutPrefTH=zeros(size(svStruct,2),1);

%% This is "frame-by-frame" plotting of spread of activity
% some of the post-processing steps need recomputing for each frame
% for tsIdx=length(tsaveIdx) % uncomment for last frame only
for tsIdx=1:length(tsaveIdx)

    % Post-processing steps repeated
    Idx=tsIdx;
    tval=tint(tsaveIdx(Idx));
    svTmp=svStruct(Idx);
    [yvsdFour,yvsdMax]=utoVSD_FCN(L,N,S,mu,svTmp,VKHatFcnLoc,VKHatFcnLat(tval),SKHat,JStructLat,JStructLoc);
    yvsdFour=yvsdFour/maxAct;
    yvsdFourSel=zeros(size(yvsdFour));
    meanysdFourSel=max(reshape(yvsdFour(:,:,:),4,N*N),[],2);
    for i=1:4
        scalefac=(1+(mean(meanysdFourSel)-meanysdFourSel(i)));
        yvsdFourSel(i,:,:)=yvsdFour(i,:,:)*scalefac;
    end
    [AngMtx,SelMtx]=AngSelFcnFCN(yvsdFourSel);
    AngMtxSel=FCNHSVAngMtx(AngMtx,SelMtx,MaxSel,brtn);
    meanAct=reshape(sum(yvsdFour,1)/4,N,N);
    AngDiffMtx=abs(AngMtx-MapMtx)/2;
    SelPrefMtx=AngDiffMtx<ADThresh&SelMtx>SelThresh;
    
    % Plots showing spread of activity Figs 5B,E; Fig 6D-E, Fig 7D etc etc
    figure(3);clf;hold on
    subplot(2,1,1);hold on
    set(gca,'position',[0.06,0.50,0.88,0.44]);
    imagesc(LinX,LinY,meanAct,[0,1.3*MeanPlateauAct])
    colormap bone
    contour(X,Y,meanAct,[ActThresh,ActThresh],'color',grey,'linewidth',1.5);
    plot(FFRadius*(cos(cc)),FFRadius*(sin(cc)),'color',red,'linewidth',1.5);
    axis off
    axis square
    set(gca,'xlim',[-xymax,xymax],'ylim',[-xymax,xymax])
    text(-5,13,['t=',num2str(tval),'ms'],'fontsize',8,'fontname','helvetica')
    subplot(2,1,2);hold on
    set(gca,'position',[0.06,0.02,0.88,0.44]);
    axis off
    axis square
    imagesc(LinX,LinY,AngMtxSel,[0,1.3*MeanPlateauSel])
    contour(X,Y,SelMtx,[SelThresh,SelThresh],'w','linewidth',1.5);
    if tsIdx==length(tsaveIdx)
        contour(X,Y,SelPrefMtx,[SelThresh,SelThresh],'y','linewidth',1);
    end
    contour(X,Y,meanAct,[ActThresh,ActThresh],'color',grey,'linewidth',1.5);
    plot(FFRadius*(cos(cc)),FFRadius*(sin(cc)),'color',red,'linewidth',1.5);
    set(gca,'xlim',[-xymax,xymax],'ylim',[-xymax,xymax])
    drawnow
    
    % Fill up vectors for plotting time histories
    AngDiffMtx=abs(AngMtx-MapMtx)/2;
    SelPrefMtx=AngDiffMtx<ADThresh&SelMtx>SelThresh;
    meanAct=reshape(sum(yvsdFour,1)/4,N,N);
    ActAreaTH(tsIdx)=numel(find(meanAct>ActThresh));
    SelAreaTH(tsIdx)=numel(find(SelMtx>SelThresh));
    SelOutTH(tsIdx)=numel(find(SelMtx>SelThresh&R>FFRadius));
    SelPrefTH(tsIdx)=numel(find(SelPrefMtx));
    SelOutPrefTH(tsIdx)=numel(find(SelPrefMtx&R>FFRadius));
    SelTH(tsIdx)=max(SelMtx(:));

end

%% Plot time course of area of spread of active and selective regions
%  e.g. Fig 5C and F
figure(4);clf;hold on
drawnow;pause(0.1)
plot([0,700],[1,1],'color',red)
l2=plot(tint(tsaveIdx),ActAreaTH/FFArea,'color',dark_green,'linewidth',1.5)
l1=plot(tint(tsaveIdx),SelAreaTH/FFArea,'color',blue,'linewidth',1.5)
plot(tint(tsaveIdx),SelOutTH/FFArea,'color',dark_grey,'linewidth',1.5)
set(gca,'xlim',[0,tint(tsaveIdx(end))])
set(findall(gcf,'-property','Fontname'),'Fontname','helvetica')
set(findall(gcf,'-property','FontSize'),'FontSize',8)
xlabel('t (ms)')
ylabel('normalised area')
legend([l1,l2],{'selective activation','general activation'});
drawnow;pause(0.1);


%% Plot the local orientation preference map
%  e.g. Fig 5A
drawnow;pause(0.1)
figure(5);clf;hold on
axis equal
drawnow;pause(0.1)
    set(gca,'position',[0.06,0.02,0.88,0.88]);
set(gca,'xlim',[-xymax,xymax],'ylim',[-xymax,xymax]);
axis off
imagesc(LinX,LinY,MapMtx)
plot(0,0,'p','markerfacecolor','w','color','k')
plot(FFRadius*(cos(cc)),FFRadius*(sin(cc)),'color',red,'linewidth',1.5);
colormap hsv
drawnow;pause(0.1);

return% comment this to also plot last figure

%% Plot values of each sub-population u in 4 panels in each figure
%  each of 4 figures is for stimulus with different orientation
%  one of these 4 figures was shown in Fig 4A for orientation 0
for ORIdx=1:4
    axlim=15;
    colormap bone
    figure(10+ORIdx);clf;hold on
%     set(gcf,'units','centimeters','position',[5 30 30 15]);
    subplot(2,2,1);hold on
    title('u_0','fontweight','normal')
    imagesc([-L L],[-L L],ufH,[-3 5]);axis equal;
    plot(ActivationRadius*(cos(cc))+IX0H,ActivationRadius*(sin(cc))+IY0H,'color',red,'linewidth',3)
    axis off
    axis([-axlim axlim -axlim axlim]);
    drawnow
    subplot(2,2,2);hold on
    title('u_{45}','fontweight','normal')
    imagesc([-L L],[-L L],ufA,[-3 5]);axis equal;
    plot(ActivationRadius*(cos(cc))+IX0A,ActivationRadius*(sin(cc))+IY0A,'color',red,'linewidth',3)
    axis off
    axis([-axlim axlim -axlim axlim]);
    drawnow
    subplot(2,2,3);hold on
    title('u_{90}','fontweight','normal')
    imagesc([-L L],[-L L],ufV,[-3 5]);axis equal;
    plot(ActivationRadius*(cos(cc))+IX0V,ActivationRadius*(sin(cc))+IY0V,'color',red,'linewidth',3)
     axis([-axlim axlim -axlim axlim]);
    drawnow
    subplot(2,2,4);hold on
    title('u_{135}','fontweight','normal')
    imagesc([-L L],[-L L],ufD,[-3 5]);axis equal;
    plot(ActivationRadius*(cos(cc))+IX0D,ActivationRadius*(sin(cc))+IY0D,'color',red,'linewidth',3)
    colorbar
    colormap bone
    axis off
    axis([-axlim axlim -axlim axlim]);
    drawnow
end


