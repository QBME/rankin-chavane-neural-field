function [SteepestIdx,PolySlope,linIdx,NRPNew]=RadFindSlope(RadialCoord,DecayOf,slopeWidth,MinRadius,NRFun)


minRadIdx=find(RadialCoord>MinRadius,1);

[~,centerIdx]=min(diff(DecayOf(minRadIdx:end)));
centerIdx=centerIdx+minRadIdx-1;


linIdx=centerIdx-slopeWidth:centerIdx+slopeWidth;
p=polyfit(RadialCoord(linIdx),DecayOf(linIdx),1);
PolySlope=p;

SteepestIdx=centerIdx;

    opts=optimset('display','off');
NRLSQFun=@(p,x) NRFun(x,p(1),p(2),p(3));
NRPar0=[5,20,0.0];
NRlb=[0.5,1,-.5];
NRub=[20,25,.5];

[NRPNew,~,~,lsqflag]=lsqcurvefit(NRLSQFun,NRPar0,RadialCoord,DecayOf,NRlb,NRub,opts);
