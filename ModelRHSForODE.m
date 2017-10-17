function [F] = RHSCFourMR2D_Test(t,u,p,wELocHat,wELatHat,wIHat...
    ,JH,JA,JV,JD,RI,X,Y) %#ok<*FNDEF>

  %% Rename parameters
  theta  = p(1);
  mu = p(2);
  IX0H  = p(3);
  IY0H  = p(4);
  % RI  = p(5); no longer used
  beta  = p(6);
  L  = p(7);
  I0 = p(8);
  IL = p(9);
%   OR = p(10); no longer used
  rho = p(11);
  IX0A  = p(13);
  IY0A  = p(14);
  IX0V  = p(15);
  IY0V  = p(16);
  IX0D  = p(17);
  IY0D  = p(18);
  Ilil  = p(19);
  N  = sqrt(size(u,1)/4);
  taup=p(20);
  beta_rec=p(24);
  % useful for centering inputs at right location
  x0=X(1,1);
  xend=X(1,end);
  xlen=2*L;
  
  %% Indices for individual sub-populations
  HIdx=1:N*N;
  AIdx=N*N+1:2*N*N;
  VIdx=2*N*N+1:3*N*N;
  DIdx=3*N*N+1:4*N*N;
  %% Get each sub-population u* from full state vector u (all populations)
  uH = reshape(u(HIdx),N,N);
  uA = reshape(u(AIdx),N,N);
  uV = reshape(u(VIdx),N,N);
  uD = reshape(u(DIdx),N,N);

  %% Inputs
  % linear ramp up of input strength starting at 2taup=20ms and reaching 1 at
  % 12taup=120ms; defined at end of section "Outline of the model" after
  % equation (7) for I_i(r)
  RampStart=2;RampFin=12;
  if t<RampStart*taup
      Ibig=0;
  elseif t>RampFin*taup
    Ibig=1;
  else      
    Ibig=(t-RampStart*taup)/(RampFin*taup);
  end
  % This is the relationship k_1=2*k_2; equivalently Ilil=0.5
  Ilil=Ilil*Ibig;
  
  % function h eq (10) in paper
  GaussRing=@(r,s,r0) exp(-(r-r0).^2/2/s^2);
  InputH=zeros(size(X));
  InputA=InputH;InputV=InputH;InputD=InputH;
  % For any orientations with zero radius (equivalently no input) the input
  % Input* remains 0
  if RI(1)~=0
      R=sqrt((mod(X-IX0H-x0,xlen)+x0).^2+(mod(Y-IY0H-x0,xlen)+x0).^2);
      % Equation (7)
      InputH = I0 * GaussRing(R,IL,RI(1));
      InputH(R<RI(1))=I0;
  end
  if RI(2)~=0
      R=sqrt((mod(X-IX0A-x0,xlen)+x0).^2+(mod(Y-IY0A-x0,xlen)+x0).^2);
      InputA = I0 * GaussRing(R,IL,RI(2));
      InputA(R<RI(2))=I0;
  end
  if  RI(3)~=0
      R=sqrt((mod(X-IX0V-x0,xlen)+x0).^2+(mod(Y-IY0V-x0,xlen)+x0).^2);
      InputV = I0 * GaussRing(R,IL,RI(3));
      InputV(R<RI(3))=I0;
  end
  if  RI(4)~=0
      R=sqrt((mod(X-IX0D-x0,xlen)+x0).^2+(mod(Y-IY0D-x0,xlen)+x0).^2);
      InputD = I0 * GaussRing(R,IL,RI(4));
      InputD(R<RI(4))=I0;
  end
  
  %% Modulate inputs but orientation map J* with strength beta 
  %  (beta_inp in paper)
  InputH=InputH.*(1+beta*JH);
  InputA=InputA.*(1+beta*JA);
  InputV=InputV.*(1+beta*JV);
  InputD=InputD.*(1+beta*JD);
  

  %% Firing rate function (and its derivative)
  % Apply firing rate function S, eq (5)
  [fH] = ComputeFiringRate(mu*uH,theta);
  [fA] = ComputeFiringRate(mu*uA,theta);
  [fV] = ComputeFiringRate(mu*uV,theta);
  [fD] = ComputeFiringRate(mu*uD,theta);
  % Preparating for applying the convolutions in Terms (3)--(4) 
  fHatHI = fft2(fH);
  fHatAI = fft2(fA);
  fHatVI = fft2(fV);
  fHatDI = fft2(fD);  
  fHatHELat = fft2(fH.*(1+beta_rec*JH));
  fHatAELat = fft2(fA.*(1+beta_rec*JA));
  fHatVELat = fft2(fV.*(1+beta_rec*JV));
  fHatDELat = fft2(fD.*(1+beta_rec*JD));
  fHatHELoc = fft2(fH);
  fHatAELoc = fft2(fA);
  fHatVELoc = fft2(fV);
  fHatDELoc = fft2(fD);

  %% RHS
  
  FH = 1/taup*(- uH + (2*L/N)^2*ifftshift(real(ifft2(fHatHELoc .* wELocHat))) ...
      - (2*L/N)^2*ifftshift(real(ifft2(fHatHI .* wIHat))) ...
      + (2*L/N)^2*ifftshift(real(ifft2(fHatHELat .* wELatHat)))...
      + Ibig*InputH + Ilil*(InputA+InputV+InputD) - rho*(uA+uV+uD));
  FA = 1/taup*(- uA + (2*L/N)^2*ifftshift(real(ifft2(fHatAELoc .* wELocHat))) ...
      - (2*L/N)^2*ifftshift(real(ifft2(fHatAI .* wIHat))) ...
      + (2*L/N)^2*ifftshift(real(ifft2(fHatAELat .* wELatHat)))...
      + Ibig*InputA + Ilil*(InputH+InputV+InputD) - rho*(uH+uV+uD));
  FV = 1/taup*(- uV + (2*L/N)^2*ifftshift(real(ifft2(fHatVELoc .* wELocHat))) ...
      - (2*L/N)^2*ifftshift(real(ifft2(fHatVI .* wIHat))) ...
      + (2*L/N)^2*ifftshift(real(ifft2(fHatVELat .* wELatHat)))...
      + Ibig*InputV + Ilil*(InputH+InputA+InputD) - rho*(uH+uA+uD));
  FD = 1/taup*(- uD + (2*L/N)^2*ifftshift(real(ifft2(fHatDELoc .* wELocHat))) ...
      - (2*L/N)^2*ifftshift(real(ifft2(fHatDI .* wIHat))) ...
      + (2*L/N)^2*ifftshift(real(ifft2(fHatDELat .* wELatHat)))...
      + Ibig*InputD + Ilil*(InputH+InputA+InputV) - rho*(uH+uA+uV));
  % returh the RHS in one big vector, same size as u that came in
  F =[ FH(:); FA(:); FV(:); FD(:)];

end
