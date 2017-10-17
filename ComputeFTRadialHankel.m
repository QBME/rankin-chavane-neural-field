function sum=ComputeFTRadialHankel(G,r,k,flag)
% G - a function handle that takes radial variable only G(r)
% r - domain over which to integrate, strictly should be r\in[0,inf)
% the waves numbers k=|k|=sqrt(kx^2+ky^2) you want to evaluate at
if nargin==4
    fbflag=flag;
else
    fbflag='ft';
end
h=(r(end)-r(1))/(length(r)-1);
sum=.5*r(1)*G(r(1))*besselj(0,k.*r(1));
for i=2:length(r)-1
    tmp=r(i)*G(r(i))*besselj(0,k.*r(i));
   sum=sum+tmp;
end
sum=sum+.5*r(end)*G(r(end))*besselj(0,k.*r(end));
if ~isreal(sum)
    warning('imaginary part in computation of hankel transform')
end
switch fbflag
    case 'ft'
        sum=2*pi*real(h*sum);
    case 'ift'
        sum=1/(2*pi)*real(h*sum);
end