
function AngMtxSel=FCNHSVAngMtx(AngMtx,SelMtx,MaxSel,brtn)
N=size(AngMtx,2);
hsvcmap=brighten(hsv,0);

AngMtxSel=zeros(N,N,3);
AngMtxHSVIdx=round(63*(AngMtx+pi)/2/pi)+1;
SelectivityMtxScl=SelMtx./MaxSel;
if max(SelectivityMtxScl(:))>1
%     warning('MaxSel too small')
    SelectivityMtxScl(SelectivityMtxScl>1)=1;
end
for i=1:N
    for j=1:N
        AngMtxSel(i,j,:)=hsvcmap(AngMtxHSVIdx(i,j),:).*SelectivityMtxScl(i,j).^brtn;
    end
end


end