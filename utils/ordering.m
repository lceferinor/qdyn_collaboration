A=[];B=[];
nxfft=6; 
nwlocal=5;
ilocal=0;
for ix=1:nxfft
    for iw=1:nwlocal
        ilocal=ilocal+1;
        A(iw,ix)=ilocal
    end
end

iloc=0;
%iw=1;
for iwx=1:nxfft*nwlocal
    iloc=iloc+1;
    iw=floor((iwx-1)/(nxfft)) +1
    idx=(iwx-iw)*nwlocal+iw-((nxfft-1)*nwlocal)*(iw-1)
    B(iloc)=idx
    %iw=floor(iwx/(nxfft)) +1
end
    


