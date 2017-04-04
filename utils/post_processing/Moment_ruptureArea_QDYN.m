%%% Specfem
 function Moment_ruptureArea_QDYN(snap)
 
 km=1000;
%snap=1015;
qdynfile=['fort.' num2str(snap)];
qdynfile0=['fort.20001'];

qdynsnap  = Qdyn_read_ox_seq(qdynfile);
qdynsnap0 = Qdyn_read_ox_seq(qdynfile0);

slip=qdynsnap.d-qdynsnap0.d;
fig=figure;scatter(qdynsnap.X/km,qdynsnap.Z/km,[],slip)
title(qdynsnap.t-qdynsnap0.t);
 
 
 A = 125*125; %m2
 
 RupArea=zeros(1,length(slip))';
 
 RupArea(find(slip>0.01))=A;
 
 RupAreaTotal = sum(RupArea)/1e6
 
 RupAreaSlip=RupArea.*slip;
 
 SumRupAreaSlip=sum(RupAreaSlip);
 
mudQDYN=40e+9; %PG warning in this case the density is 2.7*2 Vs=3km/s 

Mo =  mudQDYN.*SumRupAreaSlip
% MoQDYN = mudQDYN.*SumRupArea
% 
 Mw = (2.0/3.0)*(log10(Mo*1e7))-10.7
% MwQDYN=(2.0/3.0)*(log10(MoQDYN*1e7))-10.73
 end
 
%%%% 

