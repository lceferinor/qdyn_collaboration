%%% Specfem
 function [Mw,SumRupAreaSlip]=Moment_ruptureArea(snap)

 if nargin<1, snap=15000; end
 snaptime = snap;
 stime = snaptime*4e-3;
 snapspec = FSEM3D_snapshot(snaptime);
 slip=sqrt(snapspec.Dx.^2 + snapspec.Dz.^2);
 
 A = 125*125; %m2
 
 RupArea=zeros(1,length(slip))';
 
 RupArea(find(slip>0.01))=A;
 
 %figure,scatter(snapspec.X,snapspec.Z,[],RupArea)
 
 RupAreaTotal = sum(RupArea)
 
 RupAreaSlip=RupArea.*slip;
 
 SumRupAreaSlip=sum(RupAreaSlip);
 
 phi = 4444; %density (kg/m3);
 vs = 3000;%shear velocity (m/s)
 mud = phi*vs.^2;
% mudQDYN=40e+9; %PG warning in this case the density is 2.7*2 Vs=3km/s 
 Mo =  mud.*SumRupAreaSlip
% MoQDYN = mudQDYN.*SumRupArea
% 
 Mw = (2.0/3.0)*(log10(Mo*1e7))-10.73
% MwQDYN=(2.0/3.0)*(log10(MoQDYN*1e7))-10.73
 end

%%%% 

