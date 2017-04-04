%%% Specfem
 function Moment_ruptureArea_join()

 pth1='/Users/percy.galvez/Dropbox/qdyn-cohesion/examples/planar_faults_cohesion/QSB_Mw6/SPECFEM_out/run2';
 pth2='/Users/percy.galvez/Dropbox/qdyn-cohesion/examples/planar_faults_cohesion/QSB_Mw6/SPECFEM_out/run3';
 
 snaptime = 15000;
 stime = snaptime*4e-3;
 snapspec1 = FSEM3D_snapshot(snaptime,1,pth1);
 snapspec2 = FSEM3D_snapshot(snaptime,1,pth2);
 snapspec.Dx = snapspec1.Dx + snapspec2.Dx;
 snapspec.Dz = snapspec1.Dz + snapspec2.Dz;
 snapspec.X = snapspec1.X;
 snapspec.Z = snapspec1.Z;

 
 slip=sqrt(snapspec.Dx.^2 + snapspec.Dz.^2);
 
 A = 125*125; %m2
 
 RupArea=zeros(1,length(slip))';
 
 RupArea(find(slip>0.01))=A;
 
 figure,scatter(snapspec.X,snapspec.Z,[],slip);caxis([0 max(slip)])
 
 RupAreaTotal = sum(RupArea)
 
 RupAreaSlip=RupArea.*slip;
 
 SumRupAreaSlip=sum(RupAreaSlip);
 
 phi = 2700; %density (kg/m3);
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

