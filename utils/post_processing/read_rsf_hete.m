% Reading fault branch. 
% [X(:),Y(:),Z(:),SIGMA(:),V_0(:),TH_0(:),A(:),B(:),DC(:),V1(:),V2(:),MU_SS(:),V_SS(:),CO(:)]');
function fb=read_rsf_hete(filename)
 fb=[];
 fid = fopen(filename);
 rline=textscan(fid,'%d\n',1);
 N=rline{1};

 rlines=textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',N);
 
 fb.X = rlines{1};
 fb.Y = rlines{2};
 fb.Z = rlines{3};

 fb.SIGMA = rlines{4};
 fb.tau=rlines{5};
 fb.tau_dip=rlines{6};
 fb.V_0   = rlines{7};
 fb.MUSS  = rlines{8};

 fb.A  = rlines{9};
 fb.B  = rlines{10};
 fb.DC = rlines{11};
 fb.V_ini = rlines{12};
 fb.TH_ini = rlines{13};
 fb.C0 = rlines{14};
 
 fclose(fid);
end