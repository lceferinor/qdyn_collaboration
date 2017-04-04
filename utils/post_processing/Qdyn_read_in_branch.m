% Reading fault branch. 
% [X(:),Y(:),Z(:),SIGMA(:),V_0(:),TH_0(:),A(:),B(:),DC(:),V1(:),V2(:),MU_SS(:),V_SS(:),CO(:)]');
function fb=Qdyn_read_in_branch
 fb=[];
 fid = fopen('qdyn_branch.in');
 rline=textscan(fid,'%d\n',1);
 N=rline{1};
 rline=textscan(fid,'%15.6f\n',1);
 DIP0=rline{1};
 rline=textscan(fid,'%20.6f %20.6f\n',1);
 lam=rline{1};
 mu=rline{2};
 rlines=textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f\n',N);
 
 fb.X = rlines{1};
 fb.Y = rlines{2};
 fb.Z = rlines{3};

 fb.SIGMA = rlines{4};
 fb.V_0   = rlines{5};
 fb.TH_0  = rlines{6};

 fb.A  = rlines{7};
 fb.B  = rlines{8};
 fb.DC = rlines{9};
 fb.V1 = rlines{10};
 fb.V2 = rlines{11};
 fb.MU_SS = rlines{12};
 fb.V_SS = rlines{13};
 fb.C0 = rlines{14};
 
 fclose(fid);
end