clear;
clc;

%------------------------------
rand(1,floor(sum(100*clock)));
%------------------------------
setenv('DYLD_LIBRARY_PATH', '/usr/local/bin/')

DC_mean1 = 0.025; DC_mean2 = 0.025; % Main, Branch
dcsigma1 =  0.25; dcsigma2 = 0.25;	
col_l1 = 2.5e3; col_l2 = 2.5e3;
DC_min = 0.01;

% Constant
year = 3600*24*365;
pb= qdyn('set'); % Branch.
p = qdyn('set'); % Main.

p.NSTOP = 0;        
p.TMAX = 100000*year; % stop at v = v_th = tmax

p.OX_SEQ=1;
p.OX_DYN=1;

p.DYN_TH_ON  = 0.1;
p.DYN_TH_OFF = 0.001;

%%%p.DC=0.3;
nxout=1;       %snapshot output grid interval
p.NXOUT_DYN=1; %dynamic snapshot output grid interval

co = 4e6;        %cohesion
co_limit = 2e3;  %first X m to apply cohesion 

p.MESHDIM=2;  % Main
pb.MESHDIM=2; % Branch

p.MU=32.5e9;
pb.MU=p.MU;
p.LAM=32.5e9;
pb.LAM=p.LAM;
p.MU_SS=0.6;
pb.MU_SS=p.MU_SS;

p.V_SS=0.005/year; % 5mm/year Ventura Fault! 
pb.V_SS=p.V_SS;
p.V2=100.; 
pb.V2=p.V2;
p.V1=p.V2;
pb.V1=p.V2;
 
scl=2;
% [-18.5,-8.5] : Dip = 45.0
% [-8.5, -7.5] : Dip = 7.5
% [-7.5,0]     : Dip = 45.0
dw0=125*scl;

% Dip limits for the main fault.
Lfault=[15.0e3,10.0e3,10.5e3];
Lbfault=10.5e3;
dipb = [45.0,7.5,45.0];
NW1=Lfault(1)/dw0;
NW2=NW1+Lfault(2)/dw0;
NW3=NW2+Lfault(3)/dw0;

p.W=sum(Lfault);
pb.W = Lbfault;

pb.L=64e3;
p.L=64e3;

p.NX=p.L/dw0;
pb.NX=pb.L/dw0;
p.NW = p.W/dw0;
pb.NW = pb.W/dw0;
p.N=p.NX*p.NW;
pb.N = pb.NX*pb.NW;

p.DW(1:p.NW)=dw0;    
pb.DW(1:pb.NW)=dw0;

p.DIP_W(1:p.NW)=dipb(1);      %deep to shallow 

% Diping for Fault sections with its corresponding depth.
p.DIP_W(1:NW1)    = dipb(1);
p.DIP_W(NW1+1:NW2)= dipb(2);
p.DIP_W(NW2+1:NW3)= dipb(3);

pb.DIP_W(1:pb.NW)=135.0;

db=Lfault(1)*sind(dipb(1))+Lfault(2)*sind(dipb(2))+Lfault(3)*sind(dipb(3));    %bottom of simulation zone (depth in m)

aa0=0.01;   %p.A
sigma0=75e6;%sigma max

p.Z_CORNER=-db;
p.Y_CORNER=0;
pb.Z_CORNER=-Lfault(3)*sind(dipb(3));
pb.Y_CORNER=Lfault(1)*cosd(dipb(1))+Lfault(2)*cosd(dipb(2));

p.VS=3000.0;
pb.VS=p.VS;

r_filter1 = ceil(col_l1/dw0);
r_filter2 = ceil(col_l2/dw0);

p.IC = p.N/2;

% Inserting 
p=qdyn('set',p);
pb= qdyn('set',pb); % Branch.


ZXW = reshape(p.Z,[p.NX,p.NW]);
ZW  = ZXW(1,:);

ba0=1.5;  %b/a at seismogenic zone
bam=0.6;  %b/a at shallow/deeper part
 
%% Ventura Fault has three transitions at depth. 
% Finding layers:
dd=10e3;    %depth of L sigma change
idd=max(find(ZW+dd<=0));
d1=2e3;     %upperbound of seismogenic zone (depth in m)
id1=max(find(ZW+d1<=0));
d2=5e3;     %limit of constant b/a
id2=max(find(ZW+d2<=0));
d3=12e3;    %lowerbound of seismogenic zone (depth in m)
id3=max(find(ZW+d3<=0));
d4=15e3;    %limit of constant b/a
id4=max(find(ZW+d4<=0));
ico_limit=max(find(ZW+co_limit<=0));

temp_A(1:p.NW)   = aa0;
temp_B(1:p.NW)   = aa0*bam;
temp_DC(1:p.NW)  = 0.3;
temp_SIGMA(1:p.NW) = sigma0;
temp_CO(1:p.NW) =  co;

temp_B(1:id4)      = temp_A(1:id4).*bam;
temp_B(id4+1:id3)  = temp_A(id4+1:id3).*linspace(bam,ba0,numel(id4+1:id3));   %increasing a/b below seismogenic zone
temp_B(id3+1:id2)  = temp_A(id3+1:id2).*ba0;   %a/b < 1 in seismogenic zone
temp_B(id2+1:id1)  = temp_A(id2+1:id1).*linspace(ba0,bam,numel(id2+1:id1));
temp_B(id1+1:p.NW) = temp_A(id1+1:p.NW).*bam;  %a/b >1 at shallow part 

temp_SIGMA(1:idd)      = sigma0;
temp_SIGMA(idd+1:p.NW) = linspace(sigma0,1e6,numel(idd+1:p.NW));
temp_CO(1:ico_limit)      = 0;


tmp_B(1:pb.NW) = temp_B(NW2+1:NW3);
tmp_A(1:pb.NW) = temp_A(NW2+1:NW3);
tmp_SIGMA(1:pb.NW)=temp_SIGMA(NW2+1:NW3);
tmp_DC(1:pb.NW)=0.3;
tmp_CO(1:pb.NW)=co;

    % Main and Branching fault (pb.NW = p.NW)
for i=1:p.NW
 %  p.X((i-1)*p.NX+1:i*p.NX) = linspace(0,p.L,p.NX);
    p.A((i-1)*p.NX+1:i*p.NX)     = temp_A(i);
    p.B((i-1)*p.NX+1:i*p.NX)     = temp_B(i);
    p.SIGMA((i-1)*p.NX+1:i*p.NX) = temp_SIGMA(i);
    p.DC((i-1)*p.NX+1:i*p.NX)    = temp_DC(i);
    p.CO((i-1)*p.NX+1:i*p.NX)    = temp_CO(i);
end;
for i=1:pb.NW
    % Branching fault 
   % pb.X((i-1)*pb.NX+1:i*pb.NX) = linspace(0,pb.L,pb.NX);
    pb.A((i-1)*pb.NX+1:i*pb.NX) = tmp_A(i);
    pb.B((i-1)*pb.NX+1:i*pb.NX) = tmp_B(i);
    pb.SIGMA((i-1)*pb.NX+1:i*pb.NX) = tmp_SIGMA(i);
    pb.DC((i-1)*pb.NX+1:i*pb.NX) = tmp_DC(i);
    pb.CO((i-1)*pb.NX+1:i*pb.NX) = tmp_CO(i);
end

for i=1:1:p.N
    p.IOT(i) = 0;
end

dl=-15000;
DC1=DC_corrl_lognor(p.X,p.NX,p.NW,DC_mean1,dcsigma1,r_filter1,DC_min);
DC2=DC_corrl_lognor(pb.X,pb.NX,pb.NW,DC_mean2,dcsigma2,r_filter2,DC_min);

p.DC  = DC1;
pb.DC = DC2;

twm=100000;         %warmup time in years
ts=1000;    %simulation time in years
p.ACC = 1e-10;
Vdyn=2*mean(p.A.*p.SIGMA./p.MU.*p.VS);
%p.DYN_TH_ON=Vdyn/10.;
%p.DYN_TH_OFF=Vdyn/10.;

%------------------------------
Lc = min(p.MU.*p.DC.*p.B./((p.B-p.A).^2.*(p.SIGMA)))
Lb = min(p.MU.*p.DC./p.SIGMA./p.B)
Lnuc = 1.3774*Lb
%------------------------------

filename = ['Hete_3D_ss_twm',num2str(twm),'L',num2str(p.L/1000.),...
    'nx',num2str(p.NX),'W',num2str(p.W/1000.),'nw',num2str(p.NW),...
    'dip',num2str(dipb(1)),'DCmean1',num2str(DC_mean1),'DcV1',num2str(dcsigma1),'.mat']

p.IC=ceil(p.N/2);
dw=p.W/p.NW;
Lb_over_dw = Lb/dw;
dx=p.L/p.NX;
Lb_over_dx = Lb/dx;

p  = qdyn('set',p);
pb = qdyn('set',pb);
pb.BRANCH='.true.';
pb.STRIKE=0; 
pb.XC = 0;
pb.YC = 0;

p.V_0  = 1.01*p.V_SS ;
pb.V_0 = p.V_0;

p.NTOUT=1000;
p.NXOUT=nxout;
p.NPROCS=6;
% p.DYN_FLAG=1;
% p.DYN_M=10.^19;
% p.DYN_SKIP = 1;
[pb,ot0,ox0]  = qdyn('write',pb);
[p,ot0,ox0]  = qdyn('run',p);

%semilogy(ot0.t/year,ot0.v)
%xlabel('Time (years)');
%ylabel('Vmax');
% 
%   p.NTOUT=10;
% 
%   p.V_0 = ox0.v(:,end);
%   p.TH_0= ox0.th(:,end);
%   %p.V_0 =  (ox1.v(:,end)+ox1.v(end:-1:1,end))/2;
%   %p.TH_0=  (ox1.th(:,end)+ox1.th(end:-1:1,end))/2;
%   [p,ot1,ox1]=qdyn('run',p);
%p0=p;
%V_0 = ox1.v(:,end);
%TH_0= ox1.th(:,end);
%save(filename)  
%save('warmup_jp_s.mat','p', 'V_0', 'TH_0');

