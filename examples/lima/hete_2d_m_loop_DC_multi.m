%------------------------------
% 1.5D fault with heterogenous Dc distributions (random log normal or random log uniform) 
%------------------------------

clear;
clc;

num_run = 100;	% number of simulations
ts=100; % duration in years of each simulation
i_dc_dist=2;    %=1 for lognormal, =2 for uniform log Dc

%------------------------------
rand(1,floor(sum(100*clock)));
%------------------------------

year = 3600*24*365;

% loop over dc range, from Dc_min to DC_max=0.03+0.01*ii_dc
% for uniform log (i_dc_dist=2) simulation only 
for ii_dc=1:1:1		

% p = qdyn('set');
p.RNS_LAW=0;
p.MESHDIM=1;      %FFT enabled
p.THETA_LAW=1;
p.RNS_LAW=0;
p.MU=30e9;
p.LAM=30e9;
p.MU_SS=0.6;
p.V_SS=0.085/year;
%p.V2=100.;              %no cut off velocity
%p.V1=p.V2;
p.OX_SEQ=0;
%p.DC=0.4;

p.W=20e3;
p.L=1000e3;
p.N=256*8;
p.NX=p.N;
p.FINITE=0;

p.A=0.5e-2;
p.B=1e-2;
p.SIGMA=75e6;
  
twm=1000;         %warmup time in years
p.ACC = 1e-14;

%------------
%lognormal distribution
dis_mu=0;
dis_sigma=1.;
th_DC=0.02;     %threshhold of lower DC
filter_range=1;    %cell size of filters
co_range=1;    % co_range of cells to be set at the same value
% for uniform Log Dc
DC_min=0.03;     %threshhold of lower DC
DC_max=0.03+0.01*ii_dc;     %max DC
%------------

p = qdyn('set',p);

%Lc=Lb*(p.B/(p.B-p.A));
%disp(['  Lc=',num2str(Lc),'  L/Lc=',num2str(p.L/Lc),'  W/Lc=',num2str(p.W/Lc)]);
%Linf=2/pi*(p.B/(p.B-p.A))^2*Lb;
%disp(['  Linf=',num2str(Linf),'  L/Linf=',num2str(p.L/Linf),'  W/Linf=',num2str(p.W/Linf)]);

p.TMAX=twm*year;
p.V_0 = 1.01*p.V_SS ;

tmp_DC=p.DC;
dd=zeros(1,p.NX+filter_range*2);
p.DC=zeros(1,p.NX);

if i_dc_dist == 1
    for i=1:co_range:p.NX+filter_range*2
        dd(i:min(i+co_range-1,p.NX+filter_range*2))=mean(tmp_DC)/exp(dis_mu+.5*dis_sigma^2)*lognrnd(dis_mu,dis_sigma);
    end
    ddd=filter(ones(1,filter_range)/filter_range,1,dd);
end
if i_dc_dist == 2
    for i=1:co_range:p.NX+filter_range*2
        dd(i:min(i+co_range-1,p.NX+filter_range*2))=log10(DC_min)+(log10(DC_max)-log10(DC_min))*rand;
    end
    ddd=filter(ones(1,filter_range)/filter_range,1,dd);
    ddd=10.^ddd;
end

p.DC=ddd(filter_range+1:p.NX+filter_range);
p.DC=max(p.DC,th_DC);

%------------------------------
Lb = min(p.MU.*p.DC./p.SIGMA./p.B)
%Lnuc = 1.3774*Lb;
%------------------------------

p.IC=ceil(p.N/2);
dx=p.L/p.NX;
Lb_over_dx = Lb/dx



p.NTOUT=1000;
p.NXOUT=1;
p.NSTOP=0;

[p,ot0,ox0]  = qdyn('run',p);
semilogy(ot0.t/year,ot0.v)
xlabel('Time (years)');
ylabel('Vmax');
filename = ['Hete_2D_uni_run_2_twm',num2str(twm),'L',num2str(p.L/1000.),'nx',num2str(p.NX),'W',num2str(p.W/1000.),'DC',num2str(DC_min),'to',num2str(DC_max),'_warmup.mat']
save(filename)  

V_0=ox0.v(:,end);
TH_0=ox0.th(:,end);

clear ot0 ox0

for irun = 1:1:num_run
    filename = ['Hete_2D_uni_run_2_twm',num2str(twm),'L',num2str(p.L/1000.),'nx',num2str(p.NX),'W',num2str(p.W/1000.),...
        'DC',num2str(DC_min),'to',num2str(DC_max),'ts',num2str((irun-1)*ts),'to',num2str(irun*ts),'.mat']

    % 
       p.TMAX = ts*year;  
       p.NTOUT=10;
    % 
       p.V_0 = V_0;
       p.TH_0= TH_0;
    %   %p.V_0 =  (ox1.v(:,end)+ox1.v(end:-1:1,end))/2;
    %   %p.TH_0=  (ox1.th(:,end)+ox1.th(end:-1:1,end))/2;
       [p,ot1,ox1]=qdyn('run',p);
    %p0=p;
    V_0 = ox1.v(:,end);
    TH_0= ox1.th(:,end);
    save(filename)  

    clear ot1 ox1

end
end

