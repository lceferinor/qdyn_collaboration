%%%%%%%% File for evaluating the earthquake rupture cycles in Lima %%%%%%%%
% Adapted from "test_3dfft.m": Prof. Ampuero, Dr. Galvez
% Modified on: 06/22/17
% By: Luis Ceferino

clear;
clc;
scratch=getenv('SCRATCH');
addpath(strcat(scratch,'/qdyn_collaboration/src')); 

%------------------------------
rand(1,floor(sum(100*clock)));
%------------------------------


year = 3600*24*365; % Should we put a longer timeframe? Would this timeline
                    % allow us to see the big events (Mw>8.5)?(Luis)
                    % PG: this is only a unit.
p = qdyn('set');

%% Defining general properties
co = 4e6; % Cohesion imposed in the first kilometers depth.
p.MU = 32.5e9; % Shear Modulus.
p.LAM = p.MU;
p.VS = 3000.0;
p.MU_SS = 0.6; % Reference friction coefficient.
p.MESHDIM=2;   % 3D
p.THETA_LAW=1; % Edited on 06/21 Ageing Law
p.SIGMA=100.0e6; % PG: The normal stress. Here 100 Mpa, and later it changes
                % with depth.
p.V_SS=0.020/year; % 20mm/year;
                % PG : Using Villegas et al. (2016) GPS survey, Figure 3.
p.A=0.01; % From Tohoku. (p.B defined as function of p.B and depth later)
p.B=0.006; % From Tohoku. (Referential values, it will be redefined later)
p.DC=0.8; % From Tohoku, Tohoku is 0.8
p.SIGMA_CPL=1; % PG: This is relating with coupling of the normal stress 
                % with seismic waves.
p.RNS_LAW=0; %LC: 06/22           
p.L=800e3; % Along-strike length of the faults (in meters) extracted from 
           % the fault geometry.
WH= 200e3; % Along-slip horizontal length of the fault (in meters) extracted 
           % from the fault geometry.
p.W=WH/cos(15/180*pi); %Along-dip distance, 
           % Dip angle 15 degress according to Villegas et al.(2016). Table 1.
p.NX=2048*2/8; % So that dx is smaller than Lb
p.NW=128*8/8; % So that dw is smaller than Lb

p.Z_CORNER=-WH*tan(15/180*pi); %-60e3;
                  % PG : Should be the end of the along-dip fault. 
                  % PG: In the first 0 to 5km we will use velocity strenghtening 
                  % And from 5-30km, velocity weakening. 
p.N=p.NX*p.NW;
p.DW(1:p.NW)=p.W/p.NW;
p.DIP_W(1:p.NW)=15.0; % Extracted from the fault geometry (Luis)
p.ACC = 1e-10; %Modified LC *(06/22) Same as in Tohoku


%% Defining space varying properties
sigma0 = p.SIGMA;
p.A(1:p.NW) = p.A; % Along Dip
p.B(1:p.NW) = p.B; % Along Dip
p.SIGMA(1:p.NW) = p.SIGMA; % Along Dip
p.DC(1:p.NW) = p.DC; % Along Dip
p.Y = linspace(0,p.W,p.NW); % Along Dip
p.Z = linspace(p.Z_CORNER,0,p.NW); % Along Dip


% Initial values
p.V_0(1:p.NW) = p.V_SS; 
p.TH_0(1:p.NW) = p.DC/p.V_0; % (Along dip) From manual (TH_SS=DC/VS) Is TH_SS equal to TH_0 

% Setting temps
tmp_A=p.A;
tmp_B=p.B;
tmp_SIGMA=p.SIGMA;
tmp_DC=p.DC;
tmp_Y=p.Y;
tmp_Z=p.Z;
tmp_DW=p.DW;
tmp_DIP_W=p.DIP_W;
tmp_CO(1:p.NW)=co;
tmp_V_0=p.V_0;
tmp_TH_0=p.TH_0;

% Here changing SIGMA, A, B along depth.
ba0=1.5;  %b/a at seismogenic zone
bam=0.6;  %b/a at shallow/deeper part
ZW = p.Z; 
%% PG: Addapt the below parameters for Lima Fault has three transitions at depth. 
% LC: p.Z goes from 0 to -p.ZCorner, so the index convention has changed
% Finding layers:
% Sigma
dd = abs(0.50*p.Z_CORNER); % Ratio from Ventura (LC) 
idd=min(find(ZW+dd>=0)); %LC: Change < by >
tmp_SIGMA(1:idd)      = sigma0;
tmp_SIGMA(idd+1:p.NW) = linspace(sigma0,1e6,numel(idd+1:p.NW));

% B
d1 = abs(0.10*p.Z_CORNER); % Ratio from Ventura (LC) 
id1=min(find(ZW+d1>=0));
d2 = abs(0.25*p.Z_CORNER); % Ratio from Ventura (LC) 
id2=min(find(ZW+d2>=0));
d3 = abs(0.60*p.Z_CORNER); % Ratio from Ventura (LC) 
id3=min(find(ZW+d3>=0));
d4 = abs(0.75*p.Z_CORNER); % Ratio from Ventura (LC) 
id4=min(find(ZW+d4>=0));
tmp_B(1:id4)      = tmp_A(1:id4).*bam;
tmp_B(id4+1:id3)  = tmp_A(id4+1:id3).*linspace(bam,ba0,numel(id4+1:id3));   %increasing a/b below seismogenic zone
tmp_B(id3+1:id2)  = tmp_A(id3+1:id2).*ba0;   %a/b < 1 in seismogenic zone
tmp_B(id2+1:id1)  = tmp_A(id2+1:id1).*linspace(ba0,bam,numel(id2+1:id1));
tmp_B(id1+1:p.NW) = tmp_A(id1+1:p.NW).*bam;  %a/b >1 at shallow part 

% CO
co_limit = abs(0.1*p.Z_CORNER);
ico_limit=min(find(ZW+co_limit>=0));
tmp_CO(1:ico_limit)= 0; 


% Replicate along the strike direction
for i=1:p.NW
    p.X((i-1)*p.NX+1:i*p.NX) = linspace(0,p.L,p.NX);
    p.Y((i-1)*p.NX+1:i*p.NX) = tmp_Y(i);
    p.Z((i-1)*p.NX+1:i*p.NX) = tmp_Z(i);
    p.A((i-1)*p.NX+1:i*p.NX) = tmp_A(i);
    p.B((i-1)*p.NX+1:i*p.NX) = tmp_B(i);
    p.SIGMA((i-1)*p.NX+1:i*p.NX) = tmp_SIGMA(i);
    p.DC((i-1)*p.NX+1:i*p.NX) = tmp_DC(i);
    p.V_0((i-1)*p.NX+1:i*p.NX) = tmp_V_0(i);
    p.TH_0((i-1)*p.NX+1:i*p.NX) = tmp_TH_0(i);
    p.CO((i-1)*p.NX+1:i*p.NX)    = tmp_CO(i);
    
end

% Plot B/A ratio
%figure
%scatter3(p.X,p.Y,p.Z,3,p.A./p.B)

%% Set asperity
% Lima fault potentially 6 asperities (2007,1974,1940,1966,1970,1996)
% PG: 1970 Earhtquake is not subduction event. 1996 is too far away from Lima and
% to save some computational time, we may need to obmit this event as well. In this way
% we can use the exact dimension from Villegas fault p.L = 620Km. 
n_asperities = 5;
asp_center_x = [0.128 0.340 0.431 0.622 0.818]*p.L;
asp_center_zy = [0.735 0.296 0.875 0.315 0.253]*p.W;
asp_elip_a = [0.214 0.335 0.294 0.243 0.130]/2*p.L;% Along strike
asp_elip_b = [0.798 0.667 0.649 0.699 0.261]/2*p.W;% Along dip


%----asperity property
l_asp0=40.0*1e3;          %asperity size lower limit in m
l_asp1=40.0*1e3;          %asperity size upper limit in m
ba_asp=1.5;       %b/a of asperity
dc_asp=0.025;       % Taken from AECOM example (Luis)
dc_asp=0.025*8; % So that it runs
sigma_asp=88.2e6*2;        %sigma(asp)


twm=20000;         %warmup time in years (what is this? (Luis))
%twm=500;         

p.IOT=zeros(size(p.X));
p.IASP=zeros(size(p.X));

%-----set asperity location (For Lima, we are in customized case)
i_asp=zeros(size(p.X));
disp(['Asperities scattering: Customized ']);  


for i = 1:n_asperities
    tmp_dist_asp = (asp_center_x(i) - p.X).^2 + (asp_center_zy(i) - sqrt(p.Y.^2+p.Z.^2)).^2;
    [M,I] = min(tmp_dist_asp);
    i_asp(I) = 1;
 end


% Count asperities
asp_count_all=0;
for i=1:p.N
    if i_asp(i) == 1
        asp_count_all=asp_count_all+1;         
    end
end

asp_count=0;
%----set asperity property


% Multiple asperities (this is the case of Lima Asperity)
for i=1:p.N
    if i_asp(i) == 1
        tmp_sq_dis_asp = (asp_center_x - p.X(i)).^2 + ...
            (asp_center_zy - (p.Y(i)^2+p.Z(i)^2)).^2;
        [M,I_asp] = min(tmp_sq_dis_asp);

        p.IASP(i) = 1;
        p.IOT(i) = 1;
        asp_count=asp_count+1;
        l_asp=l_asp0+(l_asp1-l_asp0)*rand(1);
        disp(['Setting asperity: ',num2str(asp_count),'/',num2str(asp_count_all)]);
        for j=1:1:p.N
            %dd=sqrt((p.X(j)-p.X(i))^2+(p.Y(j)-p.Y(i))^2+(p.Z(j)-p.Z(i))^2);
            sq_dd_ratio=(p.X(j)-p.X(i))^2/...
                                        (asp_elip_a(I_asp)^2) +...
                         ((p.Y(j)-p.Y(i))^2+(p.Z(j)-p.Z(i))^2)/...
                                        (asp_elip_b(I_asp)^2); %Distance ratio:
                                                    % Smaller than one
                                                    % if it is inside
                                                    % the ellipe, and
                                                    % larger than one
                                                    % otherwise
            %p.B(j)=p.B(j)+(-p.B(j)+p.A(j)*ba_asp)*exp(-(dd/l_asp*2).^6);
            %p.DC(j)=p.DC(j)+(-p.DC(j)+dc_asp)*exp(-(dd/l_asp*2).^6);
            %p.SIGMA(j)=p.SIGMA(j)+(-p.SIGMA(j)+sigma_asp)*exp(-(dd/l_asp*2).^6);
            % The exponent equal 3 since the ratio is an square measure
            % of distance ratio in the ellipe. In case the asperity
            % were a circle, the result would equal distance ratio to
            % the power of 6
            p.B(j)=p.B(j)+(-p.B(j)+p.A(j)*ba_asp)*exp(-(sq_dd_ratio^3));
            p.DC(j)=p.DC(j)+(-p.DC(j)+dc_asp)*exp(-(sq_dd_ratio^3));
            p.SIGMA(j)=p.SIGMA(j)+(-p.SIGMA(j)+sigma_asp)*exp(-(sq_dd_ratio^3));            

        end
        Lc_asp=p.MU*p.DC(i)/(p.SIGMA(i)*(p.B(i)-p.A(i)));
        %disp(['  Normalized size L/Lc = ',num2str(l_asp/Lc_asp)]);
        disp(['  Normalized size along strike 2a/Lc = ',...
            num2str(2*asp_elip_a(I_asp)/Lc_asp)]);
        disp(['  Normalized size along slip 2b/Lc = ',...
            num2str(2*asp_elip_b(I_asp)/Lc_asp)]);
    end
end



% Plot B/A ratio
%figure
%scatter3(p.X,p.Y,p.Z,3,p.A./p.B)

%az = 0;
%el = 90;
%view(az, el);

%% Controlling parameters

Vdyn=2*mean(p.A.*p.SIGMA./p.MU.*p.VS);
disp(['Vdyn = ' num2str(Vdyn)]);
p.DYN_TH_ON=Vdyn/10.;
p.DYN_TH_OFF=Vdyn/10.;

%------------------------------
Lb = min(p.MU.*p.DC./p.SIGMA./p.B);
Lnuc = 1.3774*Lb;
dx = p.L/p.NX;
dw = p.W/p.NW;
disp(['dx/Lb is ', num2str(dx/Lb), '. It should be smaller than 1.']);
disp(['dw/Lb is ', num2str(dw/Lb), '. It should be smaller than 1.']);

%------------------------------


filename = ['LC_3D','L',num2str(p.L/1000.),'nx',num2str(p.NX),'W',num2str(p.W/1000.),'nw',num2str(p.NW),'.mat']
p.IC=ceil(p.N/2);


p.OX_DYN = 1;
p.OX_SEQ = 0; % fort.19 survives!

p.TMAX=twm*year;
p.NTOUT=100;
p.NXOUT=1;
p.NSTOP=0;
p.DYN_FLAG=0;
p.DYN_M=10.^19.5;
p.DYN_SKIP = 1;
p.NPROCS=16;


[p,ot1,ox1]  = qdyn('run',p);
semilogy(ot1.t/year,ot1.v)
xlabel('Time (years)');
ylabel('Vmax');
saveas(gcf,'velocity_Lima.png');
% 
%   p.TMAX = ts*year;  
%   p.NTOUT=1;
% 
%   p.V_0 = ox1.v(:,end);
%   p.TH_0= ox1.th(:,end);
%   %p.V_0 =  (ox1.v(:,end)+ox1.v(end:-1:1,end))/2;
%   %p.TH_0=  (ox1.th(:,end)+ox1.th(end:-1:1,end))/2;
%   [p,ot,ox]=qdyn('run',p);

V_0 = ox1.v(:,end);
TH_0= ox1.th(:,end);
save(filename)  
save('warmup_jp_3d.mat','p', 'V_0', 'TH_0');

