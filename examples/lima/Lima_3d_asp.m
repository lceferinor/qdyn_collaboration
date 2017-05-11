%%%%%%%% File for evaluating the earthquake rupture cycles in Lima %%%%%%%%
% Adapted from "test_3dfft.m": Prof. Ampuero, Dr. Galvez
% Modified on: 04/26/17
% By: Luis Ceferino

clear;
clc;
addpath ~/qdyn_developer/src 

%------------------------------
rand(1,floor(sum(100*clock)));
%------------------------------


year = 3600*24*365; % Should we put a longer timeframe? Would this timeline
                    % allow us to see the big events (Mw>8.5)?(Luis)
                    % PG: this is only a unit.
p = qdyn('set');

p.TMAX = 100000*year % Duration of the whole simulation.

p.DYN_TH_ON = 0.1; % Slip velocity at the nucletion phase starts.
p.DYN_TH_OFF = 0.001;

co = 4e6; % Cohesion imposed in the first kilometers depth.
co_limit = 2e3;
p.MU = 32.5e9; % Shear Modulus.
p.LAM = p.MU
p.VS = 3000.0
p.MU_SS = 0.6; % Reference friction coefficient.

p.MESHDIM=2;   % I think we can preserve the mesh size for now (Luis).
p.THETA_LAW=1; % 1) Should we change this? Would this be reasonable for the 
               % Peru subduction zone? (Luis)
               % PG: Please USE p.THETA_LAW = 1, Ageing LAW.
p.SIGMA=100.0e6; % 1) Should we change this? Would this be reasonable for the 
               % Peru subduction zone? 
               % 2) In the manual, it says that the normal stress remains 
               % constant unless it equals 1. Do we want this behavior in 
               % our model? (Luis)
               % PG: The normal stress (P.SIGMA) should change with depth and 
               %  constant value. Here 100 Mpa, but this could be changed.
%p.V_SS=1e-9;
p.V_SS=0.020/year % 20mm/year %1.97e-9;
                % PG : Using Villegas et al. (2016) GPS survey, Figure 3.
                % This is the relative velocity (m/s) of Nazca and South 
                % America plates extracted from Kendrick et al (2003) Azimut 82NE.
                % Is this OK for the steady-state slip  velocity? or do we
                % have to consider creeping of the locking slip?  (Luis)
p.A=0.003; % 1) Should we change this? Would this be reasonable for the 
           % Peru subduction zone? How does this relate to
           % different slip coupling in different locations in
           % the fault? (Luis)
           % PG: P.A and P.B should change with depth.
p.B=0.01;  % Same question than for p.A (Luis)
p.SIGMA_CPL=1; % 1) Is this related to p.SIGMA? Do you agree if we start the 
               % modeling with a coupled normal stress (p.SIGMA_CPL=1)? (Luis)
               % PG: This is relating with coupling of the normal stress with seismic waves.
%		SIGMA_CPL  normal stress coupling (only for dipping faults)
%			0 = disable
%			1 = enable
p.V2=0.01; % 1) Should we change this? Would this be reasonable for the 
           % Peru subduction zone? (Luis)
%p.L=80e3;
p.L=800e3; % Along-strike length of the faults (in meters) extracted from 
           % the fault geometry 
           % Is it the lenght of the fault? It is not clear from the manual 
           % for MESHDIM=2 (Luis)
           % PG : Yes. 
%p.W=80e3;
WH= 200e3; % Along-slip horizontal length of the fault (in meters) extracted 
           % from the fault geometry. This is the lenght of only the 
           % superficial part of the fault. 
           % 1) In the manual, it says that this parameter is distance between
           % displacment loading and fault. So, I am not sure if this
           % should be the along-slip horizontal lenght of the fault.
           % 2) In case it is, how far will the fault model cover. The 
           % initial 200 km is OK? With this horizontal 200 km, we are 
           % reaching up to a depth of 80 km. We could potentially model 
           % the initial 300 km of horizontal length and reach 110 km of 
           % depth, but the dip angle changes a bit after the 80 of depth. 
           % Would the model be able to incorporate this?
           % PG : I do not think we should go down to more than 80 km depth. 
           % At this depth, the megathurst events do not nucleate anymore
           % and normally ruptures do not go that deep.
           % 3) If we model below the 80 km of depth, we might have to
           % conseider that the properties of the fault are a bit different
           % below the 80 km, since the earthquake rates and magnitudes are
           % not as big as in the superficial part of the fault.
           % Any comments on this? (Luis)
           % PG: At this depth the code imposes creeping zone.
p.W=WH/cos(15/180*pi); %Along-dip distance, 
           % Dip angle 15 degress according to Villegas et al.(2016). Table 1.
p.NX=400;%100;
p.NW=100;%30; % Related to question on p.W. Is this related to the contact 
         % surface? (Luis)
         % PG : This is the number of elements along-dip fault. 
         % This value is too small. To get stable running the grid size (p.W/p.NW)=2km.

p.Z_CORNER=-WH/sind(15) %-60e3;
                  % PG : Should be the end of the along-dip fault. 
                  % In Km (extracted from Nocquet, 20146). A report from 
                  % Dr. Aguilar, mentions that active fault starts 30 km
                  % below the surface. Therefore in the superficial part of
                  % the fault we have: 80 km - 30 km = 50 km.
                  % Any comments on the fault starting at 30 km of depth?
                  % Should we start modeling it from depth=0 km?
                  % PG: In the first 0 to 5km we will use velocity strenghtening 
                  % And from 5-30km, velocity weakening. 
p.N=p.NX*p.NW;
p.DW(1:p.NW)=p.W/p.NW;
p.DIP_W(1:p.NW)=15.0; % Extracted from the fault geometry (Luis)
twm=0.5;
%ts=0.5;
p.ACC = 1e-14;

%% Defining space varying properties
%load('warmup_jp');
p.OX_DYN = 1;
p.OX_SEQ = 1;
%p.NX=512;
%p.RNS_LAW=0; % Is this contradicting to the V2 cutoff velocity (Delete?)
%p.N=p.NX*p.NW; % The manual only requires to define N for MESHDIM=1
%p.L=400e3;
%p.ACC=1e-10;
sigma0 = p.SIGMA;
p.A(1:p.NW) = p.A; % Along Dip
p.B(1:p.NW) = p.B; % Along Dip
p.SIGMA(1:p.NW) = p.SIGMA; % Along Dip
p.DC(1:p.NW) = p.DC; % Along Dip
p.Y = linspace(0,p.W,p.NW); % Along Dip
p.Z = linspace(0,p.Z_CORNER,p.NW); % Along Dip
% CHECK THIS
p.V_0(1:p.NW) = p.V_SS; % Along Dip (Is this OK?). What is the difference between
                        % Reference steady state and the initial V0
                        % conditions?
p.TH_0(1:p.NW) = p.DC/p.V_0; % (Along dip) From manual (TH_SS=DC/VS) Is TH_SS equal to TH_0 

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

tmp_B(1:id4)      = tmp_A(1:id4).*bam;
tmp_B(id4+1:id3)  = tmp_A(id4+1:id3).*linspace(bam,ba0,numel(id4+1:id3));   %increasing a/b below seismogenic zone
tmp_B(id3+1:id2)  = tmp_A(id3+1:id2).*ba0;   %a/b < 1 in seismogenic zone
tmp_B(id2+1:id1)  = tmp_A(id2+1:id1).*linspace(ba0,bam,numel(id2+1:id1));
tmp_B(id1+1:p.NW) = tmp_A(id1+1:p.NW).*bam;  %a/b >1 at shallow part 

tmp_SIGMA(1:idd)      = sigma0;
tmp_SIGMA(idd+1:p.NW) = linspace(sigma0,1e6,numel(idd+1:p.NW));
tmp_CO(1:ico_limit)      = 0;

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
    
    % I am not sure if we need DW per each element (Luis)
    p.DW((i-1)*p.NX+1:i*p.NX) = tmp_DW(i);
    p.DIP_W((i-1)*p.NX+1:i*p.NX) = tmp_DIP_W(i); 
end

% Lima fault potentially 6 asperities (2007,1974,1940,1966,1970,1996)
% PG: 1970 Earhtquake is not subduction event. 1996 is too far away from Lima and
% to save some computational time, we may need to obmit this event as well. In this way
% we can use the exact dimension from Villegas fault p.L = 620Km. 

n_asperities = 6;
asp_center_x = [0.128 0.340 0.431 0.622 0.758 0.818]*p.L;
asp_center_zy = [0.735 0.296 0.875 0.315 0.672 0.253]*p.W;
asp_elip_a = [0.214 0.335 0.294 0.243 0.156 0.130]/2*p.L;% Along strike
asp_elip_b = [0.798 0.667 0.649 0.699 0.425 0.261]/2*p.W;% Along dip
%----asperity spacing control
sp_int=25;   %spacing of asperities
scatter_mode=4;     % = 1 square; 
                    % = 2 triangular
                    % = 3 random
                    % = 4 Customized
bd_left=0.;       %asperities boundary
bd_right=0.;
bd_low=0.;
bd_up=0.;

asp_single = 0;

%----asperity property
l_asp0=40.0*1e3;          %asperity size lower limit in m
l_asp1=40.0*1e3;          %asperity size upper limit in m
ba_asp=1.5;       %b/a of asperity
dc_asp=0.002;       %Dc of asperity
sigma_asp=88.2e6*2;        %sigma(asp)




twm=20000;         %warmup time in years (what is this? (Luis))


p.IOT=zeros(size(p.X));
p.IASP=zeros(size(p.X));

%-----set asperity location (For Lima, we are in customized case)
i_asp=zeros(size(p.X));
switch scatter_mode
    case 1
        disp(['Asperities scattering: Square      Spacing:', num2str(sp_int)]);
        for ix=1:p.NX
            for iw=1:p.NW
                in=ix+(iw-1)*p.NX;
                if mod(ix,sp_int) == 1 && mod(iw,sp_int) == 1
                    i_asp(in)=1;
                end
            end
        end

    case 2
        sp_int_2 = round(sqrt((sp_int)^2-(sp_int/2)^2));  
        disp(['Asperities scattering: Triangular      Spacing:', num2str(sp_int)]);    
        for ix=1:p.NX
             for iw=1:p.NW
                 in=ix+(iw-1)*p.NX;
                 if mod(ix,sp_int) == 1 && mod(iw,sp_int_2*2) == 1
                     i_asp(in)=1;
                 end
                 if mod(ix,sp_int) == 1+floor(sp_int/2) && mod(iw,sp_int_2*2) == 1+sp_int_2
                     i_asp(in)=1;
                 end
             end
        end
     
    case 3
        pth = 1/(pi*sp_int^2);
        disp(['Asperities scattering: Random      Spacing:', num2str(sp_int)]);    
        for in=1:p.N
            if rand(1,1) <= pth
                i_asp(in)=1;
            end
        end

    case 4
    disp(['Asperities scattering: Customized ']);  
    

%     i_asp(63489+round(160/(400/512)))=1;
%     i_asp(64000-round(160/(400/512)))=1;

    %i_asp((63489+10):40:(64000-10))=1;
    for i = 1:n_asperities
        tmp_dist_asp = (asp_center_x(i) - p.X).^2 + (asp_center_zy(i) - sqrt(p.Y.^2+p.Z.^2)).^2;
        [M,I] = min(tmp_dist_asp);
        i_asp(I) = 1;
     end
    
end    


% Making asperities = 0 for outside asperity boundary
for ix=1:p.NX
     for iw=1:p.NW
         in=ix+(iw-1)*p.NX;
         if ix<=p.NX*bd_left || ix>=p.NX*(1-bd_right) || iw<=p.NW*bd_low || iw>=p.NW*(1-bd_up)
             i_asp(in)=0;
         end
     end
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
if asp_single == 0;
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
end



% Plot B/A ratio

scatter3(p.X,p.Y,p.Z,3,p.B./p.A>=2.0)
az = 0;
el = 90;
view(az, el);


if  asp_single == 1;
  for i=1:p.N
        if i_asp(i) == 1
            p.IASP(i) = 1;
            p.IOT(i) = 1;
            asp_count=asp_count+1;
            disp(['Setting single_cell asperity: ',num2str(asp_count),'/',num2str(asp_count_all)]);
            p.B(i)=+p.A(i)*ba_asp;
            p.DC(i)=dc_asp;
            p.SIGMA(j)=p.SIGMA(j)+(-p.SIGMA(j)+sigma_asp)*exp(-(dd/l_asp*2).^6);
            Lc_asp=p.MU*p.DC(i)/(p.SIGMA(i)*(p.B(i)-p.A(i)));
            disp(['  Normalized size Lx/Lc = ',num2str(dx/Lc_asp), '    Lw/Lc = ', num2str(dw/Lc_asp)]);
        end
  end
end

Vdyn=2*mean(p.A.*p.SIGMA./p.MU.*p.VS);
disp(['Vdyn = ' num2str(Vdyn)]);
p.DYN_TH_ON=Vdyn/10.;
p.DYN_TH_OFF=Vdyn/10.;

%------------------------------
Lb = min(p.MU.*p.DC./p.SIGMA./p.B);
Lnuc = 1.3774*Lb;
%------------------------------



filename = ['JP_3D','L',num2str(p.L/1000.),'nx',num2str(p.NX),'W',num2str(p.W/1000.),'nw',num2str(p.NW),'.mat']
p.IC=ceil(p.N/2);
dx=p.L/p.NX;
Lb_over_dx = Lb/dx
dw=p.W/p.NW;
Lb_over_dw = Lb/dw


p.TMAX=twm*year;
p.NTOUT=100;
p.NXOUT=1;
p.NSTOP=0;
p.DYN_FLAG=0;
p.DYN_M=10.^19.5;
p.DYN_SKIP = 1;

[p,ot1,ox1]  = qdyn('run',p);
semilogy(ot1.t/year,ot1.v)
xlabel('Time (years)');
ylabel('Vmax');
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

