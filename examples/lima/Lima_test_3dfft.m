%%%%%%%% File for evaluating the earthquake rupture cycles in Lima %%%%%%%%
% Adapted from "test_3dfft.m": Prof. Ampuero, Dr. Galvez
% Modified on: 04/26/17
% By: Luis Ceferino

clear;
addpath ~/qdyn_developer/src 

%------------------------------
rand(1,floor(sum(100*clock)));
%------------------------------


year = 3600*24*365; % Should we put a longer timeframe? Would this timeline
                    % allow us to see the big events (Mw>8.5)?(Luis)
p = qdyn('set');

p.MESHDIM=2;   % I think we can preserve the mesh size for now (Luis).
p.THETA_LAW=2; % 1) Should we change this? Would this be reasonable for the 
               % Peru subduction zone? (Luis)
p.SIGMA=0.5e6; % 1) Should we change this? Would this be reasonable for the 
               % Peru subduction zone? 
               % 2) In the manual, it says that the normal stress remains 
               % constant unless it equals 1. Do we want this behavior in 
               % our model? (Luis)
%p.V_SS=1e-9;
p.V_SS=1.97e-9; % This is the relative velocity (m/s) of Nazca and South 
                % America plates extracted from Kendrick et al (2003) Azimut 82NE.
                % Is this OK for the steady-state slip  velocity? or do we
                % have to consider creeping of the locking slip?  (Luis)
p.A=0.003; % 1) Should we change this? Would this be reasonable for the 
           % Peru subduction zone? How does this relate to
           % different slip coupling in different locations in
           % the fault? (Luis)
p.B=0.01;  % Same question than for p.A (Luis)
p.SIGMA_CPL=1; % 1) Is this related to p.SIGMA? Do you agree if we start the 
               % modeling with a coupled normal stress (p.SIGMA_CPL=1)? (Luis)
p.V2=0.01; % 1) Should we change this? Would this be reasonable for the 
               % Peru subduction zone? (Luis)
%p.L=80e3;
p.L=800e3; % Along-strike length of the faults (in meters) extracted from 
           % the fault geometry 
           % Is it the lenght of the fault? It is not clear from the manual 
           % for MESHDIM=2 (Luis)
%p.W=80e3;
p.W=200e3; % Along-slip horizontal length of the fault (in meters) extracted 
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
           % 3) If we model below the 80 km of depth, we might have to
           % conseider that the properties of the fault are a bit different
           % below the 80 km, since the earthquake rates and magnitudes are
           % not as big as in the superficial part of the fault.
           % Any comments on this? (Luis)

p.NX=100;
p.NW=25; % Related to question on p.W. Is this related to the contact 
         % surface? (Luis)

p.Z_CORNER=-50e3; % In Km (extracted from the fault geometry). A report from 
                  % Dr. Aguilar, mentions that active fault starts 30 km
                  % below the surface. Therefore in the superficial part of
                  % the fault we have: 80 km - 30 km = 50 km.
                  % Any comments on the fault starting at 30 km of depth?
                  % Should we start modeling it from depth=0 km?
p.N=p.NX*p.NW;
p.DW(1:p.NW)=p.W/p.NW;
p.DIP_W(1:p.NW)=13.0; % Extracted from the fault geometry (Luis)
twm=0.5;
%ts=0.5;
p.ACC = 1e-14;

%------------------------------
Lb = p.MU*p.DC/p.SIGMA/p.B;
Lnuc = 1.3774*Lb;
%------------------------------

filename = ['test_2d_fft_ab',num2str(p.A/p.B),'L',num2str(p.L/1000),'nx',num2str(p.NX),'W',num2str(p.W/1000),'nw',num2str(p.NW),'z',num2str(p.Z_CORNER/1000),'.mat']
p.IC=ceil(p.N/2);
dx=p.L/p.NX;
Lb_over_dx = Lb/dx


p = qdyn('set',p);

Lc=Lb*(p.B/(p.B-p.A));
disp(['  Lc=',num2str(Lc),'  L/Lc=',num2str(p.L/Lc),'  W/Lc=',num2str(p.W/Lc)]);
Linf=2/pi*(p.B/(p.B-p.A))^2*Lb;
disp(['  Linf=',num2str(Linf),'  L/Linf=',num2str(p.L/Linf),'  W/Linf=',num2str(p.W/Linf)]);

p.TMAX=twm*year;

%for i=1:1:floor(p.N*0.05)
%    p.V_0(i) = 1.01 *p.V_SS ;
%end
%for i=floor(p.N*0.05)+1:1:p.N
%    p.V_0(i)=p.V_SS;
%end
 p.V_0 = 1.01*p.V_SS ;
% p.V_0 = p.V_0/mean(p.V_0)*p.V_SS;
 p.V_00=p.V_0;
%   for i=1:1:p.N
%       p.V_0(i) = p.V_00(mod((i+1024),p.N)+1);
%   end
% p.V_0(1:p.N) = p.V_SS*1e-80;
% p.V_0(2) = 1;


% Would it be possible to get the geometry of the ruptured elements, their
% location and the time at which it ruptured?
p.NTOUT=10;
p.NXOUT=1;
p.NSTOP=0;

[p,ot1,ox1]  = qdyn('run',p);
%semilogy(ot1.t/year,ot1.v)
semilogy(ot1.t/year,ot1.vc)
xlabel('Time (years)');
ylabel('Vmax');
% 
%   p.TMAX = ts*year;  
%   p.NTOUT=10;
% 
%   p.V_0 = ox0.v(:,end);
%   p.TH_0= ox0.th(:,end);
%   %p.V_0 =  (ox1.v(:,end)+ox1.v(end:-1:1,end))/2;
%   %p.TH_0=  (ox1.th(:,end)+ox1.th(end:-1:1,end))/2;
%   [p,ot1,ox1]=qdyn('run',p);

save(filename)  


