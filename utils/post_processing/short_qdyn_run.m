%-----------------------------------------------
setenv('DYLD_LIBRARY_PATH', '/usr/local/bin/')
setenv('OMP_NUM_THREADS', '12')

% Loading qdyn.in files.
%-----------------------------------------------
disp('Loading qdyn_initial.in ...');
p=Qdyn_read_in();
disp('Done');

disp(['Loading QDYN output :  fort.XX ...']);
o = Qdyn_read_ox_seq('./hete_3d_ss_compute_stress_3D/fort.1001');
disp('Done');
year = 3600*24*365;

p.X=o.X;
p.Y=o.Y;
p.Z=o.Z;

p.V_0 = o.v ;
p.TH_0 = o.th;

p.NSTOP = 0;
p.TMAX = 100000*year;

p.DYN_TH_ON  = 0.1;
p.DYN_TH_OFF = 0.001;

p.N = p.NX*p.NW;

p.NTOUT=100;
p.NXOUT=1;
p.OX_SEQ=1;
% Run qdyn
[p,ot1,ox1]=qdyn('run',p);