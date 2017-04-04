disp('Loading qdyn.in ...');
p=Qdyn_read_in_old('qdyn.in');
disp('Done');

disp(['Loading QDYN output :  fort. ...']);
%o = Qdyn_read_ox_seq('fort_xyz');
%disp('Done');
XYZ=dlmread('XYZ.dat');
disp('Done');
o.X = XYZ(:,1);
o.Y = XYZ(:,2);
o.Z = XYZ(:,3);
p.X=o.X;
p.Y=o.Y;
p.Z=o.Z;

qvar=p.DC;
%qvar=p.V_0;
pf   = find(p.X==min(p.X));
Zp   = p.Z(pf);
ABp  = p.A(pf)./p.B(pf);
A_Bp = p.A(pf)-p.B(pf);
sigmap = p.SIGMA(pf);

fig=figure('Color',[1 1 1]);
subplot(1,3,1),plot(sigmap/1e6,Zp/1e3,'LineWidth',2);
ylabel('Depth (Km)','FontSize',25)
xlabel('Sigma (Mpa)','FontSize',25)

xlim([0 100]);
%ylim([-30 0]);
ylim([-20 0]);
set(gca,'FontSize',25)


subplot(1,3,2),plot(ABp,Zp/1e3,'LineWidth',2);

xlim([0.5 2]);
%ylim([-30 0]);
ylim([-20 0]);
xlabel('a/b','FontSize',25)
%ylabel('Z (km)','FontSize',20)
set(gca,'FontSize',25)
%figname=['fig_friction_profile'];
%figname=['figV0'];
%print(fig,figname,'-djpeg')



subplot(1,3,3),plot(A_Bp,Zp/1e3,'LineWidth',2);
xlim([-0.01 0.01]);
%ylim([-30 0]);
ylim([-20 0]);
xlabel('(a-b)','FontSize',25)
%ylabel('Z (km)','FontSize',20)
set(gca,'FontSize',25)
figname=['fig_friction_profile'];
%figname=['figV0'];
print(fig,figname,'-djpeg')


fig1=figure;
histogram(log(p.DC));
xlabel('log(Dc) meters','FontSize',20)
ylabel('number of fault points','FontSize',20)
xlim([log(0.011) log(max(p.DC))])
ylim([0 3500])
set(gca,'FontSize',20)
figname1=['DC_Historgram'];
%figname=['figV0'];
print(fig1,figname1,'-djpeg')





