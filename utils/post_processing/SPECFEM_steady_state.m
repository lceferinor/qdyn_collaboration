
function SPECFEM_steady_state(snaptime,cbar)

if nargin<2, cbar='off'; end

stime = snaptime*4e-3;
snapspec = FSEM3D_snapshot(snaptime);
qdynfile0=['fort.20001'];

disp('Loading qdyn.in ...');
p=Qdyn_read_in('qdyn.in');
disp('Done');

disp(['Loading QDYN output :  fort.20001 ...']);
o = Qdyn_read_ox_seq('fort.20001');
disp('Done');

p.X=o.X;
p.Y=o.Y;
p.Z=o.Z;

km = 1e3;

disp('Matching Values ...');

cc = [];
snapspec.DC = [];
snapspec.A = [];
snapspec.B = [];
snapspec.sigma = [];

snapspec.MU = p.MU;
snapspec.VS=p.VS;



for i=1:length(snapspec.X)
    
    [cc, id] = min((snapspec.X(i)*km-p.X(:)).^2 + (snapspec.Z(i)*km-p.Z(:)).^2);
    
    snapspec.DC(i)=p.DC(id);
    snapspec.A(i)=p.A(id);
    snapspec.B(i)=p.B(id);
    snapspec.sigma(i)=p.SIGMA(id);
    
    if mod(i,ceil(length(snapspec.X)/100)) == 0
        disp([num2str(i/ceil(length(snapspec.X)/100)) '% Complete'])
    end
end

Vdyninter=2*(snapspec.A.*snapspec.sigma./snapspec.MU).*snapspec.VS; 
Vdynsteady=2*((snapspec.A-snapspec.B).*snapspec.sigma./snapspec.MU).*snapspec.VS;


%qdynpar =  snapspec.Vx./Vdyninter';
qdynpar = snapspec.S;
%qdynpar = snapspec.Vx.*snapspec.S./snapspec.DC'; GAMMA

%idt1=find(qdynpar<0.9);
%qdynpar(idt1)=999;
%idt2=find(qdynpar>1.1);
%qdynpar(idt2)=999;

fig=figure;scatter(snapspec.X,snapspec.Z,[],qdynpar);
%slip=sqrt(snapspec.Dx.^2 + snapspec.Dz.^2);
%fig=figure;scatter(snapspec.X,snapspec.Z,[],slip)

title(num2str(stime));

if (strcmp(cbar,'on'))
c=colorbar;
c.Location ='northoutside';
c.Label.String = 'Slip rate (m/s)';
end
%c.Label.String = 'Slip (m)';
timesnap=num2str(stime);
%title(['Slip rate(m/s), time: ' timesnap ' s'],'FontSize',20)
%title(['time: ' timesnap ' s'],'FontSize',20)
%caxis([0 max(slip)]);
caxis([0 max(qdynpar)])
%caxis([0 10])
%caxis([0 0.001])
%caxis([0 1])


axis equal
xlim([min(snapspec.X) max(snapspec.X)]);
ylim([min(snapspec.Z) max(snapspec.Z)])
xlabel('X (km)','FontSize',20)
ylabel('Z (km)','FontSize',20)
set(gca,'FontSize',20)
figname=['fig' num2str(snaptime)];
%figname=['fig_slip' timesnap];
%print(fig,figname,'-djpeg')