% Reading rsf_hete_input...

fbh=read_rsf_hete('rsf_hete_input_file_fault2.txt');

figure(10),scatter3(fbh.X,fbh.Y,fbh.Z,[],fbh.tau);
figure(11),scatter3(fbh.X,fbh.Y,fbh.Z,[],fbh.SIGMA);
figure(12),scatter3(fbh.X,fbh.Y,fbh.Z,[],log(fbh.V_ini));
figure(13),scatter3(fbh.X,fbh.Y,fbh.Z,[],fbh.DC);

fm=read_rsf_hete('rsf_hete_input_file_fault1.txt');

figure(14),scatter(fm.Y,fm.Z,[],fm.tau);
figure(15),scatter(fm.Y,fm.Z,[],fm.SIGMA);
figure(16),scatter(fm.Y,fm.Z,[],fm.V_ini);
figure(17),scatter(fm.Y,fm.Z,[],log(fm.TH_ini));


idv=find(fbh.V_ini>1);

Fm0=FSEM3D_snapshot(42000,'./',1);
Fb0=FSEM3D_snapshot(42000,'./',2);

figure(10), subplot(2,1,1),scatter(Fm0.Y,Fm0.Z,[],abs(Fm0.Vx));
subplot(2,1,2),scatter(Fb0.Y,Fb0.Z,[],abs(Fb0.Vx))

figure(1),scatter3(Fm0.X,Fm0.Y,Fm0.Z,[],Fm0.Vx);hold on;
scatter3(Fb0.X,Fb0.Y,Fb0.Z,[],Fb0.Vx);hold on;
axis tight
axis equal

figure(2),scatter3(Fm0.X,Fm0.Y,Fm0.Z,[],Fm0.Tx);hold on;
scatter3(Fb0.X,Fb0.Y,Fb0.Z,[],Fb0.Tx);hold on;
axis tight
axis equal
