V0=load('snap_V_loc0.dat');
V1=load('snap_V_loc1.dat');
V2=load('snap_V_loc2.dat');
V3=load('snap_V_loc3.dat');
V4=load('snap_V_loc4.dat');
V5=load('snap_V_loc5.dat');


figure(10),
scatter(V0(:,1),V0(:,3),[],V0(:,3),'r');
hold on;
scatter(V1(:,1),V1(:,3),[],V1(:,3),'g');
hold on;
scatter(V2(:,1),V2(:,3),[],V2(:,3),'b');
hold on;
scatter(V3(:,1),V3(:,3),[],V3(:,3),'y');
hold on;
scatter(V4(:,1),V4(:,3),[],V4(:,3),'k');
hold on;
scatter(V5(:,1),V5(:,3),[],V5(:,3),'g');

V0g=load('snap_V_glo0.dat');
V1g=load('snap_V_glo1.dat');
V2g=load('snap_V_glo2.dat');
V3g=load('snap_V_glo3.dat');

figure(1),
scatter(V0g(:,1),V0g(:,3),[],V0g(:,4),'r');
%hold on;
figure(2)
scatter(V1g(:,1),V1g(:,3),[],V1g(:,4),'b');
%hold on;
figure(3)
scatter(V2g(:,1),V2g(:,3),[],V2g(:,4),'g');
%hold on;
figure(4)
scatter(V3g(:,1),V3g(:,3),[],V3g(:,4),'y');


V0i=load('snap_V_gli0.dat');
V1i=load('snap_V_gli1.dat');
V2i=load('snap_V_gli2.dat');
V3i=load('snap_V_gli3.dat');

figure(1),
scatter(V0i(:,1),V0i(:,3),[],V0i(:,4),'r');
%hold on;
figure(2)
scatter(V1i(:,1),V1i(:,3),[],V1i(:,4),'b');
%hold on;
figure(3)
scatter(V2i(:,1),V2i(:,3),[],V2i(:,4),'g');
%hold on;
figure(4)
scatter(V3i(:,1),V3i(:,3),[],V3i(:,4),'y');



