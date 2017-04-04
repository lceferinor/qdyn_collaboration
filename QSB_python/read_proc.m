fil2='proc000002_fault_xyz.dat'
fid2=fopen(fil2,'r');
nglob2=str2num(fgetl(fid2));
xyz2=textscan(fid2,'%f %f %f\n');
x2=xyz2{1};y2=xyz2{2};z2=xyz2{3};
fclose(fid2)

fil3='proc000003_fault_xyz.dat'
fid3=fopen(fil3,'r');
nglob3=str2num(fgetl(fid3));
xyz3=textscan(fid3,'%f %f %f\n');
x3=xyz3{1};y3=xyz3{2};z3=xyz3{3};
fclose(fid3)