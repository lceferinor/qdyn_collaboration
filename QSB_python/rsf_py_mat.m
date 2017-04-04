rsf_mat=load('rsf_hete_input_file_matlab.txt');
rsf_py=load('rsf_hete_input_file_python.txt');


figure(1),scatter(rsf_mat(:,1),rsf_mat(:,2),[],rsf_mat(:,4))
figure(2),scatter(rsf_py(:,1),rsf_py(:,2),[],rsf_py(:,4))