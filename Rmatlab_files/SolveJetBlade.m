x = dlmread('/RNAmf-Reproducibility-main/Rmatlab_files/generate_text/temp_to_matlab.txt', ',', 1, 0);
y = [];
CPU_time = [];
for i=1:size(x,1)
    a1 = x(i,1);
    a2 = x(i,2);
    t = x(i,3);
    tic; 
    [avg_u] = bladeMatlab(a1,a2,t);
    CPU_time(i) = toc;    
    y(i) = avg_u;
end
final_data = [x,y',CPU_time'];
dlmwrite('/RNAmf-Reproducibility-main/Rmatlab_files/generate_text/temp_to_r.txt',final_data,'delimiter',',','precision', 32);