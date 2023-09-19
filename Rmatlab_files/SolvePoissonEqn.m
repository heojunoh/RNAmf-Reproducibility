
x = dlmread('generate_text/temp_to_matlab.txt', ',', 1, 0);
y = [];
CPU_time = [];
for i=1:size(x,1)
    a = x(i,1);
    t = x(i,2);
    tic; 
    [avg_u] = poissonMatlab(a,t);
    CPU_time(i) = toc;    
    y(i) = avg_u;
end
final_data = [x,y',CPU_time'];
dlmwrite('generate_text/temp_to_r.txt',final_data,'delimiter',',','precision', 32);