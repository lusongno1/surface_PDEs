%turn Eig_Clu to eig_cal
eig_cal = []
path = 'NISO_k2N32';
for i=1:10
    eig_cal = [eig_cal Eig_Clu{i,2}]
end
eig_cal = eig_cal.';
dlmwrite([path '/eig_cal.txt'],eig_cal,'delimiter','\n','precision',15)