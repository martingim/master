%%makes functions for transforming from sensor data to height
%%read tables from csv files
folder = "MEK4600_G3_21_10_02/";
data0 = readtable(folder + 'calib0.csv');
datapos2 = readtable(folder + 'calibpos2.csv');
datapos4 = readtable(folder + 'calibpos4.csv');
dataneg2 = readtable(folder + 'calibneg2.csv');
dataneg4 = readtable(folder + 'calibneg4.csv');

%%make the tables into vectors
data0 = table2array(data0(:,5:8));
datapos2 = table2array(datapos2(:,5:8));
datapos4 = table2array(datapos4(:,5:8));
dataneg2 = table2array(dataneg2(:,5:8));
dataneg4 = table2array(dataneg4(:,5:8));

%%add the height to the arrays
data0 = [data0 ones(size(data0, 1), 1)*0];
datapos2 = [datapos2 ones(size(datapos2, 1), 1)*0.02];
datapos4 = [datapos4 ones(size(datapos4, 1), 1)*0.04];
dataneg2 = [dataneg2 ones(size(dataneg2, 1), 1)*-0.02];
dataneg4 = [dataneg4 ones(size(dataneg4, 1), 1)*-0.04];

%%concatenate all the data from the sensors into one array
datasensor = vertcat(data0, datapos2, datapos4, dataneg2, dataneg4);

sensor1 =  fit(datasensor(:,1), datasensor(:,5), 'poly2');
sensor2 =  fit(datasensor(:,2), datasensor(:,5), 'poly2');
sensor3 =  fit(datasensor(:,3), datasensor(:,5), 'poly2');
sensor4 =  fit(datasensor(:,4), datasensor(:,5), 'poly2');

