from pylab import *

data = loadtxt('data.txt');
intensity = 0.5;

#Get the index for 100 randomly selected samples
index = rand(200);
index = array(index*data.size,dtype = int);

data_short = data;
data_short[index] = data_short[index] + intensity*data_short[index];


#plot(data_short[index[1]-20:index[1]+20]);
#show();
savetxt('data_short.txt',data_short);
savetxt('index_short.txt',index);
