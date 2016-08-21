from pylab import *

data = loadtxt('data.txt');
length = size(data);

'''Generate Noise'''
N=length # no of data points
k=2   # no of sets of data with varying noise
# generate the data points and add noise
Y=meshgrid(data,ones(k),indexing='ij')[0] # make k copies
scl = [0.5,0.5];
n=dot(randn(N,k),diag(scl)) # generate k vectors
yy=Y+n    # add noise to signal
'''End generation of noise'''

##Separate the signal as signal from two sensors
t1 = transpose(yy[:,0]);
t2 = transpose(yy[:,1]);

'''Add the short faults to one of the sensors(t2)'''

intensity = 0.5;

#Get the index for 100 randomly selected samples
index = rand(9000);
flag = zeros(size(t1));
index = array(index*(data.size-1200),dtype = int)+1200;#Make sure there are no faults in the first 1000 samples as they are used to estimate future values
t2[index] = t2[index] + intensity*t2[index];
flag[index] = 1;


''''''
c = array([t1,t2]);

##Estimating t2 from t1.
##Consider the first M samples when calculating mean and covariance.
t2_est = [];
M = 1000; # M indicates how many previous samples to consider during estimation
u_t1 = mean(t1[0:M]);
u_t2 = mean(t2[0:M]);
C = array([c[0,0:M],c[1,0:M]]);
for i in range(M,length):
	covar = cov(C);
	#Estimate the value of t2
	t2_est.append(u_t2 + (covar[1,1]/covar[0,0])*(t1[i]-u_t1));
t2_est = array(t2_est);

diff = abs(t2_est - t2[M:]);
data_diff = zeros((size(diff),2));
data_diff[:,0] = diff;
data_diff[:,1] = flag[M:];
savetxt('diff.txt',data_diff);
print std(diff);

###Show the results
N = 100; # N is the number of samples to be clumped together during calculation of stdev
stdev = [];
i = 0;
while True:
	if i+N > size(diff):
		break;
	A = diff[i:i+N];
	stdev.append(abs(std(A)));
	i = i+N;
hist(stdev,100);
show();

plot(t2[M:M+100],'b',label = 'True Sensor Value');
#plot(t2[M:],'b',label = 'True Sensor Value');
plot(t2_est[0:100],'g',label = 'Estimated Sensor Value');
#plot(t2_est[0:],'g',label = 'Estimated Sensor Value');
legend(loc = 'upper right');
show();

'''Checking for short faults'''
threshold = 3.5;
i = 0;
short_count = 0;
ind = [];
while i<diff.size:
	if diff[i] > threshold: #Occurence of short fault
		short_count += 1;
		ind.append(i+1);
		#Go through all the samples that are part of the same short-fault
		i = i +1;
		while (i<diff.size):
			if diff[i] > threshold:
				break;
			i = i+1;
	else :
		i = i+1;
print "Number of short faults detected: {}".format(short_count);
indexes = sort(index)+1-M;
#print diff[array(index)-M];



##Find out number of true positives
correct_count = 0;
for i in range(0,len(ind)):
	if ind[i] in indexes:
		correct_count += 1;

print "Number of short faults detected correctly: {}".format(correct_count);
''''''

