# 9 states x y z vx vy vz yaw pitch roll, 6 control factors yv pv rv ax ay az

import numpy as np
import matplotlib.pyplot as plt
import sys

dti = 1.0/100.0 #rate of IMU 100HZ
dtg = 1.0 #rate of GPS 1Hz

# unpacking GPS data
data = np.loadtxt(sys.argv[1],skiprows=1);
tsg = data[:,0];
alt = data[:,3];
x = data[:,8];
y = data[:,9];
#unpacking IMU data
data = np.loadtxt(sys.argv[2],skiprows=1);
tsi = data[:,0];
ax = data[:,1]; # wrt IMU system
ay = data[:,2];
az = data[:,3] + 9.8;
yv = data[:,4];
pv = data[:,5];
rv = data[:,6];
yaw = data[:,10];
pitch = data[:,11];
roll = data[:,12];

# calculating initial velocities
vx0 = (x[2]-x[1])/dtg;
vy0 = (y[2]-y[1])/dtg;
vz0 = (alt[2]-alt[1])/dtg;

# covariance of initial state
P = np.diag([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]);

# process noise matrix
'''sGPS = 0.5*10.0*dti**2; # worst case of drift during open loop operation of IMU
sVel = 10.0*dti;
sYPR = 0.5*dti;
Q = np.diag([sGPS**2, sGPS**2, sGPS**2, sVel**2, sVel**2, sVel**2, sYPR**2, 0.0, sYPR**2, 0.0, 0.0, 0.0, 0.0]);
'''
qtrans = 0.05;
qrot = 0.05;
Q = np.diag([qtrans, qtrans, qtrans, qtrans, qtrans, qtrans, qrot, qrot, qrot]);
# measurement Jacobian
JH = np.matrix([[1,0,0,0,0,0,0,0,0],
				[0,1,0,0,0,0,0,0,0],
				[0,0,1,0,0,0,0,0,0],
				[0,0,0,1,0,0,0,0,0],
				[0,0,0,0,1,0,0,0,0],
				[0,0,0,0,0,1,0,0,0]]);

# measurement error covariance
sGPS = 0.5;
R = np.diag([sGPS**2, sGPS**2, sGPS**2, sGPS**2, sGPS**2, sGPS**2]);

# identity matrix
I = np.eye(9);
print (I);

# calculating velocities
vx = np.hstack((vx0, np.diff(x)))/dtg;
vy = np.hstack((vy0, np.diff(y)))/dtg;
vz = np.hstack((vz0, np.diff(alt)))/dtg;

# initial state
X = np.matrix([x[0], y[0], alt[0], vx[0], vy[0], vz[0], yaw[0], pitch[0], roll[0]]).T;
U = np.matrix([yv[0], pv[0], rv[0], ax[0], ay[0], az[0]]).T;
# array of states
X0 = [];
X1 = [];
X2 = [];
X3 = [];
X4 = [];
X5 = [];
X6 = [];
X7 = [];
X8 = [];
spoofx = [];
spoofy = [];
Pxx = [];

# decide if GPS reading is to be taken
diff = np.diff(np.floor(tsi));
takeGPS = np.hstack((1.0,diff));

# plotting GPS course
plt.plot(x,y,'ro');
plt.show();

# miscellaneous variables
outage = 0;
spoof = 0;

# iteration begins! Godspeed !

for i in range(0,len(ax)-1):

	U = np.matrix([yv[i], pv[i], rv[i], ax[i], ay[i], az[i]]).T;

	Ax = U[3] + 2*U[2]*X[5] - 2*U[0]*X[4];
	Ay = U[3] + 2*U[0]*X[3] - 2*U[1]*X[5];
	Az = U[4] + 2*U[1]*X[4] - 2*U[2]*X[3];

	X[0] = X[0] + X[3]*dti + Ax*0.5*dti**2;
	X[1] = X[1] + X[4]*dti + Ay*0.5*dti**2;
	X[2] = X[2] + X[5]*dti + Az*0.5*dti**2;
	X[3] = X[3] + Ax*dti;
	X[4] = X[4] + Ay*dti;
	X[5] = X[5] + Az*dti;
	X[6] = X[6] + U[0]*dti;
	if(X[6] > np.pi):
		X[6] = X[6] - 2*np.pi;
	elif(X[6] < -np.pi):
		X[6] = X[6] + 2*np.pi;
	X[7] = X[7] + U[1]*dti;
	if(X[7] > np.pi):
		X[7] = X[7] - 2*np.pi;
	elif(X[7] < -np.pi):
		X[7] = X[7] + 2*np.pi;
	X[8] = X[8] + U[2]*dti;
	if(X[8] > np.pi):
		X[8] = X[8] - 2*np.pi;
	elif(X[6] < -np.pi):
		X[8] = X[8] + 2*np.pi;
	
	# Jacobian 
	JA = np.matrix([[ 1, 0, 0, dti,   0,   0, 0, 0, 0],
					[ 0, 1, 0,   0, dti,   0, 0, 0, 0],
					[ 0, 0, 1,   0,   0, dti, 0, 0, 0],
					[ 0, 0, 0,   1,   0,   0, 0, 0, 0],
					[ 0, 0, 0,   0,   1,   0, 0, 0, 0],
					[ 0, 0, 0,   0,   0,   1, 0, 0, 0],
					[ 0, 0, 0,   0,   0,   0, 1, 0, 0],
					[ 0, 0, 0,   0,   0,   0, 0, 1, 0],
					[ 0, 0, 0,   0,   0,   0, 0, 0, 1]]);
	# error covariance prediction
	P = JA*P*JA.T + Q

	# measurement function
	hx = np.matrix([[float(X[0])],
					[float(X[1])],
					[float(X[2])],
					[float(X[3])],
					[float(X[4])],
					[float(X[5])]]);

	# GPS values to  be taken or not
	if(takeGPS[i] == 0):
		JH = np.matrix([[0,0,0,0,0,0,0,0,0],
						[0,0,0,0,0,0,0,0,0],
						[0,0,0,0,0,0,0,0,0],
						[0,0,0,0,0,0,0,0,0],
						[0,0,0,0,0,0,0,0,0],
						[0,0,0,0,0,0,0,0,0]]);
	else:
		JH = np.matrix([[1,0,0,0,0,0,0,0,0],
						[0,1,0,0,0,0,0,0,0],
						[0,0,1,0,0,0,0,0,0],
						[0,0,0,1,0,0,0,0,0],
						[0,0,0,0,1,0,0,0,0],
						[0,0,0,0,0,1,0,0,0]]);

	# calculating Kalman gain
	S = JH*P*JH.T + R;
	K = (P*JH.T) * np.linalg.inv(S);

	# updating estimate
	check = [tsg==np.floor(tsi[i])];
	if(len(x[check])!=0):
	
		x_ = float(x[check]);
		y_ = float(y[check]);
		z_ = float(alt[check]);
		if(outage == 1):
			vx_ = float(X[3]);
			vy_ = float(X[4]);
			vz_ = float(X[5]);
			outage = 0;
		else:
			vx_ = float(vx[check]);
			vy_ = float(vy[check]);
			vz_ = float(vz[check]);
		Z = np.matrix([[x_],
						[y_],
						[z_],
						[vx_],
						[vy_],
						[vz_]]);
		
		if(takeGPS[i] == 1):
			if(abs(float(X[0])-x_)>P[0,0] or abs(float(X[1])-y_)>P[1,1] or abs(float(X[2])-z_)>P[2,2]):
				spoofx.append(x_);
				spoofy.append(y_);
				spoof = 1;
				continue;

		# kalman innovation
		Y = Z - hx;
		X = X + (K*Y);
		# update error covariance
		P = (I - K*JH)*P;
				
	else:
		outage = 1;



	X0.append(float(X[0]));
	X1.append(float(X[1]));
	X2.append(float(X[2]));
	X3.append(float(X[3]));
	X4.append(float(X[4]));
	X5.append(float(X[5]));
	X6.append(float(X[6]));
	X7.append(float(X[7]));
	X8.append(float(X[8]));
	Pxx.append(P[0,0]);


plt.plot(X0,X1,'r-');
plt.plot(x,y,'go');
plt.plot(spoofx,spoofy,'yo');
plt.show();
plt.plot(Pxx,range(len(Pxx)),'r-');
plt.show();
print P;