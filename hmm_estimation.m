clear all; clc;

%% Read in the data
fileID = fopen('diff.txt','r');
data = fscanf(fileID,'%f %f',[ 2 Inf] );
data = data';
fclose(fileID);

%% Apply HMM
seq = data(:,1);
states = data(:,2);
states = states+1;
seq = uint32(seq*10)+1;

M = size(seq);
M = uint64(M(1)*6/10);

PseudoT = [10 9000;10 10];
%[T,E] = hmmestimate(seq(1:M),states(1:M));
[T,E] = hmmestimate(seq(1:M),states(1:M),'PseudoTransitions',PseudoT);

maximum = max(seq(1:M));
index_exceeding = find(seq>maximum);
seq(index_exceeding) = maximum;

likely_states = hmmviterbi(seq(M+1:end),T,E);

correct = find(likely_states == states(M+1:end)');

true_positive = find(likely_states' == 2 & states(M+1:end)==2);
false_positive = find(likely_states' == 2 & states(M+1:end)==1);
true_negative = find(likely_states' == 1 & states(M+1:end)==1);
false_negative = find(likely_states' == 1 & states(M+1:end)==2);