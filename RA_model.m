clear all
close all

% experimental condition 
Condition_antenna = [1 , 0 , 1, 1 , 0]; % 1/0 : Antennae intact or blocked
Condition_Light = ['A' ; 'B'; 'B'; 'O'; 'O']; % A/B/Off : Light from Above, from Below or Off (in the dark)

% simulation parameters 
End_sim_time=0.14; % in sec
sampleTime = 1e-5; 
time=0:sampleTime:End_sim_time;
N = 100; % neurons
Roll(:,1) = 0:360/N:360-360/N; 

% Ring attractor network
wEEk = 45.0/N; wEI = 60.0/N;
wIE = -6.0; wII = -1.0;    
gammaE = -1.5; gammaI = -7.5; 
tauE = 42.5e-3; % in sec 0.005*8.5; %7.5; %0.005;
tauI = 2.125e-3; % in sec 0.00025*8.5; %7.5; %0.00025;
wEE = zeros(N, N);
sigma = 120.0;
 for i=1:N
     for j=1:N
         tmp = abs(Roll(i) - Roll(j));
         diff = min(tmp, 360 - tmp);
         wEE(i,j) = exp(-diff^2/(2 * sigma^2));
     end
 end
wEE = wEE * wEEk;

%%%%%%%%%%%%%%%%%%%%%%%%
% Run over the different experimental conditions
%%%%%%%%%%%%%%%%%%%%%%%%
% experiment = 1 : P+, A+; Vt 
% experiment = 2 : P+, A-; Vb 
% experiment = 3 : P+, A+; Vb
% experiment = 4 : P+, A+; Vdark
% experiment = 5 : P+, A-; Vdark
for experiment = 1:5
    
Antenna = Condition_antenna(experiment);  
Light = Condition_Light(experiment);  

% Visual cue
if Light == 'A' 
    miu_V = 0; 
else    
    miu_V = 180; % conflicting information (light source at the bottom)
end
k_V = 40;
sigma_V = 10;
tmp = abs(Roll - miu_V);
diff = min(abs(tmp), 360 - abs(tmp));
X_V = k_V * exp(-diff.^2 ./ (2 * sigma_V^2)) ./ (sqrt(2 * pi) * sigma_V);

% Antennal cue  
miu_A = 0; 
k_A = 40; 
sigma_A = 5; 
tmp = abs(Roll - miu_A);
diff = min(abs(tmp), 360 - abs(tmp));
X_A = k_A * exp(-diff.^2 ./ (2 * sigma_A^2)) ./ (sqrt(2 * pi) * sigma_A);

% sensory input of the sigma-pi units
weight =100;
if Light == 'A' & Antenna == 1 % 2 inputs
    input = X_V + X_A + weight * X_V .* X_A;
end
if Light == 'B' & Antenna == 0 % 1 input
    input = X_V;
end
if Light == 'B' & Antenna == 1 % 2 inputs
    input = X_V + X_A + weight * X_V .* X_A;
end
if Light == 'O' & Antenna == 1 % 1 input
    input = X_A;
end
if Light == 'O' & Antenna == 0 % 0 input
    input = zeros(N,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%
% Ring attractor network
%%%%%%%%%%%%%%%%%%%%%%%%

% initialize to a bump shape at zero deg (mimicking prior knowledge from leg proprioception)
c = []; c(:, 1) =  wEE(:,1); 
u = []; u(1) = 1.0;

% run integration
j = 2;
for t=sampleTime:sampleTime:End_sim_time
        
    c(:, j) = c(:,j-1) + ( -c(:,j-1) + ...
        max(0, gammaE + wEE*c(:,j-1)+ wIE*u(j-1)+ input ) )*sampleTime/tauE; 
    u(j) = u(j-1) + ( -u(j-1) + max(0, gammaI+wEI*sum(c(:,j-1))+wII*u(j-1)) )*sampleTime/tauI;
    if mod(j,100) == 0
        figure(1)
        hold off, plot(Roll, input,'k')
        hold on, plot(Roll, c(:,j),'or')
        plot(Roll, c(:,j),'r')
        drawnow
    end
    j = j+1;
end

figure(2)
ci = [c; c(1,:)];
Rolli = [Roll; 360]';
polarPcolor(0:sampleTime:End_sim_time,Rolli, log(ci), 'labelR','Time (sec)', 'Nspokes',5)
shading('interp')
colorbar off

% Ring output over time 
Nt = length(c(1,:));
winner = [];
output = [];
for j=1:Nt
    output(j,1) = max(c(:,j));
    ind = find(c(:,j)==output(j,1));
    winner(j,1) = Roll(ind);    
end
figure(3)
plot(time,winner)
xlabel('Time (in sec)')
ylabel('WTA (Roll in deg)')
title('winning neuron')
axis([0 0.14 -180 180])
figure(4)
plot(time,output)
xlabel('Time (sec)')
ylabel('K (a.u.)')
title('winner activation')
drawnow


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Closed-loop control of body roll 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

roll_init = [156.9416, 164.9175, 169.0660, 170.7711, 177.7980];
data = roll_init(experiment)-winner; 
Kval = output*30; 

simin.time = time;
simin.signals.values = data;
siminK.time = time;
siminK.signals.values = Kval;

H=tf(1,[.5e-3 1]);
Hz=c2d(H,sampleTime,'tustin');

SimOut = sim('model_mouche_integrateur_V2');

colorCurve = ['b'; 'r' ; 'k'; 'y'; 'g'];
figure(5);
hold on, plot(time,roll_init(experiment)-Theta_roll, [ colorCurve(experiment) '--']);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparison with experimental data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Compare with experimental data...\n')
data_tmp=xlsread('Anna_data.xlsx');
condition = data_tmp(:,1);
t = data_tmp(:,2);
d = data_tmp(:,3);
figure(5)
ind1 = find(condition==1); % condition : Ltop_Aplus
ind2 = find(condition==2); % condition 'Lbot_Amoins'
ind3 = find(condition==3); % condition 'Lbot_Aplus'
ind4 = find(condition==4); % condition 'Dark_Aplus'
ind5 = find(condition==5); % condition 'Dark_Amoins'

hold on, plot(t(ind1), d(ind1),'b')
hold on, plot(t(ind2), d(ind2),'r')
hold on, plot(t(ind3), d(ind3),'k')
hold on, plot(t(ind4), d(ind4),'y')
hold on, plot(t(ind5), d(ind5),'g')

axis([0 0.10 -100 200])
xlabel('Time (s)')
ylabel('Roll angle (°)')
legend('M Ltop Aplus', 'M Lbot Amoins', 'M Lbot Aplus', 'M Dark Aplus', 'M Dark Amoins', 'M Ltop Aplus', 'M Lbot Amoins', 'Lbot Aplus', 'Dark Aplus', 'Dark Amoins')

