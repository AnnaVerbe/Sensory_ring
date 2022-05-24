clear all
close all

folder_name0 = uigetdir(cd,'Fichier ring attractor');
cd(folder_name0) 

Condition_antenna = [1 , 0 , 1, 1 , 0];
Condition_Light = ['A' ; 'B'; 'B'; 'O'; 'O'];

% 1 = P+, La ; A+ 
% 2 = P+, Lb ; A- 
% 3 = P+, Lb ; A+
% 4 = P+, Ldark ; A+
% 5 = P+, Ldark ; A-

for iiii = 1:5%:5
% Condition expérimentale
Antenna = Condition_antenna(iiii);
Light = Condition_Light(iiii);

% Antenna = 1/0 : presence ou absence des antennes
% Light = A/B/Off : Light Above or Below or Off (pas de lumière)

% parameters of simulation
End_sim_time=0.12;%0.14;%0.1; % in sec
sampleTime = 1e-5; 
time=0:sampleTime:End_sim_time;
%T = 0.1; % in sec
%dt = 1e-5; 
N = 100; % Number of neurons
Roll = 0:360/N:360-360/N; 


% CUES -------------------------

% En fait on utilise pas la cue 1
% On utilise la proprio uniquement pour initialiser avec un bump à zero deg
% % Cue 1: proprioception
% miu1 = 0;
% k1 = 40; 
% sigma1 = 10;

% Cue 2: Visual
if Light == 'A' 
    miu2 = 0; 
else    
    miu2 = 180; % conflicting information (light source at the bottom)
end
%k2 = 45;
%sigma2 = 4;
k2 = 40;
sigma2 = 10;

% Cue 3: Antenna  
miu3 = 0; 
%k3 = 80;
%sigma3 = 1;
k3 = 40;%100; %40
sigma3 = 5;%10; %5 ou 10

% generate cue 1  : En fait on utilise pas la cue 1
% On utilise la proprio uniquement pour initialiser avec un bump à zero deg  
% tmp = abs(Roll - miu1);
% diff = min(abs(tmp), 360 - abs(tmp));
% c1 = k1 * exp(-diff.^2 ./ (2 * sigma1^2)) ./ (sqrt(2 * pi) * sigma1);
% j=1;
% for t=0:sampleTime:End_sim_time
%     x1(:,j) = c1;
%     j = j+1;
% end
%     
% generate cue 2 -- Visual
tmp = abs(Roll - miu2);
diff = min(abs(tmp), 360 - abs(tmp));
c2 = k2 * exp(-diff.^2 ./ (2 * sigma2^2)) ./ (sqrt(2 * pi) * sigma2);
j=1;
for t=0:sampleTime:End_sim_time
    x2(:,j) = c2;
    j = j+1;
end

% %sans temps 
% % generate cue 3 
tmp = abs(Roll - miu3);
diff = min(abs(tmp), 360 - abs(tmp));
c3 = k3 * exp(-diff.^2 ./ (2 * sigma3^2)) ./ (sqrt(2 * pi) * sigma3);
j=1;
for t=0:sampleTime:End_sim_time
    x3(:,j) = c3;
    j = j+1;
end



% % avec delay
% generate cue 3
% antenna_delay = 0.04; % in sec
% tmp = abs(Roll - miu3);
% diff = min(abs(tmp), 360 - abs(tmp));
% c3 = k3 * exp(-diff.^2 ./ (2 * sigma3^2)) ./ (sqrt(2 * pi) * sigma3);
% j=1;
% for t=0:sampleTime:End_sim_time
%     x3(:,j) = zeros(size(c3));
%     if t > antenna_delay, x3(:,j) = c3; end
%     j = j+1;
% end




% Ring attractor network
wEEk = 45.0/N; wIE = 60.0/N;
wEI = -6.0; wII = -1.0;    
gammaE = -1.5; gammaI = -7.5; 

tauE = 0.005*8.5; %7.5; %0.005;
tauI = 0.00025*8.5; %7.5; %0.00025;

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
 
% Initialization 
c(:, 1) =  wEE(:,1); % initialize to a bump shape à zero deg (Proprioception)
%u(1) = 0.69;
u(1) = 1.0;

%%

        writerObj = VideoWriter(['myVideo_'  num2str(iiii) '.avi'],'Uncompressed AVI');
        writerObj.FrameRate = 32;

        open(writerObj);
 
   
        
% run integration
weight =100;
j = 2;
for t=sampleTime:sampleTime:End_sim_time
    
    if Light == 'A' & Antenna == 1 % 2 inputs
        input = x2(:,j-1) + x3(:,j-1) + weight * x2(:,j-1).*x3(:,j-1);
    end
    if Light == 'B' & Antenna == 0 % 1 input
        input = x2(:,j-1);
    end
    if Light == 'B' & Antenna == 1 % 2 inputs
        input = x2(:,j-1) + x3(:,j-1) + weight * x2(:,j-1).*x3(:,j-1);
    end
    if Light == 'O' & Antenna == 1 % 1 input
        input = x3(:,j-1);
    end
    if Light == 'O' & Antenna == 0 % 0 input
        input = zeros(N,1);
    end
    
    c(:, j) = c(:,j-1) + ( -c(:,j-1) + ...
        max(0, gammaE + wEE*c(:,j-1)+ wEI*u(1,j-1)+ input ) )*sampleTime/tauE; 
    u(1,j) = u(1,j-1) + ( -u(1,j-1) + max(0, gammaI+wIE*sum(c(:,j-1))+wII*u(1,j-1)) )*sampleTime/tauI;
   
    if mod(j,100) == 0
        fh = figure(1);

        set(gcf, 'Color', 'w')
        
        hold off, 
        plot(Roll, input,'k','LineWidth',5)
        hold on, plot(Roll, c(:,j),'or','markersize',20)
        plot(Roll, c(:,j),'r','LineWidth',5)
        xlim([0 360]);
        xticks([0 180 360])
        fh.WindowState = 'maximized';
        %t
       % drawnow
       
      drawnow
      set(gca,'FontSize',20)
      set(gca,'LineWidth',3)

      F = getframe(gcf) ;
      frame = F;    
      writeVideo(writerObj, frame);
      
      
      Roll_Tot(:,iiii) = Roll;
      Input_Tot(:,iiii) = input;
      C_Tot(j,iiii) = {c(:,j)};
      
    end
    

    j = j+1;
end
        close(writerObj);

%%
% figure(1)
% hold on 
% for jj = [1,10,50,100,300,500,1000,length(c(1,:))]
% %for jj = [1:60:length(c(1,:))]
% plot(Roll, c(:,jj),'r')
% %plot(Roll, c(:,jj),'or')
% end



%%

figure(2)
hold on 
plot(Roll, input,'k')

figure(100)
hold on
plot(Roll,wEE(:,1))
plot(Roll, c2)
plot(Roll,c3)
xlabel('Roll (deg)')
xlim([0 360])
%ylabel('Proprioception cue (a.u.)')
saveas(gcf, fullfile([ folder_name0 '\sensors_input.pdf']))
  
%%
% figure
% hold on
% %plot(c1)
% plot((1:100),c2)
% plot((1:100),c3)
% plot((1:100),c(:,length(c(1,:))))
% hold off

% figure(10)
% hold on 
% ind=find(time>=0.03);
% subplot(5,2,2*iiii-1), plot(Roll, c(:,ind(1)),'or')
% subplot(5,2,2*iiii), plot(Roll, c(:,ind(end)),'or')
% if iiii==5, subplot(5,2,2*iiii-1), xlabel('t=0.03 sec'), subplot(5,2,2*iiii), xlabel('t=0.14 sec'), end
% drawnow 

% circular variand

c_norm = [];
Nt = length(c(1,:));
for j=1:Nt
    c_norm(:,j) = c(:,j)/sum(c(:,j));
    c_var(j) = 1 - norm( [sin(deg2rad(Roll))*c_norm(:,j), cos(deg2rad(Roll))*c_norm(:,j)] ) ;
end
var_ring(iiii) = c_var(end); 

%c_norm_end = c(:,end)/sum(c(:,end));
%c_mean = rad2deg( atan2(sin(deg2rad(Roll))*c_norm, cos(deg2rad(Roll))*c_norm) );
%c_var_end = 1 - norm( [sin(deg2rad(Roll))*c_norm_end, cos(deg2rad(Roll))*c_norm_end] ) ;
%%

figure(10)
hold on 
if iiii==1, subplot(3,2,1), plot(Roll, c(:,1),'r'), hold on, plot(Roll, c(:,1),'or'), xlim([0 360]), end
subplot(3,2,iiii+1), plot(Roll, c(:,end),'r'), hold on, plot(Roll, c(:,end),'or'), xlim([0 360])
saveas(gcf, fullfile([ folder_name0 '\sensory_curves.pdf']))
  
%%


figure(11)
colorCurve = ['b'; 'r' ; 'k'; 'y'; 'g'];

hold on, plot(0:sampleTime:End_sim_time, c_var, [ colorCurve(iiii)]);
xlabel('Time (sec)')
ylabel('Circular variance')
drawnow 

%%
% miuF = mean(c(:,length(c(1,:))));%revoir
% sigmaF = std(c(:,length(c(1,:))));
% y_normF = normpdf(1:100,miuF,sigmaF);
% 
% miu1 = mean(wEE(:,1));%revoir
% sigma1 = std(wEE(:,1));
% y_norm1 = normpdf(1:100,miu1,sigma1);

% y_norm3 = normpdf(1:100,miu3,sigma3);
% y_norm2 = normpdf(1:100,miu2,sigma2);
% figure
% hold on
% plot(1:N,y_norm1)
% plot(1:N,y_norm2)
% plot(1:N,y_norm3)
% plot(1:N,y_normF)

% Evolution au cours du temps des neurones du ring
% % % % figure(2)
% % % % for i = 1:N
% % % %     hold on, plot(time,c(i,:))
% % % % end
% % % % xlabel('Time (sec')
% % % % title('Time evolution of ring neurons')


% Decodage en prenant le max d'activation
Nt = length(c(1,:));
winner = [];
output = [];
for j=1:Nt
    output(j,1) = max(c(:,j));
    ind = find(c(:,j)==output(j,1));
    winner(j,1) = Roll(ind);    
end

%%

figure(3)
hold on
plot(time,winner)
xlabel('Time')
ylabel('Roll (winner encoding )')
title('winning neuron')

figure(4)
plot(time,output)
xlabel('Time')
ylabel('winner activation')
title('winner activation')

drawnow
%%




% On lance simulink -----------

roll_init = [156.9416, 164.9175, 169.0660, 170.7711, 177.7980];

%data = 170-winner; 
data = roll_init(iiii)-winner; 
Kval = output*30; 

% if miu2 == 0
%   Kval = output*35*100; %30*output;
% else
%     Kval = 35*output; %exp(8*output);
% end

% if miu2 == 0
%  Kval = output*35*100; %30*output;
% elseif Light == 'O'
%  Kval = output*35; %25*output;   
% else
%  Kval = output*35; %30*output;
% end


% figure
% plot(output)
% hold on 
% plot(output.*output)
% plot(log2(output))
% plot(log(output))

simin.time = time;
simin.signals.values = data;
siminK.time = time;
siminK.signals.values = Kval;

H=tf(1,[.5e-3 1]);
Hz=c2d(H,sampleTime,'tustin');

SimOut = sim('model_mouche_integrateur_V2');

colorCurve = ['b'; 'r' ; 'k'; 'y'; 'g'];


figure(5);
%hold on, plot(time,170-Theta_roll, [ colorCurve(iiii) '--']);
hold on, plot(time,roll_init(iiii)-Theta_roll, [ colorCurve(iiii) '--']);
timesaved(:,iiii) = time;
PlotSaved(:,iiii) = roll_init(iiii)-Theta_roll;

%%
figure(18)
hold on
plot(0:sampleTime:End_sim_time,integrator_output,[ colorCurve(iiii)])


%figure(6);
%plot(time,Omega_roll);
%figure(7);
%plot(time,integrator_output);

%MeanCorpsTOT(:,iiii) = 170-Theta_roll;
MeanCorpsTOT(:,iiii) = roll_init(iiii)-Theta_roll;
angularspeedTOT(:,iiii) = angularspeed;

%%
figure(203)
hold on
plot(0:sampleTime:End_sim_time,-angularspeedTOT,[ colorCurve(iiii)])
xlabel('Time (s)')
ylabel('vitesse angulaire (°/s)')


% siminK_Total(iiii,:) = siminK;
% simin_Total(iiii,:) = simin;

end

%%

for trial =1:5
    fprintf('trial %d\t Circular variance ring activity = %f\n', trial, var_ring(trial)) 
end


delete([ folder_name0  '\ModelTOT.xlsx']); 
xlswrite(fullfile([ folder_name0  '\ModelTOT.xlsx']), [time' MeanCorpsTOT] )

delete([ folder_name0  '\angularspeedTOT.xlsx']); 
xlswrite(fullfile([ folder_name0  '\angularspeedTOT.xlsx']), [time' -angularspeedTOT] )

%%
close all 
colorCurve = [[0,0,245/255];[255/255,0,0];[120/255,206/255,27/255] ;[102/255,102/255,102/255]; [253/255,172/255,16/255]];

for iiii = 1:5
writerObj = VideoWriter(['curve_'  num2str(iiii) '.avi'],'Uncompressed AVI');
writerObj.FrameRate = 32;

open(writerObj);

data_tmp=xlsread('Anna_data.xlsx');
condition = data_tmp(:,1);
t = data_tmp(:,2);
d = data_tmp(:,3);

fh = figure(5);
figure(5)

axis([0 0.12 -100 200])
xlabel('Time (s)')
ylabel('Roll angle (°)')
fh.WindowState = 'maximized';
set(gcf, 'Color', 'w')
set(gca,'FontSize',20)
set(gca,'LineWidth',3)


ind1 = find(condition==1); % condition : Ltop_Aplus
ind2 = find(condition==2); % condition 'Lbot_Amoins'
ind3 = find(condition==3); % condition 'Lbot_Aplus'
ind4 = find(condition==4); % condition 'Dark_Aplus'
ind5 = find(condition==5); % condition 'Dark_Amoins'

if iiii == 1
indC = ind1;    
elseif iiii == 2
indC = ind2;  
elseif iiii == 3
indC = ind3;  
elseif iiii == 4
indC = ind4;  
elseif iiii == 5
indC = ind5;  
end

time = timesaved(:,iiii);
valroll = PlotSaved(:,iiii);
for jj = 1:length(t(t(indC)<0.12))
hold on, plot(t(indC(1:jj)), d(indC(1:jj)),'color',colorCurve(iiii,:),'LineWidth',5)
cjj = max(find(time(time < t(indC(jj)))));
hold on, plot(time(1:cjj),valroll(1:cjj),'color',colorCurve(iiii,:),'linestyle',':','LineWidth',5);
% hl = legend('trials','model');
% set(hl, 'TextColor','k', 'Color','w', 'EdgeColor','k')
drawnow
F = getframe(gcf) ;
frame = F;    
writeVideo(writerObj, frame);
end
% hold on, plot(t(ind2), d(ind2),'r')
% hold on, plot(t(ind3), d(ind3),'k')
% hold on, plot(t(ind4), d(ind4),'y')
% hold on, plot(t(ind5), d(ind5),'g')




%legend('M Ltop Aplus', 'M Lbot Amoins', 'M Lbot Aplus', 'M Dark Aplus', 'M Dark Amoins', 'M Ltop Aplus', 'M Lbot Amoins', 'Lbot Aplus', 'Dark Aplus', 'Dark Amoins')

close(writerObj);

end

%%
% 
% 
% for iiii = 1:5
%     
% ind1 = find(condition==1); % condition : Ltop_Aplus
% ind2 = find(condition==2); % condition 'Lbot_Amoins'
% ind3 = find(condition==3); % condition 'Lbot_Aplus'
% ind4 = find(condition==4); % condition 'Dark_Aplus'
% ind5 = find(condition==5); % condition 'Dark_Amoins'
% 
% if iiii == 1
% indC = ind1;    
% elseif iiii == 2
% indC = ind2;  
% elseif iiii == 3
% indC = ind3;  
% elseif iiii == 4
% indC = ind4;  
% elseif iiii == 5
% indC = ind5;  
% end
%         
%  %   
% writerObj = VideoWriter(['myVideoV2_'  num2str(iiii) '.avi'],'Uncompressed AVI');
% writerObj.FrameRate = 32;
% 
% open(writerObj);
% 
%         
% % run integration
% weight =100;
% j = 2;
% for t=sampleTime:sampleTime:End_sim_time
%     
% Roll_Tot(:,iiii) = Roll;
% Input_Tot(:,iiii) = input;
% C_Tot(j,iiii) = {c(:,j)};
%       
%     if mod(j,100) == 0
%         fh = figure(1);
% 
%         set(gcf, 'Color', 'w')
%         
%         hold off, 
%         plot(Roll_Tot(:,iiii), Input_Tot(:,iiii),'k','LineWidth',5)
%         hold on, plot(Roll_Tot(:,iiii), C_Tot(j,iiii),'or','markersize',20)
%         plot(Roll_Tot(:,iiii), C_Tot(j,iiii),'r','LineWidth',5)
%         xlim([0 360]);
%         xticks([0 180 360])
%         fh.WindowState = 'maximized';
%         %t
%        % drawnow
%        
%       drawnow
%       set(gca,'FontSize',20)
%       set(gca,'LineWidth',3)
% 
%       F = getframe(gcf) ;
%       frame = F;    
%       writeVideo(writerObj, frame);
%     end
%     
% 
%     j = j+1;
% end
%         close(writerObj);
% end
% %%
% On plotte les données expérimentales
% 1 = P+, La ; A+ 
% 2 = P+, Lb ; A- 
% 3 = P+, Lb ; A+
% 4 = P+, Ldark ; A+
% 5 = P+, Ldark ; A-

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

return





