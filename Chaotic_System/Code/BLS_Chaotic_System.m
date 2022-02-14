clear;
close all;

TT = 7;%total time
ts = 0.001;%sampling time
error = 0;%error vector
error_old = 0;%old error vector
derror = 0;%derivative of error
ss = 0;%smc controller

%%%%%% BLS parameters %%%%%%%
intput_num = 1;%input dimension
output_num = 1;%output dimension
fea_numi = 2;%the number of feature group
fea_numk = 3;%the number of feature neuron in each feature group
en_num = 1;%the number of enhancement neuron

%%%%%% initial parameters of BLS %%%%%%%
wfkli = unifrnd(-1,1,[fea_numi*fea_numk,intput_num]);
bfki = unifrnd(-1,1,[fea_numi*fea_numk,1]);
wjki = unifrnd(-1,1,[en_num,fea_numi*fea_numk]);
bj = unifrnd(-1,1,[en_num,output_num]);
wki = unifrnd(-1,1,[fea_numi*fea_numk,output_num]);
wj = unifrnd(-1,1,[en_num,output_num]);

%%%%%% learning rate %%%%%%%
ewfkli = 0.01*eye(intput_num*fea_numi*fea_numk);
ebfki = 0.01*eye(fea_numi*fea_numk);
ewjki = 0.01*eye(en_num*fea_numi*fea_numk);
ebj = 0.01*eye(en_num*output_num);
ewki = 0.01*eye(fea_numi*fea_numk);
ewj = 0.01*eye(en_num);

%%%%%% plot parameters %%%%%%%
E1 = [];
D_E1 = [];
X1 = [];
X_D1 = [];
XD1 = [];
XD_D1 =[];
U1 = [];

x = 0.2;%initial system state
x_d = 0.2;%initial system state
RR = 0.07;%attenuation coefficient

uBLS = zeros(output_num,1);%the output of BLS

%%%%%% main program %%%%%%%
for index = 0:1:TT/ts
    t = index*ts;
    xd = sin(1.1*t);
    xd_d = 1.1*cos(1.1*t);
    error = xd-x;%tracking error
    if index ~= 0
        derror = (error-error_old)/ts;
    end
    ss = 0.5*derror+10*error;
    
    %%%%%% count the output of BLS %%%%%%%
    Z = zeros(fea_numi*fea_numk,1);
    H = zeros(en_num,1);
    
    Z = wfkli*ss+bfki;
    for i = 1:fea_numi*fea_numk
        Z(i) = (exp(Z(i))-exp(-Z(i)))/(exp(Z(i))+exp(-Z(i)));
    end
    H = wjki*Z+bj;
    H = (exp(H)-exp(-H))/(exp(H)+exp(-H));
    uBLS = wki'*Z+wj'*H;
    
    uRC = ((RR^2+1)/(2*RR^2))*ss;%the output robust controller
    U_t = uRC +uBLS;%output torque
    d = sin(2.*t)+cos(2.*t);%disturbance
    % Duffing-Holmes chaotic system
    xd_dd = -0.25 * xd_d + xd - xd^3 + 0.1 * sqrt(xd^2 + xd_d^2) * sin(t) + 0.3 * cos(t);
    x_dd = -0.25 * x_d + x - x^3 + 0.1 * sqrt(x^2 + x_d^2) * sin(t) + 0.3 * cos(t) + 0.1*x + d + U_t;
    xd_d = xd_d + xd_dd * ts;
    xd = xd + xd_d * ts;
    x_d = x_d + x_dd * ts;
    x = x + x_d *ts;
    
    %%%%%% BLS parameters updating %%%%%%%
    %%%%%% wki updating %%%%%%%
    wki = wki-ewki*Z*ss';
    
    %%%%%% wj updating %%%%%%%
    wj = wj-ewj*H*ss';
    
    %%%%%% wfkli updating %%%%%%%
    wfkli_d = zeros(fea_numi*fea_numk,intput_num*fea_numi*fea_numk);
        for j = 1:fea_numi*fea_numk
            for k = 1:intput_num
                wfkli_d(j,(j-1)*intput_num+k) = (1-Z(j).^2).*ss(k);
            end
        end
    wfkli = wfkli-reshape(ewfkli*wfkli_d'*wki*ss,fea_numi*fea_numk,intput_num);
    
    %%%%%% bfki updating %%%%%%%
    bfki_d1 = zeros(1,fea_numi*fea_numk);
    bfki_d = zeros(fea_numi*fea_numk,fea_numi*fea_numk);
    for j = 1:fea_numi*fea_numk
        bfki_d1(j)=(1-Z(j).^2);
    end
    for j = 1:fea_numi*fea_numk
        bfki_d(j,:) = bfki_d1;
    end
    bfki = bfki-reshape(ebfki*bfki_d'*wki*ss,fea_numi*fea_numk,1);
    
    %%%%%% wjki updating %%%%%%%
    wjki_d = zeros(en_num,en_num*fea_numi*fea_numk);
    for j = 1:en_num
        for k = 1:fea_numi*fea_numk
            wjki_d(j,(j-1)*fea_numi*fea_numk+k) = (1-H(j).^2).*Z(k);
        end
    end
    wjki = wjki-reshape(ewjki*wjki_d'*wj*ss,en_num,fea_numi*fea_numk);
    
    %%%%%% bj updating %%%%%%%
    bj_d1 = zeros(1,en_num*output_num);
    bj_d = zeros(en_num,en_num*output_num);
    for j = 1:en_num*output_num
        bj_d1(j)=(1-H(j).^2);
    end
    for j = 1:en_num
        bj_d(j,:) = bj_d1;
    end
    bj = bj-reshape(ebj*bj_d'*wj*ss,en_num,output_num);
    
    %%%%%% plot %%%%%%%
    E1 = [E1 error];
    D_E1 = [D_E1 derror];
    X1 = [X1 x];
    X_D1 = [X_D1 x_d];
    XD1 = [XD1 xd];
    XD_D1 = [XD_D1 xd_d];
    U1 = [U1 U_t];
    
    error_old = error;
end

fprintf('The RMSE is %d\n',sqrt(sum((E1).^2)/(TT/ts)));

T = 0:0.001:7;
figure(1)
plot(XD1,XD_D1,'r',X1,X_D1,'g')
xlabel('x');
ylabel('x-dot');
title('Trajectory phase portrait of Duffing-Holmes Chaotic System');
grid on;
figure(2)
subplot(2,1,1)
plot(T,X1,'r',T,XD1,'g')
xlabel('Time(s)');
grid on;
legend(': x',': xd ');
subplot(2,1,2)
plot(T,X_D1,'r',T,XD_D1,'g')
xlabel('Time(s)');
grid on;
legend(': x-dot ',': xd-dot');
figure(3)
plot(T,E1,'r',T,D_E1,'g')
xlabel('Time(s)');
ylabel('Error');
grid on;
title('Tracking error');
legend(': e ',': e-dot');