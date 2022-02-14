clear;
close all;

TT = 7;%total time
ts = 0.001;%sampling time
error = 0;%error vector
error_old = 0;%old error vector
derror = 0;%derivative of error
ss = 0;%smc controller

%%%%%% DLRBLS parameters %%%%%%%
intput_num = 1;%input dimension
output_num = 1;%output dimension
fea_numi = 2;%the number of feature group
fea_numk = 3;%the number of feature neuron in each feature group
en_num = 1;%the number of enhancement neuron

%%%%%% initial parameters of DLRBLS %%%%%%%
wfkli = unifrnd(-1,1,[fea_numi*fea_numk,intput_num]);
bfki = unifrnd(-1,1,[fea_numi*fea_numk,1]);
wjki = unifrnd(-1,1,[en_num,fea_numi*fea_numk]);
wki = unifrnd(-1,1,[fea_numi*fea_numk,output_num]);
wj = unifrnd(-1,1,[en_num,output_num]);
wro = unifrnd(-1,1,[intput_num,1]);
wr = unifrnd(-1,1,[en_num,1]);
bb = unifrnd(-1,1,[en_num,fea_numi*fea_numk]);
cc = unifrnd(-1,1,[en_num,fea_numi*fea_numk]);

%%%%%% learning rate %%%%%%%
ewfkli = 0.01*eye(intput_num*fea_numi*fea_numk);
ebfki = 0.01*eye(fea_numi*fea_numk);
ewjki = 0.01*eye(en_num*fea_numi*fea_numk);
ewki = 0.01*eye(fea_numi*fea_numk);
ewj = 0.01*eye(en_num);
ewro = 0.01*eye(intput_num);
ewr = 0.01*eye(en_num);
ebb = 0.01*eye(en_num*fea_numi*fea_numk);
ecc = 0.01*eye(en_num*fea_numi*fea_numk);

%%%%%% plot parameters %%%%%%%
E2 = [];
D_E2 = [];
X2 = [];
X_D2 = [];
XD2 = [];
XD_D2 = [];
U2 = [];

x = 0.2;%initial system state
x_d = 0.2;%initial system state
RR = 0.07;%attenuation coefficient

uDLRBLS = zeros(output_num,1);%the output of DLRBLS
uDLRBLS_old = zeros(output_num,1);
uDLRBLS_old_old = zeros(output_num,1);
H_old = zeros(output_num,en_num);

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
    
    %%%%%% count the output of DLRBLS %%%%%%%
    Z = zeros(fea_numi*fea_numk,1);
    H = zeros(en_num,1);
    if t>0.001
        Z = wfkli*(ss.*wro.*uDLRBLS_old./uDLRBLS_old_old)+bfki;
    else
        Z = wfkli*ss+bfki;
    end
    for i=1:fea_numi*fea_numk
        Z(i) = (exp(Z(i))-exp(-Z(i)))/(exp(Z(i))+exp(-Z(i)));
    end
    for j = 1:en_num
        for k = 1:fea_numi*fea_numk
            H(j) = H(j)+(wjki(j,k)*Z(k)+wr(j).*H_old(j)-cc(j,k)).^2./(bb(j,k).^2);
        end
        H(j) = exp(-H(j));
        H_old(j) = H(j);
    end
    uDLRBLS = wki'*Z+wj'*H;
    %%%%%% record the output of DLRBLS %%%%%%%
    if t>0
        uDLRBLS_old_old = uDLRBLS_old;
    end
    uDLRBLS_old = uDLRBLS;
    
    uRC = ((RR^2+1)/(2*RR^2))*ss;%the output robust controller
    U_t = uRC +uDLRBLS;%output torque
    d = sin(2.*t)+cos(2.*t);%disturbance
    
    %Duffing-Holmes chaotic system
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
    
     
    wfkli_d = zeros(fea_numi*fea_numk,intput_num*fea_numi*fea_numk);
    if t>0.001
        for j = 1:fea_numi*fea_numk
            for k = 1:intput_num
                wfkli_d(j,(j-1)*intput_num+k) = (1-Z(j).^2).*ss(k).*wro(k).*uDLRBLS_old(k)./uDLRBLS_old_old(k);
            end
        end
    else
        for j = 1:fea_numi*fea_numk
            for k = 1:intput_num
                wfkli_d(j,(j-1)*intput_num+k) = (1-Z(j).^2).*ss(k);
            end
        end
    end
    wfkli = wfkli-reshape(ewfkli*wfkli_d'*wki*ss,fea_numi*fea_numk,intput_num);
    
    %%%%%% bfki updating %%%%%%%
    bfki_d1 = zeros(1,fea_numi*fea_numk);
    bfki_d = zeros(fea_numi*fea_numk,fea_numi*fea_numk);
    for j = 1:fea_numi*fea_numk
        bfki_d1(j)=(1-Z(j).^2);
    end
    for j=1:fea_numi*fea_numk
        bfki_d(j,:) = bfki_d1;
    end
    bfki = bfki-reshape(ebfki*bfki_d'*wki*ss,fea_numi*fea_numk,1);
    
    %%%%%% wro updating %%%%%%%
    if t>0.001
        wro_d1 = zeros(1,intput_num);
        wro_d = zeros(fea_numi*fea_numk,intput_num);
        for j = 1:intput_num
            wro_d1(j) = (1-Z(j).^2).*sum(wfkli(:,j)).*(ss(j).*uDLRBLS_old(j)./uDLRBLS_old_old(j));
        end
        for j = 1:fea_numi*fea_numk
            wro_d(j,:) = wro_d1;
        end
        wro = wro-reshape(ewro*wro_d'*wki*ss,intput_num,1);
    end
    
    %%%%%% wjki updating %%%%%%%
    wjki_d = zeros(en_num,en_num*fea_numi*fea_numk);
    for j = 1:en_num
        for k = 1:fea_numi*fea_numk
            wjki_d(j,(j-1)*fea_numi*fea_numk+k) = -2.*H(j).*Z(k).*(wjki(j,k)*Z(k)+wr(j).*H_old(j)-cc(j,k))./(bb(j,k).^2);
        end
    end
    wjki = wjki-reshape(ewjki*wjki_d'*wj*ss,en_num,fea_numi*fea_numk);
    
    %%%%%% bb updating %%%%%%%
    bb_d = zeros(en_num,en_num*fea_numi*fea_numk);
    for j = 1:en_num
        for k = 1:fea_numi*fea_numk
            bb_d(j,(j-1)*fea_numi*fea_numk+k) = 2.*H(j).*(wjki(j,k)*Z(k)+wr(j).*H_old(j)-cc(j,k)).^2./(bb(j,k).^3);
        end
    end
    bb = bb-reshape(ebb*bb_d'*wj*ss,en_num,fea_numi*fea_numk);
    
    %%%%%% cc updating %%%%%%%
    cc_d = zeros(en_num,en_num*fea_numi*fea_numk);
    for j = 1:en_num
        for k = 1:fea_numi*fea_numk
            cc_d(j,(j-1)*fea_numi*fea_numk+k) = 2.*H(j).*(wjki(j,k)*Z(k)+wr(j).*H_old(j)-cc(j,k))./(bb(j,k).^2);
        end
    end
    cc = cc-reshape(ecc*cc_d'*wj*ss,en_num,fea_numi*fea_numk);
    
    %%%%%% wr updating %%%%%%%
    wr_d1 = zeros(1,en_num);
    wr_d = zeros(en_num,en_num);
    for j = 1:en_num
        wr_d1(j) = -2.*H(j).*H_old(j).*sum(sum((wjki(j,:).*Z(:)+wr(j).*H_old(j)-cc(j,:))./(bb(j,:).^2)));
    end
    for j=1:en_num
        wr_d(j,:) = wr_d1;
    end
    wr = wr-reshape(ewr*wr_d'*wj*ss,en_num,1);
    
    %%%%%% plot %%%%%%%
    E2=[E2 error];
    D_E2=[D_E2 derror];
    X2=[X2 x];
    X_D2=[X_D2 x_d];
    XD2=[XD2 xd];
    XD_D2=[XD_D2 xd_d];
    U2=[U2 U_t];
    
    error_old = error;
end

fprintf('The RMSE is %d\n',sqrt(sum((E2).^2)/(TT/ts)));

T=0:0.001:7;
figure(1)
plot(XD2,XD_D2,'r',X2,X_D2,'g')
xlabel('x');
ylabel('x-dot');
title('Trajectory phase portrait of Duffing-Holmes Chaotic System');
grid on;
figure(2)
subplot(2,1,1)
plot(T,X2,'r',T,XD2,'g')
xlabel('Time(s)');
grid on;
legend(': x',': xd ');
subplot(2,1,2)
plot(T,X_D2,'r',T,XD_D2,'g')
xlabel('Time(s)');
grid on;
legend(': x-dot ',': xd-dot');
figure(3)
plot(T,E2,'r',T,D_E2,'g')
xlabel('Time(s)');
ylabel('Error');
grid on;
title('Tracking error');
legend(': e ',': e-dot');