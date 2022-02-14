clear;
close all;

TT = 7;%total time
ts = 0.001;%sampling time
error = 0;%error vector
error_old = 0;%old error vector
derror = 0;%derivative of error
ss = 0;%smc controller

%%%%%% SOBLS parameters %%%%%%%
intput_num = 1;%input dimension
output_num = 1;%output dimension
fea_numi = 2;%the number of feature group
fea_numk_max = 10;%the maximum number of feature neuron in each feature group
fea_numk = randi([1 10],fea_numi,1);%random the number of feature neuron in each feature group
en_num = 1;%the number of enhancement neuron

%%%%%% initial parameters of SOBLS %%%%%%%
wfkli = unifrnd(-1,1,[fea_numi*fea_numk_max,intput_num]);
bfki = unifrnd(-1,1,[fea_numi*fea_numk_max,1]);
wjki = unifrnd(-1,1,[en_num,fea_numi*fea_numk_max]);
bj = unifrnd(-1,1,[en_num,output_num]);
wki = unifrnd(-1,1,[fea_numi*fea_numk_max,output_num]);
wj = unifrnd(-1,1,[en_num,output_num]);

%%%%%% learning rate %%%%%%%
ewfkli = 0.01*eye(intput_num*fea_numi*fea_numk_max);
ebfki = 0.01*eye(fea_numi*fea_numk_max);
ewjki = 0.01*eye(en_num*fea_numi*fea_numk_max);
ebj = 0.01*eye(en_num*output_num);
ewki = 0.01*eye(fea_numi*fea_numk_max);
ewj = 0.01*eye(en_num);

ddd = 0;
dgg = 0;
dd = 0;%deleting threshold
dg = 0;%generating threshold
Pwc = 0.1;%deleting threshold coefficient
Pbc = 1;%generating threshold coefficient

%%%%%% plot parameters %%%%%%%
E3 = [];
D_E3 = [];
X3 = [];
X_D3 = [];
XD3 = [];
XD_D3 = [];
U3 = [];
NUMK3 = [];

x = 0.2;%initial system state
x_d = 0.2;%initial system state
RR = 0.07;%attenuation coefficient

%%%%%% main program %%%%%%%
for index = 0:1:TT/ts
    NUMK3 = [NUMK3 sum(fea_numk)];
    
    t = index*ts;
    xd = sin(1.1*t);
    xd_d = 1.1*cos(1.1*t);
    error = xd-x;%tracking error
    if index ~= 0
        derror = (error-error_old)/ts;
    end
    ss = 0.5*derror+10*error;
    
    %%%%%% count the output of SOBLS %%%%%%%
    Z1 = zeros(fea_numi,fea_numk_max);
    Z2 = zeros(fea_numi,fea_numk_max);
    Z = zeros(fea_numi,fea_numk_max);
    H1 = zeros(en_num,1);
    H = zeros(en_num,1);
    for i = 1:fea_numi
        for k = 1:fea_numk(i)
            Z1(i,k) = wfkli((i-1)*fea_numk_max+k)*ss+bfki((i-1)*fea_numk_max+k);
            Z2(i,k) = (exp(Z1(i,k))-exp(-Z1(i,k)))./(exp(Z1(i,k))+exp(-Z1(i,k)));
        end
    end
    for j = 1:en_num
        for i=1:fea_numi
            for k = 1:fea_numk(i)
                H1(j) = H1(j)+wjki(j,(i-1)*fea_numk_max+k)*Z2(i,k)+bj(j);
            end
        end
        H(j) = (exp(H1(j))-exp(-H1(j)))./(exp(H1(j))+exp(-H1(j)));
    end
    Z = zeros(fea_numi*fea_numk_max,1);
    for i = 1:fea_numi
        for k = 1:fea_numk_max
            Z((i-1)*fea_numk_max+k) = Z2(i,k);
        end
    end
    uSOBLS = wki'*Z+wj'*H;
    
    uRC = ((RR^2+1)/(2*RR^2))*ss;%the output robust controller
    U_t = uRC +uSOBLS;%output torque
    d = sin(2.*t)+cos(2.*t);%disturbance
    %Duffing-Holmes chaotic system
    xd_dd = -0.25 * xd_d + xd - xd^3 + 0.1 * sqrt(xd^2 + xd_d^2) * sin(t) + 0.3 * cos(t);
    x_dd = -0.25 * x_d + x - x^3 + 0.1 * sqrt(x^2 + x_d^2) * sin(t) + 0.3 * cos(t) + 0.1*x + d + U_t;
    xd_d = xd_d + xd_dd * ts;
    xd = xd + xd_d * ts;
    x_d = x_d + x_dd * ts;
    x = x + x_d *ts;
    
    %%%%%% TOPSIS %%%%%%%
    for i = 1:fea_numi
        for k = 1:fea_numk_max
            if sum(Z2(:,k))~=0
                ZZZ(i,k) = Z2(i,k)/sqrt(sum(Z2(:,k).^2));
            else
                ZZZ(i,k)=0;
            end
        end
    end
    
    for k = 1:fea_numk_max
        ma = ZZZ(:,k);
        summ = 0;
        for idx = 1:length(ma)
            if ma(idx) > 0
                summ = summ+ma(idx)*log(ma(idx));
            end
        end
        E(k) = -1/log(fea_numi)*summ;
        D(k) = 1-E(k);
    end
    
    for k = 1:fea_numk_max
        fai(k)=D(k)./sum(D);
    end
    
    for i = 1:fea_numi
        for k = 1:fea_numk_max
            v(i,k)=ZZZ(i,k).*fai(k);
        end
    end
    
    for k = 1:fea_numk_max
        vmax(k) = max(v(:,k));
        vmin(k) = min(v(:,k));
    end
    
    for i = 1:fea_numi
        smax(i) = sqrt(sum((v(i,:)'-vmax(:)).^2));
        smin(i) = sqrt(sum((v(i,:)'-vmin(:)).^2));
        wc(i) = smax(i)/(smax(i)+smin(i));
        bc(i) = smin(i)/(smax(i)+smin(i));
    end
    
    maa = wc(:);
    maxx = -inf;
    for idx = 1:length(maa)
        if maa(idx) < 1 && maa(idx) > 0.00001
            if maa(idx) > maxx
                maxx = maa(idx);
            end
        end
    end
    if t > 0
        if ddd > maxx && maxx~=-inf
            ddd = maxx;
        end
    else
        ddd = maxx;
    end
    dd = Pwc*ddd;
    
    for i=1:fea_numi
        for k = 1:fea_numk(i)
            if Z2(i,k)<dd
                if fea_numk(i)>1
                    fea_numk(i) = fea_numk(i)-1;
                    for k1 = k:fea_numk(i)
                        for l = 1:intput_num
                            wfkli((i-1)*fea_numk_max+k1,l) = wfkli((i-1)*fea_numk_max+k1+1,l);
                        end
                        bfki((i-1)*fea_numk_max+k1) = bfki((i-1)*fea_numk_max+k1+1);
                        for j = 1:en_num
                            wjki(j,(i-1)*fea_numk_max+k1) = wjki(j,(i-1)*fea_numk_max+k1+1);
                        end
                    end
                end
            end
        end
    end
    
    maa = bc(:);
    minn = inf;
    for idx = 1:length(maa)
        if maa(idx) > 0.00001 && maa(idx) < 1
            if maa(idx) < minn
                minn = maa(idx);
            end
        end
    end
    if t>0
        if dgg < minn && minn ~= inf
            dgg = minn;
        end
    else
        dgg = minn;
    end
    dg = Pbc*dgg;
    Z_max = max(Z);
    if Z_max<dg
        for i = 1:fea_numi
            if fea_numk(i)<fea_numk_max
                fea_numk(i) = fea_numk(i)+1;
            end
        end
    end
    
    %%%%%% SOBLS parameters updating %%%%%%%
    %%%%%% wki updating %%%%%%%
    wki = wki-ewki*Z*ss';
    
    %%%%%% wj updating %%%%%%%
    wj = wj-ewj*H*ss';
    
    %%%%%% wfkli updating %%%%%%%
    wfkli_d = zeros(fea_numi*fea_numk_max,intput_num*fea_numi*fea_numk_max);
    for j = 1:fea_numi*fea_numk_max
        for k = 1:intput_num
            wfkli_d(j,(j-1)*intput_num+k) = (1-Z(j).^2).*ss(k);
        end
    end
    wfkli = wfkli-reshape(ewfkli*wfkli_d'*wki*ss,fea_numi*fea_numk_max,intput_num);
    
    %%%%%% bfki updating %%%%%%%
    bfki_d1 = zeros(1,fea_numi*fea_numk_max);
    bfki_d = zeros(fea_numi*fea_numk_max,fea_numi*fea_numk_max);
    for j = 1:fea_numi*fea_numk_max
        bfki_d1(j)=(1-Z(j).^2);
    end
    for j=1:fea_numi*fea_numk_max
        bfki_d(j,:) = bfki_d1;
    end
    bfki = bfki-reshape(ebfki*bfki_d'*wki*ss,fea_numi*fea_numk_max,1);
    
    %%%%%% wjki updating %%%%%%%
    wjki_d = zeros(en_num,en_num*fea_numi*fea_numk_max);
    for j = 1:en_num
        for k = 1:fea_numi*fea_numk_max
            wjki_d(j,(j-1)*fea_numi*fea_numk_max+k) = (1-H(j).^2).*Z(k);
        end
    end
    wjki = wjki-reshape(ewjki*wjki_d'*wj*ss,en_num,fea_numi*fea_numk_max);
    
    %%%%%% bj updating %%%%%%%
    bj_d1 = zeros(1,en_num*output_num);
    bj_d = zeros(en_num,en_num*output_num);
    for j = 1:en_num*output_num
        bj_d1(j)=(1-H(j).^2);
    end
    for j=1:en_num
        bj_d(j,:) = bj_d1;
    end
    bj = bj-reshape(ebj*bj_d'*wj*ss,en_num,output_num);
    
    %%%%%% plot %%%%%%%
    E3 = [E3 error];
    D_E3 = [D_E3 derror];
    X3 = [X3 x];
    X_D3 = [X_D3 x_d];
    XD3 = [XD3 xd];
    XD_D3 = [XD_D3 xd_d];
    U3 = [U3 U_t];
    
    error_old = error;
end

fprintf('The RMSE is %d\n',sqrt(sum((E3).^2)/(TT/ts)));

T = 0:0.001:7;
figure(1)
plot(XD3,XD_D3,'r',X3,X_D3,'g')
xlabel('x');
ylabel('x-dot');
title('Trajectory phase portrait of Duffing-Holmes Chaotic System');
grid on;
figure(2)
subplot(2,1,1)
plot(T,X3,'r',T,XD3,'g')
xlabel('Time(s)');
grid on;
legend(': x',': xd ');
subplot(2,1,2)
plot(T,X_D3,'r',T,XD_D3,'g')
xlabel('Time(s)');
grid on;
legend(': x-dot ',': xd-dot');
figure(3)
plot(T,E3,'r',T,D_E3,'g')
xlabel('Time(s)');
ylabel('Error');
grid on;
title('Tracking error');
legend(': e ',': e-dot');
figure(4)
plot(T,NUMK3);
axis([-1 8 0 21]);
xlabel('Time(s)');
ylabel('Number of the remaining feature node');