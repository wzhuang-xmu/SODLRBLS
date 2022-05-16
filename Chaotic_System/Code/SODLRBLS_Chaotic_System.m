clear;
close all;

TT = 7;%total time
ts = 0.001;%sampling time
error = 0;%error vector
error_old = 0;%old error vector
derror = 0;%derivative of error
ss = 0;%smc controller

%%%%%%%%%%%%%%%%%%SODLRBLS parameters%%%%%%%%%%%%%%%%
intput_num = 1;%input dimension
output_num = 1;%output dimension
fea_numi = 2;%the number of feature group
fea_numk_max = 10;%the maximum number of feature neuron in each feature group
fea_numk = randi([1 10],fea_numi,1);%random the number of feature neuron in each feature group
en_num = 1;%the number of enhancement neuron

%%%%%% SODLRBLS parameters %%%%%%%
wfkli = unifrnd(-1,1,[fea_numi*fea_numk_max,intput_num]);
bfki = unifrnd(-1,1,[fea_numi*fea_numk_max,1]);
wjki = unifrnd(-1,1,[en_num,fea_numi*fea_numk_max]);
wki = unifrnd(-1,1,[fea_numi*fea_numk_max,output_num]);
wj = unifrnd(-1,1,[en_num,output_num]);
wro = unifrnd(-1,1,[intput_num,1]);
wr = unifrnd(-1,1,[en_num,1]);
bb = unifrnd(-1,1,[en_num,fea_numi*fea_numk_max]);
cc = unifrnd(-1,1,[en_num,fea_numi*fea_numk_max]);

%%%%%% learning rate %%%%%%%
ewfkli = 0.01*eye(intput_num*fea_numi*fea_numk_max);
ebfki = 0.01*eye(fea_numi*fea_numk_max);
ewjki = 0.01*eye(en_num*fea_numi*fea_numk_max);
ewki = 0.01*eye(fea_numi*fea_numk_max);
ewj = 0.01*eye(en_num);
ewro = 0.01*eye(intput_num);
ewr = 0.01*eye(en_num);
ebb = 0.01*eye(en_num*fea_numi*fea_numk_max);
ecc = 0.01*eye(en_num*fea_numi*fea_numk_max);

ddd = 0;
dgg = 0;
dd = 0;%deleting threshold
dg = 0;%generating threshold
Pwc = 0.1;%deleting threshold coefficient
Pbc = 1;%generating threshold coefficient

%%%%%% plot parameters %%%%%%%
E4 = [];
D_E4 = [];
X4 = [];
X_D4 = [];
XD4 = [];
XD_D4 = [];
U4 =  [];
NUMK4 = [];

x = 0.2;%initial system state
x_d = 0.2;%initial system state
RR = 0.07;%attenuation coefficient

uSODLRBLS = zeros(output_num,1);%the output of SODLRBLS
uSODLRBLS_old = zeros(output_num,1);
uSODLRBLS_old_old = zeros(output_num,1);
H_old = zeros(output_num,en_num);

%%%%%% main program %%%%%%%
for index = 0:1:TT/ts
    NUMK4 = [NUMK4 sum(fea_numk)];
    
    t = index*ts;
    xd = sin(1.1*t);
    xd_d = 1.1*cos(1.1*t);
    error = xd-x;%tracking error
    if index ~= 0
        derror = (error-error_old)/ts;
    end
    ss = 0.5*derror+10*error;
    
    %%%%%% count the output of SODLRBLS %%%%%%%
    Z1 = zeros(fea_numi,fea_numk_max);
    Z2 = zeros(fea_numi,fea_numk_max);
    Z = zeros(fea_numi,fea_numk_max);
    H = zeros(en_num,1);
    for i = 1:fea_numi
        for k = 1:fea_numk(i)
            if t>0.001
                Z1(i,k) = wfkli((i-1)*fea_numk_max+k,:)*(ss.*wro.*uSODLRBLS_old./uSODLRBLS_old_old)+bfki((i-1)*fea_numk_max+k);
            else
                Z1(i,k) = wfkli((i-1)*fea_numk_max+k,:)*ss+bfki((i-1)*fea_numk_max+k);
            end
            Z2(i,k) = (exp(Z1(i,k))-exp(-Z1(i,k)))./(exp(Z1(i,k))+exp(-Z1(i,k)));
        end
    end
    Z = zeros(fea_numi*fea_numk_max,1);
    for i = 1:fea_numi
        for k = 1:fea_numk_max
            Z((i-1)*fea_numk_max+k) = Z2(i,k);
        end
    end
    for j = 1:en_num
        for i = 1:fea_numi
            for k = 1:fea_numk(i)
                H(j) = H(j)+(wjki(j,(i-1)*fea_numk_max+k)*Z2(i,k)+wr(j).*H_old(j)-cc(j,(i-1)*fea_numk_max+k)).^2./(bb(j,(i-1)*fea_numk_max+k).^2);
            end
        end
        H(j) = exp(-H(j));
        H_old(j) = H(j);
    end
    uSODLRBLS = wki'*Z+wj'*H;
    
    %%%%%% record the output of SODLRBLS %%%%%%%
    if t>0
        uSODLRBLS_old_old = uSODLRBLS_old;
    end
    uSODLRBLS_old = uSODLRBLS;
    
    uRC = ((RR^2+1)/(2*RR^2))*ss;%the output robust controller
    U_t = uRC +uSODLRBLS;%output torque
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
                ZZZ(i,k) = Z2(i,k)./sqrt(sum(Z2(:,k).^2));%进行规范化处理，值限定在[0,1]之间
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
    
    %%%%%% SODLRBLS parameters updating %%%%%%%
    %%%%%% wki updating %%%%%%%
    wki = wki-ewki*Z*ss';
    
    %%%%%% wj updating %%%%%%%
    wj = wj-ewj*H*ss';
    
    %%%%%% wfkli updating %%%%%%%
    wfkli_d = zeros(fea_numi*fea_numk_max,intput_num*fea_numi*fea_numk_max);
    if t>0.001
        for j = 1:fea_numi*fea_numk_max
            for k = 1:intput_num
                wfkli_d(j,(j-1)*intput_num+k) = (1-Z(j).^2).*ss(k).*wro(k).*uSODLRBLS_old(k)./uSODLRBLS_old_old(k);
            end
        end
    else
        for j = 1:fea_numi*fea_numk_max
            for k = 1:intput_num
                wfkli_d(j,(j-1)*intput_num+k) = (1-Z(j).^2).*ss(k);
            end
        end
    end
    wfkli = wfkli-reshape(ewfkli*wfkli_d'*wki*ss,fea_numi*fea_numk_max,intput_num);
    
    %%%%%% bfki updating %%%%%%%
    bfki_d1 = zeros(1,fea_numi*fea_numk_max);
    bfki_d = zeros(fea_numi*fea_numk_max,fea_numi*fea_numk_max);
    for j = 1:fea_numi*fea_numk_max
        bfki_d1(j) = (1-Z(j).^2);
    end
    for j = 1:fea_numi*fea_numk_max
        bfki_d(j,:) = bfki_d1;
    end
    bfki = bfki-reshape(ebfki*bfki_d'*wki*ss,fea_numi*fea_numk_max,1);
    
    %%%%%% wro updating %%%%%%%
    if t>0.001
        wro_d1 = zeros(1,intput_num);
        wro_d = zeros(fea_numi*fea_numk_max,intput_num);
        for j = 1:intput_num
            wro_d1(j) = (1-Z(j).^2).*sum(wfkli(:,j)).*(ss(j).*uSODLRBLS_old(j)./uSODLRBLS_old_old(j));
        end
        for j = 1:fea_numi*fea_numk_max
            wro_d(j,:) = wro_d1;
        end
        wro = wro-reshape(ewro*wro_d'*wki*ss,intput_num,1);
    end
    
    %%%%%% wjki updating %%%%%%%
    wjki_d = zeros(en_num,en_num*fea_numi*fea_numk_max);
    for j = 1:en_num
        for k = 1:fea_numi*fea_numk_max
            wjki_d(j,(j-1)*fea_numi*fea_numk_max+k) = -2.*H(j).*Z(k).*(wjki(j,k)*Z(k)+wr(j).*H_old(j)-cc(j,k))./(bb(j,k).^2);
        end
    end
    wjki = wjki-reshape(ewjki*wjki_d'*wj*ss,en_num,fea_numi*fea_numk_max);
    
    %%%%%% bb updating %%%%%%%
    bb_d = zeros(en_num,en_num*fea_numi*fea_numk_max);
    for j = 1:en_num
        for k = 1:fea_numi*fea_numk_max
            bb_d(j,(j-1)*fea_numi*fea_numk_max+k) = 2.*H(j).*(wjki(j,k)*Z(k)+wr(j).*H_old(j)-cc(j,k)).^2./(bb(j,k).^3);
        end
    end
    bb = bb-reshape(ebb*bb_d'*wj*ss,en_num,fea_numi*fea_numk_max);
    
    %%%%%% cc updating %%%%%%%
    cc_d = zeros(en_num,en_num*fea_numi*fea_numk_max);
    for j = 1:en_num
        for k = 1:fea_numi*fea_numk_max
            cc_d(j,(j-1)*fea_numi*fea_numk_max+k) = 2.*H(j).*(wjki(j,k)*Z(k)+wr(j).*H_old(j)-cc(j,k))./(bb(j,k).^2);
        end
    end
    cc = cc-reshape(ecc*cc_d'*wj*ss,en_num,fea_numi*fea_numk_max);
    
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
    E4=[E4 error];
    D_E4=[D_E4 derror];
    X4=[X4 x];
    X_D4=[X_D4 x_d];
    XD4=[XD4 xd];
    XD_D4=[XD_D4 xd_d];
    U4=[U4 U_t];
    
    error_old = error;
end

fprintf('The RMSE is %d\n',sqrt(sum((E4).^2)/(TT/ts)));

T=0:0.001:7;
figure(1)
plot(XD4,XD_D4,'r',X4,X_D4,'g')
xlabel('x');
ylabel('x-dot');
title('Trajectory phase portrait of Duffing-Holmes Chaotic System');
grid on;
figure(2)
subplot(2,1,1)
plot(T,X4,'r',T,XD4,'g')
xlabel('Time(s)');
grid on;
legend(': x',': xd ');
subplot(2,1,2)
plot(T,X_D4,'r',T,XD_D4,'g')
xlabel('Time(s)');
grid on;
legend(': x-dot ',': xd-dot');
figure(3)
plot(T,E4,'r',T,D_E4,'g')
xlabel('Time(s)');
ylabel('Error');
grid on;
title('Tracking error');
legend(': e ',': e-dot');
figure(4)
plot(T,NUMK4);
axis([-1 8 0 21]);
xlabel('Time(s)');
ylabel('Number of the remaining feature neuron');
