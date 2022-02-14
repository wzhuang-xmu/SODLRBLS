clear;
close all;

%%%%%% manipulator parameters %%%%%%%
a1 = 0.6; a2 = 0.5; a3 = 0.4;%the length of each joint of the manipulator
l1 = 1/2*a1; l2 = 1/2*a2; l3 = 1/2*a3;%the position of the center of mass of each joint of the manipulator
ml1 = 3; ml2 = 1.8; ml3 = 1.5;%the weight of each joint of the manipulator
mm1 = 0.3; mm2 = 0.3; mm3 = 0.3;%the weight of each rotor of the manipulator
Il1 = 50.45e-3; Il2 = 32.68e-3; Il3 = 30.47e-3;%the moment of inertia of each joint of the manipulator
Im1 = 12.0e-3; Im2 = 12.0e-3; Im3 = 12.0e-3; %the moment of inertia of each rotor of the manipulator
kr1 = 1; kr2 = 1; kr3 = 1;%motor reduction ratio
grv = 9.8;%gravitational acceleration
teda1 = -0.3;%initial angle of joint 1
teda2 = 0.1;%initial angle of joint 2
teda3 = -0.4;%initial angle of joint 3
q_old = [teda1 teda2 teda3]';%angular displacement vector for joint position
qd_old = [0 0 0]';%angular velocity vector
qdd_old = [0 0 0]';%angular acceleration vector
xd = [0 0 0]';%position vector
xd_d = [0 0 0]';%velocity vector
xd_dd = [0 0 0]';%acceleration vector
error = [0 0 0]';%error vector
error_old = [0 0 0]';%old error vector
derror = [0 0 0]';%derivative of error
ss = [0 0 0]';%smc controller
Tao = [0 0 0]';%output torque
Q = q_old;%angular displacement
qd = qd_old;%angular velocity
qdd = qdd_old;%angular acceleration

%%%%%% SODLRBLS parameters %%%%%%%
intput_num = 3;%input dimension
output_num = 3;%output dimension
fea_numi = 3;%the number of feature group
fea_numk_max = 10;%the maximum number of feature neuron in each feature group
fea_numk = randi([1 10],intput_num,1);%random the number of feature neuron in each feature group
en_num = 3;%the number of enhancement neuron

%%%%%% initial parameters of SODLRBLS %%%%%%%
wfkli = unifrnd(-0.01,0.01,[fea_numi*fea_numk_max,intput_num]);
bfki = unifrnd(-0.01,0.01,[fea_numi*fea_numk_max,1]);
wjki = unifrnd(-0.01,0.01,[en_num,fea_numi*fea_numk_max]);
wki = unifrnd(-0.01,0.01,[fea_numi*fea_numk_max,output_num]);
wj = unifrnd(-0.01,0.01,[en_num,output_num]);
wro = unifrnd(-0.01,0.01,[intput_num,1]);
wr = unifrnd(-0.01,0.01,[en_num,1]);
bb = unifrnd(-0.01,0.01,[en_num,fea_numi*fea_numk_max]);
cc = unifrnd(-0.01,0.01,[en_num,fea_numi*fea_numk_max]);

%%%%%% learning rate %%%%%%%
ewfkli = 0.0001*eye(intput_num*fea_numi*fea_numk_max);
ebfki = 0.0001*eye(fea_numi*fea_numk_max);
ewjki = 0.0001*eye(en_num*fea_numi*fea_numk_max);
ewki = 0.0001*eye(fea_numi*fea_numk_max);
ewj = 0.0001*eye(en_num);
ewro = 0.0001*eye(intput_num);
ewr = 0.0001*eye(en_num);
ebb = 0.0001*eye(en_num*fea_numi*fea_numk_max);
ecc = 0.0001*eye(en_num*fea_numi*fea_numk_max);

ddd = 0;
dgg = 0;
dd = 0;%deleting threshold
dg = 0;%generating threshold
Pwc = 0.1;%deleting threshold coefficient
Pbc = 0.1;%generating threshold coefficient

%%%%%% plot parameters %%%%%%%
REF11 = [ ];
REF12 = [ ];
REF13 = [ ];
Q11 = [ ];
Q12 = [ ];
Q13 = [ ];
E11 = [ ];
E12 = [ ];
E13 = [ ];
URC11 = [];
URC12 = [];
URC13 = [];
NUMK = [];

RMSE = zeros(output_num,1);%root mean square error
RR = 0.066*eye(output_num);%attenuation coefficient
Allts = 25;%total time
ts = 0.001;%sampling time
uSODLRBLS = zeros(output_num,1);%the output of SODLRBLS
uSODLRBLS_old = zeros(output_num,1);
uSODLRBLS_old_old = zeros(output_num,1);
H_old = zeros(en_num,1);

%%%%%% main program %%%%%%%
for t = 0:ts:Allts
    NUMK = [NUMK sum(fea_numk)];
    
    Q = q_old;%update Q
    qd = qd_old;%update qd
    qdd = qdd_old;%update qdd
    s1 = sin(Q(1));s2 = sin(Q(2));s3 = sin(Q(3));%joint 1,2,3 Y-axis
    c1 = cos(Q(1));c2 = cos(Q(2));c3 = cos(Q(3));%joint 1,2,3 X-axis
    s23 = sin(Q(2)+Q(3));c23 = cos(Q(2)+Q(3));%horizontal and vertical coordinates relative to a13
    %when the manipulator starts to move, the manipulator needs to track the first reference pattern; 
    %after 15 seconds, the manipulator needs to track the second reference pattern
    if t >15
        ref = 0.5*[0.5*(sin(2*t)+cos(t+1)) sin(2*t)*cos(t+1) 0.25-0.1*sin(t+1)-0.1*sin(2*t)]';
    else
        ref = 0.5*[0.5*(sin(t+2.5)+0.75*cos(2*t+1.5)) (sin(t)+sin(2*t)) 0.2*cos(2*t)-0.2*sin(t)]';
    end
    error = ref-Q;
    RMSE = RMSE + error.^2;
    if t ~= 0
        derror = (error-error_old)/ts;
    end
    ss = 0.55*derror+[10;10;10].*error;
    
    %%%%%% count the output of SODLRBLS %%%%%%%
    Z1 = zeros(fea_numi,fea_numk_max);
    Z2 = zeros(fea_numi,fea_numk_max);
    Z = zeros(fea_numi*fea_numk_max,1);
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
    
    uRC = ((RR^2+eye(3))/(2*RR^2))*ss;%the output robust controller
    Tao = uRC +uSODLRBLS;
    
    %%%%%% dynamic function %%%%%%%
    %%%%%% inertia matrix %%%%%%%
    bq11 = Il1 + Im2 + Il2 + Il3 + Im3 + Im1*kr1^2 + ml3*(a2*c2 + l3*c23)^2 + a2^2*c2^2*mm3 + c2^2*l2^2*ml2;
    bq12 = 0;
    bq13 = 0;
    bq21 = 0;
    bq22 = ml3*(a2^2+l3^2+2*a2*l3*c3) + a2^2*mm3 + l2^2*ml2 + Il2 + Il3 + Im3 + Im2*kr2^2;
    bq23 = Il3 + Im3*kr3 + a2*l3*ml3*c3 + l3^2*ml3;
    bq31 = 0;
    bq32 = Im3*kr3 + Il3 + a2*l3*ml3*c3 + l3^2*ml3;
    bq33 = l3^2*ml3 + Il3 + Im3*kr3^2;
    Mq = [bq11 bq12 bq13; bq21 bq22 bq23; bq31 bq32 bq33];
    
    %%%%%% centripetal matrix %%%%%%%
    c111 = 0;
    c112 = 1/2*(2*ml3*(a2*c2 + l3*c23)*(-a2*s2-l3*s23)-2*(a2^2*mm3+l2^2*ml2)*s2);
    c113 = 1/2*(2*ml3*(a2*c2 + l3*c23)*l3*(-s23));
    c121 = c112;
    c122 = 0;
    c123 = 0;
    c131 = c113;
    c132 = c123;
    c133 = 0;
    c211 = -c112;
    c212 = -c122;
    c213 = -c123;
    c221 = -c122;
    c222 = 0;
    c223 = 1/2*ml3*(-2*a2*l3*s3);
    c231 = c213;
    c232 = c223;
    c233 = -a2*l3*ml3*s3;
    c311 = -c113;
    c312 = c213;
    c313 = 0;
    c321 = c312;
    c322 = -c223;
    c323 = 0;
    c331 = 0;
    c332 = 0;
    c333 = 0;
    cq11 = c111*qd(1) + c112*qd(2) + c113*qd(3);
    cq12 = c121*qd(1) + c122*qd(2) + c123*qd(3);
    cq13 = c131*qd(1) + c132*qd(2) + c133*qd(3);
    cq21 = c211*qd(1) + c212*qd(2) + c213*qd(3);
    cq22 = c221*qd(1) + c222*qd(2) + c223*qd(3);
    cq23 = c231*qd(1) + c232*qd(2) + c233*qd(3);
    cq31 = c311*qd(1) + c312*qd(2) + c313*qd(3);
    cq32 = c321*qd(1) + c322*qd(2) + c323*qd(3);
    cq33 = c331*qd(1) + c332*qd(2) + c333*qd(3);
    Cq_qd = [cq11 cq12 cq13; cq21 cq22 cq23; cq31 cq32 cq33];
    
    %%%%%% gravity matrix %%%%%%%
    Gq = [0;
        grv*ml3*(a2*c2 + c23*l3) + c2*grv*l2*ml2 + a2*c2*grv*mm3;
        c23*grv*l3*ml3];
    
    Tao_d = 2*[0.2*sin(2*t);0.1*cos(2*t);0.1*sin(t)];%external disturbance
    QDD = Mq\(Tao-Tao_d-Cq_qd*qd-Gq);
    qdd_old = QDD;
    qd_old = qd_old+qdd_old*ts;
    q_old = q_old+qd_old*ts;
    
    %%%%%% TOPSIS %%%%%%%
    for i = 1:fea_numi
        for k = 1:fea_numk_max
            if sum(Z2(:,k)) ~= 0
                ZZZ(i,k) = Z2(i,k)/sqrt(sum(Z2(:,k).^2));
            else
                ZZZ(i,k) = 0;
            end
        end
    end
    
    for k = 1:fea_numk_max
        ma = ZZZ(:,k);
        summ = 0;
        for x = 1:length(ma)
            if ma(x) > 0
                summ = summ+ma(x)*log(ma(x));
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
    for x = 1:length(maa)
        if maa(x) < 1 && maa(x) > 0.00001
            if maa(x) > maxx
                maxx = maa(x);
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
    for x = 1:length(maa)
        if maa(x) > 0.00001 && maa(x) < 1
            if maa(x) < minn
                minn = maa(x);
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
        for i=1:fea_numi
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
        bfki_d1(j)=(1-Z(j).^2);
    end
    for j = 1:fea_numi*fea_numk_max
        bfki_d(j,:) = bfki_d1;
    end
    bfki = bfki-reshape(ebfki*bfki_d'*wki*ss,fea_numi*fea_numk_max,1);
    
    %%%%%% wro updating %%%%%%%
    if t>0.001
        wro_d1 = zeros(1,intput_num);
        wro_d = zeros(intput_num*fea_numk_max,intput_num);
        for j = 1:intput_num
            wro_d1(j) = (1-Z(j).^2).*sum(wfkli(:,j)).*(ss(j).*uSODLRBLS_old(j)./uSODLRBLS_old_old(j));
        end
        for j = 1:intput_num*fea_numk_max
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
    for j = 1:en_num
        wr_d(j,:) = wr_d1;
    end
    wr = wr-reshape(ewr*wr_d'*wj*ss,en_num,1);
    
    %%%%%% plot %%%%%%%
    REF11 = [REF11 ref(1)];
    REF12 = [REF12 ref(2)];
    REF13 = [REF13 ref(3)];
    Q11 = [Q11 Q(1)];
    Q12 = [Q12 Q(2)];
    Q13 = [Q13 Q(3)];
    E11 = [E11 error(1)];
    E12 = [E12 error(2)];
    E13 = [E13 error(3)];
    URC11 = [URC11 uRC(1)];
    URC12 = [URC12 uRC(2)];
    URC13 = [URC13 uRC(3)];
    
    error_old = error;
end

T = 0:0.001:25;
figure(1)
plot(T,REF11,'r',T,Q11,'g');
xlabel('Time(s)');
ylabel('Angle 1(rad)');
figure(2)
plot(T,REF12,'r',T,Q12,'g');
xlabel('Time(s)');
ylabel('Angle 2(rad)');
figure(3)
plot(T,REF13,'r',T,Q13,'g');
xlabel('Time(s)');
ylabel('Angle 3(rad)');
figure(4)
plot(T,E11,'r');
xlabel('Time(s)');
ylabel('Error(rad)');
figure(5)
plot(T,E12,'r');
xlabel('Time(s)');
ylabel('Error(rad)');
figure(6)
plot(T,E13,'r');
xlabel('Time(s)');
ylabel('Error(rad)');
figure(7)
plot(T,NUMK);
xlabel('Time(s)');
ylabel('Number of the remaining feature neuron');
axis([-2 27 0 31]);

fprintf('The last error1 is %d;The RMSE1 is %d\n',error(1),sqrt(sum(E11.^2)/(t/ts)));
fprintf('The last error2 is %d;The RMSE2 is %d\n',error(2),sqrt(sum(E12.^2)/(t/ts)));
fprintf('The last error3 is %d;The RMSE3 is %d\n',error(3),sqrt(sum(E13.^2)/(t/ts)));