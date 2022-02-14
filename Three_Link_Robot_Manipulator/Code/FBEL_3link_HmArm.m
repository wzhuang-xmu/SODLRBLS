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

%%%%%% FBEL parameters %%%%%%%
mns = [-2.1:0.6:2.1; -2.1:0.6:2.1; -2.1:0.6:2.1];
v_std = 1.2*ones(3,8);
v_w = unifrnd(-1,1,3,8,3);
v_w_d = zeros(3,8,3);
wght = unifrnd(-1,1,3,8,3);
wght_d = zeros(3,8,3);
a_w = zeros(3,1);
o_w = zeros(3,1);
u_exp = zeros(3,8);
Rd_p = zeros(3,1);
output = zeros(3,1);

%%%%%% learning rate %%%%%%%
alpha = 0.0001.*[0.05,0.05,0.05];
beta = 0.0001.*[0.05,0.05,0.05];
b_gains = 0.0001*ones(1,3);
c_gains = 0.0001.*[1.00 1.00 1.00];

%%%%%% plot parameters %%%%%%%
REF31 = [ ];
REF32 = [ ];
REF33 = [ ];
Q31 = [ ];
Q32 = [ ];
Q33 = [ ];
E31 = [ ];
E32 = [ ];
E33 = [ ];
URC31 = [];
URC32 = [];
URC33 = [];

RMSE = zeros(3,1);%root mean square error
Allts = 25;%total time
ts = 0.001;%sampling time
RR=0.066*eye(3);%attenuation coefficient
uFBEL_old=zeros(3,1);

%%%%%% main program %%%%%%%
for t=0:ts:Allts
    Q=q_old;%update Q
    qd=qd_old;%update qd
    qdd=qdd_old;%update qdd
    if t>0
        uFBEL_old=uFBEL;
    end
    s1 = sin(Q(1)); s2 = sin(Q(2)); s3 = sin(Q(3));%joint 1,2,3 Y-axis
    c1 = cos(Q(1)); c2 = cos(Q(2)); c3 = cos(Q(3));%joint 1,2,3 X-axis
    s23 = sin(Q(2) + Q(3)); c23 = cos(Q(2) + Q(3));%horizontal and vertical coordinates relative to a13
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
    
    %%%%%% count the output of FBEL %%%%%%%
    for i = 1:3    
        for j = 1:8
            u_exp(i,j) = exp(-(ss(i)-mns(i,j))^2/(v_std(i,j)^2));
        end
    end
    for o_d = 1:3
        a_w(o_d) = sum(diag(u_exp*v_w(:,:,o_d)'));
        o_w(o_d) = sum(diag(u_exp*wght(:,:,o_d)'));
    end
    output = a_w - o_w;
    uFBEL = output;
    
    uRC = ((RR^2+eye(3))/(2*RR^2))*ss;%the output robust controller    
    Tao = uRC + uFBEL;    
    
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
    QDD = inv(Mq)*(Tao-Tao_d-Cq_qd*qd-Gq);
    qdd_old=QDD;
    qd_old=qd_old+qdd_old*ts;
    q_old=q_old+qd_old*ts;
    
    %%%%%% FBEL parameters updating %%%%%%%
    for i = 1:3
        Rd_p(i) = b_gains(i)*ss(i) + c_gains(i)*uFBEL(i);
    end
    for i = 1:3
        for j = 1:8
            for o_d = 1:3
                v_w_d(i,j,o_d) = alpha(o_d)*(u_exp(i,j)*(max(0,Rd_p(o_d)-a_w(o_d))));
                wght_d(i,j,o_d) = beta(o_d)*(u_exp(i,j)*(uFBEL(o_d)-Rd_p(o_d)));
            end
        end
    end
    v_w = v_w + v_w_d;
    wght = wght + wght_d;
    
    %%%%%% plot %%%%%%%
    REF31=[REF31 ref(1)];
    REF32=[REF32 ref(2)];
    REF33=[REF33 ref(3)];
    Q31=[Q31 Q(1)];
    Q32=[Q32 Q(2)];
    Q33=[Q33 Q(3)];
    E31=[E31 error(1)];
    E32=[E32 error(2)];
    E33=[E33 error(3)];
    URC31=[URC31 uRC(1)];
    URC32=[URC32 uRC(2)];
    URC33=[URC33 uRC(3)];
    
    error_old = error;
end

T = 0:0.001:25;
figure(1)
plot(T,REF31,'r',T,Q31,'g');
xlabel('Time(s)');
ylabel('Angle 1(rad)');
figure(2)
plot(T,REF32,'r',T,Q32,'g');
xlabel('Time(s)');
ylabel('Angle 2(rad)');
figure(3)
plot(T,REF33,'r',T,Q33,'g');
xlabel('Time(s)');
ylabel('Angle 3(rad)');
figure(4)
plot(T,E31,'r');
xlabel('Time(s)');
ylabel('Error(rad)');
figure(5)
plot(T,E32,'r');
xlabel('Time(s)');
ylabel('Error(rad)');
figure(6)
plot(T,E33,'r');
xlabel('Time(s)');
ylabel('Error(rad)');

fprintf('The last error1 is %d;The RMSE1 is %d\n',error(1),sqrt(sum(E31.^2)/(t/ts)));
fprintf('The last error2 is %d;The RMSE2 is %d\n',error(2),sqrt(sum(E32.^2)/(t/ts)));
fprintf('The last error3 is %d;The RMSE3 is %d\n',error(3),sqrt(sum(E33.^2)/(t/ts)));