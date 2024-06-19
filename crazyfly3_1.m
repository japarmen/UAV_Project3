%% PARÂMETROS
clc
clear all;
initPlots;

% Simulation parameters
nx = 6 ;
ny = 6 ;
nu = 3 ;
x0 = [0; 0; 0; 0 ; 0 ; 0] ;
Dt = 0.001 ;
t = 0:Dt:5 ;
g = 9.81 ;
rho = 1.225 ;
Tmax = .7848 ;
lx = 0.65 ;
ly = 0.065 ;
lz = 0.029 ;
l = 0.046 ;
Ax = ly * lz ;
Ay = lx * lz ;
Az = lx * ly ;
m = 0.028 ;
cd = 2.3 ;
betax = .5*cd*rho*Ax ;
betay = .5*cd*rho*Ay ;
betaz = .5*cd*rho*Az ;
D = diag([ betax , betay , betaz ]) ;
zI = [0 ; 0 ; -1] ;

%% SIMULAÇÃO NÃO LINEAR
T = 0.4*Tmax ;
fp = [0;0;T] ;

phi = deg2rad(0.5);
theta = deg2rad(0.25);
psi = deg2rad(0);

u_NL = [  T  ;  phi;   theta;  psi]*ones(size(t));


% simulate nonlinear system
Nsim = length(t);
x = zeros(nx,Nsim);
y = zeros(ny,Nsim);
x(:,1) = x0;

for k = 1:Nsim
    % prepare variables:
    T = u_NL(1,k);
    lbd = u_NL(2:4,k);
    R = Euler2R(lbd);
    p = x(1:3,k);
    v = x(4:6,k);
   
  
    
    % compute state derivative:
    p_dot = Euler2R(lbd)*v;
    ua = g*R'*zI + (1/m)*fp ;
    v_dot = ua-(1/m)*R*D*R'* v .* abs(v) ; %V_dot= fa + ua 
    

    %parte3
    ua = g*R'*zI + (1/m)*fp ;
    v_dot = ua-(1/m)*R*D*R'* v .* abs(v) ; %V_dot= fa + ua 

    x_dot = [p_dot;v_dot];
    % integrate state
    x(:,k+1) = x(:,k) + Dt*x_dot;

    % compute current output:
    y(:,k) = x(:,k);
end


figure;
plot(t,y(1,:),t,y(2,:),t,y(3,:),t,y(4,:),'-.',t,y(5,:),'-.',t,y(6,:), '-.');
grid on;
xlabel('$$t$$ [s]');
legend('$$p_x$$','$$p_y$$','$$p_z$$','$$v_x$$','$$v_y$$','$$v_z$$');


figure;
plot(t,u_NL(1,:),t,u_NL(2,:)*180/pi,t,u_NL(3,:)*180/pi,t,u_NL(4,:)*180/pi);
grid on;
xlabel('$$t$$ [s]');
legend('$$T$$','$$\phi$$','$$\theta$$','$$\psi$$');


%% DERIVADAS A E B

clear T px py pz phi theta psi vx vy vz ua ax ay az
syms T px py pz phi theta psi vx vy vz ua ax ay az

ua = [ax;ay;az] ;
v = [vx;vy;vz] ;
lbd = [phi;theta;psi] ;


rotx = [	1	,	0			,	 0         ;
			0	,	cos(phi)	,	-sin(phi)  ;
			0	,	sin(phi)	,	 cos(phi)  ];

roty = [	 cos(theta)	,	0	,	sin(theta)  ;
			0				,	1	,	0       ;
			-sin(theta)	,	0	,	cos(theta)	];

rotz = [	cos(psi)	,	-sin(psi)	,	0   ;
			sin(psi)	,	 cos(psi)	,	0   ;
			0			,	 0		,	1	];

R = rotz*roty*rotx;

p_dot = R*v;
v_dot = ua-(1/m)*R*D*R'* v .* abs(v) ;

a1 = diff(p_dot,px) ; a2 = diff(p_dot,py) ; a3 = diff(p_dot,pz) ;
a4 = diff(p_dot,vx) ; a5 = diff(p_dot,vy) ; a6 = diff(p_dot,vz) ;
a7 = diff(v_dot,px) ; a8 = diff(v_dot,py) ; a9 = diff(v_dot,pz) ;
a10 = diff(v_dot,vx) ; a11 = diff(v_dot,vy) ; a12 = diff(v_dot,vz) ;

A = [ a1(1) , a2(1) , a3(1) , a4(1) , a5(1) , a6(1) ;
      a1(2) , a2(2) , a3(2) , a4(2) , a5(2) , a6(2) ;
      a1(3) , a2(3) , a3(3) , a4(3) , a5(3) , a6(3) ;
      a7(1) , a8(1) , a9(1) , a10(1) , a11(1) , a12(1) ;
      a7(2) , a8(2) , a9(2) , a10(2) , a11(2) , a12(2) ;
      a7(3) , a8(3) , a9(3) , a10(3) , a11(3) , a12(3) ] ;

b1 = diff(p_dot,ax) ; b2 = diff(p_dot,ay) ; b3 = diff(p_dot,az) ; 
b4 = diff(v_dot,ax) ; b5 = diff(v_dot,ay) ; b6 = diff(v_dot,az) ; 

B = [ b1(1) , b2(1) , b3(1)  ;
      b1(2) , b2(2) , b3(2)  ; 
      b1(3) , b2(3) , b3(3)  ;
      b4(1) , b5(1) , b6(1)  ;
      b4(2) , b5(2) , b6(2)  ;
      b4(3) , b5(3) , b6(3)  ] ;


valores = {T, px, py, pz, phi, theta, psi, vx, vy, vz, ax, ay, az};
num_valores = {0.4*Tmax, 0, 0, 5, deg2rad(0.5), deg2rad(0.25), 0, 5, 5, 5, 1, .5, 2};

A = subs(A, valores, num_valores);
B = subs(B, valores, num_valores);

A = double(A)
B = double(B)

C = eye(nx);

D = zeros(nx,nu);

%% LINEAR SYSTEM


%simulate linear system:
sys = ss(A,B,C,D);
T = 0.4*Tmax ;
phi = deg2rad(0.5) ;
theta = deg2rad(0.25) ;
psi = 0 ;
px = 0 ;
py = 0 ;
pz = 5 ;
vx = 5 ;
vy = 5 ;
vz = 5 ;
ax = 1 ;
ay = .5 ;
az = 2 ;

u_L = [ ax ; ay ; az ]*(t>=0);
xop = [px ; py ; pz ; vx ; vy ; vz ] ;
y_L = lsim(sys,u_L,t,xop)';

figure;
plot(t,y_L(1,:),t,y_L(2,:),t,y_L(3,:),t,y_L(4,:),'-.',t,y_L(5,:),'-.',t,y_L(6,:), '-.');
grid on;
xlabel('$$t$$ [s]');
legend('$$p_x$$','$$p_y$$','$$p_z$$','$$v_x$$','$$v_y$$','$$v_z$$');



figure;
plot(t,u_L(1,:),t,u_L(2,:),t,u_L(3,:));
grid on;
xlabel('$$t$$ [s]');
legend('$$a_x$$','$$a_y$$','$$a_z$$');


%% LQR CONTROLLER

% Weighting matrices Q and R
Q = eye(nx); 
R = eye(nu); 

[K,P] = lqr(A, B, Q, R)

% Define the closed-loop system
A_cl = A - B * K;
B_cl = B;
C_cl = eye(nx); 
D_cl = zeros(nx,nu);

% Simulate the closed-loop ans step response 
sys_cl = ss(A_cl, B_cl, C_cl, D_cl);

x0 = [1; 0; 0; 0 ; 0 ; 0] ;

[y, t, x] = initial(sys_cl, x0, t);

figure;
plot(t,x(:,1),t,x(:,2),t,x(:,3),t,x(:,4),'-.',t,x(:,5),'-.',t,x(:,6), '-.');
xlabel('$$t$$ [s]');
legend('$$p_x$$','$$p_y$$','$$p_z$$','$$v_x$$','$$v_y$$','$$v_z$$');
grid on;

figure
step(sys_cl)
