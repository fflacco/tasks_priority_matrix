%Simulating ROMEO
%      Copyright 2016 Fabrizio Flacco
%  
%      Licensed under the Apache License, Version 2.0 (the "License");
%      you may not use this file except in compliance with the License.
%      You may obtain a copy of the License at
%      http://www.apache.org/licenses/LICENSE-2.0
%  
%      Unless required by applicable law or agreed to in writing, software
%      distributed under the License is distributed on an "AS IS" BASIS,
%      WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%      See the License for the specific language governing permissions and
%      limitations under the License.

T=0.001; %sampling time
Ttot=5; % total simulation time

q0=deg2rad([0,10,20,10,0,0,0,0,-45,0,0,0,0,-10,20,-10,0,0,45,0,0])';%initial configuration

n=21;
I=eye(n);


%initial position
p0_lfoot=get_RLEG_LLEG_POSITION(q0);
p0_lhand=get_RLEG_LHAND_POSITION(q0);
p0_rhand=get_RLEG_RHAND_POSITION(q0);
p0_com=get_RLEG_COM_POSITION(q0);

%current configuration
q=q0;

%controller gains
Kp0=100;
Kp1=100;
Kp2=100;
Kp3=100;

Tcircle=2;
r2=0.1; %circle radius
r3=0.1; %circle radius

desiredVelocity2=zeros(3,1);
desiredPosition2=p0_lhand;

desiredVelocity3=zeros(3,1);
desiredPosition3=p0_rhand;

%saving
QV=zeros(Ttot/T+1,21);
T0d=zeros(Ttot/T+1,3);
T1d=zeros(Ttot/T+1,3);
T2d=zeros(Ttot/T+1,3);
T3d=zeros(Ttot/T+1,3);
extTime=zeros(Ttot/T+1,1);

i=1;
for t=0:T:Ttot
   
    %save current configuration
    QV(i,:)=q';
    T0d(i,:)=p0_com';
    T1d(i,:)=p0_lfoot';
    T2d(i,:)=desiredPosition2';
    T3d(i,:)=desiredPosition3';

    %TASK 0: keep com  to the initial XY position
    p_com=get_RLEG_COM_POSITION(q);
    Acom=get_RLEG_COM_JACOBIAN(q);
    A0=Acom(1:2,:);
    
    b0= Kp0 * (p0_com(1:2) - p_com(1:2));
    
    %TASK 1: keep left foot to the initial position
    p_lfoot=get_RLEG_LLEG_POSITION(q);
    A1=get_RLEG_LLEG_JACOBIAN(q);

    b1= Kp1 * (p0_lfoot - p_lfoot);
    
  
    %TASK 2: circle with the left hand+
    p_lhand=get_RLEG_LHAND_POSITION(q);
    A2=get_RLEG_LHAND_JACOBIAN(q);
 
	desiredVelocity2(1)=0;
	desiredVelocity2(2)=-2*(2*pi*r2*cos((2*pi*(t))/Tcircle))/Tcircle;
	desiredVelocity2(3)=-(2*pi*r2*sin((2*pi*(t))/Tcircle))/Tcircle;

	b2= Kp2*(desiredPosition2-p_lhand);
   	b2= b2 + desiredVelocity2;

	desiredPosition2(1)=p0_lhand(1);
	desiredPosition2(2)=  -2*r2*sin(2*pi*(t)/Tcircle)+p0_lhand(2);
	desiredPosition2(3)= (r2*cos(2*pi*(t)/Tcircle)-r2)  + p0_lhand(3);

    %TASK 3: circel with the right hand
    p_rhand=get_RLEG_RHAND_POSITION(q);
    A3=get_RLEG_RHAND_JACOBIAN(q);
    
    desiredVelocity3(1)=0;
	desiredVelocity3(2)=2*(2*pi*r3*cos((2*pi*t)/Tcircle))/Tcircle;
	desiredVelocity3(3)=-(2*pi*r3*sin((2*pi*t)/Tcircle))/Tcircle;


	b3= Kp3*(desiredPosition3-p_rhand);
	b3= b3 + desiredVelocity3;

	desiredPosition3(1)=p0_rhand(1);
	desiredPosition3(2)= 2*r3*sin(2*pi*t/Tcircle)+ p0_rhand(2);
	desiredPosition3(3)= r3*cos(2*pi*t/Tcircle)-r3  + p0_rhand(3);    
    
    
    
    % solve tasks with priority
    eps=1e-1;
    lambda=0.5e-1;
    
    
    A=[A0;A1;A2;A3];
    idx=[2,3,3,3];
    b=[b0;b1;b2;b3];
    tic    
        [Q,R]=qr(A',0);
        X=tasksPriorityMatrix(R,idx,lambda,eps);
        dq_prref=Q*pinv(R')*X*b;    
    extTime(i)=toc;
  %  TIME1=[TIME1;time];

    q= q + dq_prref*T;
    %t=t+T;
    i=i+1;

    if (mod(t,1)==0)
        t
    end
end


%saving with different naem
QV_2=QV;
T0d_2=T0d;
T1d_2=T1d;
T2d_2=T2d;
T3d_2=T3d;
extTime_2=extTime;



