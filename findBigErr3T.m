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

clear 


err=0;
minerr=1e-2;

a=-pi*1;
b=pi*1;
i=0;

la=0.1;
lb=1;

bestQ=zeros(6,1);
bestL=zeros(6,1);
maxErr=0;
Err2=0;
n=10000;

dxa=-100;
dxb=100;

E11=[0];
E21=[0];
E31=[0];
TIME1=[];


E12=[0];
E22=[0];
E32=[0];
TIME2=[];

E1C=[0];
E2C=[0];
E3C=[0];
TIMEC=[];

E1prref=[0];
E2prref=[0];
E3prref=[0];
TIMEprref=[];


E1NC=[];
E2NC=[];

delta=1e-10;


max_err11=0;



while ((err<minerr)&&(i<n))
    q= a + (b-a).*rand(6,1);
    l= la + (lb-la).*rand(6,1);
    J1=getJacobianNRL(q,l);
    J2a=getJacobianNRL(q(1:4),l(1:4));
    J2=[J2a,[0 0;0 0]];
    J3a=getJacobianNRL(q(1:2),l(1:2));
    J3=[J3a,[0 0 0 0;0 0 0 0]];

   
    
    dx1=dxa + (dxb-dxa).*rand(2,1);
    dx2=dxa + (dxb-dxa).*rand(2,1);
    dx3=dxa + (dxb-dxa).*rand(2,1);




    I=eye(6);
    tic
        [pJ1,Q1]=damped_pinv_P8(J1);
        dq1=pJ1*dx1;
        P1=I-Q1;
        [pJ12,Q12]=damped_pinv_P8(J2*P1);    
        dq12=dq1 + pJ12*(dx2 -J2*dq1);
        P12=P1 - Q12;
        [pJ123,Q123]=damped_pinv_P8(J3*P12);    
        dq123=dq12 + pJ123*(dx3-J3*dq12);
    time=toc;
    TIME1=[TIME1;time];


    
 

    
    
    tic
        [pJ1,Q1]=damped_pinv_P8(J1);
        P1=I-Q1;
        [pJ12,Q12]=damped_pinv_P8(J2*P1); 
        P12=P1 - Q12;
        dqC=pJ1*dx1+P1*damped_pinv_P8(J2)*dx2 + P12*damped_pinv_P8(J3)*dx3;
    time=toc;
    TIMEC=[TIMEC;time]; 
    
    dq=zeros(6,1);
    dqi=zeros(6,1);
    pJ=[];
    P=[];
    tic
        [dq,pJ,P]=RP_transitionC(J3,dx3,1,dq,dqi,pJ,P);
        [dq,pJ,P]=RP_transitionC(J2,dx2,1,dq,dqi,pJ,P);
        [dq,pJ,P]=RP_transitionC(J1,dx1,1,dq,dqi,pJ,P);
    time=toc;
    TIME2=[TIME2;time];    
    dq321a=dq;
    
    A=[J1;J2;J3];
    tic
        [Q,R]=qr(A',0);
        idx=[2,2,2];
        X=tasksPriorityMatrix(R,idx);
        dq_prref=Q*(damped_pinv_P8(R')*(X*[dx1;dx2;dx3]));    
    time=toc;
    TIMEprref=[TIMEprref;time];    
    
    
    err11=100*norm(J1*dq123-dx1)/norm(dx1);
    err21=100*norm(J2*dq123-dx2)/norm(dx2);
    err31=100*norm(J3*dq123-dx3)/norm(dx3);
    err12=100*norm(J1*dq321a-dx1)/norm(dx1);
    err22=100*norm(J2*dq321a-dx2)/norm(dx2);
    err32=100*norm(J3*dq321a-dx3)/norm(dx3);
    err1C=100*norm(J1*dqC-dx1)/norm(dx1);
    err2C=100*norm(J2*dqC-dx2)/norm(dx2);
    err3C=100*norm(J3*dqC-dx3)/norm(dx3);   
    err1prref=100*norm(J1*dq_prref-dx1)/norm(dx1);
    err2prref=100*norm(J2*dq_prref-dx2)/norm(dx2);
    err3prref=100*norm(J3*dq_prref-dx3)/norm(dx3);
    
   E11=[E11;err11];
   E21=[E21;err21];
   E31=[E31;err31];


   
     E12=[E12;err12];
    E22=[E22;err22];
     E32=[E32;err32];

    E1C=[E1C;err1C];
    E2C=[E2C;err2C];
    E3C=[E3C;err3C];


    E1prref=[E1prref,err1prref];
    E2prref=[E2prref,err2prref];
    E3prref=[E3prref,err3prref];
    
    
    
    i=i+1;
    
    if (mod(i,10*floor(n/100))==0)
        perc=(i/n)*100
    end
    
end


clc
disp('STANDARD')
avg_std=[mean(E11) mean(E21) mean(E31)]
std_std=[std(E11) std(E21) std(E31)]
max_std=[max(E11) max(E21) max(E31)]
maxTime=max(TIME1(2:end))


disp('RP')
avg_RP=[mean(E12) mean(E22) mean(E32)]
std_RP=[std(E12) std(E22) std(E32)]
max_RP=[max(E12) max(E22) max(E32)]
maxTime=max(TIME2(2:end))

 
%  disp('CHIAVERINI')
%  [max(E1C) max(E2C) max(E3C)]
%  maxTime=max(TIMEC(2:end))

disp('PRREF')
avg_proposed=[mean(E1prref) mean(E2prref) mean(E3prref)]
std_proposed=[std(E1prref) std(E2prref) std(E3prref)]
max_proposed=[max(E1prref) max(E2prref) max(E3prref)]
maxTime=max(TIMEprref(2:end))

