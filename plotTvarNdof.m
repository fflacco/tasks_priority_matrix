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

n=30;
m=3;

l=5;




diffMAX=0;

n_tests=1000;

maxTstd=[];
maxTprref=[];
meanTstd=[];
meanTprref=[];

maxTk=[];
meanTk=[];

for n=20:20:300;
    Tstd=[];
    Tprref=[];
    Tk=[];
    
    I=eye(n);
    idx=ones(1,l)*m;

    for i=1:n_tests


        A=rand(m*l,n);
        b=rand(m*l,1);

        P=I;
        x=zeros(n,1);

        tic

        for i=1:l
         pA=damped_pinv(A(1+(i-1)*m:i*m,:)*P);
         x=x + pA*(b(1+(i-1)*m:i*m,1) - A(1+(i-1)*m:i*m,:)*x);
         P= P - pA*A(1+(i-1)*m:i*m,:)*P;
        end
        time=toc;
        Tstd=[Tstd,time];
        xstd=x;


        


        tic
        [Q,R]=qr(A',0);
        X=tasksPriorityMatrix(R,idx);
        xprref=Q*damped_pinv(R')*X*b;
        time=toc;
        Tprref=[Tprref,time];
        

        
    end
perc=0.99; %we used the percentile to remove outliers due to the operative system
    
maxTstd=[maxTstd,quantile(Tstd,perc)];
maxTprref=[maxTprref,quantile(Tprref,perc)];
meanTstd=[meanTstd,mean(Tstd)];
meanTprref=[meanTprref,mean(Tprref)];
maxTk=[maxTk,quantile(Tk,perc)];
meanTk=[meanTk,mean(Tk)];

n
end



T=20:20:300;





set(0, 'defaultTextInterpreter', 'latex'); 

hfig=figure;
cmToInch=0.393701;
resolution=200;
w=4*cmToInch*resolution
h=2*cmToInch*resolution

set(hfig,'Position',[0 0 w h])




plot(T,maxTstd*1000,'LineWidth',2);
grid on
hold on
plot(T,meanTstd*1000,'b--','LineWidth',2);
plot(T,maxTprref*1000,'r','LineWidth',2);
plot(T,meanTprref*1000,'r--','LineWidth',2);
xlabel('number od DOF')
ylabel('execution time [ms]')
xlim([20,300]);
legend('max_{std}','mean_{std}','max_{tpm}','mean_{tpm}','Location','Best')
%legend('$e^{\mathcal A}$','$e^{\mathcal B}$','Location','Best')
set(gca,'LineWidth',1.2)



set(findall(gcf,'type','text'),'fontSize',10);









    