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

close all


set(0, 'defaultTextInterpreter', 'none'); 
hfig=figure
set(gcf,'Position',[0 0 800 1000])




q=QV;
samples=size(q,1);


T=0.001;
tt=0:T:samples*T-T;


for i=1:samples
    rhand(i,:)=get_RLEG_RHAND_POSITION(q(i,:))';
    lhand(i,:)=get_RLEG_LHAND_POSITION(q(i,:))';
    lleg(i,:)=get_RLEG_LLEG_POSITION(q(i,:))';
    com(i,:)=get_RLEG_COM_POSITION(q(i,:))';
    
    error_com(i)=norm(T0d(i,1:2)-com(i,1:2));
    error_lleg(i)=norm(T1d(i,:)-lleg(i,:));
    error_lhand(i)=norm(T2d(i,:)-lhand(i,:));
    error_rhand(i)=norm(T3d(i,:)-rhand(i,:));
%    error_posture(i)=norm(q(1,:)-q(i,:));
end




subplot(3,1,1:2)
[RL_HIP,HIP_LL,HIP_RH,HIP_LH,COM]=get_RomeoPoints(q(1,:));

RL_HIP_X=RL_HIP(1,:);
RL_HIP_Y=RL_HIP(2,:);
RL_HIP_Z=RL_HIP(3,:);
hPlotRobot1=plot3(RL_HIP_X,RL_HIP_Y,RL_HIP_Z,'o-k','LineWidth',3);
set(hPlotRobot1,'XDataSource','RL_HIP_X')
set(hPlotRobot1,'YDataSource','RL_HIP_Y')
set(hPlotRobot1,'ZDataSource','RL_HIP_Z')

hold on

HIP_LL_X=HIP_LL(1,:);
HIP_LL_Y=HIP_LL(2,:);
HIP_LL_Z=HIP_LL(3,:);
hPlotRobot2=plot3(HIP_LL_X,HIP_LL_Y,HIP_LL_Z,'o-b','LineWidth',3);
set(hPlotRobot2,'XDataSource','HIP_LL_X')
set(hPlotRobot2,'YDataSource','HIP_LL_Y')
set(hPlotRobot2,'ZDataSource','HIP_LL_Z')

HIP_RH_X=HIP_RH(1,:);
HIP_RH_Y=HIP_RH(2,:);
HIP_RH_Z=HIP_RH(3,:);
hPlotRobot3=plot3(HIP_RH_X,HIP_RH_Y,HIP_RH_Z,'o-r','LineWidth',3);
set(hPlotRobot3,'XDataSource','HIP_RH_X')
set(hPlotRobot3,'YDataSource','HIP_RH_Y')
set(hPlotRobot3,'ZDataSource','HIP_RH_Z')

HIP_LH_X=HIP_LH(1,:);
HIP_LH_Y=HIP_LH(2,:);
HIP_LH_Z=HIP_LH(3,:);
hPlotRobot4=plot3(HIP_LH_X,HIP_LH_Y,HIP_LH_Z,'o-m','LineWidth',3);
set(hPlotRobot4,'XDataSource','HIP_LH_X')
set(hPlotRobot4,'YDataSource','HIP_LH_Y')
set(hPlotRobot4,'ZDataSource','HIP_LH_Z')

CoM_X=COM(1,:);
CoM_Y=COM(2,:);
CoM_Z=COM(3,:);
hPlotRobot5=plot3(CoM_X,CoM_Y,CoM_Z,'xk','MarkerSize',8,'LineWidth',3);
set(hPlotRobot5,'XDataSource','CoM_X')
set(hPlotRobot5,'YDataSource','CoM_Y')
set(hPlotRobot5,'ZDataSource','CoM_Z')


RH_X=rhand(1,1);
RH_Y=rhand(1,2);
RH_Z=rhand(1,3);
hPlotRH=plot3(RH_X,RH_Y,RH_Z,'r','LineWidth',1);
set(hPlotRH,'XDataSource','RH_X')
set(hPlotRH,'YDataSource','RH_Y')
set(hPlotRH,'ZDataSource','RH_Z')

 LH_X=lhand(1,1);
 LH_Y=lhand(1,2);
 LH_Z=lhand(1,3);
 hPlotLH=plot3(LH_X,LH_Y,LH_Z,'m','LineWidth',1);
 set(hPlotLH,'XDataSource','LH_X')
 set(hPlotLH,'YDataSource','LH_Y')
 set(hPlotLH,'ZDataSource','LH_Z')

LL_X=lleg(1,1);
LL_Y=lleg(1,2);
LL_Z=lleg(1,3);
hPlotLL=plot3(LL_X,LL_Y,LL_Z,'b','LineWidth',1);
set(hPlotLL,'XDataSource','LL_X')
set(hPlotLL,'YDataSource','LL_Y')
set(hPlotLL,'ZDataSource','LL_Z')

% H_X=hip(1,1);
% H_Y=hip(1,2);
% H_Z=hip(1,3);
% hPlotH=plot3(H_X,H_Y,H_Z,'k');
% set(hPlotH,'XDataSource','H_X')
% set(hPlotH,'YDataSource','H_Y')
% set(hPlotH,'ZDataSource','H_Z')

C_X=com(1,1);
C_Y=com(1,2);
C_Z=com(1,3);
hPlotC=plot3(C_X,C_Y,C_Z,'k','LineWidth',1);
set(hPlotC,'XDataSource','C_X')
set(hPlotC,'YDataSource','C_Y')
set(hPlotC,'ZDataSource','C_Z')

grid on
xlabel('X','FontSize',10);
ylabel('Y','FontSize',10);
zlabel('Z','FontSize',10);

axis equal

az=130;
el=40;
view([az,el]);
xlim([0,0.8])


subplot(3,1,3)
TT=tt(1);
E_COM=error_com(1);
E_LLEG=error_lleg(1);
E_RHAND=error_rhand(1);
E_LHAND=error_lhand(1);
hplotEL=plot(TT,E_LLEG,'b','LineWidth',1.2);
hold on
set(hplotEL,'XDataSource','TT');
set(hplotEL,'YDataSource','E_LLEG')
hplotEH=plot(TT,E_RHAND,'r','LineWidth',1.2);
set(hplotEH,'XDataSource','TT');
set(hplotEH,'YDataSource','E_RHAND')
hplotELH=plot(TT,E_LHAND,'m','LineWidth',1.2);
set(hplotELH,'XDataSource','TT');
set(hplotELH,'YDataSource','E_LHAND')
hplotEC=plot(TT,E_COM,'k','LineWidth',1.2);
set(hplotEC,'XDataSource','TT');
set(hplotEC,'YDataSource','E_COM')

grid on
xlabel('time [s]')
s=sprintf('task errors \n[m,10*rad]');
ylabel(s)
s=sprintf('right hand ');
legend('left leg',s,'left hand')
set(gca,'LineWidth',1.2,'FontSize',10)
%ylim([0,0.3])
xlim([0,Ttot]);

f=1;
for i=1:samples
     
    if (mod(i,20)==1)  
        
        [RL_HIP,HIP_LL,HIP_RH,HIP_LH,COM]=get_RomeoPoints(q(i,:));
        RL_HIP_X=RL_HIP(1,:);
        RL_HIP_Y=RL_HIP(2,:);
        RL_HIP_Z=RL_HIP(3,:);
        refreshdata(hPlotRobot1)

        HIP_LL_X=HIP_LL(1,:);
        HIP_LL_Y=HIP_LL(2,:);
        HIP_LL_Z=HIP_LL(3,:);
        refreshdata(hPlotRobot2)

        HIP_RH_X=HIP_RH(1,:);
        HIP_RH_Y=HIP_RH(2,:);
        HIP_RH_Z=HIP_RH(3,:);
        refreshdata(hPlotRobot3)

        HIP_LH_X=HIP_LH(1,:);
        HIP_LH_Y=HIP_LH(2,:);
        HIP_LH_Z=HIP_LH(3,:);
        refreshdata(hPlotRobot4)
        
        CoM_X=COM(1,:);
        CoM_Y=COM(2,:);
        CoM_Z=COM(3,:);
        refreshdata(hPlotRobot5)

        RH_X=rhand(1:i,1);
        RH_Y=rhand(1:i,2);
        RH_Z=rhand(1:i,3);
        refreshdata(hPlotRH)
        
        

         LH_X=lhand(1:i,1);
         LH_Y=lhand(1:i,2);
         LH_Z=lhand(1:i,3);
         refreshdata(hPlotLH)

        LL_X=lleg(1:i,1);
        LL_Y=lleg(1:i,2);
        LL_Z=lleg(1:i,3);
        refreshdata(hPlotLL)

%         H_X=hip(1:i,1);
%         H_Y=hip(1:i,2);
%         H_Z=hip(1:i,3);
%         refreshdata(hPlotH)   
        
        C_X=com(1:i,1);
        C_Y=com(1:i,2);
        C_Z=com(1:i,3);
        refreshdata(hPlotC) 
       
        subplot(3,1,1:2)
        axis equal
       
        TT=tt(1:i);

 
        E_COM=error_com(1:i);
        E_LLEG=error_lleg(1:i);
        E_RHAND=error_rhand(1:i);
        E_LHAND=error_lhand(1:i);
        refreshdata(hplotEL)        
        refreshdata(hplotEH)        
        refreshdata(hplotELH)        
        refreshdata(hplotEC)        
        
        
       % az=az+0.5;
       % view([az,el]);
      %  pause(1e-3);
        M(f)=getframe(hfig);
       %pause(1e-5);
        f=f+1;
    end

end
disp('VIDEO PERFORMED')
%movie2avi(M, 'videoRomeoTest1.avi','fps',25, 'compression', 'none');
%close all


