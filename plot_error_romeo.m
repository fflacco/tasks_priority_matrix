set(0, 'defaultTextInterpreter', 'latex'); 

hfig=figure(1);
clf
cmToInch=0.393701;
resolution=200;
w=4*cmToInch*resolution
h=2*cmToInch*resolution

set(hfig,'Position',[0 0 w h])

q=QV;
samples=size(q,1);


T=0.001;
tt=0:T:samples*T-T;

error_com=[];
error_lleg=[];
error_rhand=[];
error_lhand=[];

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




plot(tt,error_com,'k','LineWidth',2);
hold on
plot(tt,error_lleg,'b','LineWidth',2);


plot(tt,error_lhand,'m','LineWidth',2);
plot(tt,error_rhand,'r','LineWidth',2);

grid on
xlabel('time [s]')
s=sprintf('task errors [m]');
ylabel(s)
s=sprintf('right hand ');
legend('COM','left foot','left hand',s)
set(gca,'LineWidth',1.2,'FontSize',10)
%ylim([0,0.3])
xlim([0,Ttot]);