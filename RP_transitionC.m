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

function [dq,pJ,P,fullRank,fullRank1]=RP_transitionC(J,dx,h,dq_prev,dq_old,pJ_prev,P_prev)
    
    [m,n]=size(J);
    Im=eye(m);
    I=eye(n);
    
    
    
    if (h<1e-50)
        dq=dq_prev;
        P=P_prev;
        pJ=pJ_prev;
        fullRank1=1;
        fullRank=1;
        return;
    end
    
   %pJ_prev=pJ_prev*1e2;
    
    %Ja=h*J;
    
    if isempty(pJ_prev)
        [pJ,Q,fullRank1]=damped_pinv_P8(J);
        P=I-Q;
        dq=dq_prev + h*pJ*(dx- J*dq_prev);  
        fullRank=1;
    else
        E=J*P_prev;
        [pE,Q,fullRank]=damped_pinv_P8(E);
        P=P_prev-Q;
        %P=P_prev-pE*E;
       
        
        
      %  fullRank=0; %avoid to use P8 directly
         if fullRank
               T=pE;
               pJ=[(pJ_prev-T*J*pJ_prev), T];
   
               dq=dq_prev + h*P_prev*T*(dx- J*dq_prev);
               
               fullRank1=1;
           else
            K_inv=Im+(Im-E*pE)*J*(pJ_prev*pJ_prev')*J'*(Im-E*pE);
            T=pE + (I-pE*J)*(pJ_prev*pJ_prev')*J'*(K_inv\(Im-E*pE));

            pJ=[(pJ_prev-T*J*pJ_prev), T];
            fullRank=0;
            fullRank1=1;
            [pJT,Qq,fullRank1]=damped_pinv_P8(J*T);
        %    pJT=pinv(J*T);
            dq=dq_prev + h*T*pJT*(dx - J*dq_prev); %
        end
    end
    
    
end