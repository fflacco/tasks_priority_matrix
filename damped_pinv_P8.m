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



function [invJ,Q,fullRank]=damped_pinv_P8(J)
    lambda_max=1e-12;
    epsilon=1e-8;
    epsilonQ=1e-20;
    
    r=0;
    m=size(J,1);
    n=size(J,2);
    m=min(n,m);
    if (m==0)
        Q=zeros(n);
        invJ=J'*0;
    else
    
        [U,S,V]=svd(J',0);
        sigma=diag(S);

        if (sigma(m)>epsilon)
            iS=inv(S);     
            invJ=U*iS*V';
            Q=U*U';
            fullRank=1;        
        else
            %disp('d');
            lambda2=(1-(sigma(m)/epsilon)*(sigma(m)/epsilon))*lambda_max*lambda_max;
            for i=1:m
                if (sigma(i)> epsilonQ) 
                    r=r+1;
                end
                sigma(i)=(sigma(i)/(sigma(i)*sigma(i)+lambda2));
            end
            invJ=U(:,1:r)*diag(sigma(1:r))*V(:,1:r)';
            if (r>0)
                subU=U(:,1:r);
            else
                subU=zeros(n,1);
            end
            Q=subU*subU';
            fullRank=0;
        end
    end

    %invJ=pinv(J,1e-8);

end



