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

function [J]=getJacobianNR(q,l)

    n=length(q);

    Q=0;
    for i=1:n
        Q=Q+q(i);
    end
    
    prevS=0;
    prevC=0;
    
    for j=0:n-1;
        J(1,n-j)=-l(n-j)*sin(Q) + prevS;
        J(2,n-j)=l(n-j)*cos(Q) + prevC;
        
        prevS=J(1,n-j);
        prevC=J(2,n-j);
        
        Q=Q-q(n-j);
        
    end
    




end