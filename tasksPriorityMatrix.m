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
 

function F=tasksPriorityMatrix(bF,idx,lambda,eps)

	if nargin < 4	
    	if nargin < 3	
            if nargin < 2
                error('Not enough inputs')
            end
            lambda=1e-12; 
        end
        eps=1e-8; 
	end
		
	m=size(bF,1);
	
	i=1;
	i_idx=1;
	while (i<=m)
		j=i+idx(i_idx)-1;
		
		pR=damped_pinv(bF(i:j,i:j),lambda,eps);        
        bF(i:j,i:end)=pR*bF(i:j,i:end); %first step
				
        
		bF(1:i-1,i:end)=bF(1:i-1,i:end)...
		  -bF(1:i-1,i:j)*bF(i:j,i:end); %second step
				
		i=j+1;
		i_idx=i_idx+1;
	end
	F=bF';
end

