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


function X=get_RLEG_HIP_POSITION(q)


	HIP_DHTable =[ pi/2, 0, 0.0673, pi/2;...
		   	pi/2, 0, 0, q(1)+pi/2;...
		   	0, 0.29, 0, q(2);...
			0, 0.32, 0, -q(3)];
		   		

  

      
        A0 = DH(HIP_DHTable(1,1), HIP_DHTable(1,2), HIP_DHTable(1,3), HIP_DHTable(1,4));
    	A1 = DH(HIP_DHTable(2,1), HIP_DHTable(2,2), HIP_DHTable(2,3), HIP_DHTable(2,4));
    	A2 = DH(HIP_DHTable(3,1), HIP_DHTable(3,2), HIP_DHTable(3,3), HIP_DHTable(3,4));
    	A3 = DH(HIP_DHTable(4,1), HIP_DHTable(4,2), HIP_DHTable(4,3), HIP_DHTable(4,4));
  
    	A_EE = A0*A1*A2*A3;

	X = A_EE(1:3,4);

            
end