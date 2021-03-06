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


function X=get_RLEG_LLEG_POSITION(q)

	RLEG_LLEG_DHTable=[ pi/2, 0, 0.0673, pi/2;...
		   		pi/2, 0, 0, q(1)+pi/2;...
		   		0, 0.29, 0, q(2);...
		   		0, 0.32, 0, -q(3);...
		   		-pi/2, 0, 0, q(4);...
		   		-pi/2, 0, 0, q(5)-pi/2;...
		   		pi, 0.192, 0, q(6);...
		   		-pi/2, 0, 0, q(12);...
		   		-pi/2, 0, 0, q(13)-pi/2;...
		   		0, 0.32, 0, q(14);...
		   		0, 0.29, 0, q(15);...
		   		pi/2, 0, 0, q(16);...
		   		0, 0.0673, 0, q(17)];


  

      
        A0 = DH(RLEG_LLEG_DHTable(1,1), RLEG_LLEG_DHTable(1,2), RLEG_LLEG_DHTable(1,3), RLEG_LLEG_DHTable(1,4));
    	A1 = DH(RLEG_LLEG_DHTable(2,1), RLEG_LLEG_DHTable(2,2), RLEG_LLEG_DHTable(2,3), RLEG_LLEG_DHTable(2,4));
    	A2 = DH(RLEG_LLEG_DHTable(3,1), RLEG_LLEG_DHTable(3,2), RLEG_LLEG_DHTable(3,3), RLEG_LLEG_DHTable(3,4));
    	A3 = DH(RLEG_LLEG_DHTable(4,1), RLEG_LLEG_DHTable(4,2), RLEG_LLEG_DHTable(4,3), RLEG_LLEG_DHTable(4,4));
    	A4 = DH(RLEG_LLEG_DHTable(5,1), RLEG_LLEG_DHTable(5,2), RLEG_LLEG_DHTable(5,3), RLEG_LLEG_DHTable(5,4));
    	A5 = DH(RLEG_LLEG_DHTable(6,1), RLEG_LLEG_DHTable(6,2), RLEG_LLEG_DHTable(6,3), RLEG_LLEG_DHTable(6,4));
    	A6 = DH(RLEG_LLEG_DHTable(7,1), RLEG_LLEG_DHTable(7,2), RLEG_LLEG_DHTable(7,3), RLEG_LLEG_DHTable(7,4));
    	A7 = DH(RLEG_LLEG_DHTable(8,1), RLEG_LLEG_DHTable(8,2), RLEG_LLEG_DHTable(8,3), RLEG_LLEG_DHTable(8,4));
    	A8 = DH(RLEG_LLEG_DHTable(9,1), RLEG_LLEG_DHTable(9,2), RLEG_LLEG_DHTable(9,3), RLEG_LLEG_DHTable(9,4));
    	A9 = DH(RLEG_LLEG_DHTable(10,1), RLEG_LLEG_DHTable(10,2), RLEG_LLEG_DHTable(10,3), RLEG_LLEG_DHTable(10,4));
    	A10 = DH(RLEG_LLEG_DHTable(11,1), RLEG_LLEG_DHTable(11,2), RLEG_LLEG_DHTable(11,3), RLEG_LLEG_DHTable(11,4));
    	A11 = DH(RLEG_LLEG_DHTable(12,1), RLEG_LLEG_DHTable(12,2), RLEG_LLEG_DHTable(12,3), RLEG_LLEG_DHTable(12,4));
    	AEE = DH(RLEG_LLEG_DHTable(13,1), RLEG_LLEG_DHTable(13,2), RLEG_LLEG_DHTable(13,3), RLEG_LLEG_DHTable(13,4));
  
    	A_EE = A0*A1*A2*A3*A4*A5*A6*A7*A8*A9*A10*A11*AEE;

	X = A_EE(1:3,4);

            
end