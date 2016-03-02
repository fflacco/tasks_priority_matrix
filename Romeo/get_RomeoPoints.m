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


function [RL_HIP,HIP_LL,HIP_RH,HIP_LH,COM]=get_RomeoPoints(q)

	HIP_DHTable =[ pi/2, 0, 0.0673, pi/2;...
		   	pi/2, 0, 0, q(1)+pi/2;...
		   	0, 0.29, 0, q(2);...
			0, 0.32, 0, -q(3)];
		   		
        A0 = DH(HIP_DHTable(1,1), HIP_DHTable(1,2), HIP_DHTable(1,3), HIP_DHTable(1,4));
    	A1 = DH(HIP_DHTable(2,1), HIP_DHTable(2,2), HIP_DHTable(2,3), HIP_DHTable(2,4));
    	A2 = DH(HIP_DHTable(3,1), HIP_DHTable(3,2), HIP_DHTable(3,3), HIP_DHTable(3,4));
    	A3 = DH(HIP_DHTable(4,1), HIP_DHTable(4,2), HIP_DHTable(4,3), HIP_DHTable(4,4));
  

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
        	
        A4 = DH(RLEG_LLEG_DHTable(5,1), RLEG_LLEG_DHTable(5,2), RLEG_LLEG_DHTable(5,3), RLEG_LLEG_DHTable(5,4));
    	A5 = DH(RLEG_LLEG_DHTable(6,1), RLEG_LLEG_DHTable(6,2), RLEG_LLEG_DHTable(6,3), RLEG_LLEG_DHTable(6,4));
    	A6 = DH(RLEG_LLEG_DHTable(7,1), RLEG_LLEG_DHTable(7,2), RLEG_LLEG_DHTable(7,3), RLEG_LLEG_DHTable(7,4));
    	A7 = DH(RLEG_LLEG_DHTable(8,1), RLEG_LLEG_DHTable(8,2), RLEG_LLEG_DHTable(8,3), RLEG_LLEG_DHTable(8,4));
    	A8 = DH(RLEG_LLEG_DHTable(9,1), RLEG_LLEG_DHTable(9,2), RLEG_LLEG_DHTable(9,3), RLEG_LLEG_DHTable(9,4));
    	A9 = DH(RLEG_LLEG_DHTable(10,1), RLEG_LLEG_DHTable(10,2), RLEG_LLEG_DHTable(10,3), RLEG_LLEG_DHTable(10,4));
    	A10 = DH(RLEG_LLEG_DHTable(11,1), RLEG_LLEG_DHTable(11,2), RLEG_LLEG_DHTable(11,3), RLEG_LLEG_DHTable(11,4));
    	A11 = DH(RLEG_LLEG_DHTable(12,1), RLEG_LLEG_DHTable(12,2), RLEG_LLEG_DHTable(12,3), RLEG_LLEG_DHTable(12,4));
    	A12 = DH(RLEG_LLEG_DHTable(13,1), RLEG_LLEG_DHTable(13,2), RLEG_LLEG_DHTable(13,3), RLEG_LLEG_DHTable(13,4));
        
        J1=A0;
        J2=J1*A1;
        J3=J2*A2;
        J4=J3*A3;
        
        RL_HIP=[[0;0;0],J1(1:3,4),J2(1:3,4),J3(1:3,4),J4(1:3,4)];
        
        J5=J4*A4;
        J6=J5*A5;
        J7=J6*A6;
        J8=J7*A7;
        J9=J8*A8;
        J10=J9*A9;
        J11=J10*A10;
        J12=J11*A11;
        J13=J12*A12;
        
        HIP_LL=[J4(1:3,4),J5(1:3,4),J6(1:3,4),J7(1:3,4),J8(1:3,4),J9(1:3,4),J10(1:3,4),J11(1:3,4),J12(1:3,4),J13(1:3,4)];
        
        J_bridge=J7;
        
        	RLEG_RHAND_DHTable=[pi/2, 0, 0.0673, pi/2;...
		   		pi/2, 0, 0, q(1)+pi/2;...
		   		0, 0.29, 0, q(2);...
		   		0, 0.32, 0, -q(3);...
		   		-pi/2, 0, 0, q(4);...
		   		-pi/2, 0, 0, q(5)-pi/2;...
		   		0, 0.096, 0.4, q(6);...
		   		pi/2, 0, 0, q(7)-pi/2;...
		   		pi, 0.005, 0.19, 0;...
		   		pi/2, 0, 0, q(8);...
		   		-pi/2, 0, 0, q(9)-pi/2;...
		   		0, 0, 0.205, 0;...
		   		pi/2, 0, 0, q(10);...
		   		-pi/2, 0, 0, q(11);...
		   		0, 0, 0.1823, pi];

      
    	A4 = DH(RLEG_RHAND_DHTable(5,1), RLEG_RHAND_DHTable(5,2), RLEG_RHAND_DHTable(5,3), RLEG_RHAND_DHTable(5,4));
    	A5 = DH(RLEG_RHAND_DHTable(6,1), RLEG_RHAND_DHTable(6,2), RLEG_RHAND_DHTable(6,3), RLEG_RHAND_DHTable(6,4));
    	A6 = DH(RLEG_RHAND_DHTable(7,1), RLEG_RHAND_DHTable(7,2), RLEG_RHAND_DHTable(7,3), RLEG_RHAND_DHTable(7,4));
    	A7 = DH(RLEG_RHAND_DHTable(8,1), RLEG_RHAND_DHTable(8,2), RLEG_RHAND_DHTable(8,3), RLEG_RHAND_DHTable(8,4));
    	A8 = DH(RLEG_RHAND_DHTable(9,1), RLEG_RHAND_DHTable(9,2), RLEG_RHAND_DHTable(9,3), RLEG_RHAND_DHTable(9,4));
    	A9 = DH(RLEG_RHAND_DHTable(10,1), RLEG_RHAND_DHTable(10,2), RLEG_RHAND_DHTable(10,3), RLEG_RHAND_DHTable(10,4));
    	A10 = DH(RLEG_RHAND_DHTable(11,1), RLEG_RHAND_DHTable(11,2), RLEG_RHAND_DHTable(11,3), RLEG_RHAND_DHTable(11,4));
    	A11 = DH(RLEG_RHAND_DHTable(12,1), RLEG_RHAND_DHTable(12,2), RLEG_RHAND_DHTable(12,3), RLEG_RHAND_DHTable(12,4));
    	A12 = DH(RLEG_RHAND_DHTable(13,1), RLEG_RHAND_DHTable(13,2), RLEG_RHAND_DHTable(13,3), RLEG_RHAND_DHTable(13,4));
    	A13 = DH(RLEG_RHAND_DHTable(14,1), RLEG_RHAND_DHTable(14,2), RLEG_RHAND_DHTable(14,3), RLEG_RHAND_DHTable(14,4));
    	A14 = DH(RLEG_RHAND_DHTable(15,1), RLEG_RHAND_DHTable(15,2), RLEG_RHAND_DHTable(15,3), RLEG_RHAND_DHTable(15,4));
        
        J5=J4*A4;
        J6=J5*A5;
        J7=J6*A6;
        J8=J7*A7;
        J9=J8*A8;
        J10=J9*A9;
        J11=J10*A10;
        J12=J11*A11;
        J13=J12*A12;
        J14=J13*A13;
        J15=J14*A14;
        
        HIP_RH=[J4(1:3,4),J5(1:3,4),J6(1:3,4),J7(1:3,4),J8(1:3,4),J9(1:3,4),J10(1:3,4),J11(1:3,4),J12(1:3,4),J13(1:3,4),J14(1:3,4),J15(1:3,4)];
        
        
        	RLEG_LHAND_DHTable=[ pi/2, 0, 0.0673, pi/2;...
		   		pi/2, 0, 0, q(1)+pi/2;...
		   		0, 0.29, 0, q(2);...
		   		0, 0.32, 0, -q(3);...
		   		-pi/2, 0, 0, q(4);...
		   		-pi/2, 0, 0, q(5)-pi/2;...
		   		0, 0.096, 0.4, q(6);...
		   		-pi/2, 0, 0, q(7)-pi/2;...
		   		0, 0.005, 0.19, 0;...
		   		pi/2, 0, 0, q(18);...
		   		-pi/2, 0, 0, q(19)-pi/2;...
		   		0, 0, 0.205, 0;...
		   		pi/2, 0, 0, q(20);...
		   		-pi/2, 0, 0, q(21);...
		   		0, 0, 0.1823, pi];


      
    	A5 = DH(RLEG_LHAND_DHTable(6,1), RLEG_LHAND_DHTable(6,2), RLEG_LHAND_DHTable(6,3), RLEG_LHAND_DHTable(6,4));
    	A6 = DH(RLEG_LHAND_DHTable(7,1), RLEG_LHAND_DHTable(7,2), RLEG_LHAND_DHTable(7,3), RLEG_LHAND_DHTable(7,4));
    	A7 = DH(RLEG_LHAND_DHTable(8,1), RLEG_LHAND_DHTable(8,2), RLEG_LHAND_DHTable(8,3), RLEG_LHAND_DHTable(8,4));
    	A8 = DH(RLEG_LHAND_DHTable(9,1), RLEG_LHAND_DHTable(9,2), RLEG_LHAND_DHTable(9,3), RLEG_LHAND_DHTable(9,4));
    	A9 = DH(RLEG_LHAND_DHTable(10,1), RLEG_LHAND_DHTable(10,2), RLEG_LHAND_DHTable(10,3), RLEG_LHAND_DHTable(10,4));
    	A10 = DH(RLEG_LHAND_DHTable(11,1), RLEG_LHAND_DHTable(11,2), RLEG_LHAND_DHTable(11,3), RLEG_LHAND_DHTable(11,4));
    	A11 = DH(RLEG_LHAND_DHTable(12,1), RLEG_LHAND_DHTable(12,2), RLEG_LHAND_DHTable(12,3), RLEG_LHAND_DHTable(12,4));
    	A12 = DH(RLEG_LHAND_DHTable(13,1), RLEG_LHAND_DHTable(13,2), RLEG_LHAND_DHTable(13,3), RLEG_LHAND_DHTable(13,4));
    	A13 = DH(RLEG_LHAND_DHTable(14,1), RLEG_LHAND_DHTable(14,2), RLEG_LHAND_DHTable(14,3), RLEG_LHAND_DHTable(14,4));
    	A14 = DH(RLEG_LHAND_DHTable(15,1), RLEG_LHAND_DHTable(15,2), RLEG_LHAND_DHTable(15,3), RLEG_LHAND_DHTable(15,4));

        %J6=J5*A5;
        J7=J6*A6;
        J8=J7*A7;
        J9=J8*A8;
        J10=J9*A9;
        J11=J10*A10;
        J12=J11*A11;
        J13=J12*A12;
        J14=J13*A13;
        J15=J14*A14;
        
        HIP_LH=[J_bridge(1:3,4),J7(1:3,4),J8(1:3,4),J9(1:3,4),J10(1:3,4),J11(1:3,4),J12(1:3,4),J13(1:3,4),J14(1:3,4),J15(1:3,4)];
    % HIP_LH=[J_bridge(1:3,4),J9(1:3,4),J10(1:3,4),J11(1:3,4),J12(1:3,4),J13(1:3,4),J14(1:3,4),J15(1:3,4)];
    
        NECK_COM_DHTable=[ 0, 0, -0.25, 0];
    	ACoM = DH(NECK_COM_DHTable(1,1), NECK_COM_DHTable(1,2), NECK_COM_DHTable(1,3), NECK_COM_DHTable(1,4));
        JCoM=J7*ACoM;
        
        
        COM=JCoM(1:3,4);

end