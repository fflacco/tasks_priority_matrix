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


function jacobian=get_RLEG_LHAND_JACOBIAN(q)


    RLEG_LHAND_DHTable =[ pi/2, 0, 0.0673, pi/2;...
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

	% Transformation matrices
    	A0 = DH(RLEG_LHAND_DHTable(1,1), RLEG_LHAND_DHTable(1,2), RLEG_LHAND_DHTable(1,3), RLEG_LHAND_DHTable(1,4));
    	A1 = DH(RLEG_LHAND_DHTable(2,1), RLEG_LHAND_DHTable(2,2), RLEG_LHAND_DHTable(2,3), RLEG_LHAND_DHTable(2,4));
    	A2 = DH(RLEG_LHAND_DHTable(3,1), RLEG_LHAND_DHTable(3,2), RLEG_LHAND_DHTable(3,3), RLEG_LHAND_DHTable(3,4));
    	A3 = DH(RLEG_LHAND_DHTable(4,1), RLEG_LHAND_DHTable(4,2), RLEG_LHAND_DHTable(4,3), RLEG_LHAND_DHTable(4,4));
    	A4 = DH(RLEG_LHAND_DHTable(5,1), RLEG_LHAND_DHTable(5,2), RLEG_LHAND_DHTable(5,3), RLEG_LHAND_DHTable(5,4));
    	A5 = DH(RLEG_LHAND_DHTable(6,1), RLEG_LHAND_DHTable(6,2), RLEG_LHAND_DHTable(6,3), RLEG_LHAND_DHTable(6,4));
    	A6 = DH(RLEG_LHAND_DHTable(7,1), RLEG_LHAND_DHTable(7,2), RLEG_LHAND_DHTable(7,3), RLEG_LHAND_DHTable(7,4));
    	A7 = DH(RLEG_LHAND_DHTable(8,1), RLEG_LHAND_DHTable(8,2), RLEG_LHAND_DHTable(8,3), RLEG_LHAND_DHTable(8,4));
    	A8 = DH(RLEG_LHAND_DHTable(9,1), RLEG_LHAND_DHTable(9,2), RLEG_LHAND_DHTable(9,3), RLEG_LHAND_DHTable(9,4));
    	A9 = DH(RLEG_LHAND_DHTable(10,1), RLEG_LHAND_DHTable(10,2), RLEG_LHAND_DHTable(10,3), RLEG_LHAND_DHTable(10,4));
    	A10 = DH(RLEG_LHAND_DHTable(11,1), RLEG_LHAND_DHTable(11,2), RLEG_LHAND_DHTable(11,3), RLEG_LHAND_DHTable(11,4));
    	A11 = DH(RLEG_LHAND_DHTable(12,1), RLEG_LHAND_DHTable(12,2), RLEG_LHAND_DHTable(12,3), RLEG_LHAND_DHTable(12,4));
    	A12 = DH(RLEG_LHAND_DHTable(13,1), RLEG_LHAND_DHTable(13,2), RLEG_LHAND_DHTable(13,3), RLEG_LHAND_DHTable(13,4));
    	A13 = DH(RLEG_LHAND_DHTable(14,1), RLEG_LHAND_DHTable(14,2), RLEG_LHAND_DHTable(14,3), RLEG_LHAND_DHTable(14,4));
    	AEE = DH(RLEG_LHAND_DHTable(15,1), RLEG_LHAND_DHTable(15,2), RLEG_LHAND_DHTable(15,3), RLEG_LHAND_DHTable(15,4));


        A_0 = A0;
    	A_1 = A0*A1;
    	A_2 = A0*A1*A2;
    	A_3 = A0*A1*A2*A3;
    	A_4 = A0*A1*A2*A3*A4;
    	A_5 = A0*A1*A2*A3*A4*A5;
    	A_6 = A0*A1*A2*A3*A4*A5*A6;
    	A_7 = A0*A1*A2*A3*A4*A5*A6*A7;
    	A_8 = A0*A1*A2*A3*A4*A5*A6*A7*A8;
    	A_9 = A0*A1*A2*A3*A4*A5*A6*A7*A8*A9;
    	A_10 = A0*A1*A2*A3*A4*A5*A6*A7*A8*A9*A10;
    	A_11 = A0*A1*A2*A3*A4*A5*A6*A7*A8*A9*A10*A11;
    	A_12 = A0*A1*A2*A3*A4*A5*A6*A7*A8*A9*A10*A11*A12;
    	A_13 = A0*A1*A2*A3*A4*A5*A6*A7*A8*A9*A10*A11*A12*A13;
    	A_EE = A0*A1*A2*A3*A4*A5*A6*A7*A8*A9*A10*A11*A12*A13*AEE;


       % Rotations matrices
        R_0 = A_0(1:3,1:3);
        R_1 = A_1(1:3,1:3);
        R_2 = A_2(1:3,1:3);
        R_3 = A_3(1:3,1:3);
        R_4 = A_4(1:3,1:3);
        R_5 = A_5(1:3,1:3);
        R_6 = A_6(1:3,1:3);
        R_7 = A_7(1:3,1:3);
        R_8 = A_8(1:3,1:3);
        R_9 = A_9(1:3,1:3);
        R_10 = A_10(1:3,1:3);
        R_11 = A_11(1:3,1:3);
        R_12 = A_12(1:3,1:3);
        R_13 = A_13(1:3,1:3);
        R_EE = A_EE(1:3,1:3);

        z0 = R_0(1:3,3);
    	z1 = R_1(1:3,3);
        z2 = R_2(1:3,3);
        z3 = R_3(1:3,3);
        z4 = R_4(1:3,3);
        z5 = R_5(1:3,3);
        z6 = R_6(1:3,3);
        z7 = R_7(1:3,3);
        z8 = R_8(1:3,3);
        z9 = R_9(1:3,3);
        z10 = R_10(1:3,3);
        z11 = R_11(1:3,3);
        z12 = R_12(1:3,3);
        z13 = R_13(1:3,3);
        zEE = R_EE(1:3,3);


        p0 = A_0(1:3,4);
        p1 = A_1(1:3,4);
        p2 = A_2(1:3,4);
        p3 = A_3(1:3,4);
        p4 = A_4(1:3,4);
        p5 = A_5(1:3,4);
        p6 = A_6(1:3,4);
        p7 = A_7(1:3,4);
        p8 = A_8(1:3,4);
        p9 = A_9(1:3,4);
        p10 = A_10(1:3,4);
        p11 = A_11(1:3,4);
        p12 = A_12(1:3,4);
        p13 = A_13(1:3,4);
        pEE = A_EE(1:3,4);


    zero=zeros(3,1);
	jacobian=[cross(z0,pEE-p0), cross(z1,pEE-p1), -cross(z2,pEE-p2), cross(z3,pEE-p3), cross(z4,pEE-p4), cross(z5,pEE-p5), cross(z6,pEE-p6), zero, zero, zero, zero, zero, zero, zero, zero, zero, zero,  cross(z8,pEE-p8), cross(z9,pEE-p9), cross(z11,pEE-p11), cross(z12,pEE-p12)];
           
end	

	%jacobian << z0.cross(pEE-p0), z1.cross(pEE-p1), z2.cross(pEE-p2), z3.cross(pEE-p3), z4.cross(pEE-p4), z5.cross(pEE-p5), z6.cross(pEE-p6), zero, zero, zero, zero, zero, zero, zero, zero, zero, zero, z8.cross(pEE-p8), z9.cross(pEE-p9), z11.cross(pEE-p11), z12.cross(pEE-p12);

