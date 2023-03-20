% calc_Rd_2_0.m

close all; clear all; clc

syms q2 qd2

R_1_0 = [ 0 1 0
         -1 0 0
          0 0 1]
      
R_2_1 = [cos(q2) -sin(q2) 0
         sin(q2)  cos(q2) 0
               0        0 1]
           
R_2_0 = R_1_0 * R_2_1

diff_R_2_0 = diff(R_2_0,q2)

Rd_2_0 = diff_R_2_0 * qd2

S_w_2 = simplify(Rd_2_0 * transpose(R_2_0))