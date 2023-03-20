% inv_3_by_3_matrix.m

close all; clear all; clc

syms Pxx Pxy Pxz Pyx Pyy Pyz Pzx Pzy Pzz
syms Cx Cy Cz

P = [Pxx Pxy Pxz
     Pyx Pyy Pyz
     Pzx Pzy Pzz]
 
 C = [Cx  0  0
       0 Cy  0
       0  0 Cz]
   
L = P * C

S = inv(eye(3) + L)