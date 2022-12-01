%% Function to give insight in singular value decomposition
% The svd command is there, but where is the fun in using that if you can
% build it yourself?
clear all; close all; clc

% Generate a matrix
X = magic(5)

[Ucheck,Scheck,Vcheck] = svd(X);


n = size(X,1);      % Rows
p = size(X,2);      % Columns 

% [U,S,V] = svd(X)

