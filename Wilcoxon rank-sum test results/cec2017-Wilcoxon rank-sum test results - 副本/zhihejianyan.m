clear all
clc
A=importdata('MTV-SCA.txt');
B=importdata('IGKSOme.txt');
for i=1:30
   zhihe(i,1)=ranksum(A(i,:),B(i,:)); 
end