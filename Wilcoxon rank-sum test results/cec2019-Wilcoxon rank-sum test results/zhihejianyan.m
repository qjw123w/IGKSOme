clear all
clc
A=importdata('IGKSO2.txt');
B=importdata('IGKSOme.txt');
for i=1:10
   zhihe(i,1)=ranksum(A(i,:),B(i,:)); 
end