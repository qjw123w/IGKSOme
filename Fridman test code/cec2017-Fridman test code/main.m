clc
clear
% A=importdata('TTHHO.txt');
% % A=A';
% B=importdata('GKSO.txt');
% % B=B';
% C=importdata('IGKSOme.txt');
% % C=C';
% D=importdata('IGKSO.txt');
% % D=D';
% E=importdata('IGKSO2.txt');
% % E=E';
% F=importdata('MGKSO.txt');
% % F=F';
% H=importdata('MDBO.txt');
% % H=H';
% I=importdata('ALA.txt');
% % I=I';
% J=importdata('MTV-SCA.txt');
% % J=J';
% G=[C B I D E F J H A];
G=importdata('zong.txt');
for i=1:29
% % 按从小到大排序并获取索引
[~, sorted_indices(i,:)] = sort(G(i,:),'ascend');
end
rank_order = zeros(size(G));
for i=1:29
for j = 1:9
    rank_order(i,sorted_indices(i,j)) = j;
end
end
for j=1:9
    mean_index(j)=mean(rank_order(:,j));
end
[~, original_indices1] = sort(mean_index);
