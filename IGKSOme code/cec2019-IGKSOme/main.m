clc
clear
tic
N=50;
runs=20;
hanshu.pj=[];
hanshu.dd=[];
IIGKSOme_huatu1= repmat(hanshu, 10, 1);
for i =1:10
    F=i;
    for lun=1:runs
        [lb,ub,dim,fobj] = Get_Functions_cec2019(F);       
%         [best_fit,GKSO_curve,X]=cheng14(N,Max_iter,lb,ub,dim,fobj);
         [bestfit,Bestcost1,Bestcost2,pop,tk]=cheng42(fobj,lb,ub,dim,N);
         fit_1(lun)=bestfit;
        hanshupj(lun,:)=Bestcost1;
        hanshudd(lun,:)=Bestcost2;
        diedaihui_IIGKSOme1(i,:)=Bestcost1;
    end
    IIGKSOme_huatu1(i).pj=hanshupj;
    IIGKSOme_huatu1(i).dd=hanshudd;
%     bbb(t_b).hanshu=ceshihanshu(i1);bbb(t_b).jilu=biaoji;
%     t_b=t_b+1;
    jilulunci(i,:)=fit_1;    
    fit_ave(i)=mean(fit_1)-1;
    fit_min(i)=min(fit_1);
    fit_std(i)=std(fit_1);
    fprintf('%d       %d\n',fit_ave(i),fit_std(i));
 end 
t1=toc 
% diedaihui_MTVSCA=huatu.pj;
save('IIGKSOme_huatu1.mat','diedaihui_IIGKSOme1')