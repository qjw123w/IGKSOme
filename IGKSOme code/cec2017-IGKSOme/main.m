clc
clear
tic
% ceshihanshu=[1;2;3;4;6;9;12;15;16;17;18;19;20];
% ceshihanshu=[1;2;3;4;5;6;9;15;16;17;18];
% ceshihanshu=[1;2;3;4;6;9;15;16;17;18;19;20];
% ceshihanshu=[1 ;2 ;3 ;4;5;6 ;7;8,;9;10;12;16];
% ceshihanshu=[1;2 ;3;4;6;9;12;15;19;20];
ceshihanshu=[1:30];
run=20;
c=1;
t_b=1;
search_agent=50;
duoyu=[100 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900 ,...
    2000 2100 2200 2300 2400 2500 2600 2700 2800 2900 3000];


hanshu.pj=[];
hanshu.dd=[];
IIGKSOme_huatu= repmat(hanshu, 30, 1);
%%%%%%%%%%%%%%%%%%%%%%测试
if c==1 
    fhd=str2func('cec17_func') ;
end
for i1=1    :    length(ceshihanshu)
% for i1=6
        if i1~=2
%% % %              FF=duoyu(i);
    FF=ceshihanshu(i1)*100;
    for j1=1:run
        %             fobj=func_name(i);
        fobj=ceshihanshu(i1);
        %             [lb,ub,dim]=Get_Functions_detailsCEC(func_name(i));
        [lb,ub,dim]=Get_Functions_detailsCEC(ceshihanshu(i1));
        %         [bestfit,pop,Bestcost1,Bestcost2,FE,it]=CBSO11(fobj,fhd,lb,ub,dim,search_agent,FF);
        % [bestfit,pop,Bestcost1,Bestcost2,FE,it]=CBSO(fobj,fhd,lb,ub,dim,search_agent,FF);
        [bestfit,Bestcost1,Bestcost2,pop,t]=cheng42(fobj,fhd,lb,ub,dim,search_agent,FF);
%  [bestfit,Bestcost1,Bestcost2,pop,t]=dbo11(fobj,fhd,lb,ub,dim,search_agent,FF);
        % [bestfit,pop,Bestcost1,Bestcost2,FE,it,zongjie,zs]=CBSO_0401gai(fobj,fhd,lb,ub,dim,search_agent,FF);
        % [bestfit,pop,Bestcost1,Bestcost2,FE,it]=shishi19(fobj,fhd,lb,ub,dim,search_agent,FF);
%         biaoji(j1).zong=zongjie;
        fit_1(j1)=bestfit;
        pop2(j1).lunci=pop;%pop;
        pop2(j1).fitmin=bestfit;
        jiluit(i1,j1)=t;%it;
        hanshupj(j1,:)=Bestcost1;
        hanshudd(j1,:)=Bestcost2;
        diedaihui_IIGKSOme(i1,:)=Bestcost1;
    end
    IIGKSOme_huatu(i1).pj=hanshupj;
    IIGKSOme_huatu(i1).dd=hanshudd;
%     bbb(t_b).hanshu=ceshihanshu(i1);bbb(t_b).jilu=biaoji;
%     t_b=t_b+1;
    jilulunci(i1,:)=fit_1;    
    [~,nn]=min(fit_1);
    POP=pop2(nn).lunci;
    meilunceshizhi(1,:)=fit_1;
    fit_ave(i1)=mean(fit_1);
    fit_min(i1)=min(fit_1);
    fit_std(i1)=std(fit_1);
    jk=ceshihanshu(i1);
    fprintf('%d       %d\n',fit_ave(i1),fit_std(i1));
    biaoge(i1).name=ceshihanshu(i1);
    biaoge(i1).pinjun=fit_ave(i1);
    biaoge(i1).zuixiao=fit_min(i1);
    pophuizhong(ceshihanshu(i1)).zhongti=POP;
    pophuizhong(ceshihanshu(i1)).shiyindu=fit_min(i1);
        end 
end
t1=toc 
% diedaihui_MTVSCA=huatu.pj;
save('IIGKSOme_huatu.mat','diedaihui_IIGKSOme')



pop1.lunci=[];
pop1.fitmin=[];
pop2 = repmat(pop1, 30, 1);



% _________________________________________________________________________%
%  Mantis Search Algorithm (MSA) source codes demo 1.0               %
%                                                                         %
%  Developed in MATLAB R2019A                                      %
%                                                                         %
%  Author and programmer: Mohamed Abdel-Basset (E-mail: mohamedbasset@ieee.org) & Reda Mohamed (E-mail: redamoh@zu.edu.eg)                              %
%                                                                         %
%   Main paper: Abdel-Basset, M., Mohamed, R.                                    %
%               Mantis Search Algorithm: A novel bio-inspired algorithm for global optimization and engineering design problems,                         %
%               Computer Methods in Applied Mechanics and Engineering, in press              %
%                                                                         %
% _________________________________________________________________________%

