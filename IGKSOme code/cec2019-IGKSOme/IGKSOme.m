function [bestfit,Bestcost1,Bestcost2,pop,tk]=cheng42(fobj,lb,ub,dim,N)
X=lb+(ub-lb).*rand(N,dim);
evaluation_count = 0;
 bestfit1=inf;
 MaxFES=100000;
for i=1:N
    fit(i,1)=feval(fobj,X(i,:));
    evaluation_count = evaluation_count + 1;
            if bestfit1> fit(i,1)
        bestfit1= fit(i,1);
    end    
    Bestcost1(evaluation_count)=bestfit1;
    %     ch(i+1)=mod(ch(i)+0.2-0.5*sin(2*pi*ch(i))/(2*pi),1);
end
[bestfit1k,best_index]=min(fit);
xbest=X(best_index,:);
Max_iter=1000;
% curve=zeros(Max_iter,1); %每一代的最优值基于代数画一个曲线
w(1)=0.1;
A=[];
% t=1;
t=0;
while evaluation_count<100000
% for t=1:Max_iter
    t = t +1;
      Bestcost2(t)=bestfit1;
    %     B=[];
    %     B1=[];
    w(t+1)=1-2*(w(t))^4;
    %         p(t)=1/(1+1.5*exp((10*t/Max_iter)-5))-0.1*rand;
    p(t)=(1-(t/Max_iter)^0.25+abs(w(t+1))*((t/Max_iter)^0.25-(t/Max_iter)^3));  %equation 27
    %    p(t)=2*(1-(t/Max_iter)^(t/Max_iter));
    %     p(t)=(1/(1+1.5*exp((10*t/Max_iter)-5)))+0.1*abs(w(t+1));
    %     p(t)=1/(1+1.5*exp((10*t/Max_iter)-5))-0.1*rand;
    %       p(t)=exp((Max_iter-t)/Max_iter)-fit(i)/(sum(fit)+eps);
    beta_min=0.2;
    beta_max=1.2;
    beta=beta_min+( beta_max-beta_min)*(1-(t/Max_iter)^3)^2;
    alpha=abs(beta.*sin(3*pi/2+sin(3*pi*beta/2)));
    %         fit1=mean(fit);
    [fit1, fit_indice] = sort(fit, 'ascend');
    %     N2=20-round(0.1*N*t/Max_iter);
    N1=10-round(0.1*N*((t/Max_iter)^2)); %equation 17
    X1=X(fit_indice(1:15),:);
    %         X3=X(fit_indice(1:N1),:);
    %     sim=zeros(size(X1,1),1);
    XU=[];XU_fit=[];
    %     for i=1:size(X1,1)
    %         sumx=0;sum1=0;sum2=0;
    %         for j=1:dim
    %             sumx=sumx+(xbest(1,j)*X1(i,j));
    %             sum1=sum1+(xbest(1,j)^2);
    %             sum2=sum2+(X1(i,j)^2);
    %         end
    %         sim(i,1)=sumx/(sqrt(sum1)*sqrt(sum2));
    %     end
    %
    %     [sort_sim,sim_index]=sort(sim);
    %     sim=[];
    %     sim_indexs=sim_index(1:N1-(round(0.6*N1-0.2*N1*t/Max_iter)));
    %     remind= setdiff(1:15,sim_indexs);
    %     sim_indexs1=remind(1:round(0.6*N1-0.2*N1*t/Max_iter));
    %     X3=[X1(sim_indexs,:);X1( sim_indexs1,:)];
    %     fit3=[fit1(sim_indexs);fit1(sim_indexs1)];
    [X3, fit3,sim_selected_idx] = selectPopulation(X, fit, t, Max_iter, N1,dim);%Fitness-Cosine Similarity Hybrid Population Partitioning
    % 找出 X 中不属于 X3 的行的索引
    [~, idxxx] = ismember(X, X3, 'rows');
    non_zero_mask = (idxxx ~= 0);
    A_non_zero = idxxx(non_zero_mask);
    original_indices = find(non_zero_mask); % 记录非零元素的原始位置
    
    % 步骤2：找到非零元素中不相邻的重复值后索引
    [~, ~, ic] = unique(A_non_zero, 'stable');
    counts = accumarray(ic, 1);
    duplicate_values = find(counts > 1);
    
    back_indices_relative = [];
    for val = duplicate_values'
        idx = find(ic == val);
        if length(idx) >= 2
            % 取每组重复值的第二个及之后的索引（相对A_non_zero的位置）
            back_indices_relative = [back_indices_relative; idx(2:end)];
        end
    end
    
    % 步骤3：映射回原始数组的索引
    back_indices_original = original_indices(back_indices_relative);
    remain_indices = find(idxxx == 0);  % 找到所有未匹配的行索引
    remain_indices=[remain_indices;back_indices_original];
    % 提取剩余个体
    X2 = X(remain_indices, :);
    fit2 = fit(remain_indices);
    X=[X3;X2];
    fit=[fit3;fit2];
    
    %     remind= setdiff(1:15,sim_index(1:round(0.5*N1)));
    %     sim_indexs=sim_index(1:round(0.5*N1));
    %     X3=[X1(sim_indexs,:);X1(remind(1:N1-round(0.5*N1)),:)];
    %     fit3=[fit1(sim_indexs);fit1(remind(1:N1-round(0.5*N1)))];
%         XU=[XU;X3];
%         XU_fit=[XU_fit;fit3];
    for i=1:N
          if mod(t,100)==0
                 new_X(i,:)=(min(X)+max(X)-X(i,:));%equation 15
               else
             
        %         if ismember(X(i,:), X3, 'rows')
        if i<=size(X3,1)
            %             lamda=2*rand(1,dim)-1;
            rk=randperm(N1,2);
            q=((max(fit)-fit(i))/(max(fit)-min(fit))).^((t/Max_iter)^2);   %%equation 20
            yita=1;
            if rand<=0.5
                betau=(rand(1,dim).*2).^(1/(1+yita));
            else
                betau=(1./(2-2.*rand(1,dim))).^(1/(1+yita));%equation 24
            end
            xbest21= 0.5*((1-betau).*X1(1,:)+(1+betau).*X1(min(sim_selected_idx),:));%equation 22
            xbest22= 0.5*((1+betau).*X1(1,:)+(1-betau).*X1(min(sim_selected_idx),:));%equation 23
            if i==fit_indice(1)
                xbest2=X(i,:);
            else
                X0=[xbest21;xbest22;X(i,:)];
                xbest2=X0(randperm(3,1),:); %equation 21
            end
            new_X(i,:)=xbest2+q.*(X3(rk(1),:)-X3(rk(2),:));  %%equation 19
%             for j=1:dim
%                 if new_X(i,j)>ub(:,j)
%                     new_X(i,j)=ub(:,j)-mod(new_X(i,j)-ub(:,j),ub(:,j)-lb(:,j));
%                 elseif new_X(i,j)<lb(:,j)
%                     new_X(i,j)=lb(:,j)+mod(lb(:,j)-new_X(i,j),ub(:,j)-lb(:,j));
%                 end
%             end
%             new_fit(i,1)=feval(fobj,new_X(i,:));
%             evaluation_count = evaluation_count + 1;
%             if new_fit(i,1)<fit(i,1)
%                 XU=[XU;new_X(i,:)];
%                 XU_fit=[XU_fit;new_fit(i)];
%                 [~,min_indexXU]=max( XU_fit);
%                 if size(XU,1)>10-round(0.1*N*(t/Max_iter)^2)
%                     XU(min_indexXU, :) = [];
%                     %                             XU(randi(size(XU, 1)), :) = [];
%                     XU_fit(min_indexXU) = [];
%                 end
%                 fit(i,1)=new_fit(i,1);
%                 X(i,:)=new_X(i,:);
%             end
        else
%             if mod(t,100)~=0
                % if rand<0.5
                r3=randperm(N,2);
                lamda=randsrc(1,1,[-1,1]);
                while i==r3(1)||i==r3(2)
                    r3=randperm(N,2);
                end
                %                                 rr=(2*rand(1,dim)-1).*0.25*(1-(t/Max_iter)^2);
                %                         rr=(0.1+0.8*(cos(pi*t/(2*Max_iter))));
                if t/Max_iter<0.5
                    rr=0.2+(1-t/Max_iter)^0.25;
                else
                    rr=0.2+0.3*(t/Max_iter)^2;%equation 26
                end
                rr1=randperm(N1,1);
                z=zeros(1,dim);
                z1=randi([1,dim],1);
                %                                                 z1=round(0.8*dim-round((t/Max_iter)*0.5*dim));
                positions = randperm(dim,z1);
                z(positions) = 1;
                xpool=X(1:N1,:);
                X3_suiji=xpool(randperm(N1,1),:);
%                 if isempty(XU)
%                      X3_suiji=X3(rr1,:);
%                 else
%                     X3_suiji=XU(randperm(size(XU,1),1),:);
%                 end
                %                 if isequal(X(i,:),X3_suiji)&&ss==1
                %                          u3=setdiff(size(XU,1),X3_rand);
                %                          X3_rand=u3(randperm(length(u3),1));
                %                          X3_suiji=XU(X3_rand,:);
                %                 elseif isequal(X(i,:),X3_suiji)&&ss==0
                %                          u2=setdiff(1:N1,rr1);
                %                          rr1=u2(randperm(length(u2),1));
                %                          X3_suiji=X3(rr1,:);
                %                 end
                new_X(i,:)=X(i,:)+z.*lamda*p(t).^2.*(X(r3(1),:)-X(r3(2),:))+z.*(rr).*(X3_suiji-X(i,:));%equation 25
           
%             for j=1:dim
%                 if new_X(i,j)>ub(:,j)
%                     new_X(i,j)=ub(:,j)-mod(new_X(i,j)-ub(:,j),ub(:,j)-lb(:,j));
%                 elseif new_X(i,j)<lb(:,j)
%                     new_X(i,j)=lb(:,j)+mod(lb(:,j)-new_X(i,j),ub(:,j)-lb(:,j));
%                 end
%             end
%             new_fit(i,1)=feval(fobj,new_X(i,:));
%             evaluation_count = evaluation_count + 1;
%             if new_fit(i,1)<fit(i,1)
%                 %                 XU=[XU;new_X(i,:)];
%                 %                 XU=[XU;new_X(i,:)];
%                 %                 XU_fit=[XU_fit;new_fit(i)];
%                 %                 [~,min_indexXU]=max( XU_fit);
%                 %                 if size(XU,1)>15-round(0.1*N*(t/Max_iter)^2)zz
%                 %                     XU(min_indexXU, :) = [];
%                 %                     %                             XU(randi(size(XU, 1)), :) = [];
%                 %                     XU_fit(min_indexXU) = [];
%                 %                 end
%                 fit(i,1)=new_fit(i,1);
%                 X(i,:)=new_X(i,:);
%             end
        end
          end
        for j=1:dim
            if new_X(i,j)>ub(:,j)
                new_X(i,j)=ub(:,j)-mod(new_X(i,j)-ub(:,j),ub(:,j)-lb(:,j));
            elseif new_X(i,j)<lb(:,j)
                new_X(i,j)=lb(:,j)+mod(lb(:,j)-new_X(i,j),ub(:,j)-lb(:,j));
            end
        end
        new_fit(i,1)=feval(fobj,new_X(i,:));
        evaluation_count = evaluation_count + 1;
        if bestfit1> new_fit(i,1)
            bestfit1= new_fit(i,1);
        end
         Bestcost1( evaluation_count)=bestfit1;
        if new_fit(i,1)<fit(i,1)
            %                 XU=[XU;new_X(i,:)];
            %                 XU=[XU;new_X(i,:)];
            %                 XU_fit=[XU_fit;new_fit(i)];
            %                 [~,min_indexXU]=max( XU_fit);
            %                 if size(XU,1)>15-round(0.1*N*(t/Max_iter)^2)zz
            %                     XU(min_indexXU, :) = [];
            %                     %                             XU(randi(size(XU, 1)), :) = [];
            %                     XU_fit(min_indexXU) = [];
            %                 end
            fit(i,1)=new_fit(i,1);
            X(i,:)=new_X(i,:);
        end
    end
      if(min(fit)<bestfit1k)
    [bestfit1k,best_index]=min(fit);
        xbest=X(best_index,:);
    end  
    %             for j=1:dim
    %                 if new_X(i,j)>ub(:,j)
    %                     new_X(i,j)=ub(:,j)-mod(new_X(i,j)-ub(:,j),ub(:,j)-lb(:,j));
    %                 elseif new_X(i,j)<lb(:,j)
    %                     new_X(i,j)=lb(:,j)+mod(lb(:,j)-new_X(i,j),ub(:,j)-lb(:,j));
    %                 end
    %             end
    %             new_fit(i,1)=feval(fobj,new_X(i,:));
    %             evaluation_count = evaluation_count + 1;
    %             if new_fit(i,1)<fita(i,1)
    %                 %                 XU=[XU;new_X(i,:)];
    % %                 XU=[XU;new_X(i,:)];
    % %                 XU_fit=[XU_fit;new_fit(i)];
    % %                 [~,min_indexXU]=max( XU_fit);
    % %                 if size(XU,1)>15-round(0.1*N*(t/Max_iter)^2)zz
    % %                     XU(min_indexXU, :) = [];
    % %                     %                             XU(randi(size(XU, 1)), :) = [];
    % %                     XU_fit(min_indexXU) = [];
    % %                 end
    %                 fita(i,1)=new_fit(i,1);
    %                 Xa(i,:)=new_X(i,:);
    %             end
    %     end
    %     xpool=Xa(1:size(X3,1),:);
    %         else
    %
    %                                                             z=zeros(1,dim);
    % %                                                             z1=randi([1,dim],1);
    %                                                             z1=round(0.8*dim-round((t/Max_iter)*0.3*dim));
    %                                                             positions = randperm(dim,z1);
    %                                                             z(positions) = 1;
    %      t
    %             for i=size(X3,1)+1:size(X,1)
    %             % %
    %             if mod(t,100)~=0
    %                 % if rand<0.5
    %                 r3=randperm(N,2);
    %                 lamda=randsrc(1,1,[-1,1]);
    %                 while i==r3(1)||i==r3(2)
    %                     r3=randperm(N,2);
    %                 end
    %                 %                                 rr=(2*rand(1,dim)-1).*0.25*(1-(t/Max_iter)^2);
    %                 %                         rr=(0.1+0.8*(cos(pi*t/(2*Max_iter))));
    %                 if t/Max_iter<0.5
    %                     rr=0.2+(1-t/Max_iter)^0.25;
    %                 else
    %                     rr=0.2+0.3*(t/Max_iter)^2;
    %                 end
    %                 rr1=randperm(N1,1);
    %                 z=zeros(1,dim);
    %                 z1=randi([1,dim],1);
    %                 %                                                 z1=round(0.8*dim-round((t/Max_iter)*0.5*dim));
    %                 positions = randperm(dim,z1);
    %                 z(positions) = 1;
    % %                 if isempty(XU)
    %                     X3_suiji=xpool(randperm(size(X3,1),1),:);
    %                     %                     ss=0;
    % %                 else
    % %                     X3_suiji=XU(randperm(size(XU,1),1),:);
    % %                     %                     ss=1;
    % %                 end
    %                 %                 if isequal(X(i,:),X3_suiji)&&ss==1
    %                 %                          u3=setdiff(size(XU,1),X3_rand);
    %                 %                          X3_rand=u3(randperm(length(u3),1));
    %                 %                          X3_suiji=XU(X3_rand,:);
    %                 %                 elseif isequal(X(i,:),X3_suiji)&&ss==0
    %                 %                          u2=setdiff(1:N1,rr1);
    %                 %                          rr1=u2(randperm(length(u2),1));
    %                 %                          X3_suiji=X3(rr1,:);
    %                 %                 end
    %                 new_X(i,:)=Xa(i,:)+z.*lamda*p(t).^2.*(X(r3(1),:)-X(r3(2),:))+z.*(rr).*(X3_suiji-Xa(i,:));
    %             else
    %                 new_X(i,:)= (min(X)+max(X)-Xa(i,:));
    %             end
    %             for j=1:dim
    %                 if new_X(i,j)>ub(:,j)
    %                     new_X(i,j)=ub(:,j)-mod(new_X(i,j)-ub(:,j),ub(:,j)-lb(:,j));
    %                 elseif new_X(i,j)<lb(:,j)
    %                     new_X(i,j)=lb(:,j)+mod(lb(:,j)-new_X(i,j),ub(:,j)-lb(:,j));
    %                 end
    %             end
    %             new_fit(i,1)=feval(fobj,new_X(i,:));
    %             evaluation_count = evaluation_count + 1;
    %             if new_fit(i,1)<fita(i,1)
    %                 fita(i,1)=new_fit(i,1);
    %                 Xa(i,:)=new_X(i,:);
    %             end
    %        end
    %     end
%     if(min(fit)<best_fit)
%         [best_fit,best_index]=min(fit);
%         xbest=X(best_index,:);
%     end
    U=0;
%                 [~,fit_indexs]=sort(fit);
%                 Xs=X(fit_indexs,:);
    for i=1:N
        
        %                                 PA=[X;A];
        %                                    r1=randi(size(PA,1));
        %         l1=randsrc(1,1,[0,1]);
        %         l2=randsrc(1,1,[0,1]);
        %         l3=randsrc(1,1,[-1,1]);
        l1 = randi([0,1]);
        l2 = randi([0,1]);
        l3 = randi([-1,1]);
        %                 a1= l1*2*(0.5+(1-cos(pi*t/(2*Max_iter)))/2)+(1-l1);
        
        %                 a3= l1*2*(0.5+(1-sin(pi*t/Max_iter))/2)+(1-l1);
        a1=l1*2*(0.5+sin(pi*t/Max_iter)/2)+(1-l1);%equation 33
        a4=l1*2*(0.5+(1-sin(pi*t/Max_iter))/2)+(1-l1);%equation 36
        a2= l1*(0.5+(1-cos(pi*t/(2*Max_iter)))/2)+(1-l1);%equation 34
        a3= l1*(0.5+cos(pi*t/Max_iter)/2)+(1-l1);%equation 35
        a5= l1*2*(0.5+(cos(pi*t/(Max_iter)))/2)+(1-l1);
        a6=l1*2*(0.5+sin(pi*t/Max_iter)/2)+(1-l1);
        %         a4=l1*2*(0.5+cos(pi*t/Max_iter)/2)+(1-l1);
        %         a4= l1*(0.5+sin(pi*t/Max_iter)/2)+(1-l1);
        %         a2= l1*(0.5+(cos(pi*t/(2*Max_iter)))/2)+(1-l1);
        %         a1= l1*2*(0.5+(1-cos(pi*t/(2*Max_iter)))/2)+(1-l1);
        %         a3= l1*2*(0.5+(cos(pi*t/(2*Max_iter)))/2)+(1-l1);
        %         a4= l1*(0.5+(1-cos(pi*t/(2*Max_iter)))/2)+(1-l1);
        rho=alpha.*(2.*rand-1);
        available = setdiff(1:N, i);
        u = available(randperm(length(available), 3));
        %         u1 = available(randperm(length(available), 1));
        %         u =randperm(N, 3);
        %         if i==u(3)
        %             candidates=setdiff(1:50,i);
        %             u(3)=randi(length(candidates));
        %         end
        k1=(-1+2*rand)*(1-(t/Max_iter)^0.25);%equation 32
        k2=normrnd(0,1-(t/Max_iter)^2);
        while abs(k2)>1
            k2=normrnd(0,1-(t/Max_iter)^2);
        end
        %         Xk=l2*(X(i,:)-X(u1,:))+X(u1,:);
         Xk=l2*(X(i,:)-X(u(3),:))+X(u(3),:);
%          Xk1=l2*(X(i,:)-Xs(aa1,:))+Xs(aa1,:);
        [~,fit_indexs]=sort(fit);
        Xs=X(fit_indexs,:);
        aa1=randperm(5,1);
        if isequal(Xk,Xs(aa1,:))
            l2=~l2;
            %             Xk=l2*(X(i,:)-X(u1,:))+X(u1,:);
            Xk=l2*(X(i,:)-X(u(3),:))+X(u(3),:);%equation 29
        end
        X5=min(X)+rand(1,dim).*(max(X)-min(X));
        X6=min(X)+rand(1,dim).*(max(X)-min(X));
        %         P_ar= calculatePopulationDiversity(X);
        nVOL = calculateHypervolumeDiversity3(X, lb, ub);%equation 37
        %                                              D = calculatePopulationDiversity(X);
        %                                             [h,h1] = calculatePopulationDiversity(X,ub,lb,xbest);% 计算多样性测度 h(k)
        %                                              Diver = calculatePopulationDiversity(X,ub,lb);
        for j=1:dim
            %             if rand<0.5
            if nVOL<0.15
                %               if P_ar<0.2*(t/Max_iter)^0.25
                new_X(i,j)=X(i,j)+(k1).*a1*(Xs(aa1,j)-Xk(1,j)) +a2*(X(u(1),j)-X(u(2),j))+k2*rho.*a3*( X5(:,j)- X6(:,j))*U;%+k2*rho.*a3*( X5(:,j)- X6(:,j))*U;%equation 28
                % %                 if abs(X(i,j) - new_X(i,j)) < 1e-10
                % %                     U=1;
                % %                 else
                % %                     U=0;
                % %                 end
            else
                % %                 %                 %                 %                if t/Max_iter<0.5
                % %                 %                 %                 f_values = fit - min(fit);
                % %                 %                 %                 F =  (1 - exp(-f_values / max(f_values + eps)));
                % %                 %                 %                 new_X(i,j)=X(i,j)+F(i).*(Xs(aa1,j)-X(i,j));
                % %                 %                 % %             else
                % %                 %                 % %                 new_X(i,j)=xbest(:,j)+exp(-t/Max_iter)*randn;
                % %                 %                 % %             end
%                 new_X(i,j)=Xk(1,j)+(k1).*a6*(Xs(aa1,j)-Xk(1,j))+(a2)*(X(u(1),j)-X(u(2),j))+k2*rho.*a3*( X5(:,j)- X6(:,j))*U;
                  new_X(i,j)=Xs(aa1,j)+(k1).*a4*(Xs(aa1,j)-Xk(1,j))+(a2)*(X(u(1),j)-X(u(2),j))/2+k2*rho.*a3*( X5(:,j)- X6(:,j))*U;%equation 28
                % %                 %                   [new_X(i,:), ~] = gradient_local_search(xbest, min(fit), grad_obj, lb, ub, 1000, 1e-6);
            end
            if abs(X(i,j) - new_X(i,j)) < 1e-10%equation 30
                U=1;
            else
                U=0;
            end
            if new_X(i,j)>ub(:,j)
                new_X(i,j)=ub(:,j)-min(new_X(i,j)-ub(:,j),ub(:,j)-lb(:,j)).*rand;
            elseif new_X(i,j)<lb(:,j)
                new_X(i,j)=lb(:,j)+min(lb(:,j)-new_X(i,j),ub(:,j)-lb(:,j)).*rand;
            end
            %                         if new_X(i,j)<lb(:,j)
            %                                         new_X(i,j)=min(X(:,j))+rand.*(max(X(:,j))-min(X(:,j)));
            %                             new_X(i,j)=lb(:,j)+(ub(:,j)-lb(:,j))*abs(new_X(i,j)-lb(:,j))/abs(new_X(i,j)-ub(:,j));
            %                         elseif new_X(i,j)>ub(:,j)
            %                             new_X(i,j)=lb(:,j)+(ub(:,j)-lb(:,j))*abs(new_X(i,j)-ub(:,j))/abs(new_X(i,j)-lb(:,j));
            %                         end
        end
        new_fit(i,1)=feval(fobj,new_X(i,:));
        evaluation_count = evaluation_count + 1;
                if bestfit1> new_fit(i,1)
            bestfit1= new_fit(i,1);
        end
         Bestcost1( evaluation_count)=bestfit1;

        if new_fit(i,1) < fit(i,1)    % 小于原有值就更新
            %                                                         A = [A; X(i,:)];
            %                                                         [~,min_indexA]=max(fit);
            %                                                         if size(A, 1) > 50
            % %                                                             A(randi(size(A, 1)), :) = [];    % 保持A的数目不超过popsize
            %                                                               A(min_indexA, :) = [];
            %                                                         end
            fit(i,1)=new_fit(i,1);
            X(i,:)=new_X(i,:);
            %             pre_X(i,:)=X(i,:);
            %                         memory_k2 = [memory_k2; k2];
            %                         if size(memory_k2, 1) > 50
            %                             memory_k2(randi(size(memory_k2, 1)), :) = [];    % 保持A的数目不超过popsize
            %                         end
        end
        
    end
     if(min(fit)<bestfit1k)
    [bestfit1k,best_index]=min(fit);
        xbest=X(best_index,:);
    end  
        pop=X;
        bestfit=Bestcost1( evaluation_count);
        Bestcost1(MaxFES+1: evaluation_count)=[];
        Bestcost2(Max_iter+1:t-1)=[];
        tk=t-1;  
end
end
% end
