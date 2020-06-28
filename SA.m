function [Fit_and_p,FVr_bestmemit, fitMaxVector] = ...
    SA(deParameters,caseStudyData,otherParameters,low_habitat_limit,up_habitat_limit,iRuns)
%%
%-----This is just for notational convenience and to keep the code uncluttered.--------
I_NP         = deParameters.I_NP; % 种群数 = 10
I_D          = numel(up_habitat_limit); %Number of variables or dimension 3408
deParameters.nVariables=I_D;
popmin = low_habitat_limit; %变量下限 1*3048
popmax = up_habitat_limit; %变量上限  1*3048
I_itermax    = deParameters.I_itermax; %number of max iterations/gen
fnc=  otherParameters.fnc; %选择fitness function
fitMaxVector = nan(1,I_itermax); %1*500的NaN矩阵
T=2000;%初始温度
alpha=0.99;%退火系数
de_save(1,1)=1;
%%
%参数初始化
maxg=deParameters.I_itermax; %迭代次数500次
sizepop=I_NP; %总群规模 10个小球
%初始速度和总群上下边界值
popminMatrix=repmat(popmin,I_NP,1); %产生10*3408的矩阵，每列相同，每行为变量下限
popmaxMatrix=repmat(popmax,I_NP,1); %产生10*3408的矩阵，每列相同，每行为变量上限
deParameters.minPositionsMatrix=popminMatrix; %存入deParameters中
deParameters.maxPositionsMatrix=popmaxMatrix; %存入deParameters中

%%产生初始代
% generate initial population.随机生成的第一代10*3408的子代（10为子代数量，3408为变量数）
rand('state',otherParameters.iRuns) %Guarantee same initial population
pop=genpop(I_NP,I_D,popminMatrix,popmaxMatrix); %随机生成的第一代10*3408的子代（10为子代数量，3408为变量数）
[fitness, ~]=feval(fnc,pop,caseStudyData,otherParameters);%10个子代的fitness function值，10*1
%找最好的小球
[bestfitness,bestindex]=min(fitness);
fitnesszbest=bestfitness; %全局最佳适应度1*1
fitMaxVector(1,1) = fitnesszbest;%将初始最佳子代的fitness值存入fitMaxVector中，1*500
fitMaxVector_alt(1,1) = fitnesszbest;%将初始最佳子代的fitness值存入fitMaxVector中，1*500
FVr_bestmemit = pop(bestindex,:);%将初始最佳子代存入fitMaxVector中，1*3408
zbest=pop(bestindex,:); %全局最佳个体1*3048
%%迭代寻优，跳499次
for i=1:maxg-1 %跳499次
    %产生随机扰动
    pop_new=disturb(pop,popmin,popmax);
    fitness_old=fitness;%更新fitness_old
    %计算
    [fitness, ~]=feval(fnc,pop_new,caseStudyData,otherParameters);%10个子代的fitness function值，10*1
    %选择新一代种族
    for j=1:sizepop
        %退火概率
        prob=exp(-abs((fitness(j)-fitness_old(j))/fitness_old(j))*4*i/(maxg+1));%2,退火概率
        %prob=exp(-(fitness(j)-fitness_old(j))/T);%4,退火概率
        if fitness(j)<=fitness_old(j)
            pop(j,:)=pop_new(j,:);
        %若新的fitness大于旧的fitness，则以rand为条件视情况是否接受新的    
        else 
            if rand<=prob  
                pop(j,:)=pop_new(j,:);
            else %不更新pop，选取旧的pop及其对应的fitness
                fitness(j,:)=fitness_old(j,:);
            end
        end
        de_save(j,i)=prob;
    end
    T=T*alpha;
    %群体最优更新
    [bestfitness,bestindex]=min(fitness);
    fitMaxVector_alt(1,i+1) = bestfitness;
    if bestfitness<=fitnesszbest
        zbest=pop(bestindex,:); %全局最佳个体1*3048
        fitnesszbest=bestfitness; %全局最佳适应度1*1
    end
    fitMaxVector(1,i+1) = fitnesszbest;%将初始最佳子代的fitness值存入fitMaxVector中，1*500 
    FVr_bestmemit = zbest; % best member of current iteration
end
Fit_and_p=[fitMaxVector(1,i+1) 0];
%结果分析，10个球到了同一个地方
figure(1)
plot(fitMaxVector,'Linewidth',2)%
hold on
plot(fitMaxVector_alt,'Linewidth',2)%
title(['适应度曲线 ' '终止代数=' num2str(maxg)]);
grid on
xlabel('进化代数' ); ylabel('适应度');
end

% VECTORIZED THE CODE INSTEAD OF USING FOR 创建初始族
function pop=genpop(a,b,lowMatrix,upMatrix) %用于产生初始族，a=10,b=3408
    pop=unifrnd(lowMatrix,upMatrix,a,b);
end

%产生随机扰动
function pop_new=disturb(pop,lowMatrix,upMatrix)
    for i=1:size(pop,1)
        dist_step=(rand-0.5).*(upMatrix-lowMatrix)/10;
        pop_new(i,:)=pop(i,:)+dist_step;
        [idx_low] = find(pop_new(i,:)<lowMatrix);
        [idx_up] = find(pop_new(i,:)>upMatrix);
        while ~isempty(idx_low)
            step_up=abs(rand*dist_step(1,idx_low)/10);
            pop_new(i,idx_low)=pop_new(i,idx_low)+step_up;
            [idx_low] = find(pop_new(i,:)<lowMatrix);
        end
        while ~isempty(idx_up)
            step_back=abs(rand*dist_step(1,idx_up)/10);
            pop_new(i,idx_up)=pop_new(i,idx_up)-step_back;
            [idx_up] = find(pop_new(i,:)>upMatrix);
        end
    end
end
