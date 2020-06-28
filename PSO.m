function [Fit_and_p,FVr_bestmemit, fitMaxVector] = ...
    PSO(deParameters,caseStudyData,otherParameters,low_habitat_limit,up_habitat_limit,iRuns)
%%
%-----This is just for notational convenience and to keep the code uncluttered.--------
I_NP         = deParameters.I_NP; % 种群数 = 10
F_weight     = deParameters.F_weight; %Mutation factor 0.5
F_CR         = deParameters.F_CR; %Recombination constant 0.9
I_D          = numel(up_habitat_limit); %Number of variables or dimension 3408
deParameters.nVariables=I_D;
popmin = low_habitat_limit; %变量下限 1*3048
popmax = up_habitat_limit; %变量上限  1*3048
I_itermax    = deParameters.I_itermax; %number of max iterations/gen
fnc=  otherParameters.fnc; %选择fitness function

fitMaxVector = nan(1,I_itermax); %1*500的NaN矩阵
%%
%参数初始化
%粒子群算法中的俩个参数

c1 = 1.49445; %惯量因子，向个体的历史最佳方向移动
c2 = 1.49445; %惯量因子，向所有10个个体的历史最佳方向移动
w_max=1.4;
w_min=0.4;

maxg=deParameters.I_itermax; %进化次数 移动500回
sizepop=I_NP; %总群规模 10只鸟
%初始速度和总群上下边界值
popminMatrix=repmat(popmin,I_NP,1); %产生10*3408的矩阵，每列相同，每行为变量下限
popmaxMatrix=repmat(popmax,I_NP,1); %产生10*3408的矩阵，每列相同，每行为变量上限
deParameters.minPositionsMatrix=popminMatrix; %存入deParameters中
deParameters.maxPositionsMatrix=popmaxMatrix; %存入deParameters中
Vgap=(popmax-popmin)/3; %速度范围,通过变量上下限设定，上帝让我飞这么快啊【速度上下限需要优化】
VminMatrix=repmat(-Vgap,I_NP,1); %产生10*3408的矩阵，每列相同，每行为变量下限
VmaxMatrix=repmat(Vgap,I_NP,1); %产生10*3408的矩阵，每列相同，每行为变量上限

%%产生初始粒子和速度，上帝说要有鸟，就有了鸟
% generate initial population.随机生成的第一代10*3408的子代（10为子代数量，3408为变量数）
rand('state',otherParameters.iRuns) %Guarantee same initial population
pop=genpop(I_NP,I_D,popminMatrix,popmaxMatrix); %随机生成的第一代10*3408的子代（10为子代数量，3408为变量数）
V=genpop(I_NP,I_D,VminMatrix,VmaxMatrix); %随机生成的第一代10*3408子代的速度（10为子代数量，3408为变量数）
[fitness, ~]=feval(fnc,pop,caseStudyData,otherParameters);%10个子代的fitness function值，10*1

%%找最好的染色体，最大的小肉块啊
[bestfitness,bestindex]=min(fitness);
zbest=pop(bestindex,:); %全局最佳
gbest=pop; %个体最佳
fitnessgbest=fitness; %个体最佳适应度
fitnesszbest=bestfitness; %全局最佳适应度

fitMaxVector(1,1) = fitnesszbest;%将初始最佳子代的fitness值存入fitMaxVector中，1*500
FVr_bestmemit = pop(bestindex,:);
%%迭代寻优，我要找499啊，老婆在家孵蛋啊
for i=1:maxg-1 %移动499次
    w=w_max-i*(w_max-w_min)/maxg;%权重线性递减的PSO算法
    for j = 1:sizepop %10个个体
        %速度更新，要找最好吃的肉啊
        V(j,:) = w*V(j,:) + c1*rand*(gbest(j,:)-pop(j,:))+c2*rand*(zbest-pop(j,:));
        V(j,: )= update(V(j,: ),VminMatrix(j,:),VmaxMatrix(j,:));
        %移动
        pop(j,:)=pop(j,:)+V(j,:);
        pop(j,: )= update(pop(j,: ), popminMatrix(j,:), popmaxMatrix(j,:));
        %自适应变异，洒鸟分了心啊
        if rand>0.9
            k=ceil(I_D*rand);
            pop(j,k)=popmin(k)+(popmax(k)-popmin(k))*rand;
        end
    end   
    %适应度值，函数最小的点（鸟）
    [fitness, ~]=feval(fnc,pop,caseStudyData,otherParameters);
    %个体最优更新，我找到最好吃的肉要实时追踪啊
    for j=1:sizepop %10个个体
        if fitness(j) < fitnessgbest(j)
            gbest(j,:) = pop(j,:);
            fitnessgbest(j) = fitness(j);
        end
    end
    %群体最优更新，我是一只善于观察和沟通并成功让酋长漂亮女儿当我老婆的鸟
    [smallest_fit,smallest_index]=min(fitness);
    if smallest_fit< fitnesszbest
        zbest = pop(smallest_index,:);
        fitnesszbest = fitness(smallest_index);
    end
    fitMaxVector(1,i+1) = fitnesszbest;%将初始最佳子代的fitness值存入fitMaxVector中，1*500 
    FVr_bestmemit = zbest; % best member of current iteration
end
Fit_and_p=[fitMaxVector(1,i+1) 0];
%结果分析，20个鸟到了同一个地方，来啊互相伤害啊
figure(1)
plot(fitMaxVector,'Linewidth',2)%
title(['适应度曲线 ' '终止代数=' num2str(maxg)]);
grid on
xlabel('进化代数' ); ylabel('适应度');
%结果输出
%zbest; %最佳个体
%fitnesszbest; %最优值
end

% VECTORIZED THE CODE INSTEAD OF USING FOR 创建初始族
function pop=genpop(a,b,lowMatrix,upMatrix) %用于产生初始族，a=10,b=3408
    pop=unifrnd(lowMatrix,upMatrix,a,b);
end

%判断新产的个体或速度是否各变量满足上下限要求，将不满足的替换成上下限的值
function p=update(p,lowMatrix,upMatrix)
        %[popsize,dim]=size(p);
        [idx] = find(p<lowMatrix);
        p(idx)=lowMatrix(idx);
        [idx] = find(p>upMatrix);
        p(idx)=upMatrix(idx);
end