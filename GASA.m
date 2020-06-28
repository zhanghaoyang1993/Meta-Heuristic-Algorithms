function [Fit_and_p,FVr_bestmemit, fitMaxVector] = ...
    GASA(deParameters,caseStudyData,otherParameters,low_habitat_limit,up_habitat_limit,iRuns)
%%
%-----This is just for notational convenience and to keep the code uncluttered.--------
I_NP         = deParameters.I_NP; % 种群数 = 10
sizepop1 = floor(I_NP/2); %GA 种群数
sizepop2 = I_NP-sizepop1; %SA种群数
elite_num = 1; %精英数量
mutation_rate_norm     = 0.008; %变异概率 0.008
mutation_rate_eli     = 0; %精英变异概率 0
w=2; %非线性变换系数
elitism = 1; % elitism: 输入是否精英选择
T=2000;%初始温度
alpha=0.99;%退火系数
I_D = numel(up_habitat_limit); %Number of variables or dimension 3408
deParameters.nVariables=I_D;
popmin = low_habitat_limit; %变量下限 1*3048
popmax = up_habitat_limit; %变量上限  1*3048
I_itermax    = deParameters.I_itermax; %number of max iterations/gen
fnc=  otherParameters.fnc; %选择fitness function
fitMaxVector = nan(1,I_itermax); %1*500的NaN矩阵
fitMaxVector_draw = nan(1,I_itermax); %1*500的NaN矩阵
%%
%参数初始化
maxg=deParameters.I_itermax; %进化次数 进化500回
sizepop=I_NP; %总群规模 10个染色体
%初始速度和总群上下边界值
popminMatrix=repmat(popmin,I_NP,1); %产生10*3408的矩阵，每列相同，每行为变量下限
popmaxMatrix=repmat(popmax,I_NP,1); %产生10*3408的矩阵，每列相同，每行为变量上限
deParameters.minPositionsMatrix=popminMatrix; %存入deParameters中
deParameters.maxPositionsMatrix=popmaxMatrix; %存入deParameters中

%%产生初始代
% generate initial population.随机生成的第一代10*3408的子代（10为子代数量，3408为变量数）
rand('state',otherParameters.iRuns) %Guarantee same initial population
pop=genpop(I_NP,I_D,popminMatrix,popmaxMatrix); %随机生成的第一代10*3408的子代（10为子代数量，3408为变量数）
pop1 = pop(1:sizepop1,:);
pop2 = pop(sizepop1+1:I_NP,:);
[fitness, ~]=feval(fnc,pop,caseStudyData,otherParameters);%10个子代的fitness function值，10*1
fitness1=fitness(1:sizepop1,:);
fitness2=fitness(sizepop1+1:I_NP,:);
%找最好的染色体
[bestfitness,bestindex]=min(fitness);
fitnesszbest=bestfitness; %全局最佳适应度1*1
fitMaxVector(1,1) = fitnesszbest;%将初始最佳子代的fitness值存入fitMaxVector中，1*500
fitMaxVector_draw(1,1) = fitnesszbest;%将初始最佳子代的fitness值存入fitMaxVector_draw中，用于画图 1*500
FVr_bestmemit = pop(bestindex,:);%将初始最佳子代存入fitMaxVector中，1*3408
zbest=pop(bestindex,:); %全局最佳个体1*3048
%%迭代寻优，遗传499代
for i=1:maxg-1 %遗传499代
        %排序
        [fitness_top,pop_top]=rank(fitness,sizepop,pop);%10个子代一起排
        fitness1_temp=[fitness_top(1:4,:);fitness_top(sizepop-sizepop1+5:sizepop,:)];%取前4个子代和倒数sizepop1-4个子代
        pop1_temp=[pop_top(1:4,:);pop_top(sizepop-sizepop1+5:sizepop,:)];%%取前2个子代和倒数sizepop1-2个子代
        fitness1=fitness1_temp;
        pop1=pop1_temp;
        fitness_norm=(fitness1-fitness1(1))./(fitness1(sizepop1)-fitness1(1)); %fitness归一化
        fitness_norm_rev=1-fitness_norm;
        fitness_adj=(exp(w.*fitness_norm_rev)-1)./(exp(w)-1);%将fitness_norm非线性变换，因为此处fitness越小，被选中的概率越大
        for j=1:sizepop1
            if j==1
                fitness_adj_sum(1)=fitness_adj(1);
            else
                fitness_adj_sum(j)=fitness_adj_sum(j-1)+fitness_adj(j);
            end
        end
        %selection和crossover，选则父母并产生新子代
        pop1 = sele_and_cross(fitness_adj,fitness_adj_sum,pop1,sizepop1,elitism,elite_num);  
        %mutation 变异
        pop1=mutation(pop1,mutation_rate_norm,mutation_rate_eli,sizepop1,popmin,popmax,elite_num);
        %计算
        %[fitness1, ~]=feval(fnc,pop1,caseStudyData,otherParameters);%10个子代的fitness function值，10*1   
%%
%SA
        %产生随机扰动
        %pop_new_2=disturb(pop_top5,popmin,popmax);
        pop_new_2=disturb(pop1_temp,popmin,popmax);
        pop_new=[pop1;pop_new_2];
        fitness_old=fitness2;%更新fitness_old
        %计算
        [fitness, ~]=feval(fnc,pop_new,caseStudyData,otherParameters);%10个子代的fitness function值，10*1
        fitness1=fitness(1:sizepop1,:);
        fitness2=fitness(sizepop1+1:I_NP,:);
        %选择新一代种族
        for j=1:sizepop2
            %退火概率
            prob=exp(-abs((fitness2(j)-fitness_old(j))/fitness_old(j))*4*i/(maxg+1));%2,退火概率:https://www.sciencedirect.com/science/article/pii/S0959652619330707
            %prob=exp(-(fitness2(j)-fitness_old(j))/T);%4,退火概率
            if fitness2(j)<=fitness_old(j)
                pop2(j,:)=pop_new_2(j,:);
            %若新的fitness大于旧的fitness，则以rand为条件视情况是否接受新的    
            else 
                if rand<=prob  
                    pop2(j,:)=pop_new_2(j,:);
                else %不更新pop，选取旧的pop及其对应的fitness
                    fitness2(j,:)=fitness_old(j,:);
                end
            end
        end
        T=T*alpha;
    fitness=[fitness1;fitness2];
    pop=[pop1;pop2];
    %群体最优更新
    [bestfitness,bestindex]=min(fitness);
    fitMaxVector_draw(1,i+1) = bestfitness;%将初始最佳子代的fitness值存入fitMaxVector中，用于画图，1*500
    if bestfitness<=fitnesszbest
        zbest=pop(bestindex,:); %全局最佳个体1*3048
        fitnesszbest=bestfitness; %全局最佳适应度1*1
    end
    fitMaxVector(1,i+1) = fitnesszbest;%将初始最佳子代的fitness值存入fitMaxVector中，1*500 
    FVr_bestmemit = zbest; % best member of current iteration
end
Fit_and_p=[fitMaxVector(1,i+1) 0];
%结果分析
figure(1)
plot(fitMaxVector_draw,'Linewidth',2)%
title(['适应度曲线 ' '终止代数=' num2str(maxg)]);
grid on
xlabel('进化代数' ); ylabel('适应度');
end

% VECTORIZED THE CODE INSTEAD OF USING FOR 创建初始族
function pop=genpop(a,b,lowMatrix,upMatrix) %用于产生初始族，a=10,b=3408
    pop=unifrnd(lowMatrix,upMatrix,a,b);
end

%判断新产的个体是否各变量满足上下限要求，将不满足的替换成上下限的值
function p=update(p,lowMatrix,upMatrix)
    %[popsize,dim]=size(p);
    [idx] = find(p<lowMatrix);
    p(idx)=lowMatrix(idx);
    [idx] = find(p>upMatrix);
    p(idx)=upMatrix(idx);
end


%对个体按适应度大小进行排序 fitness（10*1）从小到大排序；
function [fitness,pop]=rank(fitness,sizepop,pop)
    % 冒泡排序
    for i=1:sizepop
        min_index = i;
        for j = i+1:sizepop
            if fitness(j) < fitness(min_index)
                min_index = j;
            end
        end
        if min_index ~= i
            % 交换 fitness(i) 和 fitness(min_index) 的值
            temp = fitness(i);
            fitness(i) = fitness(min_index);
            fitness(min_index) = temp;
            % 此时 fitness_value(i) 的适应度在[i,population_size]上最小

            % 交换 population(i) 和 population(min_index) 的染色体串
            temp_chromosome = pop(i,:);
            pop(i,:) = pop(min_index,:);
            pop(min_index,:) = temp_chromosome;
        end    
    end
end

%选择操作，从pop中选择，得到新的population（10*3408），如果为精英操作，则保留第一项
function pop = sele_and_cross(fitness_adj,fitness_adj_sum,pop,sizepop,elitism,elite_num)

% 是否精英选择
if elitism==1
    p = sizepop-elite_num;
else
    p = sizepop;
end

for i=1:p
    r1 = rand * fitness_adj_sum(sizepop);  % 生成一个随机数，在[0,总适应度]之间
    r2 = rand * fitness_adj_sum(sizepop);  % 生成一个随机数，在[0,总适应度]之间
    while r1==r2
        r2 = rand * fitness_adj_sum(sizepop);  % 确保母代父代不同
    end    
    p1 = min(find(fitness_adj_sum > r1));  % 母代序号
    p2 = min(find(fitness_adj_sum > r2));  % 父代序号
    pop1=pop(p1,:); %母代
    pop2=pop(p2,:); %父代
    fitness_cross=fitness_adj(p1)+fitness_adj(p2);
    for j=1:size(pop,2)
        if rand * fitness_cross<=fitness_adj(p1)
            pop_new(i,j)=pop1(1,j);
        else
            pop_new(i,j)=pop2(1,j);
        end
    end 
end
% 是否精英选择
if elitism==1
    pop(elite_num+1:sizepop,:) = pop_new;
else
    pop = pop_new;
end

end

%mutation，变异
function pop=mutation(pop,mutation_rate_norm,mutation_rate_eli,sizepop,popmin,popmax,elite_num)
    for i=1:elite_num % 精英变异
        for j=1:size(pop,2)
            if rand < mutation_rate_eli
                pop(i,j) = popmin(1,j)+rand*(popmax(1,j)-popmin(1,j));
            end
        end
    end
    for i=elite_num+1:sizepop % 民众变异
        for j=1:size(pop,2)
            if rand < mutation_rate_norm
                pop(i,j) = popmin(1,j)+rand*(popmax(1,j)-popmin(1,j));
            end
        end
    end
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