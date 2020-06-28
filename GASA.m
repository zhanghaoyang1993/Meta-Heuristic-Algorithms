function [Fit_and_p,FVr_bestmemit, fitMaxVector] = ...
    GASA(deParameters,caseStudyData,otherParameters,low_habitat_limit,up_habitat_limit,iRuns)
%%
%-----This is just for notational convenience and to keep the code uncluttered.--------
I_NP         = deParameters.I_NP; % ��Ⱥ�� = 10
sizepop1 = floor(I_NP/2); %GA ��Ⱥ��
sizepop2 = I_NP-sizepop1; %SA��Ⱥ��
elite_num = 1; %��Ӣ����
mutation_rate_norm     = 0.008; %������� 0.008
mutation_rate_eli     = 0; %��Ӣ������� 0
w=2; %�����Ա任ϵ��
elitism = 1; % elitism: �����Ƿ�Ӣѡ��
T=2000;%��ʼ�¶�
alpha=0.99;%�˻�ϵ��
I_D = numel(up_habitat_limit); %Number of variables or dimension 3408
deParameters.nVariables=I_D;
popmin = low_habitat_limit; %�������� 1*3048
popmax = up_habitat_limit; %��������  1*3048
I_itermax    = deParameters.I_itermax; %number of max iterations/gen
fnc=  otherParameters.fnc; %ѡ��fitness function
fitMaxVector = nan(1,I_itermax); %1*500��NaN����
fitMaxVector_draw = nan(1,I_itermax); %1*500��NaN����
%%
%������ʼ��
maxg=deParameters.I_itermax; %�������� ����500��
sizepop=I_NP; %��Ⱥ��ģ 10��Ⱦɫ��
%��ʼ�ٶȺ���Ⱥ���±߽�ֵ
popminMatrix=repmat(popmin,I_NP,1); %����10*3408�ľ���ÿ����ͬ��ÿ��Ϊ��������
popmaxMatrix=repmat(popmax,I_NP,1); %����10*3408�ľ���ÿ����ͬ��ÿ��Ϊ��������
deParameters.minPositionsMatrix=popminMatrix; %����deParameters��
deParameters.maxPositionsMatrix=popmaxMatrix; %����deParameters��

%%������ʼ��
% generate initial population.������ɵĵ�һ��10*3408���Ӵ���10Ϊ�Ӵ�������3408Ϊ��������
rand('state',otherParameters.iRuns) %Guarantee same initial population
pop=genpop(I_NP,I_D,popminMatrix,popmaxMatrix); %������ɵĵ�һ��10*3408���Ӵ���10Ϊ�Ӵ�������3408Ϊ��������
pop1 = pop(1:sizepop1,:);
pop2 = pop(sizepop1+1:I_NP,:);
[fitness, ~]=feval(fnc,pop,caseStudyData,otherParameters);%10���Ӵ���fitness functionֵ��10*1
fitness1=fitness(1:sizepop1,:);
fitness2=fitness(sizepop1+1:I_NP,:);
%����õ�Ⱦɫ��
[bestfitness,bestindex]=min(fitness);
fitnesszbest=bestfitness; %ȫ�������Ӧ��1*1
fitMaxVector(1,1) = fitnesszbest;%����ʼ����Ӵ���fitnessֵ����fitMaxVector�У�1*500
fitMaxVector_draw(1,1) = fitnesszbest;%����ʼ����Ӵ���fitnessֵ����fitMaxVector_draw�У����ڻ�ͼ 1*500
FVr_bestmemit = pop(bestindex,:);%����ʼ����Ӵ�����fitMaxVector�У�1*3408
zbest=pop(bestindex,:); %ȫ����Ѹ���1*3048
%%����Ѱ�ţ��Ŵ�499��
for i=1:maxg-1 %�Ŵ�499��
        %����
        [fitness_top,pop_top]=rank(fitness,sizepop,pop);%10���Ӵ�һ����
        fitness1_temp=[fitness_top(1:4,:);fitness_top(sizepop-sizepop1+5:sizepop,:)];%ȡǰ4���Ӵ��͵���sizepop1-4���Ӵ�
        pop1_temp=[pop_top(1:4,:);pop_top(sizepop-sizepop1+5:sizepop,:)];%%ȡǰ2���Ӵ��͵���sizepop1-2���Ӵ�
        fitness1=fitness1_temp;
        pop1=pop1_temp;
        fitness_norm=(fitness1-fitness1(1))./(fitness1(sizepop1)-fitness1(1)); %fitness��һ��
        fitness_norm_rev=1-fitness_norm;
        fitness_adj=(exp(w.*fitness_norm_rev)-1)./(exp(w)-1);%��fitness_norm�����Ա任����Ϊ�˴�fitnessԽС����ѡ�еĸ���Խ��
        for j=1:sizepop1
            if j==1
                fitness_adj_sum(1)=fitness_adj(1);
            else
                fitness_adj_sum(j)=fitness_adj_sum(j-1)+fitness_adj(j);
            end
        end
        %selection��crossover��ѡ��ĸ���������Ӵ�
        pop1 = sele_and_cross(fitness_adj,fitness_adj_sum,pop1,sizepop1,elitism,elite_num);  
        %mutation ����
        pop1=mutation(pop1,mutation_rate_norm,mutation_rate_eli,sizepop1,popmin,popmax,elite_num);
        %����
        %[fitness1, ~]=feval(fnc,pop1,caseStudyData,otherParameters);%10���Ӵ���fitness functionֵ��10*1   
%%
%SA
        %��������Ŷ�
        %pop_new_2=disturb(pop_top5,popmin,popmax);
        pop_new_2=disturb(pop1_temp,popmin,popmax);
        pop_new=[pop1;pop_new_2];
        fitness_old=fitness2;%����fitness_old
        %����
        [fitness, ~]=feval(fnc,pop_new,caseStudyData,otherParameters);%10���Ӵ���fitness functionֵ��10*1
        fitness1=fitness(1:sizepop1,:);
        fitness2=fitness(sizepop1+1:I_NP,:);
        %ѡ����һ������
        for j=1:sizepop2
            %�˻����
            prob=exp(-abs((fitness2(j)-fitness_old(j))/fitness_old(j))*4*i/(maxg+1));%2,�˻����:https://www.sciencedirect.com/science/article/pii/S0959652619330707
            %prob=exp(-(fitness2(j)-fitness_old(j))/T);%4,�˻����
            if fitness2(j)<=fitness_old(j)
                pop2(j,:)=pop_new_2(j,:);
            %���µ�fitness���ھɵ�fitness������randΪ����������Ƿ�����µ�    
            else 
                if rand<=prob  
                    pop2(j,:)=pop_new_2(j,:);
                else %������pop��ѡȡ�ɵ�pop�����Ӧ��fitness
                    fitness2(j,:)=fitness_old(j,:);
                end
            end
        end
        T=T*alpha;
    fitness=[fitness1;fitness2];
    pop=[pop1;pop2];
    %Ⱥ�����Ÿ���
    [bestfitness,bestindex]=min(fitness);
    fitMaxVector_draw(1,i+1) = bestfitness;%����ʼ����Ӵ���fitnessֵ����fitMaxVector�У����ڻ�ͼ��1*500
    if bestfitness<=fitnesszbest
        zbest=pop(bestindex,:); %ȫ����Ѹ���1*3048
        fitnesszbest=bestfitness; %ȫ�������Ӧ��1*1
    end
    fitMaxVector(1,i+1) = fitnesszbest;%����ʼ����Ӵ���fitnessֵ����fitMaxVector�У�1*500 
    FVr_bestmemit = zbest; % best member of current iteration
end
Fit_and_p=[fitMaxVector(1,i+1) 0];
%�������
figure(1)
plot(fitMaxVector_draw,'Linewidth',2)%
title(['��Ӧ������ ' '��ֹ����=' num2str(maxg)]);
grid on
xlabel('��������' ); ylabel('��Ӧ��');
end

% VECTORIZED THE CODE INSTEAD OF USING FOR ������ʼ��
function pop=genpop(a,b,lowMatrix,upMatrix) %���ڲ�����ʼ�壬a=10,b=3408
    pop=unifrnd(lowMatrix,upMatrix,a,b);
end

%�ж��²��ĸ����Ƿ����������������Ҫ�󣬽���������滻�������޵�ֵ
function p=update(p,lowMatrix,upMatrix)
    %[popsize,dim]=size(p);
    [idx] = find(p<lowMatrix);
    p(idx)=lowMatrix(idx);
    [idx] = find(p>upMatrix);
    p(idx)=upMatrix(idx);
end


%�Ը��尴��Ӧ�ȴ�С�������� fitness��10*1����С��������
function [fitness,pop]=rank(fitness,sizepop,pop)
    % ð������
    for i=1:sizepop
        min_index = i;
        for j = i+1:sizepop
            if fitness(j) < fitness(min_index)
                min_index = j;
            end
        end
        if min_index ~= i
            % ���� fitness(i) �� fitness(min_index) ��ֵ
            temp = fitness(i);
            fitness(i) = fitness(min_index);
            fitness(min_index) = temp;
            % ��ʱ fitness_value(i) ����Ӧ����[i,population_size]����С

            % ���� population(i) �� population(min_index) ��Ⱦɫ�崮
            temp_chromosome = pop(i,:);
            pop(i,:) = pop(min_index,:);
            pop(min_index,:) = temp_chromosome;
        end    
    end
end

%ѡ���������pop��ѡ�񣬵õ��µ�population��10*3408�������Ϊ��Ӣ������������һ��
function pop = sele_and_cross(fitness_adj,fitness_adj_sum,pop,sizepop,elitism,elite_num)

% �Ƿ�Ӣѡ��
if elitism==1
    p = sizepop-elite_num;
else
    p = sizepop;
end

for i=1:p
    r1 = rand * fitness_adj_sum(sizepop);  % ����һ�����������[0,����Ӧ��]֮��
    r2 = rand * fitness_adj_sum(sizepop);  % ����һ�����������[0,����Ӧ��]֮��
    while r1==r2
        r2 = rand * fitness_adj_sum(sizepop);  % ȷ��ĸ��������ͬ
    end    
    p1 = min(find(fitness_adj_sum > r1));  % ĸ�����
    p2 = min(find(fitness_adj_sum > r2));  % �������
    pop1=pop(p1,:); %ĸ��
    pop2=pop(p2,:); %����
    fitness_cross=fitness_adj(p1)+fitness_adj(p2);
    for j=1:size(pop,2)
        if rand * fitness_cross<=fitness_adj(p1)
            pop_new(i,j)=pop1(1,j);
        else
            pop_new(i,j)=pop2(1,j);
        end
    end 
end
% �Ƿ�Ӣѡ��
if elitism==1
    pop(elite_num+1:sizepop,:) = pop_new;
else
    pop = pop_new;
end

end

%mutation������
function pop=mutation(pop,mutation_rate_norm,mutation_rate_eli,sizepop,popmin,popmax,elite_num)
    for i=1:elite_num % ��Ӣ����
        for j=1:size(pop,2)
            if rand < mutation_rate_eli
                pop(i,j) = popmin(1,j)+rand*(popmax(1,j)-popmin(1,j));
            end
        end
    end
    for i=elite_num+1:sizepop % ���ڱ���
        for j=1:size(pop,2)
            if rand < mutation_rate_norm
                pop(i,j) = popmin(1,j)+rand*(popmax(1,j)-popmin(1,j));
            end
        end
    end
end

%��������Ŷ�
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