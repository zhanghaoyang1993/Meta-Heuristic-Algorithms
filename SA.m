function [Fit_and_p,FVr_bestmemit, fitMaxVector] = ...
    SA(deParameters,caseStudyData,otherParameters,low_habitat_limit,up_habitat_limit,iRuns)
%%
%-----This is just for notational convenience and to keep the code uncluttered.--------
I_NP         = deParameters.I_NP; % ��Ⱥ�� = 10
I_D          = numel(up_habitat_limit); %Number of variables or dimension 3408
deParameters.nVariables=I_D;
popmin = low_habitat_limit; %�������� 1*3048
popmax = up_habitat_limit; %��������  1*3048
I_itermax    = deParameters.I_itermax; %number of max iterations/gen
fnc=  otherParameters.fnc; %ѡ��fitness function
fitMaxVector = nan(1,I_itermax); %1*500��NaN����
T=2000;%��ʼ�¶�
alpha=0.99;%�˻�ϵ��
de_save(1,1)=1;
%%
%������ʼ��
maxg=deParameters.I_itermax; %��������500��
sizepop=I_NP; %��Ⱥ��ģ 10��С��
%��ʼ�ٶȺ���Ⱥ���±߽�ֵ
popminMatrix=repmat(popmin,I_NP,1); %����10*3408�ľ���ÿ����ͬ��ÿ��Ϊ��������
popmaxMatrix=repmat(popmax,I_NP,1); %����10*3408�ľ���ÿ����ͬ��ÿ��Ϊ��������
deParameters.minPositionsMatrix=popminMatrix; %����deParameters��
deParameters.maxPositionsMatrix=popmaxMatrix; %����deParameters��

%%������ʼ��
% generate initial population.������ɵĵ�һ��10*3408���Ӵ���10Ϊ�Ӵ�������3408Ϊ��������
rand('state',otherParameters.iRuns) %Guarantee same initial population
pop=genpop(I_NP,I_D,popminMatrix,popmaxMatrix); %������ɵĵ�һ��10*3408���Ӵ���10Ϊ�Ӵ�������3408Ϊ��������
[fitness, ~]=feval(fnc,pop,caseStudyData,otherParameters);%10���Ӵ���fitness functionֵ��10*1
%����õ�С��
[bestfitness,bestindex]=min(fitness);
fitnesszbest=bestfitness; %ȫ�������Ӧ��1*1
fitMaxVector(1,1) = fitnesszbest;%����ʼ����Ӵ���fitnessֵ����fitMaxVector�У�1*500
fitMaxVector_alt(1,1) = fitnesszbest;%����ʼ����Ӵ���fitnessֵ����fitMaxVector�У�1*500
FVr_bestmemit = pop(bestindex,:);%����ʼ����Ӵ�����fitMaxVector�У�1*3408
zbest=pop(bestindex,:); %ȫ����Ѹ���1*3048
%%����Ѱ�ţ���499��
for i=1:maxg-1 %��499��
    %��������Ŷ�
    pop_new=disturb(pop,popmin,popmax);
    fitness_old=fitness;%����fitness_old
    %����
    [fitness, ~]=feval(fnc,pop_new,caseStudyData,otherParameters);%10���Ӵ���fitness functionֵ��10*1
    %ѡ����һ������
    for j=1:sizepop
        %�˻����
        prob=exp(-abs((fitness(j)-fitness_old(j))/fitness_old(j))*4*i/(maxg+1));%2,�˻����
        %prob=exp(-(fitness(j)-fitness_old(j))/T);%4,�˻����
        if fitness(j)<=fitness_old(j)
            pop(j,:)=pop_new(j,:);
        %���µ�fitness���ھɵ�fitness������randΪ����������Ƿ�����µ�    
        else 
            if rand<=prob  
                pop(j,:)=pop_new(j,:);
            else %������pop��ѡȡ�ɵ�pop�����Ӧ��fitness
                fitness(j,:)=fitness_old(j,:);
            end
        end
        de_save(j,i)=prob;
    end
    T=T*alpha;
    %Ⱥ�����Ÿ���
    [bestfitness,bestindex]=min(fitness);
    fitMaxVector_alt(1,i+1) = bestfitness;
    if bestfitness<=fitnesszbest
        zbest=pop(bestindex,:); %ȫ����Ѹ���1*3048
        fitnesszbest=bestfitness; %ȫ�������Ӧ��1*1
    end
    fitMaxVector(1,i+1) = fitnesszbest;%����ʼ����Ӵ���fitnessֵ����fitMaxVector�У�1*500 
    FVr_bestmemit = zbest; % best member of current iteration
end
Fit_and_p=[fitMaxVector(1,i+1) 0];
%���������10������ͬһ���ط�
figure(1)
plot(fitMaxVector,'Linewidth',2)%
hold on
plot(fitMaxVector_alt,'Linewidth',2)%
title(['��Ӧ������ ' '��ֹ����=' num2str(maxg)]);
grid on
xlabel('��������' ); ylabel('��Ӧ��');
end

% VECTORIZED THE CODE INSTEAD OF USING FOR ������ʼ��
function pop=genpop(a,b,lowMatrix,upMatrix) %���ڲ�����ʼ�壬a=10,b=3408
    pop=unifrnd(lowMatrix,upMatrix,a,b);
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
