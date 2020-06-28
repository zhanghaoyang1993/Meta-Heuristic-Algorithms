function [Fit_and_p,FVr_bestmemit, fitMaxVector] = ...
    PSO(deParameters,caseStudyData,otherParameters,low_habitat_limit,up_habitat_limit,iRuns)
%%
%-----This is just for notational convenience and to keep the code uncluttered.--------
I_NP         = deParameters.I_NP; % ��Ⱥ�� = 10
F_weight     = deParameters.F_weight; %Mutation factor 0.5
F_CR         = deParameters.F_CR; %Recombination constant 0.9
I_D          = numel(up_habitat_limit); %Number of variables or dimension 3408
deParameters.nVariables=I_D;
popmin = low_habitat_limit; %�������� 1*3048
popmax = up_habitat_limit; %��������  1*3048
I_itermax    = deParameters.I_itermax; %number of max iterations/gen
fnc=  otherParameters.fnc; %ѡ��fitness function

fitMaxVector = nan(1,I_itermax); %1*500��NaN����
%%
%������ʼ��
%����Ⱥ�㷨�е���������

c1 = 1.49445; %�������ӣ���������ʷ��ѷ����ƶ�
c2 = 1.49445; %�������ӣ�������10���������ʷ��ѷ����ƶ�
w_max=1.4;
w_min=0.4;

maxg=deParameters.I_itermax; %�������� �ƶ�500��
sizepop=I_NP; %��Ⱥ��ģ 10ֻ��
%��ʼ�ٶȺ���Ⱥ���±߽�ֵ
popminMatrix=repmat(popmin,I_NP,1); %����10*3408�ľ���ÿ����ͬ��ÿ��Ϊ��������
popmaxMatrix=repmat(popmax,I_NP,1); %����10*3408�ľ���ÿ����ͬ��ÿ��Ϊ��������
deParameters.minPositionsMatrix=popminMatrix; %����deParameters��
deParameters.maxPositionsMatrix=popmaxMatrix; %����deParameters��
Vgap=(popmax-popmin)/3; %�ٶȷ�Χ,ͨ�������������趨���ϵ����ҷ���ô�찡���ٶ���������Ҫ�Ż���
VminMatrix=repmat(-Vgap,I_NP,1); %����10*3408�ľ���ÿ����ͬ��ÿ��Ϊ��������
VmaxMatrix=repmat(Vgap,I_NP,1); %����10*3408�ľ���ÿ����ͬ��ÿ��Ϊ��������

%%������ʼ���Ӻ��ٶȣ��ϵ�˵Ҫ���񣬾�������
% generate initial population.������ɵĵ�һ��10*3408���Ӵ���10Ϊ�Ӵ�������3408Ϊ��������
rand('state',otherParameters.iRuns) %Guarantee same initial population
pop=genpop(I_NP,I_D,popminMatrix,popmaxMatrix); %������ɵĵ�һ��10*3408���Ӵ���10Ϊ�Ӵ�������3408Ϊ��������
V=genpop(I_NP,I_D,VminMatrix,VmaxMatrix); %������ɵĵ�һ��10*3408�Ӵ����ٶȣ�10Ϊ�Ӵ�������3408Ϊ��������
[fitness, ~]=feval(fnc,pop,caseStudyData,otherParameters);%10���Ӵ���fitness functionֵ��10*1

%%����õ�Ⱦɫ�壬����С��鰡
[bestfitness,bestindex]=min(fitness);
zbest=pop(bestindex,:); %ȫ�����
gbest=pop; %�������
fitnessgbest=fitness; %���������Ӧ��
fitnesszbest=bestfitness; %ȫ�������Ӧ��

fitMaxVector(1,1) = fitnesszbest;%����ʼ����Ӵ���fitnessֵ����fitMaxVector�У�1*500
FVr_bestmemit = pop(bestindex,:);
%%����Ѱ�ţ���Ҫ��499���������ڼҷ�����
for i=1:maxg-1 %�ƶ�499��
    w=w_max-i*(w_max-w_min)/maxg;%Ȩ�����Եݼ���PSO�㷨
    for j = 1:sizepop %10������
        %�ٶȸ��£�Ҫ����óԵ��Ⱑ
        V(j,:) = w*V(j,:) + c1*rand*(gbest(j,:)-pop(j,:))+c2*rand*(zbest-pop(j,:));
        V(j,: )= update(V(j,: ),VminMatrix(j,:),VmaxMatrix(j,:));
        %�ƶ�
        pop(j,:)=pop(j,:)+V(j,:);
        pop(j,: )= update(pop(j,: ), popminMatrix(j,:), popmaxMatrix(j,:));
        %����Ӧ���죬��������İ�
        if rand>0.9
            k=ceil(I_D*rand);
            pop(j,k)=popmin(k)+(popmax(k)-popmin(k))*rand;
        end
    end   
    %��Ӧ��ֵ��������С�ĵ㣨��
    [fitness, ~]=feval(fnc,pop,caseStudyData,otherParameters);
    %�������Ÿ��£����ҵ���óԵ���Ҫʵʱ׷�ٰ�
    for j=1:sizepop %10������
        if fitness(j) < fitnessgbest(j)
            gbest(j,:) = pop(j,:);
            fitnessgbest(j) = fitness(j);
        end
    end
    %Ⱥ�����Ÿ��£�����һֻ���ڹ۲�͹�ͨ���ɹ�������Ư��Ů���������ŵ���
    [smallest_fit,smallest_index]=min(fitness);
    if smallest_fit< fitnesszbest
        zbest = pop(smallest_index,:);
        fitnesszbest = fitness(smallest_index);
    end
    fitMaxVector(1,i+1) = fitnesszbest;%����ʼ����Ӵ���fitnessֵ����fitMaxVector�У�1*500 
    FVr_bestmemit = zbest; % best member of current iteration
end
Fit_and_p=[fitMaxVector(1,i+1) 0];
%���������20������ͬһ���ط������������˺���
figure(1)
plot(fitMaxVector,'Linewidth',2)%
title(['��Ӧ������ ' '��ֹ����=' num2str(maxg)]);
grid on
xlabel('��������' ); ylabel('��Ӧ��');
%������
%zbest; %��Ѹ���
%fitnesszbest; %����ֵ
end

% VECTORIZED THE CODE INSTEAD OF USING FOR ������ʼ��
function pop=genpop(a,b,lowMatrix,upMatrix) %���ڲ�����ʼ�壬a=10,b=3408
    pop=unifrnd(lowMatrix,upMatrix,a,b);
end

%�ж��²��ĸ�����ٶ��Ƿ����������������Ҫ�󣬽���������滻�������޵�ֵ
function p=update(p,lowMatrix,upMatrix)
        %[popsize,dim]=size(p);
        [idx] = find(p<lowMatrix);
        p(idx)=lowMatrix(idx);
        [idx] = find(p>upMatrix);
        p(idx)=upMatrix(idx);
end