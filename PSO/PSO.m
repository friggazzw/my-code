%% ��ʼ����Ⱥ  
clear
%% Sphere
clear
f= @(x) x .* sin(x) .* cos(2 * x) - 2 * x .* sin(3 * x) +3 * x .* sin(4 * x); % �������ʽ    % �������������Сֵ  


N = 20;                         % ��ʼ��Ⱥ����  
d = 1;                          % ���н�ά��    ֵ��y��ֵ
ger = 100;                      % ����������       
limit = [0, 50];               % ����λ�ò�������  
vlimit = [-10, 10];               % �����ٶ�����  
w = 0.8;                        % ����Ȩ��  
c1 = 0.5;                       % ����ѧϰ����  
c2 = 0.5;                       % Ⱥ��ѧϰ����   
figure(1);
ezplot(f,[0,limit(2)]);   %����

x = limit(1) + (  limit( 2 ) -  limit( 1)  ) .* rand(N, d);     %��ʼ��Ⱥ��λ�� rand(N, d)���Ǹ�����  

v = rand(N, d);                  % ��ʼ��Ⱥ���ٶ�  
xm = x;                          % ÿ���������ʷ���λ��   20*1�ľ���������
ym = zeros(1, d);                % ��Ⱥ����ʷ���λ��       1*1������� 0
fxm = ones(N, 1)*inf;               % ÿ���������ʷ�����Ӧ��   20*1  ��ΪҪ����Сֵ����Ϊ��֪���������ֵ�Ƕ���
                                     %ֱ������󣬺���ֻҪ��ֵ����С���Ͱ��Ǹ���С��ֵ������
                                     %�漰�������ʱ���Ȼ��20*1��N*1��Ⱥ��1*1
fym = inf;  % ��Ⱥ��ʷ�����Ӧ��
hold on  
plot(xm, f(xm), 'ro');title('��ʼ״̬ͼ');       %��ɫԲȦ���
%% Ⱥ�����  
iter = 1;  
record = zeros(ger, 1);          % ��¼��  
while iter <= ger  
     fx = f(x) ;               % ���嵱ǰ��Ӧ��    ��ֵ����fx������fxm���� 
     for i = 1:N                  %��20�����ӣ�����ѭ��20��
        if fx(i)  <fxm(i)             %iter = 1 ʱ�϶��ǲ�����forѭ����
            fxm(i) = fx(i);     % ���¸�����ʷ�����Ӧ��  
            xm(i,:) = x(i,:);   % ���¸�����ʷ���λ��(ȡֵ)  �ѵ�i�е�Ԫ�ظ���ǰ�淢��   (i,:)��i������Ԫ����ɵ�������
        end   
     end 
    if  min(fxm)  < fym 
        [fym, nmin] = min(fxm);   % ����Ⱥ����ʷ�����Ӧ��    fym��¼��Сֵ��nmin��¼�к� 
        ym = xm(nmin, :);         % ����Ⱥ����ʷ���λ��      xm�ĵ�nmin��Ԫ�� 
    end 
    v = v * w + c1 * rand * (xm - x) + c2 * rand * (repmat(ym, N, 1) - x);% �ٶȸ���  
    % �߽��ٶȴ���  
    v(v > vlimit(2)) = vlimit(2);  
    v(v < vlimit(1)) = vlimit(1);  
    x = x + v;% λ�ø���  
    % �߽�λ�ô���  
    x(x > limit(2)) = limit(2);  
    x(x < limit(1)) = limit(1);  
    record(iter) = fym;%��Сֵ��¼      
    x0 = 0 : 0.01 : limit(2);                %���ԽСͼԽ����,������
    figure(2);
    subplot(1,2,1)
    plot(x0, f(x0), 'b-', x, f(x), 'ro');title('״̬λ�ñ仯')
    subplot(1,2,2);plot(record);title('������Ӧ�Ƚ�������')  
    pause(0.1)                         %��ͣ���Ը��õĿ���Ѱ�ŵĹ���
    iter = iter+1;  

end  

figure(3);plot(x0, f(x0), 'b-', x, f(x), 'ro');title('����״̬λ��')      %����ѭ��������ִ�У�����ֵ�����յ�״̬
disp(['��Сֵ��',num2str(fym)]);  
disp(['����ȡֵ��',num2str(ym)]);  

%% ������ʼ���Ӻ��ٶ�
for i=1:sizepop
    %�������һ����Ⱥ
    pop(i,:)=5*rands(1,2);    %��ʼ��Ⱥ
    V(i,:)=rands(1,2);  %��ʼ���ٶ�
    %������Ӧ��
    fitness(i)=fun(pop(i,:));   %Ⱦɫ�����Ӧ��
end

%% ���弫ֵ��Ⱥ�弫ֵ
[bestfitness bestindex]=min(fitness);
zbest=pop(bestindex,:);   %ȫ�����
gbest=pop;    %�������
fitnessgbest=fitness;   %���������Ӧ��ֵ
fitnesszbest=bestfitness;   %ȫ�������Ӧ��ֵ

%% ����Ѱ��
for i=1:maxgen  %��������  
    
    for j=1:sizepop  %��Ⱥ��ģ
        
        %�ٶȸ���
        V(j,:) = V(j,:) + c1*rand*(gbest(j,:) - pop(j,:)) + c2*rand*(zbest - pop(j,:));
        V(j,find(V(j,:)>Vmax))=Vmax;
        V(j,find(V(j,:)<Vmin))=Vmin;
        
        %��Ⱥ����
        pop(j,:)=pop(j,:)+0.5*V(j,:);
        pop(j,find(pop(j,:)>popmax))=popmax;
        pop(j,find(pop(j,:)<popmin))=popmin;
        
        %��Ӧ��ֵ
        fitness(j)=fun(pop(j,:)); 
   
    end
    
    for j=1:sizepop
        
        %�������Ÿ���
        if fitness(j) < fitnessgbest(j)
            gbest(j,:) = pop(j,:);
            fitnessgbest(j) = fitness(j);
        end
        
        %Ⱥ�����Ÿ���
        if fitness(j) < fitnesszbest
            zbest = pop(j,:);
            fitnesszbest = fitness(j);
        end
    end 
    yy(i)=fitnesszbest;    
        
end
%% �������
plot(yy)
title('���Ÿ�����Ӧ��','fontsize',12);
xlabel('��������','fontsize',12);ylabel('��Ӧ��','fontsize',12);