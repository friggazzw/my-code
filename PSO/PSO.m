%% 初始化种群  
clear
%% Sphere
clear
f= @(x) x .* sin(x) .* cos(2 * x) - 2 * x .* sin(3 * x) +3 * x .* sin(4 * x); % 函数表达式    % 求这个函数的最小值  


N = 20;                         % 初始种群个数  
d = 1;                          % 可行解维数    值求y的值
ger = 100;                      % 最大迭代次数       
limit = [0, 50];               % 设置位置参数限制  
vlimit = [-10, 10];               % 设置速度限制  
w = 0.8;                        % 惯性权重  
c1 = 0.5;                       % 自我学习因子  
c2 = 0.5;                       % 群体学习因子   
figure(1);
ezplot(f,[0,limit(2)]);   %曲线

x = limit(1) + (  limit( 2 ) -  limit( 1)  ) .* rand(N, d);     %初始种群的位置 rand(N, d)这是个矩阵  

v = rand(N, d);                  % 初始种群的速度  
xm = x;                          % 每个个体的历史最佳位置   20*1的矩阵即列向量
ym = zeros(1, d);                % 种群的历史最佳位置       1*1矩阵标量 0
fxm = ones(N, 1)*inf;               % 每个个体的历史最佳适应度   20*1  因为要求最小值，因为不知道它的最大值是多少
                                     %直接无穷大，后面只要有值比他小，就把那个较小的值赋给它
                                     %涉及到个体的时候必然是20*1即N*1，群体1*1
fym = inf;  % 种群历史最佳适应度
hold on  
plot(xm, f(xm), 'ro');title('初始状态图');       %红色圆圈标记
%% 群体更新  
iter = 1;  
record = zeros(ger, 1);          % 记录器  
while iter <= ger  
     fx = f(x) ;               % 个体当前适应度    将值赋给fx就能与fxm比了 
     for i = 1:N                  %有20个粒子，所以循环20次
        if fx(i)  <fxm(i)             %iter = 1 时肯定是不进入for循环的
            fxm(i) = fx(i);     % 更新个体历史最佳适应度  
            xm(i,:) = x(i,:);   % 更新个体历史最佳位置(取值)  把第i行的元素赋给前面发的   (i,:)第i行所有元素组成的行向量
        end   
     end 
    if  min(fxm)  < fym 
        [fym, nmin] = min(fxm);   % 更新群体历史最佳适应度    fym记录最小值，nmin记录行号 
        ym = xm(nmin, :);         % 更新群体历史最佳位置      xm的第nmin行元素 
    end 
    v = v * w + c1 * rand * (xm - x) + c2 * rand * (repmat(ym, N, 1) - x);% 速度更新  
    % 边界速度处理  
    v(v > vlimit(2)) = vlimit(2);  
    v(v < vlimit(1)) = vlimit(1);  
    x = x + v;% 位置更新  
    % 边界位置处理  
    x(x > limit(2)) = limit(2);  
    x(x < limit(1)) = limit(1);  
    record(iter) = fym;%最小值记录      
    x0 = 0 : 0.01 : limit(2);                %间隔越小图越精密,有上限
    figure(2);
    subplot(1,2,1)
    plot(x0, f(x0), 'b-', x, f(x), 'ro');title('状态位置变化')
    subplot(1,2,2);plot(record);title('最优适应度进化过程')  
    pause(0.1)                         %暂停可以更好的看到寻优的过程
    iter = iter+1;  

end  

figure(3);plot(x0, f(x0), 'b-', x, f(x), 'ro');title('最终状态位置')      %整体循环结束才执行，所以值是最终的状态
disp(['最小值：',num2str(fym)]);  
disp(['变量取值：',num2str(ym)]);  

%% 产生初始粒子和速度
for i=1:sizepop
    %随机产生一个种群
    pop(i,:)=5*rands(1,2);    %初始种群
    V(i,:)=rands(1,2);  %初始化速度
    %计算适应度
    fitness(i)=fun(pop(i,:));   %染色体的适应度
end

%% 个体极值和群体极值
[bestfitness bestindex]=min(fitness);
zbest=pop(bestindex,:);   %全局最佳
gbest=pop;    %个体最佳
fitnessgbest=fitness;   %个体最佳适应度值
fitnesszbest=bestfitness;   %全局最佳适应度值

%% 迭代寻优
for i=1:maxgen  %进化次数  
    
    for j=1:sizepop  %种群规模
        
        %速度更新
        V(j,:) = V(j,:) + c1*rand*(gbest(j,:) - pop(j,:)) + c2*rand*(zbest - pop(j,:));
        V(j,find(V(j,:)>Vmax))=Vmax;
        V(j,find(V(j,:)<Vmin))=Vmin;
        
        %种群更新
        pop(j,:)=pop(j,:)+0.5*V(j,:);
        pop(j,find(pop(j,:)>popmax))=popmax;
        pop(j,find(pop(j,:)<popmin))=popmin;
        
        %适应度值
        fitness(j)=fun(pop(j,:)); 
   
    end
    
    for j=1:sizepop
        
        %个体最优更新
        if fitness(j) < fitnessgbest(j)
            gbest(j,:) = pop(j,:);
            fitnessgbest(j) = fitness(j);
        end
        
        %群体最优更新
        if fitness(j) < fitnesszbest
            zbest = pop(j,:);
            fitnesszbest = fitness(j);
        end
    end 
    yy(i)=fitnesszbest;    
        
end
%% 结果分析
plot(yy)
title('最优个体适应度','fontsize',12);
xlabel('进化代数','fontsize',12);ylabel('适应度','fontsize',12);