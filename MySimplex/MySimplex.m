 %MySimplex ver2.2
%求解最大化问题，支持等式、大于等于、小于等于约束，标准形式为
%    Maximize   cx
%       s.t.    Ax<= b
%               Dx = e
%输入为c, A, b, D, e
%输出为最优解x0，取最优解时的目标函数值optimal_result
%该程序有两处不足在程序注释中有说明
%@作者：周泽宇 @学号：18221031 
%@日期：2020/4/2
function [x0,optimal_result]=MySimplex(c,A,b,D,e)
    %保证c是行向量，b是列向量
    if size(c,2)==1
        c=c';
    end
    if size(b,1)==1
        b=b';
    end
    if size(e,1)==1
        e=e';
    end
    [m1,m2,n1,n2,b2,cB,cB2,new_A]=GetNewA(c,A,b,D,e);
    %第一阶段
    [x1,counter]=Solve(new_A,b2,cB);
    %判断是否陷入死循环
    
    %此程序寻找初始基变量的方法是列出所有可能的基变量的全排列，然后遍历该排列找到第一个可行的基矩阵
    %所以当约束条件体量庞大时该程序也失去意义，为此程序不足之一
    if counter>nchoosek(n1+m1+m2+n2,m1+m2) 
        x0='No optimal solution';
        optimal_result='None';
        return
    else
        optimal_result1=cB*x1;
    end
    %第二阶段
    if optimal_result1~=0
        disp('没有基本可行解');
        optimal_result='None';
        return
    else
        A_without_art=new_A(:,1:n1+m1);
        [x0,counter]=Solve(A_without_art,b2,cB2);
        if counter>nchoosek(n1+m1+m2+n2,m1+m2)
            x0='No optimal solution';
            optimal_result='None';
            return
        else
            optimal_result=cB2*x0;
        end
    end
end

%核心程序，输入为扩展形式的等式约束和价值系数，输出最优解与迭代次数
function [x0,counter]=Solve(A,b,c)
    m=size(A,1);
    n=size(A,2);
    %找到起始可行解
    all_possible_basic=nchoosek(1:n,m);
    for i=1:nchoosek(n,m)
        B=A(:,all_possible_basic(i,:));
       if det(B)<=0
           if i==nchoosek(n,m)
               disp('没有找到可行的B，请检查模型');
           end
           continue
       else
           if sum(B\b<0)>0
               if i==nchoosek(n,m)
                    disp('没有找到可行的B，请检查模型');
               end
               continue
           else
               basic=all_possible_basic(i,:);%将基变量的序号储存在basic中
               break
           end
       end
    end
    cb=c(all_possible_basic(i,:));
    x0=zeros(n,1);
    x0(basic)=B\b;
    counter=0;
    while(1)
        counter=counter+1;
        if counter>nchoosek(n,m)
            disp('程序陷入了循环，请检查模型');
            break
        end
        %找到最大的检验数确定进基变量
        z=cb/B*A-c;
        %矩阵中零元素过多或其他情况会导致计算精度降低，用AvoidUnknownError来避免这个玄学错误
        %此处为本程序的不足之二
        z=AvoidUnknownError(z);
        b_bar=B\b;
        [zmin,k]=FindFirstMin(z);
        %disp(zmin);
        if zmin>=0
            break
        else
            y=B\A(:,k);
            %判断问题是否无界，但如果是最后一次迭代产生了最优解，而此时若y也小于0，仍会输出这条语句
            if sum(y>0)==0
                disp('该问题无界');
                disp('（但如果最后一次迭代产生了最优解，而y也小于0，仍会输出这条语句但存在可行解）');
                return
            else
                temp=ones(1,m);
                for j=1:m
                    if y(j)>0
                        temp(j)=b_bar(j)./y(j);
                    else
                        temp(j)=inf;
                    end
                end
                %找到最小比值确定离基变量
                [~,r]=FindFirstMin(temp);
                basic(r)=k;
                B=A(:,basic);
                cb=c(basic);
                if det(B)==0
                    disp('B又突然不可逆了');%Debug,现实应该不会发生
                end
                x0=zeros(n,1);
                x0(basic)=B\b;
                x0=AvoidUnknownError(x0);
            end
        end
    end
end

%找到向量第一个最小的元素，避免循环
function [min,loc]=FindFirstMin(y)
    n=size(y,2);
    min=y(1);loc=1;
    for i=1:n
        if y(i)<min
            min=y(i);loc=i;
        end
    end
end

%得到引入松弛变量与人工变量后的扩展形式的系数矩阵
function [m1,m2,n1,n2,b2,cB,cB2,new_A]=GetNewA(c,A,b,D,e)
    m1=size(A,1);
    m2=size(D,1);
    if size(A,2)>0
        n1=size(A,2);
    else
        n1=size(D,2);
    end
    slack_A=[[A;D],[diag(ones(m1,1));zeros(m2,m1)]];
    greater_than=find(b<0);
    b(greater_than)=-b(greater_than);
    slack_A(greater_than,:)=-slack_A(greater_than,:);
    n2=size(greater_than,1);
    if n2>0
        art_A=[slack_A,zeros(m1+m2,n2)];
        for i=greater_than
            j=1;
            art_A(i,m1+n1+j)=1;
            j=j+1;
        end
    else
        art_A=slack_A;
    end
    
    if m2>0
        new_A=[art_A,[zeros(m1,m2);diag(ones(m2,1))]];
    else
        new_A=art_A;
    end
    b2=[b;e];
    cB2=[c,zeros(1,m1)];
    %第一阶段
    cB=[zeros(1,n1+m1),-1*ones(1,m2+n2)];
end

function y=AvoidUnknownError(x)
    flag=0;
    if size(x,2)==1
        flag=1;
        x=x';
    end
    for i =1:size(x,2)
        if abs(x(i))<1e-6
            x(i)=0;
        end
        y=x;
    end
    if flag==1
        y=y';
    end
end

        