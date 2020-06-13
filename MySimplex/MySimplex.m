 %MySimplex ver2.2
%���������⣬֧�ֵ�ʽ�����ڵ��ڡ�С�ڵ���Լ������׼��ʽΪ
%    Maximize   cx
%       s.t.    Ax<= b
%               Dx = e
%����Ϊc, A, b, D, e
%���Ϊ���Ž�x0��ȡ���Ž�ʱ��Ŀ�꺯��ֵoptimal_result
%�ó��������������ڳ���ע������˵��
%@���ߣ������� @ѧ�ţ�18221031 
%@���ڣ�2020/4/2
function [x0,optimal_result]=MySimplex(c,A,b,D,e)
    %��֤c����������b��������
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
    %��һ�׶�
    [x1,counter]=Solve(new_A,b2,cB);
    %�ж��Ƿ�������ѭ��
    
    %�˳���Ѱ�ҳ�ʼ�������ķ������г����п��ܵĻ�������ȫ���У�Ȼ������������ҵ���һ�����еĻ�����
    %���Ե�Լ�����������Ӵ�ʱ�ó���Ҳʧȥ���壬Ϊ�˳�����֮һ
    if counter>nchoosek(n1+m1+m2+n2,m1+m2) 
        x0='No optimal solution';
        optimal_result='None';
        return
    else
        optimal_result1=cB*x1;
    end
    %�ڶ��׶�
    if optimal_result1~=0
        disp('û�л������н�');
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

%���ĳ�������Ϊ��չ��ʽ�ĵ�ʽԼ���ͼ�ֵϵ����������Ž����������
function [x0,counter]=Solve(A,b,c)
    m=size(A,1);
    n=size(A,2);
    %�ҵ���ʼ���н�
    all_possible_basic=nchoosek(1:n,m);
    for i=1:nchoosek(n,m)
        B=A(:,all_possible_basic(i,:));
       if det(B)<=0
           if i==nchoosek(n,m)
               disp('û���ҵ����е�B������ģ��');
           end
           continue
       else
           if sum(B\b<0)>0
               if i==nchoosek(n,m)
                    disp('û���ҵ����е�B������ģ��');
               end
               continue
           else
               basic=all_possible_basic(i,:);%������������Ŵ�����basic��
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
            disp('����������ѭ��������ģ��');
            break
        end
        %�ҵ����ļ�����ȷ����������
        z=cb/B*A-c;
        %��������Ԫ�ع������������ᵼ�¼��㾫�Ƚ��ͣ���AvoidUnknownError�����������ѧ����
        %�˴�Ϊ������Ĳ���֮��
        z=AvoidUnknownError(z);
        b_bar=B\b;
        [zmin,k]=FindFirstMin(z);
        %disp(zmin);
        if zmin>=0
            break
        else
            y=B\A(:,k);
            %�ж������Ƿ��޽磬����������һ�ε������������Ž⣬����ʱ��yҲС��0���Ի�����������
            if sum(y>0)==0
                disp('�������޽�');
                disp('����������һ�ε������������Ž⣬��yҲС��0���Ի����������䵫���ڿ��н⣩');
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
                %�ҵ���С��ֵȷ���������
                [~,r]=FindFirstMin(temp);
                basic(r)=k;
                B=A(:,basic);
                cb=c(basic);
                if det(B)==0
                    disp('B��ͻȻ��������');%Debug,��ʵӦ�ò��ᷢ��
                end
                x0=zeros(n,1);
                x0(basic)=B\b;
                x0=AvoidUnknownError(x0);
            end
        end
    end
end

%�ҵ�������һ����С��Ԫ�أ�����ѭ��
function [min,loc]=FindFirstMin(y)
    n=size(y,2);
    min=y(1);loc=1;
    for i=1:n
        if y(i)<min
            min=y(i);loc=i;
        end
    end
end

%�õ������ɳڱ������˹����������չ��ʽ��ϵ������
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
    %��һ�׶�
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

        