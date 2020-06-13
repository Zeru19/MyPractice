%求k条最短路径算法
%输出：precedessors 每个结点的前点，distance 起始点到每个结点的最短距离，
%       min_path 起始点到每个结点的最短路径;
%输入：graph 图的带权邻接矩阵，source 起始点，terminal 终点
function [precedessors,distance,min_path]=kdijkstra(graph,source,terminal)
    m=size(graph,2);
    unsolved = 1:m;
    solved=zeros(1,m);
    precedessors=repelem('-',m-1,m);
    pre_num=repelem(0,m);
    min_path_num=zeros(1,m);
    min_path=repelem({cell(1,1)},m);
    distance=repelem(inf,m);
    distance(source)=0;
    pre_num(source)=1;
    min_path_num(source)=1;
    min_path{source}{1}=num2str(source);
    while (size(find(unsolved>0),2)&&sum(distance(find(unsolved>0))<inf))
        tmp = distance;
        tmp(find(solved>0))=0;
        u = find(tmp==min(distance(find(unsolved>0))));
        if length(u)>1
            u=u(1);
        end
        solved(u)=u;
        unsolved(u)=0;
        res_dist=graph(u,:);
        for j=find(unsolved>0)
            if distance(j)>distance(u)+res_dist(j)
                distance(j)=distance(u)+res_dist(j);
                pre_num(j)=1;
                precedessors(:,j)=repelem('-',m-1)';
                precedessors(1,j)=num2str(u);
                min_path_num(j)=min_path_num(u);
                min_path{j}={[min_path{u}{1},'->',num2str(j)]};
                for index = 2:length(min_path{u})
                    min_path{j}{index}=[min_path{u}{index},'->',num2str(j)];
                end
            elseif distance(j)==distance(u)+res_dist(j)
                pre_num(j)=pre_num(j)+1;
                precedessors(pre_num(j),j)=num2str(u);
                min_path_num(j)=min_path_num(u)+min_path_num(j);
                for index = 1:length(min_path{u})
                    min_path{j}{length(min_path{j})+1}=[min_path{u}{index},'->',num2str(j)];
                end
            end
        end
    end
    s=distance(terminal);
    %show_shortest_path2terminal(min_path,terminal,s);
    show_shortest_path2all(min_path,distance);
end

function show_shortest_path2terminal(min_path,terminal,s)
    disp(['最短路长度为',num2str(s)]);
    disp('最短路有');
    disp(min_path{terminal});
end

function show_shortest_path2all(min_path,distance)
    for i=1:length(min_path)
        disp(['起始点到结点',num2str(i),'的最短长度为',num2str(distance(i)),'，最短路有:']);
        for j=1:length(min_path{i})
            disp(min_path{i}{j});
        end
    end
end