function [precedessors,distance]=dijkstra(graph,source,terminal)
    m=size(graph,2);
    unsolved = 1:m;
    solved=zeros(1,m);
    precedessors=repelem('-',m);
    distance=repelem(inf,m);
    distance(source)=0;
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
                precedessors(j)=num2str(u);
            end
        end
    end
    s=distance(terminal);
    show_path(source,terminal,precedessors,s);
end

function show_path(source,terminal,precedessors,s)
    path = num2str(terminal);
    node = terminal;
    while precedessors(node)~='-'
        path=[(precedessors(node)),'->', path];
        if precedessors(node)==num2str(source)
            break
        end
        node = str2num(precedessors(node));
    end
    disp('最短路为');
    disp(path);
    disp(['最短路长度为',num2str(s)]);
end