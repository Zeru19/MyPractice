graph=[0 1 2 3 4 5 inf
       inf 0 inf inf inf inf 6
       inf inf 0 inf inf inf 5
       inf inf inf 0 inf inf 3
       inf inf inf inf 0 inf 2
       inf inf inf inf inf 0 1
       inf inf inf inf inf inf 0];

seervada = [0 2 5 4 inf inf inf
            2 0 2 inf 7 inf inf
            5 2 0 1 4 3 inf
            4 inf 1 0 inf 4 inf
            inf 7 4 inf 0 1 5
            inf inf 3 4 1 0 7
            inf inf inf inf 5 7 0];
        
test = [ 0   4   3   6  inf inf inf inf inf inf inf 
         4   0  inf  5   3  inf inf inf inf inf inf 
         3  inf  0   4  inf  6  inf inf inf inf inf 
         6   5   4   0   2   5   2  inf inf inf inf 
        inf  3  inf  2   0  inf  2   4  inf inf inf 
        inf inf  6   5  inf  0   1  inf  2   5  inf 
        inf inf inf  2   2   1   0   2   5  inf inf 
        inf inf inf inf  4  inf  2   0   2  inf  7  
        inf inf inf inf inf  2   5   2   0   3   8  
        inf inf inf inf inf  5  inf inf  3   0   4  
        inf inf inf inf inf inf inf  7   8   4   0 ];
%[precedessors,distance]=dijkstra(graph,1,7);
[~,distance,min_path]=kdijkstra(test,1,11);
