clear;clc;
c=[320,320,320,400,400,400,360,360,360,290,290,290];
A=[20,0,0,16,0,0,25,0,0,13,0,0;
    0,20,0,0,16,0,0,25,0,0,13,0;
    0,0,20,0,0,16,0,0,25,0,0,13;
    500,0,0,700,0,0,600,0,0,400,0,0;
    0,500,0,0,700,0,0,600,0,0,400,0;
    0,0,500,0,0,700,0,0,600,0,0,400];
b=[12,18,10,700,900,5000];
D=[-3,2,0,-3,2,0,-3,2,0,-3,2,0;
    0,5,-9,0,5,-9,0,5,-9,0,5,-9];
e=[0,0];
[x0,op]=MySimplex2(c,A,b,D,e);
disp('最优解x*=');
disp(x0);
disp(['此时的目标函数值为',num2str(op)]);