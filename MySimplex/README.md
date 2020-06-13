# 运筹大作业：MATLAB实现单纯形法

## 基本信息

​		此代码由知行1801班**周泽宇**（学号**18221031**）参照教材与老师发放的PPT独立完成，除了测试老师给的模型，还测试了书上大部分例子，结果均正确且精确。但由于水平限制，除了有两处缺陷，在代码优化与性能方面也可能有所欠缺。

​		第一代`MySimplex ver1`只支持扩展形式（已引入松弛变量改写）的线性规划问题，且不支持大于等于约束；

​		第二代`MySimplex ver2`则在第一代的基础上大大改善，三种约束都支持，使用两阶段法实现。下面只说明第二代。

## 输入与输出

​		本程序支持的线性规划模型的标准形式为：
$$
\begin{aligned}
\text{Maximize}\quad &Z=\textbf{cx}\\
\text{subject to}\quad &A\textbf{x}\leq \textbf{b}\\
&D\textbf{x}=\textbf{e}
\end{aligned}
$$
其中大于等于约束取负与小于等于约束合并。

**输入**为$(\textbf{c},A,\textbf{b},D,\textbf{e})$

**输出**为$(\textbf{x}^{(0)},Z_{max})$，其中$\textbf{x}^{(0)}$为最优BF解，$Z_{max}$为取最优解时的目标函数值。

MATLAB中使用方式为

`[x0,optimal_result]=MySimplex(c,A,b,D,e)`

没有等式约束或不等式约束用空矩阵代替即可。

## 基本思路

代码的主要框架为：
1. 先将输入的$\textbf{c}$向量转换为行向量，$\textbf{b},\textbf{e}$向量转换为列向量，以便后面计算；
2. 通过函数`GetNewA`引入松弛变量、人工变量，输出改写后的扩展形式的系数矩阵与其右端向量、不带人工变量的系数矩阵、第一、二阶段的价值系数向量；
3. 再通过函数`Solve`求解第一阶段，`Solve`函数会判断模型是否有界、是否有基本可行解。最终输出最优解与取最优解时的目标函数值，若值为0则进入下一阶段，若不为0，则输出没有基本可行解；
4. 第二阶段，再用`Solve`函数解剔除掉人工变量的线性规划模型，得到最优解或其他返回值（无最优解时）。
## 不足

- 此代码寻找初始可行解的算法是，列出所有变量的组合然后选取第一个可行的基变量组合，所以在模型体量庞大时，组合数也相应很大，程序会运行很慢甚至出错；
- 在矩阵接近奇异值时算法不够稳定，使用了`AvoidUnknownError`函数来强制令足够小的数为0，可能会造成精度损失，但目前未出现错误。