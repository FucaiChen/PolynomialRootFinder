# Design doc for Spline based Reference Line

## Basics

### Definition

#### Spline

- A variable $\boldsymbol{t} \in R$

- A vector $\boldsymbol{T}$:
$$
\boldsymbol{T} = \begin{bmatrix}
1 \\ t \\ \vdots \\ t^{n-1}
\end{bmatrix}
\in R^n
$$

- Vector $\boldsymbol{P}$ and $\boldsymbol{Q}$ as a 2d spline coefficients:
  $$
  \boldsymbol{P} = 
  \begin{bmatrix}
  p_{0} \\ \vdots \\ p_{n-1}
  \end{bmatrix}
  $$
  $$
  \boldsymbol{Q} = 
  \begin{bmatrix}
  q_{0} \\ \vdots \\ q_{n-1}
  \end{bmatrix}
  $$

- A spline segment can be defined as:

$$
\boldsymbol{p} = 
\begin{bmatrix}
x \\ y
\end{bmatrix} = 
\begin{bmatrix}
TP \\ TQ
\end{bmatrix}= 
\begin{bmatrix}
T && 0 \\
0 && T
\end{bmatrix}
\begin{bmatrix}
P \\ Q
\end{bmatrix}.
\quad
\boldsymbol{p} \in R^2
$$


## reference
- https://zh.wikipedia.org/zh-hans/%E4%B8%89%E6%AC%A1%E6%96%B9%E7%A8%8B#%E6%B1%82%E6%A0%B9%E5%85%AC%E5%BC%8F%E6%B3%95
- https://zh.wikipedia.org/zh-hans/%E5%9B%9B%E6%AC%A1%E6%96%B9%E7%A8%8B#%E8%A7%A3%E5%B5%8C%E5%A5%97%E7%9A%84%E9%99%8D%E4%BD%8E%E6%AC%A1%E6%95%B0%E7%9A%84%E4%B8%89%E6%AC%A1%E6%96%B9%E7%A8%8B