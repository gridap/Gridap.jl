Let us consider the transient problem $A(u,\dot{u}) = 0$. We can consider that this is a ODE, since the spatial part is not important here. The most standard situation is the one in which
$$
A(u,d_t{u}) = Md_t{u} + K u,
$$
i.e., $A$ is linear with respect to $\dot{u}$, but we can consider the more general case here. 

Our motivation here is to split the time domain into time steps, and for each time step $[t^n,t^{n+1}]$, create an approximation of the map $u^n \mapsto u^{n+1}$ such that $A(R(u^{n+1},u^n),\Delta_t(u^{n+1},u^n)) = 0$. The operator $\Delta_t$ is an approximation of the time derivative and the operator $R$ is some time approximation of $u$ in $[t^n,t^{n+1}]$.  

Let us consider the $θ$-method to fix ideas. In this case, we consider at each time step the initial value $u^n$, we approximate the time derivative using finite differences as $\Delta_t(u^n,u^{n+1}) = \frac{u^{n+1}-u^{n}}{\delta t}$. For Backward-Euler, we compute $R(u^n,u^{n+1}) = u^{n+1}$, for Forward Euler $R(u^n,u^{n+1}) = u^{n}$, and for Crank-Nicolson $R(u^n,u^{n+1}) = \frac{u^{n+1}+u^{n}}{2}$ (or more precisely, $R(u^n,u^{n+1}) = u^{n+1}(t-t^n) + u^n(t^{n+1}-t)$).

Now, we want to approximate the operator using any of these methods. We can readily check that we can write the Newton linearisation of any of these problems as: 
$$
[\frac{∂A}{∂u}\frac{\partial R}{\partial u^{n+1}}  + \frac{∂A}{∂\dot{u}}\frac{\partial \Delta}{\partial u^{n+1}} ] \delta u^{n+1}  = - A(R(u^{n+1},u^n),\Delta_t(u^{n+1},u^n)).
$$
We can denote $J_0 \doteq \frac{∂A}{∂u}$ and $J_1 \doteq \frac{∂A}{∂\dot{u}}$.

E.g., for BE we have 
$$
\frac{\partial R}{\partial u^{n+1}} = 1, \qquad \frac{\partial \Delta}{\partial u^{n+1}} = 1/δt,
$$
for FE we have
$$
\frac{\partial R}{\partial u^{n+1}} = 0, \qquad \frac{\partial \Delta}{\partial u^{n+1}} = 1/δt,
$$
and for CN we have 
$$
\frac{\partial R}{\partial u^{n+1}} = 1/2, \qquad \frac{\partial \Delta}{\partial u^{n+1}} = 1/δt.
$$
Analogously, we can define $\gamma_0 \doteq \frac{\partial R}{\partial u^{n+1}}$ and $\gamma_1 \doteq \frac{\partial \Delta}{\partial u^{n+1}}$. Note that the standard Jacobian is pre-multiplied by $\gamma_0$. Thus, it makes no sense to compute $J_0$ when $\gamma_0 = 0$. The same happens for the case of $\gamma_1 = 0$ and $J_1$. But this is not the case for ODEs. $J_1$ is the mass matrix.      

We have decided to write the problem in terms of $u^{n+1}$, but we could do something different. E.g., for the $\theta$-method, we can write the problem in terms of $u^{n+\theta}$ for $\theta > 0$. In this case, we have to define the $R$ and $Δ_t$ operators in terms of $u^{n+\theta}$ and $u^n$ and perform exactly as above. The only difference is a final update step. 

From the discussion above, it seems quite clear that FE should work with the current machinery, since we would be computing $M + 0*L$. That is the reason why I say that the current machinery should work for FE but we need to avoid computing $L$ by using a if statement in the code, and only compute it for $\gamma_0 > 0$.

When FE works, we can start the discussion about the general RK solver.
