### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# ╔═╡ f17103ea-06bf-11f1-a2b0-79e68ed152eb
md"""# Project_02 - Multibody kinematic modeling

![Dual slider kinematics project](https://raw.githubusercontent.com/cooperrc/me5180-project_02/refs/heads/main/dual-slider.svg)

In this project, a rigid bar is connected to two sliding pistons along
the diagonal tracks. As the pistons move along the tracks, the rigid bar rotates at a constant rate, $\dot{\theta}_3 = 2~rad/s$. The figure above has three _relative_ ccoordinate systems that move with the bodies:

1.  $x_1-y_1-$ describes piston 1 position and orientation, $\theta_1$
2.  $x_2-y_2-$ describes piston 2 position and orientation, $\theta_2$
3.  $x_3-y_3-$ describes the rigid bar position and orientation, $\theta_3$

Each of the pistons are on tracks at $\pm 45^o$ and the rotating rigid
bar is 10 cm. The hinges are mounted to the center of the pistons
connecting the ends of the rigid bar. 
 
In this project, you need to 

1. determine constraint equations $C(\mathbf{q},~t)$
2. solve for the velocities, $\dot{q}$ and accelerations, $\ddot{q}$
3. visualize the motion of the system as the rigid bar goes through at least one full rotation
"""

# ╔═╡ 0db4e504-bf02-4b66-928e-d08661f3b9f2
md"""
## Global Reference Frame
We can initialize a global reference frame $[x_G, y_G]$ with directions $[\hat{i}_0, \hat{j}_0]$, with origin at the crossing point of the two tracks.

The piston locations and the rigid bar location can be described using their respective local coordinates and transformed using $\mathbf{A}_i$ to the global coordinates for each local frame $i$

$$\mathbf{A}_i = \begin{bmatrix} \cos\theta_i & -\sin\theta_i \\ \sin\theta_i & \cos\theta_i \end{bmatrix}$$


**Piston 1** is in the $[x_1, y_1]$ reference frame which applies a fixed $\theta_1 = 45°$ rotation from $x_G$ to $x_1$. The motion of the piston is constrained along $x_1$ so that $y_1 = 0$. Therefore, the motion of piston 1 can be described as $[s_1, 0]$ in the $[x_1, y_1]$ reference frame and
$$\mathbf{r}_1 = \mathbf{A}_1 \begin{bmatrix} s_1 \\ 0 \end{bmatrix}$$
in the $[x_G, y_G]$ reference frame, where $s_1$ is the displacement along the track.

**Piston 2** is in the $[x_2, y_2]$ reference frame which applies a fixed $\theta_2 = 45°$ rotation from $x_G$ to $x_2$. The motion of the piston is constrained along $x_2$ so that $y_2 = 0$. Therefore, the motion of piston 2 can be described as $[s_2, 0]$ in the $[x_2, y_2]$ reference frame and
$$\mathbf{r}_2 = \mathbf{A}_2 \begin{bmatrix} s_2 \\ 0 \end{bmatrix}$$ in the $[x_G, y_G]$ reference frame, where $s_2$ is the displacement along the track.


The **rigid bar's** center of mass is located in the $[x_3, y_3]$ reference frame which applies a $\theta_3(t)$ rotation from $x_G$ to $x_3$. The bar has a length of $L$. The rigid bar constrains the motion of the pistons so that the pistons are always a length $L$ apart. In the bar's own reference frame piston 1 sits at $[-L/2, 0]$ and piston 2 sits at $[+L/2, 0]$, so in the global frame they are 

$$\mathbf{r}_{B_1} = \mathbf{r}_3 + \mathbf{A}_3 \begin{bmatrix} -L/2 \\ 0 \end{bmatrix}, \qquad \mathbf{r}_{B_2} = \mathbf{r}_3 + \mathbf{A}_3 \begin{bmatrix} +L/2 \\ 0 \end{bmatrix}$$

This applies the additional constraints $\mathbf{r}_{B_1} = \mathbf{r}_1$ and $\mathbf{r}_{B_2} = \mathbf{r}_2$, which fix both piston positions once $\theta_3$ is known, and equivalently enforce $$|\mathbf{r}_2 - \mathbf{r}_1| = L$$
"""
#// Brooke Parker Thibodeau 4/3/26 //

# ╔═╡ df624f7f-8b4c-42bc-9fa7-70329f585f48
md"""
## Constraint Equations
The constraints of the system are written as $\mathbf{C}(\mathbf{q}, t) = \mathbf{0}$. The hinges connected the pistons to the rigid bar apply constraints:

$$\mathbf{r}_{B_1} - \mathbf{r}_1 = \mathbf{0} \qquad \rightarrow \qquad 
\begin{bmatrix} x_3 \\ y_3 \end{bmatrix} + \mathbf{A}_3 \begin{bmatrix} -L/2 \\ 0 \end{bmatrix} - \mathbf{A}_1 \begin{bmatrix} s_1 \\ 0 \end{bmatrix} = \begin{bmatrix} 0 \\ 0 \end{bmatrix}$$

$$\mathbf{r}_{B_2} - \mathbf{r}_2 = \mathbf{0} \qquad \rightarrow \qquad \begin{bmatrix} x_3 \\ y_3 \end{bmatrix} + \mathbf{A}_3 \begin{bmatrix} +L/2 \\ 0 \end{bmatrix} - \mathbf{A}_2 \begin{bmatrix} s_2 \\ 0 \end{bmatrix}= \begin{bmatrix} 0 \\ 0 \end{bmatrix}$$

The problem also specifies that the rigid bar rotation is fixed at $\dot{\theta}_3 = 2rad/s$, which applies another constraint:

$$\theta_3 - \dot{\theta}_3 t = 0$$

This can be summarized in one constraint matrix where $\mathbf{q} = [s_1,\ s_2,\ x_3,\ y_3,\ \theta_3]^T$ and $\theta_1$, $\theta_2$ are fixed parameters.

$$\mathbf{C}(\mathbf{q}, t) = \begin{bmatrix} x_3 - \dfrac{L}{2}\cos\theta_3 - s_1\cos\theta_1 \\ y_3 - \dfrac{L}{2}\sin\theta_3 - s_1\sin\theta_1 \\ x_3 + \dfrac{L}{2}\cos\theta_3 - s_2\cos\theta_2 \\ y_3 + \dfrac{L}{2}\sin\theta_3 - s_2\sin\theta_2 \\ \theta_3 - \dot{\theta}_3\, t \end{bmatrix} = \mathbf{0}$$
"""
#// Brooke Parker Thibodeau 4/3/26 //

# ╔═╡ 9179eed4-fc3c-40b2-970c-2f3f4ff8cc53
begin
	L  = 0.10;        # bar length [m]
	ω  = 2.0;         # prescribed angular velocity [rad/s]
	θ1 = π/4;         # piston 1 track angle [rad]
	θ2 = -π/4;        # piston 2 track angle [rad]
end
#// Brooke Parker Thibodeau 4/3/26 //

# ╔═╡ 8a6c526f-80f9-435d-93c3-0b802d886b7d
function C(q, t)
	s1, s2, x3, y3, θ3 = q
	return [
		x3 - (L/2)*cos(θ3) - s1*cos(θ1),
		y3 - (L/2)*sin(θ3) - s1*sin(θ1),
		x3 + (L/2)*cos(θ3) - s2*cos(θ2),
		y3 + (L/2)*sin(θ3) - s2*sin(θ2),
		θ3 - ω*t
	]
end
#// Brooke Parker Thibodeau 4/3/26 //

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.12.4"
manifest_format = "2.0"
project_hash = "71853c6197a6a7f222db0f1978c7cb232b87c5ee"

[deps]
"""

# ╔═╡ Cell order:
# ╠═f17103ea-06bf-11f1-a2b0-79e68ed152eb
# ╠═0db4e504-bf02-4b66-928e-d08661f3b9f2
# ╠═df624f7f-8b4c-42bc-9fa7-70329f585f48
# ╠═9179eed4-fc3c-40b2-970c-2f3f4ff8cc53
# ╠═8a6c526f-80f9-435d-93c3-0b802d886b7d
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
