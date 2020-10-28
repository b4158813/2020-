$$
f(\vec r+d\vec r,\vec \xi + \vec adt,t+dt)d\vec r d\vec \xi = f(\vec r,\vec \xi, t)d\vec r d\vec \xi
$$


$$
(\frac{\partial f}{\partial t})_{运动} = -\vec \xi \cdot \frac{\partial f}{\partial \vec r} - \vec a \cdot \frac{\partial f}{\partial \vec \xi}
$$

$$
\frac{\partial f}{\partial t} = (\frac{\partial f}{\partial t})_{运动} + (\frac{\partial f}{\partial t})_{碰撞}
$$

$$
\frac{\partial f}{\partial t} + \vec \xi \cdot \frac{\partial f}{\partial \vec r} + \vec a \cdot \frac{\partial f}{\partial \vec \xi} = (\frac{\partial f}{\partial t})_{碰撞}
$$

$$
\frac{\partial f}{\partial t} + \vec \xi \cdot \frac{\partial f}{\partial \vec r} + \vec a \cdot \frac{\partial f}{\partial \vec \xi} = \iint(f'f_1'-ff_1)d_D|\vec g| \cos \theta d\Omega d\vec {\xi_1}
$$

$$
f = n\frac{1}{(2\pi R_gT)^{\frac{D}{2}}}\exp{[-\frac{(\vec \xi-\vec u)^2}{2R_gT}]}
$$

$$
\frac{\partial(n\overline \phi)}{\partial t}+\frac{\partial}{\partial \vec r}\cdot (n\overline{\phi\vec{\xi}}) - n(\frac{\partial \overline \phi}{\partial t} + \overline{\vec\xi \cdot \frac{\partial \phi}{\partial \vec r}} + \vec a\cdot \overline{\frac{\partial \phi}{\partial \vec\xi}}) = \int\phi J(ff_1)d\xi
$$

$$
\frac{\partial f}{\partial t} + \vec \xi \cdot \frac{\partial f}{\partial \vec r} + \vec a \cdot \frac{\partial f}{\partial \vec \xi} = \Omega_f = \nu(f^{eq}-f)
$$

$$
\tau_0 = \frac{1}{\nu}
$$

$$
\frac{\partial f_\alpha}{\partial t} + \vec e\cdot \nabla f_\alpha = -\frac{1}{\tau_0}(f_\alpha-f_\alpha^{eq}) + F_{\alpha}
$$

$$
f_\alpha^{eq} = \rho \omega_\alpha\left[1 + \frac{\vec e_\alpha \cdot \vec u}{2R_gT} + \frac{(\vec e_\alpha\cdot \vec u)^2}{2R_g^2T^2} - \frac{u^2}{2R_gT}\right] + O(u^3)
$$

$$
f_\alpha(\vec r+\vec e_\alpha \delta t,t + \delta t) - f_\alpha(\vec r,t) = -\frac{1}{\tau}\left[f_\alpha(\vec r,t) - f_\alpha^{eq}(\vec r,t)\right] + \delta F_\alpha(\vec r,t)
$$

$$
\frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \vec u) =0
$$

$$
\frac{\partial \vec u}{\partial t} + \nabla \cdot ( \vec u\vec u) = - \frac{1}{\rho_0}\nabla p + \nu\nabla\cdot(\nabla \vec u + (\nabla \vec u)^T)
$$

