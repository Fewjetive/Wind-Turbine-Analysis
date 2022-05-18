# Final Report

## Data Files

1. NACA 6409
   * [Profile](/data/NACA%206409/n6409.dat "Title")
   * angle of attack vs. drag coefficient
     1. [Re=10000, Ncrit=5](/data/NACA%206409/xf-n6409-il-100000-n5.txt)
     2. [Re=10000, Ncrit=9](/data/NACA%206409/xf-n6409-il-1000000.txt)
     3. [Re=50000, Ncrit=5](/data/NACA%206409/xf-n6409-il-50000-n5.txt)
     4. [Re=50000, Ncrit=9](/data/NACA%206409/xf-n6409-il-50000.txt)
     5. [Re=200000, Ncrit=5](/data/NACA%206409/xf-n6409-il-200000-n5.txt)
     6. [Re=200000, Ncrit=9](/data/NACA%206409/xf-n6409-il-200000-n9.txt)
     7. [Re=500000, Ncrit=5](/data/NACA%206409/xf-n6409-il-500000-n5.txt)
     8. [Re=500000, Ncrit=9](/data/NACA%206409/xf-n6409-il-500000.txt)
     9. [Re=1000000, Ncrit=5](/data/NACA%206409/xf-n6409-il-1000000-n5.txt)
     10. [Re=1000000, Ncrit=9](/data/NACA%206409/xf-n6409-il-1000000.txt)
2. SG 6043
* [Profile](/data/SG%206043/sg6043.dat.txt "Title")
   * angle of attack vs. drag coefficient
     1. [Re=10000, Ncrit=5](/data/NACA%206409/xf-n6409-il-10000-n5.txt)
     2. [Re=10000, Ncrit=9](/data/SG%206043/xf-sg6043-il-10000.txt)
     3. [Re=50000, Ncrit=5](/data/SG%206043/xf-sg6043-il-50000-n5.txt)
     4. [Re=50000, Ncrit=9](/data/SG%206043/xf-sg6043-il-50000.txt)
     5. [Re=200000, Ncrit=5](/data/SG%206043/xf-sg6043-il-100000-n5.txt)
     6. [Re=200000, Ncrit=9](/data/SG%206043/xf-sg6043-il-200000.txt)
     7. [Re=500000, Ncrit=5](/data/SG%206043/xf-sg6043-il-500000-n5.txt)
     8. [Re=500000, Ncrit=9](/data/SG%206043/xf-sg6043-il-500000.txt)
     9. [Re=1000000, Ncrit=5](/data/SG%206043/xf-sg6043-il-1000000-n5.txt)
     10. [Re=1000000, Ncrit=9](/data/SG%206043/xf-sg6043-il-1000000.txt)
## Variables

### 1. Given or Specified

#### (a) Constants

* $B$: #blades, set to $B=3$ $\checkmark$
* $U_1$: wind speed, set $U_1=20m/s$ $\checkmark$
* $R$: maximum radius, set to $R=1$ $\checkmark$
* $n$: # roatations per seconds, assigned arbitrarily. Let $n=1$ $\checkmark$
* $D=2R$, diameter of turbine

#### (b) Ranging in Numerical Interval: To be Investigated

* $r$: current displacement relative to hub, ranging from $0$ to $R$ $\checkmark$
* $\Theta$: pitch angle, specified by ourselve. Ranging from $0$ to $90$$\checkmark$
* $\alpha$: angle of attack, ranging from $-8.75$ to $19.25$ with incremental of $0.25$ $\checkmark$

#### (c) Function of $r$

* $c$: camber thickness, given by turbine profile and also function of $r$. Need to be specified. $\checkmark$
* $\sigma=\frac{Bc}{2\pi r}$, solidity

#### (d) From Database

* $C_L$: lift coefficient, given by database with respect to $R_e$, Ncrit and $\alpha$ $\checkmark$
* $C_D$: lift coefficient, given by database with respect to $R_e$, Ncrit and $\alpha$ $\checkmark$

### 2. Easy to Calculate

* $J=\frac{U_1}{nD}$, advance ratio
* $\Omega=2\pi n$, angular speed

### 3. Hard to Find: Iterative Method

* $a=\frac{U_1-U_2}{U_1}$ air speed reduction ratio.
* $b=\frac{U_\Theta}{\Omega r}-1$ fix ratio
* $\phi=\Theta + \alpha$, flow angle

They can be find using the following equation repeatedly
$$a=\frac{\sigma(1-a)}{4\sin ^2\phi}(C_L\cos \phi+C_D\sin \phi)$$
$$b=\frac{\sigma(1-a)}{4\sin ^2\phi}\times\frac{J}{2\pi r/D}(C_L \sin \phi - C_D \cos \phi)$$
$$\phi = \tan^{-1}\Big(\frac{U_1(1-a)}{\Omega r(1+b)}\Big)$$

## Airfoil Design

