---
id: ab-benchmarks-param-search
area: absorbing_boundary
title: "Absorbing-boundary feasibility + MFA reflectivity sweep"
status: done
hypothesis: "A sin^2 CAP / mask-function absorber in inq-study yields clean reflectivity curves eps(E,L) usable to lock production absorber parameters."
handover: docs/handovers/absorbing-boundary.md
tasks:
  - { name: "CAB feasibility study in inq-study", done: true }
  - { name: "MFA implementation + tests", done: true }
  - { name: "reflectivity eps(E,L) sweep + free-prop runs", done: true }
  - { name: "eps curves + study notebook", done: true }
blocked_reason: ""
---

# Benchmarking and Parameter searching for the Absorbing boundary conditions. 

<identity>
You are a scientific computing researches working on first principles simulations. You have a good understanding of first principles domain and are excellent at writing scientific standard code. You adhere to the rules, principles, workflows established in this repository. 
</identity>

<description>
In this task, I want to plan an overnight task. There are two tasks fo this project. Firstly, I want to conduct a feasibility analysis of using complex absorbing boundary in this project. You should use inq-study repository (and not inq/ at all; it is to be left as is, unchanged). I want you to evaluate if a certain plan I've come up with would work in this setting. What are the main challenges? How can they be overcome have to be critically thought out and reasoned.  

The second task is to implement mask function absorbers using inq-stack library. The theory regarding the same can be found in the paper ResearchProject/literature/tddft-quantum-projectile/resources/modeling-electron-dynamics-coupled to-continuum-states-in-finite-volume.pdf. In this study, we are also going to produce reflectivity curves using epsilon as a function of E (energy) and buffer width L. Now, these curves will be used to make a decision at a later point about the parameter to choose for my actual line of research.  

For both the tasks, I want you to make a ipynb file in docs/ in an appropriately named folder, showcasing all the results, and then presenting the results. 

I would want to run these tasks overnight, without my involvement so that I have the results in the morning. So, I would like to use /grill-with-docs plugin to ensure that ambiguities are removed and a clear, crystallised plan is formed.  
</description>




## Tasks

<task>
<name>Feasibility study:</name>

Here are my initial thoughts about implementing a CAB in inq-study. Choose a AB scheme - consider sin^2 potential.
Because I chose CAB, i would need to modify the inq-study source code to have this CAB in here. However, the CAB is restricted to a certain domain in real space. So, I need to figure out a way to encode this. Also another complication might be the fact that INQ works in momentum space. However, it implements the kinetic energy opreator in momentum space. The potentials are applied in real space. This has to be checked. Consider, the previous step has been completed correctly. Then, I run the same setup as the paper, at the said time tau, I integrate over all space, and find only the reflected wavepacket. This expressed as a percentage of the initial norm of the wavepacket (assume 1), will give me error epsilon. This is representative of the S.

This is the workflow I want to implement for these boundary conditions. I want you to carefully look at each step mentioned here, critically think if each step is achievable, and if so how? If not easily achievable, what are the challenges? 

</task>


<task>
<name>Implementation of mask function absorbers, benchmarking and parameter study</name>

To implement this mask function absorbers, this is my suggetsion - atleast to my mind of masking functions absorbers, I just would need to change the unitary evolution propagator of the INQ code so that at each iteration, a new psi(t+dt) = M(x)*U*psi(t), is this right? So this must be done each iteration at the propagator level. Considering that psi(t+dt) is automatically produced by the Unitary operator, then, I can perhaps use my wrapper library to change the wavefunction each timestep. So, what i would do is 

get the result of psi(t+dt) = U*psi(t)
psi'(t+dt) = M(x)*psi(t+dt) Would this be equivalent?

So this work can be done entirely in inq-stack wrapper library and runs can be made to test this. 

There are specific steps in this plan. 
1. Implement the mask function absorbers carefully. Write tests to check its working. Fix any errors that might come up.  
2. Replicate the study performed in the paper using a single masked function. Use the function that was suggestd in the text as equation 13. Specifically, we are going to produce reflectivity curves of using the simulation runs (and analytical runs). We should also have a few free propagation runs mixed in here, only a few, randomly selected to enhance the dataset we're producing. 
3. Plot the epsilon curves with different Ls at different wavepacket energy. 
4. Visualise it in the ipynb alongside all the other important information.  
</task>

<calculating_epsilon>
In this task, epsilon as a function of E (energy of the wavepacket) and L (buffer width of the absorber) is a measure of the reflection this potential produces. This is going to be calculated by comparing the simulation runs at the time t = tau. The selection of time tau is described in the paper. In simple words, it is the time that the freely propagating wavepacket (analytical) would be beyond the box. In the simulation, as mentioned in the paper, the reflected wavepacket reaches its initialised position of -3sigma after completing a return journey. 
</calculating_epsilon>

<paper_method>
I have tried extracting the main information about the paper's method in this prompt. You should validate this against the actual paper, and work in a loop to ensure you have loaded the right parameters for everything. 

The Geometry Is Correct With One Small Refinement
The total box spans [−X,L][-X, L]
[−X,L] where X=6σX = 6\sigma
X=6σ and LL
L is the absorber width. So the total length is 6σ+L6\sigma + L
6σ+L, not X+LX + L
X+L as a separate thing — XX
X *is* 6σ6\sigma
6σ. The wavepacket starts at x0=−3σx_0 = -3\sigma
x0​=−3σ, which places it exactly halfway between the left wall at −6σ-6\sigma
−6σ and the start of the absorber at x=0x = 0
x=0. This gives it room to propagate freely before hitting the absorber, and room behind it so the reflected wave can develop without immediately hitting the left wall.

The Propagation Time and Why
The total propagation time is:
τ=2(3σ+L)k0\tau = \frac{2(3\sigma + L)}{k_0}τ=k0​2(3σ+L)​
The numerator 2(3σ+L)2(3\sigma + L)
2(3σ+L) is the total distance the wavepacket must travel, broken into two legs:

Forward leg: from x0=−3σx_0 = -3\sigma
x0​=−3σ to the far end of the absorber at x=Lx = L
x=L, a distance of 3σ+L3\sigma + L
3σ+L
Return leg: the reflected wave travels back from x=Lx = L
x=L to x=−3σx = -3\sigma
x=−3σ, another 3σ+L3\sigma + L
3σ+L

The denominator k0k_0
k0​ is the group velocity. So τ\tau
τ is precisely the time for the wavepacket to make one full round trip. The simulation is stopped *before* the reflected wave reaches the left wall at −6σ-6\sigma
−6σ and bounces again — that second bounce would contaminate the measurement.

The σ\sigma
σ Choice and Its Physical Meaning
With σ=42/k0\sigma = 4\sqrt{2}/k_0
σ=42​/k0​, the uncertainty principle gives a momentum-space width (standard deviation):
σk=1σ2=k08\sigma_k = \frac{1}{\sigma\sqrt{2}} = \frac{k_0}{8}σk​=σ2​1​=8k0​​
The spreading velocity of the wavepacket (how fast the Gaussian envelope broadens in real space) is:
vσ=σk=k08v_\sigma = \sigma_k = \frac{k_0}{8}vσ​=σk​=8k0​​
Since the translational velocity is k0k_0
k0​, you have vσ=k0/8≪k0v_\sigma = k_0/8 \ll k_0
vσ​=k0​/8≪k0​. The packet moves much faster than it spreads — it stays coherent and well-localized throughout the propagation. This is why the wavepacket has a well-defined energy even near E≈0E \approx 0
E≈0, and why you can trust ϵ(k0)\epsilon(k_0)
ϵ(k0​) as a good approximation to the survival probability at momentum k0k_0
k0​.

Grid Spacing to Energy Cutoff
The grid spacing Δx\Delta x
Δx sets a maximum representable momentum through the Nyquist criterion. The smallest wavelength the grid can represent is λmin=2Δx\lambda_{min} = 2\Delta x
λmin​=2Δx, giving:
kmax=πΔxk_{max} = \frac{\pi}{\Delta x}kmax​=Δxπ​
The kinetic energy cutoff is then Ecut=kmax2/2E_{cut} = k_{max}^2/2
Ecut​=kmax2​/2 in atomic units:
For Δx=0.1\Delta x = 0.1
Δx=0.1 a.u. (standard case):
kmax=π0.1=31.4 a.u.k_{max} = \frac{\pi}{0.1} = 31.4 \text{ a.u.}kmax​=0.1π​=31.4 a.u.
Ecut=(31.4)22=494 Hartree=494×27.21≈13,440 eV≈13.4 keVE_{cut} = \frac{(31.4)^2}{2} = 494 \text{ Hartree} = 494 \times 27.21 \approx \mathbf{13{,}440 \text{ eV} \approx 13.4 \text{ keV}}Ecut​=2(31.4)2​=494 Hartree=494×27.21≈13,440 eV≈13.4 keV
For Δx=0.05\Delta x = 0.05
Δx=0.05 a.u. (high-energy case, E>2E > 2
E>2 keV):
kmax=π0.05=62.8 a.u.k_{max} = \frac{\pi}{0.05} = 62.8 \text{ a.u.}kmax​=0.05π​=62.8 a.u.
Ecut=(62.8)22=1974 Hartree≈53,700 eV≈53.7 keVE_{cut} = \frac{(62.8)^2}{2} = 1974 \text{ Hartree} \approx \mathbf{53{,}700 \text{ eV} \approx 53.7 \text{ keV}}Ecut​=2(62.8)2​=1974 Hartree≈53,700 eV≈53.7 keV
The reason the paper switches to Δx=0.05\Delta x = 0.05
Δx=0.05 a.u. above 2 keV is exactly this: at 2 keV the de Broglie wavelength is about 0.52 a.u., which means only about 5 grid points per wavelength at Δx=0.1\Delta x = 0.1
Δx=0.1 a.u. — dangerously coarse. Halving the spacing doubles the grid points per wavelength and restores accuracy.

The E=5k02/4E = 5k_0^2/4
E=5k02​/4 Formula
The standard non-relativistic kinetic energy in atomic units is:
E=k022 Hartree=13.6×k02 eVE = \frac{k_0^2}{2} \text{ Hartree} = 13.6 \times k_0^2 \text{ eV}E=2k02​​ Hartree=13.6×k02​ eV
The paper's formula E=5k02/4E = 5k_0^2/4
E=5k02​/4 (in eV) differs from this standard expression. The 5/4 prefactor is not derivable from simple unit conversion or from the standard kinetic energy of the Gaussian, and the paper does not explicitly explain the origin. The most likely interpretation is that it encodes a specific property of how the peak energy of the wavepacket maps onto the reflection curve's dominant feature, possibly folding in the momentum spread. What you can be confident about for replication is that:

The energy axis scales as E∝k02E \propto k_0^2
E∝k02​, as expected for a free non-relativistic particle
The formula maps one-to-one with k0k_0
k0​, so it is just a relabeling of the horizontal axis
For your own implementation, using the standard formula E=k02/2E = k_0^2/2
E=k02​/2 Hartree =13.6k02= 13.6k_0^2
=13.6k02​ eV will give you the same physics — you would just rescale the horizontal axis accordingly


How the Ideal Free Propagation Is Calculated
The free-propagation reference ψ0(t,x)\psi_0(t,x)
ψ0​(t,x) — the "what would have happened without the absorber" — can be computed one of two ways:
Analytically: The free propagation of a Gaussian wavepacket has a closed-form solution:
ψ0(t,x)=1(πσ(t)2)1/4exp⁡(−(x−x0−k0t)22σ(t)2+ik0x−ik02t2+iϕ(t))\psi_0(t,x) = \frac{1}{(\pi\sigma(t)^2)^{1/4}} \exp\left(-\frac{(x - x_0 - k_0 t)^2}{2\sigma(t)^2} + ik_0 x - i\frac{k_0^2 t}{2} + i\phi(t)\right)ψ0​(t,x)=(πσ(t)2)1/41​exp(−2σ(t)2(x−x0​−k0​t)2​+ik0​x−i2k02​t​+iϕ(t))
where σ(t)=σ1+t2/σ4\sigma(t) = \sigma\sqrt{1 + t^2/\sigma^4}
σ(t)=σ1+t2/σ4​ is the spreading width and ϕ(t)\phi(t)
ϕ(t) is a known phase. This requires no numerical simulation at all.
Numerically: Run the same propagation again, with exactly the same initial state and time-stepping scheme, but with VCAP=0V_{CAP} = 0
VCAP​=0. This is probably what the Octopus/QuantumPy codes do in practice, since it ensures that any numerical artifacts from the finite difference scheme and time-stepping cancel identically when you compute ψ−ψ0\psi - \psi_0
ψ−ψ0​. You are measuring only the effect of the CAP, not the discretization error.
The reflection error at time τ\tau
τ is then just the integrated squared difference between the two:
ϵ(k0)=∫−XLdx ∣ψ(τ,x)−ψ0(τ,x)∣2\epsilon(k_0) = \int_{-X}^{L} dx \, |\psi(\tau, x) - \psi_0(\tau, x)|^2ϵ(k0​)=∫−XL​dx∣ψ(τ,x)−ψ0​(τ,x)∣2
Since ψ0(τ,x)≈0\psi_0(\tau, x) \approx 0
ψ0​(τ,x)≈0 in [−X,0][-X, 0]
[−X,0] at late times (the free wavepacket has moved past x=0x=0
x=0), this simplifies to measuring whatever is left in the inner region — which is purely the reflected wave.

The Main Experiment Summarized for Replication
The logic of the experiment is:

Fix LL
L — the absorber width
Sweep k0k_0
k0​ across many values (giving different energies EE
E)
For each k0k_0
k0​: set σ=42/k0\sigma = 4\sqrt{2}/k_0
σ=42​/k0​, place the wavepacket at x0=−3σx_0 = -3\sigma
x0​=−3σ, propagate for time τ=2(3σ+L)/k0\tau = 2(3\sigma + L)/k_0
τ=2(3σ+L)/k0​, compute ϵ\epsilon
ϵ
Repeat for many values of LL
L
Plot ϵ(E,L)\epsilon(E, L)
ϵ(E,L) as the reflection curve — this is Figure 3 through 10 in the paper
</paper_method>


<rules>
1. Never change or modify inq/ library. Any analysis that is to be done with the knowledge of the inq library, can be done using the inq-study/ library which is an exact replica. 
2. 
</rules>