Before the mass adaptability was added to the simulations, we only had sigma and energy to be the parameters in the simulation. 

By adding mass to the mix, we have three parameters that can be used to make a simulation. Let's think critically what these parameters do and how they are linked
1. Energy -> defines k_0
2. sigma_real -> defines the spread of k_0, sigma_p
3. mass of the particles -> (don't know)


Consider the original spreading equation of the free wavepacket. In the limit of the gaussian width is much smaller than (h_bar*t/(2m)), the growth in the size of the wavepacket is linear in time. This means that the inverse of the initial wavepacket width determines the growth rate of the wavepacket in time. Now, perhaps, a better way of looking at this is the characteristic length at which this expansion occurs. We can look at it this way. 

In this limit, we also observe that, for two identical particles, apart from the mass, find that the the time evolved gaussian spreading is also inversely proportional to the mass. This means that, for a simulation run that we are happy with, we can tune the mass so that the gaussian broadening is acceptable to us. This makes it feasible for us to conduct our simulations. 

However, changing the mass comes with some costs of its own. Think in terms of the momentum of a wavepacket. We find that that the momentum of a wavepacket is given by - 
$h_bar*k_0$ where k_0 is the mean. The K.E of the wavepacket is given by p^2/2m, which turns out to be $h_bar^2*k_0^2/2*m$. To understand what's happening, let's compare the 100 eV electron to the muon. 

The electron, expressed in plane waves would have a mean momentum. Consider this mean momentum k_e. Consider the mean momentum of the muon (or other particle) to be k_m. The ration of (k_m/k_e)^2 is given by the m_m, where m_m is the mass of the particle in terms of electron mass. plugging the numbers in for an electron and a muon, we see that, to describe a muon with the same energy as an electron, we would need about 14 times more momentum (and hence higher e_cutoff, as the we consider plane waves with m=1). 

To summarise, the mass knob has the following effects - 
1. Benefit: Increasing the mass reduces the spreading (almost inversely proportionally in a certain liimt)
2. Cost: Increasing the mass also increases the momentum of the projectile and hence demands a higher e_cutoff to make the particular simulation. 
3. Cost: Increasing the mass also increases the traversal time also increases with the same energy
4. 