In this notes, I will document all of my thoughts and understanding I have about the enregy decomposition of the classical and wavepacket runs. I would want to turn this into a skill. This skill is going to have all the understanding that is required for an agent to use it, and run systematic analysis between classical and wavepacket runs. Now, after running experiments, I have concluded that the difference between the perturbation potential (gaussian potential) representing a classical projectile and the wavepacket; more of the eneryg difference is understood physically. Hence, we are going to use this method for the classical projecitle. 

### Energy decomposition in the wavepacket run
E_total = E_kinetic + E_hartree + E_xc + E_external
E_kinetic - the summation of the kinetic energies of all the orbitals
E_hartree - the electrostatic repulsion between the electrons (including the wavepacket)
E_external - the electrostatic attaction between the electornic system and the background potential. This is given by INQ. 
E_xc - exchange correlation energy of the entire system. Important to note that, as this interaction energy is local, hence, there is no Wavepacket - jellium slab exchange-correlation effect happening (perhaps correlation might happen. need to verify this)


### Energy decomposition in the classical run
E_total = E_kinetic + E_hartree + E_xc + E_external + E_projectile_background + E_projectile_KE
E_kinetic - the summation of the kinetic energies of all the orbitals
E_hartree - the electrostatic repulsion between the electrons. However, in this case, there is no wavepacket. 
E_external - the electrostatic energy between the electornic system and the background potential and the pseudpotential of the classical projectile. This is where the electrostatic repulsion of the classical particle with the electornic system is stored. 
E_xc - exchange correlation energy of the entire system. 
E_projectile_KE - kinetic enregy of the projectile $(1/2)*m*v^2$ due to the velocity of the projectile v
E_projectile_background - this is the coulombic energy between the classical projectile and the positive background. This not tracked internally by INQ. However, we have built a inq-kit method to keep a track of this. 

### 
Through rigorous experimentations, we have found the following differences between a wavpacket and a classical run, of the exact same configuration (just a classical pseudopotential is used instead of a wavepacket. )

1. E_kinetic - In this kinetic enregy term, we find that the the wavepacket run has a higher energy. This energy is understood as the <p>^2/(2*m) + wavepacket localisation energy. (Add the formula for the wavepacket localisation energy here.) So, essentially, the added kinetic energy is only due to that of the wavepacket; its localisation and its energy due to its momentum. 
2. E_hartree - The key difference between the wavepacket and the classical case is that the coulombic interaction between the wavepacket and the electronic system, and the wavepacket and the wavepacket is present in this term for the wavepacket run. However, the classical case does not have this interaction here. 
3. E_xc - The difference in E_xc between the wavepacket and the classical runs is r independent. Here, r is the distance of the projectile from the system. This remains same, and can be interpretted as the added exchange correlation due to the wavepacket in this case. 
4. E_external - The difference here is that, in the wavepacket case, the total electronic system's coulombic interaction with the positive background potential is present here. However, for the classical case, this term encodes the coulombic repulsion energy between the electronic system and the radial potential of the pseudopotential. So, this is a negative charge (projectile) and negative charge electronic system interaction. 
5. E_projectile_background - This term only exists in the classical case, and encodes the repulsion energy between the positive background potential and the classical projectile. 


Now, we have identified that d(E_Hartree + E_external) - E_projectile_background is an important observable. Here, "d" represents the difference between the wavepacket and the classical cases. So, d(E_Hartree + E_external) = (E_Hartree + E_external)_WP - (E_Hartree + E_external)_classical. In essence, this experession captures the coulombic energy between the projectile and the positive background and the electronic system. The difference between these, would be zero, if we were to say that the wavepacket exactly approximates the classical projectile. Now, there is going to be some difference that we measure here. 


### Presenting results from the experiments

In the simple experiment runs, we found the difference d(E_Hartree + E_external) - E_projectile_background to be about 7.4 eV. I still need to interpret this difference. (there is an argument that a +1 charged cell, which is the wavepacket cell adn the 0 charged cell, which is the classical cell have a tiny difference in energies. It is called the gauge. I am not entirely sure this is right, but the current estimation is that it is on the order of 1 eV. This has to be verified. Even this leaves behind 5.4 eV that is not accounted for.) Perhaps, one other effect that is missing in the self interaction of the wavepacket here. This has to be analytically calculated and then, the remaining energy must be understood. 

We found the E_xc_WP - E_xc_GS  = -16.47 eV. This can be interpretted as the wavepacket's xc, or the xc change due to the waveapcket alone. This is because, the wavepacket is at a distance from the electronic system, and becase E_xc is local, this energy never disappears and can be attributed to the WP alone. 

The difference in the E_kinetic of the wavepacket and the classical runs has been foudn to be 81.7 eV. As the momentum of the projectile in both the cases was taken to be 0, we can say that 81.7 eV is localisation energy of the wavepacket alone. 

Perhaps, what is interprettable in this analysis is the (E_xc_WP - E_xc_GS) + d(E_Hartree + E_external) - E_projectile_background. Conceptually, if the wavepacket were identical to the classical projectile, then this would be 0. Now, by setting the gaussian widths of the wavepacket and the classical projectiles so that they are as similar as possible, I would expect this to be as close to a constant (E_H[WP-WP]) as possible. So, in a sense, we can think of
d(E_Hartree + E_external) - E_projectile_background 
as being a constant for any simulation setup with the same configuration and is r independent. Unless I am missing something. (I want you to check all the arguments carefully, make this computation step by step, and then we are going to analyse the difference, and hypothesise if there is a physical reason why this difference exists.)

Now, let's combine the analytical expectation of what the E_hartree[WP-WP] is, which is essentially wavepacket self interaction is 21.7 eV. This self hartree term has been confirmed by using the actual wavepacket gaussian. It should also be mentioned that there is negligible distortion in the jellium slab density. So, it contributes no change in hartree term. This difference is higher than the measured difference of 7.4 eV. So the difference of 21.7-7.4 = 14.3 eV.


The E_xc_WP is -16.47 eV. the analytical self interaction energy is 21.5 eV. Hence, the true self interaction eneryg is 5 eV. Meaning, this is the allowed difference betweeh the classical and wavepacket runs. 

Correlation is also local. So there, is no long range correlation value. 


Now, even after extensive experimentation, I couldn't come up with a physical explanation for 14.3 eV which is the missing d(E_Hartree + E_external) - E_projectile_background energy. However, serendepitously, it seems like, replacing the pseudopotential with a pertubation of a classical pseudopotential explains the missing energy (almost exactly). This is explained and the key results are verified here - docs/notes/gaussian-pertubation-for-classical-simulation.md


## Gaussian Perturbation
In the case of the gaussian pertubation, I am able to explain more of the energy. For specific deconstruction of what energies are calculated, what differences are physical, and how to interpret them, and the specific values from the experiment are given in the file gaussian-pertubation-for-classical-simulation.md. 



## Building and Testing the skill
We are going to build and test the skill that is able to log the energy differences between the classical and the wavepacket cases. As explained above, it should keep track of all the energies for both the runs. At each timestep, it must compute all of the energies. Then, going from timestep t=0, the skill must look at what energies are changing and provide a physical interpretation for the same. For example, 

- if the residual d(E_Hartree + E_external) - E_projectile_background increases, that means that the wavepacket run's coulombic interaction of the wavepacket and the electronic system + the background is greater than that of the classical projectile. Now, this idea must be combined with information about centroid of the wavepacket and the position of the classical projecitle to suggest, if there is some velocity changing or etc driving this. Or, perhaps, the gaussian broadening of the wavepacket. 

Now, we should expect there to be some butterfly effect. I mean, due to early induced changes, there are going to be further downstream changes. After a certain time, these changes become very hard to interpret. So, we are going to try to explicitly interpret the changes in the first few timesteps, and then suggest a general trend for the rest. We are going to claim these are the "quantum" effects of the simulation. 

To test these, we are going to make pair simulation, one classical and one wavepacket with exact configuration (apart from the description of the projectile). then, we are going to simulate this for some timesteps, and then try to examine the differences. I want you to start by analysing one simulation. Make a run notebook for both the simulations, and analysis notebook for the task give here. I am going to suggest corrections if there are any. 

Then, we are going to make a few more runs that I want to analyse. For example, a wider gaussian sigma. We could tune the parameters so that we can change the gaussian broadening (so that there is negligible) in the wavepacket run. Then, we get one type of results.



## Note: 
We do not have to change any internal workings of inq. We can add modules in inq-kit though. The aim has to be that, we keep track of all the energies for both classical and wavepacket runs at all the timesteps. The simulations need to be outputting all of these energies so that they can be analysed in conjucntion in post processing. 