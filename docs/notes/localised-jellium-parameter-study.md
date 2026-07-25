In this file, I am going to record my thoughts as I review the jellium parameter study campaign.

## Classical vs. WP initial energy at different positions of the projectile from the slab. 

In this task, the main question I was asking was - "What is the difference in the total energy of the system between the WP GS and the classical projectile GS when the projectile is initialised at different positions in the simulation cell?"

To answer this question, I had a stationary WP and classical projectile at different distances r from the slab. The total enegy of the system was measured for each simlation and it was plotted. 

- E_total of only jellium slab: -108.53 Ha. 

For both the classical and WP cases, the simulation was run for two steps and the total energy was calculated by considering the energy at the timestep 0. The difference in energy E_total(0)-E_GS was calculated both the runs. This was identified as the energy of inserting the projectile. 

Questions: 
1. I understand that in the WP insertion case, the Hartree term, the kinetic term and LDA self interaction energy term increase. However, in the classical case, I understand that the V_ext term increases (changes). How is this accounted for in INQ? Is the external potential also recorded? 

Suggested: Apparently, INQ stores E_external in a variable that is never really streamed outside to the data. A small change in inqkit would have be fine to get all these components. With this, we also need need p2 periodicity (for the localised jellium slab) but not periodic in the z direction. The new results will be more meaningful in this respect. 

2. Also, in this H0, is the periodicity 3? Answer: Yes


## Energy decomposition of the total energy as a function of distance r from the slab

I am now trying to interpret the results. Check my interpretation and answer my questions.

**Energy components to document**

So, you should now be documenting all the components of energy. What are all the components of energy — Hartree, external potential, kinetic energy and exchange correlation. Perhaps, while using a pseudopotential you might get a local and non-local term for the external potential. I guess there is also the ionic energy included here.

I want you to arrange the bars in a logical order so that I can understand what's happening easily. However, the others are the same right?

Now, I also want to look at an individual run, with a bar chart displaying the total energy decomposed into all the segments that are documented. The total segments should add up to the total energy of the system.

**Interpretation — periodicity 2 re-run with full measured energy decomposition**

*Classical case:* I understand that in the classical case, as the distance increases between the projectile and the slab, the total energy decreases. This makes sense as the external potential energy term must be decreasing. Now, in the classical plot, you have not shown decomposition of energies. Why so? It would be beneficial to compare these here.

*Wavepacket case:* Secondly, looking at the wavepacket run, I observe that the total energy remains the same. However, the external potential and the U_H terms increase and decrease respectively. Let's understand this better. The U_H term decreases as the distance increases and comes to zero at a distance close to 40 bohr. So, this is the electrostatic repulsion between the jellium slab and the wavepacket. Now, this should be positive energy right, as the wavepacket and the jellium are both negative and hence the potential energy is positive. Now, the external potential energy is the coulombic interaction between the background positive charge and the wavepacket. Hence, this energy is negative and hence tends to 0.

Given that the distribution of the jellium system and the positive potential is almost the same except some differences, the plot argues that the difference between the positive and negative charge contribution to the system almost cancels out and is negligible.

Coming to the classical plot, we clearly see that the total energy in the system falls as r increases. The total energy here is positive. However, we have to remember that E_ext is probably overlapping with the E_total energy (confirm this).

But, both the graphs suggest that the self interaction energy is very small and can be effectively ignored in the simulation. Secondly, it argues that if we are to start from a system where the initial coulombic interaction is 0, then about 40 bohr is the distance we should be looking to go to. I also am skeptical that the total energy is 0 here. I think we need a few more runs with longer r than 40 to conclude this.

**Questions:**
1. Does the positive charge and negative electronic systems has identical distributions? Can we plot the anisotropy between the two and visualise them in xz, yx and xy plots for the ground state? This way, I can understand if they are the same or not.
2. Is this interpretation right?

**Requested runs / next steps:**
Make new runs with the jellium slab further down the z axis so as to give more distance between the projectile and the localised jellium slab. Do this for both classical and wavepacket runs, and update the plot.


### Single run energy decomposition
I am now looking ath teh H0_p2_interpretation ipynb. Here, the energy of a single run is deconstructred. However, I have a few questions here. I feel like the sign of the numbers is flipped. So, I am a bit skeptical if this is right. 

Questions: 
1. in the wavepacket case, the wavepacket is a negative charge and so is the electronic localised jellium system. So, why is U_H negative?

ANSWER: Apparently, in the periodicity p3 (periodic BC), we find that the G=0 term is dropped completely. But, in the open z p2 (periodic only in x and y axes), we fidn that the G=0 term contributes 0.5*r_c^2, wher r_c is approximately the box size. This actually blows up (is a constant though), and alters the signs of all the components. However, this term apparently has the opposite sign. How does this explain the E_ext being positive though? 

2.  Similarly, why is E_ext positive when it represents teh attractive energy?
3. In this case the kinetic enrgy about 160 eV. Is this the total kinetic energy of the electronic system? That is, including the jellium slab system electrons?


## Theoretical simulation of the jellium slab and expectations

In the theoretical modelling of this system, we considered the electronic system as a perfect slab of jellium that is periodic in the x and y directions and have a thickness in the z directions. The effect of the slab was found by modelling it as a continuum of infinite plates. Here are some important results

Infinite plate
1. The total potential and hence the potential energy due to infinite plate is proportional to the distance r from the plate. Meaning, as r increases, the distance to the slab from the projectile, the potential energy also increases. This can be rationalised by thining about the electric field generated. The electric field of an infinite plate has a constant magnitude and direction at any position in space. Hence, in a uniform field, we have this linear relationship between the potential energy and the distance. This is the case for capacitors (if I remember this correctly). 

Slab (Localised Jellium)
The slab is a cumulation of infinite plates. In this case, it can be shown that the total coulombic energy is the same as if placing an infinite plate (with all the charge), right at the center of the slab. The potential enrgy would be the same. Hence, the infinite plate should be a good model to approximate the coulombic interations between - 
a. Projectile and the electronic system
b. Projectile and the uniform localised constant positive background

Now, it is apparent that we don't see that the total energy of the system increasing with distance r. This is because, we have two slabs, one positive and one negative overlapped on one another. Hence, the simplistic modelling suggests that, the physical effects of coulombic repulsion are mainly due to anisotropies in the system. 

Now, consider the potential energy inside the localised jellium system. By Gauss' law, the electric field right at the center is 0. As we move towards one of the ends, the field one one side decreases, while on the other side increases. One way to look at this is that, for a section of the slab on one side of the point we're considering, by symmetry, there's another the other side of the jellium which applied a force in the opposite direction. Hence, the only force that counts is from the part of the jellium slab that is away from the cancelled part. This makes it so that the total field start increasing as we move out. So, the electric field, contrary to the outside the slab is linear. Now, the work done is field integrated with distance (times charge), which would be parabolic.  

This also suggests that the component of hartree having WP-slab component must could tend to 0 if there are random anisotropies.  

Also, this explains the tending to 0 of the U_coulombic in the wavepacket case. However, this does not explain the decay in the case of a classical projectile. What really prompts this? I think this is an important question to answer. Could the cutoff radius for the gaussian coulombic potential a reason? could this be it?

In fact, it can be argued that the difference between the E_ext and the U_Hartree's 

Questions
1. Wait, how was it that in H0_base_difference ipynb, we found that the E_hartree tended to 0? This should be impossible right? This is because, the jellium slab itself must have repulsion amongst its electrons (which I would expect to be far larger in magnitude to the WP-slab interaction). Why is this not found?

2. The simple model of stacked infinte plates explains the U_hartee = 0 in the WP case. But, in the case of the classical projectile, this does not explain the decay. I wonder if the cutoff length of the gaussian potential of the classical projectile might be an issue. Can we test it somehow? 

- Further deliberation
I want another plot in the ipynb where only yhe net real density affect is shown. This must show
  that at any distance r, the coulombic attraction is almost 0. This would validate the WP run. But,
  we need to think carefully about what's happening with the classical projectile run. You have 
  potentially found a link that might help understand why the eneryg is 0 beyond r = 50 bohr or so.
  I now, want to make another test, where you make a bunch of these pseudopotentials, each with its
  own cutoff radius ranging from 10 bohr to the 40 bohr (let's say 4 such pseudopotentials). For
  each, make the delta E_total graph. For each pseudopotential, annotate where the cutoff is. This
  would help me understand if the key difference is here. However, that does not answer the central
  question I have of whether this has an appreciable impact on the simulation and if this affects
  the simulation at all. Now, the analogue of E_external for the wave packet run would be placed in
  the E_ion-ion for the classical projectile case. This, I would expect to have an opposite trend to
  the U_external. However, we've actually ploted E_total - E_GS here if I am right. This raises a
  question, is the E_ion-ion added to the total energy calculation. If so, get it in the above runs
  to confirm it



Reason why the E_hartree of the classical and the wavepacket runs cannot be directly compared
It is suggested that, the total energy in fourier space has a component G=0 which is chosen as a reference. For periodicity 3, it is 0. However, for periodicity 2, it is $0.5*rc^2$. rc here is the length of the box, in the case of the long z boxes being used, it is 120 bohr. Hence, the V(G=0) term is $0.5*120^2$, which is a huge number. This coupled with the fact that there are 82 electrons in the classical run, and the WP has 83 electrons, direct comparison would fail as the constant shift would be the $(0.5*120^2)*(charge)$. However, there is a huge constant shift between the hartree and the external potential terms, but in the opposite directions. Hence, the sum of these two ensures that the constant term cancels off. Hence, the only physical signal is the E_hartree + E_external and the difference between the clasical and the wavepacket runs. 


Potential Learning from these results: 
1. If the radial cutoff of the classical projectile's radial potential affects the contributing energy, that tells us that the classical simulations either have or don't have longer range effects due to the cutoff. Now, the wavepacket, although with the same gaussian width, might not be causing long range effects at all. There are few ways this might play out in my mind. Firsly, there might be some actual physics going on here, such as some screening effect, or correlation effect that might be lowering the coulombic impact of the wavepacket's charge on the jellium system. To look into if this is the case, I would carefully consider a classical and a wavepacket run that are comparable, and with all the components of decomposed energy in it. then, I would carefully look for difference in the E_xc term. Then, perhaps, it is wise to conduct a holistic study of the total energy difference between the classical and the wavepacket cases. Although the energy exists in different stores between these two runs, its onlt the U_ion-ion component that I belive is missing, and hence classical projecile simulations must have a lower total energy. Now, the difference can be approximated to a good degree by the analytical model we've developed of the infinite slabs put together jellium slab potential energy. I want to check, if the difference in energy agrees to the prediction by this equation. Up until now, we were scouting for physics, and some sanity checks. However, the second possibility is that there is a computationa
2. The above thoughts prompt a study of the coulombic potential produced due to the wavepacket and how it compares to the classical (defined to be a gaussian radial potential). Is there something we can do to answer these questions?


### Theoretical model using dipole and other fancy stuff

The aim here is to explain the linear decrease in the U_hartree term as the distance is increased that is observed in the plot. It is important to do this as we are deriving physically interpretable outcomes from here. 

Using what we've learnt from the crude model, 

### Examining dynamics carefully for where the differences between the runs exist
1. To understand the differences between the classical and wavepacket runs: Once I undersand the energy decomposition of the classical and wavepacket runs completely, and understand what fundamental differences exist between the two (in GS), if any exist. This allows me to go timestep by timestep, and then, using my understand of where the energy is stored in each simulation, I can try to extract specific differences between the runs. This way, I find something that is different between the two runs. This difference can be claimed (after some verification) to be the "quantum effects" that classical projectiles don't display. 
2. 

While thinking about the analysis I've done here, today, it occured to me that, the more physical thing, that causes the coulombic dyamics is the net charge in the simulation. It is only the net charge that affects the coulombic dynamics of the simulation. 