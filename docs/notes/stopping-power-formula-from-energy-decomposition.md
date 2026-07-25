In this notes, I write down my thoughts about formulae that I can come up with, that can physically represent the energy absorbed by the slab or the energy lost by the projectile. Using this, we can find the stopping power of the system-projectile combination. We are going to follow the campaign named "formula-sp-from-energy-decomp". 


## Formula suggestions

### C1

This definition is based on the energy gained by the slab. This should be applicable for both the classical and wavepacket runs. 

S = delta E_target/ delta s
Now, here delta E_target(t) is defined as the difference between the E_electronic(t) - E_electronic(0). 

E_electronic here is defined as sum of kinetic energy of the slab, sum of the coulombic energy of slab-slab, the sum of coulombic energy of slab-background and finally the exchange correlation. 

E_electronic = T_slab + E_ss + E_sb + E_xc

Notice, there is projectile interaction energy at all here. 

Let's consider the application of this formula in classical case. IN classical case, the T_slab is easy to find (the default T of the orbitals). E_ss term is essentially the hartree term in the simuation and again easy to calculate. Now, E_sb has to be calculated individually as E_ext has both E_sb and E_sp in it. So, we need to calculate it individually in the simulation (as we already are) and use it instead. Finally, E_xc corrects for any exchange and correlation effects. So, in classical case, the formula can be applied rather easily. 

Now, let's think about its accuract here. When the projectile enters the box, and is traversing through it, the slab is energised. The total energy it has is precisely as described above. Is there any other term we are missing here? What about the E_sp? The slab projectile interaction energy, why are we not including it in the energy of the slab? Is there a reason? I can understand if the ideas presented here is that at t=0, E_sp is negligible and at t=t_final, it is negligible as well. Under this case, the formula we've mentioned above is good. However, in the simulations we have run, is this the case? I don't think so. Usually, there is a substantial E_sp initially and finally. How should we treat this.  

Now, let's think about the wavepacket run. IN the wavepacket run, T_slab definition might be vauge. This can only be obtained when we assume the orbitals of the slab really represent the slab under excitation. Meaning, when the projectile is actively interacting with the slab system, we might have it that, the KS orbitals might not represent WP and slab systems separately. So, this might fall here. Similar, E_ss, E_sb terms face the same problem. However, if we work under the assumption that KS orbitals physically repsent the slab and wavepacket clearly, then, under this big assumption, we can compute the same energy target. Let's do this, but I need the caveat. 


### C2
THis is essentially the classical definition of stopping power using a pseudopotential (or a classical projectile). I have no objections to this apart from the fact that, I need to see the steps in calculating this value clearly in the notebook. Meaning, I need to see delta E(t) plot, then delta E_kinetic_projectile(t) plot, then, mention clearly which region was considered, then show what best fit line was used, then give the value. 


It would be intersting to look at C1 and C2 definitions for a classical run. Ideally, these two should agree with each other. 

### C3
IN this definition, we consider the energy lost by the projectile. 

We define, in the wavpacket run, 
E_proj_WP is T_WP + E_pp + E_ps + E_pb. Now, this again comes with the caveat that we assume that the projectile is well represented by the KS orbital when it is interacting with the system. Now, this makes this a weak definition. But, I would want you to calculate the stopping power using this. 

Now, as you've pointed, the zero point energy is defined by the localistion energy. So, we can essentially subtract this at the start to ensure that the others represent the true energy of the wavepacket (with the same kinetic energy due to momentum). 

### A1
This definition has already been used. We need to refine this. Up until now, I have seen that the calculated energies using this method are too huge. Meaning, when compared to lindhard bulk jellium prediction, these enegies are much higher. this can be found in one of the campaigns where we had a few well behaved energy plateaued runs for wavepacket projectile in localised jellium. using the plateau, and E_absorpbed, I found s(v). With the knowledge you have, you can find these runs.  

Now a few things are possible- 
a. bulk jellium does not represent the localised jellium system with the same density. hence lindhard s(v) are wrong. so, we need to make comparisons with classical projectiles. 
b. there is some mistake in the calculations of calculating s(v). Perhaps, we needed to use sum of decomposition enregy terms instead of e_total. Meaning, somehow, e_total(t_final) - e_total(0) somehow had other energy terms which we should not be tracking
c. CAP causes some anamoly in the total energy curve hence altering the total energy of the system, and hence producing the wrong stopping power. 


So, we need to understand why this problem arises. I need claude to examine these carefully, use fable 5 to think about these three cases, run tests and experiments to verify them, and the give a conclusion. 

### Conclusion: Take 1 
Overall, there are some suggestions, which are not entirely novel, but worth implementing. However, I would also like for you to think if there are any other setup modifications, which allow us to simplify the energy structure so that we can effectively calculate the stopping power. (What I mean by setup modifications is, for example, the reason i had a long simulation for localised jellium run to run is that E_total(t_final) would have the projectile and any overflow from initial slab to be absorbed. Hence whatever remains can be thought of as jellium. Hence, we are able to say that E_total(t_final) is E_slab(t_final). So you can come up with such like suggestions. Doesn't have to be the exact same, but some setups help get a simple structure out).


