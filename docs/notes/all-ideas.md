In this notebook, I am going to catalogue all of the ideas and tests I want to perform. This is a quick review before I start writing my report. As I plan my report and start planning it, I am going to look at all the possible areas that the research can end up in. Look for ideas I've forgotten to implement.

## Places to look at
- Journal
- Dissertation journal: 
    - Meeting Notes
    - Investigations
    - Explorations
- Talks

Where else could I have my ideas?


## Ideas

What's the best way to categorise these ideas? Firstly, there might be a lot of different ideas. Some I would want to test, some I can discard, some might be future exploration if required. Also, I might have questions that I have written down for myself. So, I can turn questions, or ideas into potential experiments that I can conduct. 


## Questions
1. Why do a lot of the WP runs in localised jellium have an oscillation in time of energy? Perhaps, this is to do with the fact that the hamiltonian has become time dependent?
   - Do higher density of slabs give cleaner energy absorbed profile? Need to check
2. Localised jellium system:
    - Why is there some leaking density out of the prescribed region and more consequentially, why is there an accumulation at the edges? 
3. Bulk Jellium: 
    - Perhaps, the idea was that, if the definition of stopping power in terms of its component is determined, then we can use bulk jellium to define stopping power. What we learn from there can be used in this system.  
4. CAPs
    - Localised Jellium: Does having a CAP, and the fact that the hamiltonian is not Hermitian, make the simulaiton un usable? Is this the reason behind the energy oscillation in the simulations?
    - 
5. Stopping Power Definitions
    - There are three definitions we are working on for stopping power at the moment. These are - 
        - Using energy decomposition, we come up with a formula for stopping power in a medium
        - Stopping power definition using localised jellium simulation
        - 
6. Classical Projectile
    - Verifying Nazarov Gross using a classical projectile
## Potential Tasks
1. Comparison with Nazarov Gross: Build runs with different mass of the projectile to examine if, even qualitatively, this phenomena is observed. 
2. Can we build a force field for a given projectile in the simulation box, then use it to propagate the projectile in ehrenfest way? (Can use linear interpolation where possible.)



## Tasks
1. Localised Jellium Quantum Stopping Power: Carefully thinking about the previous results S(v) using the localised jellium setup, we found that the stopping power was too much. One reason for this might be that, in the initial setup, E_sp is not negligle (as we wanted this to be). However, at t = t_final, E_sp is negligivle as the excess projectile is absorbed. Now, I want to make a rough calculation in these runs, if the subtracting this energy of E_sp from the run, and then measuring E_absorbed, and updating the stopping power calculation, if this would give a result that is close to the stopping power using lindhard bulk?
- The other thing to try would be to remove the localisation energy and the absorbed electron energy .