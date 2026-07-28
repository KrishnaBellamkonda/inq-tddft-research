In this plan, I am going to register all of my thoughts for the second report that I am making. Then, I would also have to change the first report so that everything flows as expected. This notes file can be thought of as me thinking aloud. So, things I say here, only when are affirmed with "LOCKED", are fixed decisions. Others are thoughts I have/. 

## Previous Report Structure

Below are the sections I've had in my previous report. In here, 
1. Introduction
2. Literature Review
3. Methods
4. Initial results and discussion
5. Summary and conclusions
6. Limitations and Outlook

I think, I need to work on the feedback. Secondly, I need to rewrite section 1, atleast make one pass on it, to ensure that the it can be read together with report 2. So, I need a pass here.  


## Report 2 Requriements
### Project Requirements
Here are are the tasks I need to in order to get marked accordingly. 
1. Refinement of project direction and subsequent development (20%)
   I believe this is to do with, how I've changed the direction of my project, how I've handled it, and how I have attempted to answer the core questions that were asked of me. I need to justify the direction I've taken, and the decisions that were considered. The justification must be with data.  
2. Resutls, analysis and discussion (40%):
   Here  is where we dwell directly into the crux of what was done. Present quantitative results, with uncertainty, and a critical discussion vs. literature. 
3. Conclusion and Future directions (20%):
   Here is where, I summarise the key findings of the project, provide insights, and propose justified next steps. 
4. Marks assgned for the quality of the write up
   Polised narrative, figure quality and correct referencing. 


7500 words

## Collated Data
In this section, I want to summarise for myself the campaigns I've run and the kind of runs I've been making. Then, I intend to use the key results from here, presentation i've made amongst others to make a presentation that is worthy of the final report. 


Content
1. Absorbing boundary conditions reflectivity: To motivate paramter choices for CAP 
2. Jellium CAP baselines: Consider 
12. simga 0.5 classical sweep: To motivate the paramter of wavepacket chosen
20. Localised jellium GS study: To motivate that the GS is well understood. Insights from here can help me add nuggets of insight into the reports

(I would want to pursue a few of these)
Campaigns done internally to fix results 
3. Sopping power and fourier training

4. Debugging quantum SP
7. Jellium sigma sweep
8. Energy oscillatio ndiagnosis
9. CAP wrap around fix
10. Localised jellium + scattering
14. Annular jellim S(v) vs wall r_s
15. Energy book keeping analysis
18. KS orbital vs HF orbital
19. Graphene and cornen CAP scattering
21. Pairwise comparisoiion of Quantum vs classical WP stopping 
22. WP and classical twin
24. ML patetern finding 
27. Nazarov Gross mass sweep
28. Li multi k w_peak draft
29. sigma =- 0.5 production absorber baselines
30. High density classical S(v)
33. Stopping from decomposed ledger
35. Local jellium dynamics 
35. Effective mass tuned bands (what's this? I think, this is to mention that different types of materials, metals, semi condictions etc can be modelled using jellium)
37. Sigma effect on stopping power
38. Grazing wavepacket off graphene or coronene (with CAP)

Additional Tangential Work
1. Quantum Kick Extension



## Strucutre

### Suggested Structure in the Project Briefing
1. Refined project direction
2. Code validation and verification
3. Results
4. post processing and resutls
5. Discussion and conclusions
5. References

### Running Thoughts About Sections Needed
I think, we'd need a mini-introduction section for report 2 which will be focused on the 

1. Introduction (Does this have the theory requried?)
(Literature Review requried?)
2. Methods (introduce the new methods that I've applied and what has changed from before)
3. 


### My thoughts on how the narrative of the report must be like
Here are my thoughts given all of the information I've collated here. What do I think I should add to this report. What are my real contributions? 

In the second leg of the project, I've mainly focused on trying to come up with a definition of stopping power for classical and wavepacket cases in jellium. The core idea that I've pursued is localised jellium. Now, localised jellium has differnet electrostatics to bulk jellium, due to surface effects, however, we've chosen this target as this provides a way to define stopping power in an intuitive way. The idea is simple, we shoot the wavepacket into the jellium. There are parts of the projectile that are reflected and the transmitted. We use absorbing boundary potentials to capture the lost potentials.  Now, in such systems, I've been trying to use different ways to defiuning stopping power 
1. Stopping power using the energy absorbed by the localised jellium system
2. Stopping power using the KS orbital of the wavepacket
3. Stopping power using energy decomposition that I've done

I would want to plot the resutls and then examine the results critically. This is probably the end game in terms of what results I've obtained. here, we can then extend the definitions into different velocity and other regimes. Perhaps, I can also complete the different things such as Nazarov Gross. 

For this narrative, I need to come up with a reasonable enough structure that can use. 

At some points, perhpas, I would also require to make workflow diagrams showcasing the method changes I've done to have heavy electrons, new classical projectiles, CAPs etc. #


## Draft 1 Skeleton: 

1. Introduction:

- Connection to the previous report: Mention that wave-packet propagation validated in vacuum and on coronene. Then, some initial analysis of the wavepacket propagation in bulk jellium were considered. Then, loss function was suggested as a method to arrive at stopping power. 

- The resolution failure of the loss-function route: The energy resolution of the loss function scales inversely with total propagation time. The resolution achieved was insufficient to trust the extracted stopping power. The computational cost to achieve required resolution is very high. 

Figure 1. Loss function computed from the achievable propagation time, with the target resolution marked. Shows directly that the peak structure needed for stopping power extraction is not resolved.

- The periodic boundary condition problem: Under periodic boundary conditions the projectile re-enters the cell and re-interacts with the target contaminating the energy readings. This makes a single traversal impossible to isolate.


- Objectives of this report: Construct stopping power definitions that do not depend on the loss function. Build a geometry in which those definitions are simultaneously computable.
Compare wave-packet and classical projectiles within that geometry. Establish where the quantum description departs from the classical one.


2. Methods

Opening line of this section should signpost explicitly that Sections 2.1 and 2.2 constitute the refinement of project direction and subsequent development, so that it is findable against the marking scheme.

2.1 Complex absorbing potentials (Boundaries)
- Introduce how CAPs are and how they are implemented in INQ (with my modification)

Figure 2: Workflow schematic diagram explaining how Complex Absorbing Potentials (CAPs), Jellium Slabs and classical projectiles were implemented in the code. Mention clearly the modifications that were made. 


2.2 Jellium Slab: 
- Move to a finite jellium slab with open boundaries along the propagation direction. 
- Absorb outgoing amplitude with a complex absorbing potential.
- Mention and justify grid spacing, time step, cell size, slab thickness (and others). 
- Observables include decomposed energy components E_hartree, E_kinetic, E_external and E_ss, E_ps, E_pp, E_sb etc. (where s - slab, p - projectile, b - background)
- Mention that the aim is to obtain clear energy observable. As the CAPs absorb density that flows out, we would have a case where the total energy would plateau. Hence, we can arrive at an estimate of the absorbed energy by the jellium slab. 
- Mention the parameters of the jellium slab used, justify



2.3 Classical Projectile 
- Using a pseudopotential (gaussian potential)
- Using native INQ perturbation method
- Justify the gaussian width of the radial potential (in the pseudopotential). 

Figure 3: S(v) plot of classical projectile modelled using pseudopotential with different widths. Overlay contains the literature expectation of the stopping power. By comparison, sigma = 0.5 bohr was chosen as the ideal. 

2.4 Mass of electrons
-  Mention clearly how the mass of electrons was changes, and specifically, what equations were altered for this to be changed. 

3. Results and Discussion


  - 3.1 Wave-packet through the jellium slab
      3.1.1 The slab and its ground state
        Keep minimal. Ground-state density profile shown only as confirmation the setup is as intended; refer back to the setup figure in Section 2.2 rather than repeating it.
        Quantify the self-interaction energy for a single electron in this setup, and flag that it will reappear in Section 3.3.
      
      Figure 5. Total energy against time, before and after the absorbing potential correction.
      
      3.1.2 Anatomy of a single run
      One representative run dissected in full: reflected and transmitted fractions, wave-packet spreading, energy channels against time.
      Demonstration that the energy ledger closes.
   
     Figure 6. Stacked decomposed energies (E_pp, E_ss etc. ) against time for a single wave-packet run, with the residual on a second axis showing closure. Annotate which slice each of the three definitions reads.
   
      3.2 Classical versus wave-packet stopping power
      3.2.1 Comparison of the three definitions
      Definition one: from the Kohn-Sham orbital of the projectile.
      Definition two: from the energy absorbed by the slab.
      Definition three: from the decomposed energy from simulation runs.

      For each: what it physically measures, what it assumes, what it cannot separate.
      The closure relation between the three, following from energy conservation.
      Treatment of the reflected component: which definitions can isolate the transmitted part and which cannot, and how this is handled.
      Treatment of the traversal length: the operational choice, and its sensitivity.
      
      [OPEN] Do they agree once the ledger is closed? Agreement and disagreement need different framing and a different figure.
      If they agree, the quantity is robust and the agreement is the result.
      If they disagree, the size and sign of the disagreement is itself the result, and the explanation is the physics.

      Figure 7. Stopping power S(v) as a function of velocity for wave-packet and classical projectiles, with the three wave-packet definitions shown separately.


      3.3 Wave-packet against classical
      The matched comparison, both projectiles through the identical slab.
      Report the magnitude of the excess of the wave-packet result over the classical reference.
     
      Figure 7. Stopping power for wave-packet and classical projectiles, with the three wave-packet definitions shown separately.

      3.3.1 Where the excess comes from
     
      [OPEN]

         Figure 8. Decomposition of the wave-packet minus classical gap into contributing terms.
      
      
      
      3.4 Approach to the classical limit
      Two available control parameters: projectile velocity and projectile mass.
      Using a combination of these to find the classical limit for the wavepacket. Simulating it, and finally, finding the difference. 
      
      
      Mass is the cleaner knob: increasing mass shortens the de Broglie wavelength and suppresses spreading without changing the excitation spectrum being probed.
      Identify the regime in which the quantum result converges towards the classical one, and quantify the deviation where it does not.
      [OPEN] Number of converged velocities determines whether this is a curve or a small set of points.
      Figure 9. Stopping power against projectile mass, or against velocity, with the classical reference marked as the limiting value.
      
      3.4.1 Structural signature in the induced density
      Motivation: the energy-based comparison is compromised by a known deficiency of the functional, so seek a signature independent of absolute energies.
      Mode decomposition of the induced density for both projectile types; state the rank in each case.
      Argue this is a robust quantum fingerprint that survives the self-interaction problem.
      Figure 10. Leading modes of the induced density for both cases, with the mode energy spectrum showing the difference in rank.




4. Conclusion

[OPEN] To be decieded based on the results

5. Future Directions
5.1 Layered jellium for real materials: Successive layers of differing density of jellium to emulate the response of real materials (PENN paper).
Follows directly from having a working finite slab.
5.2 Tuned band structure: Changing effective mass of electrons to emulate semiconductors and other band structures (Changing ).
5.3 Comparison with TD-HF: time-dependent Hartree-Fock comparison, to test whether the projectile orbital is well approximated.
5.4 Extension of the definitions: Wider velocity range, alternative target geometries, comparison against the loss-function route once resolution is affordable.
5.5 Nazarov Gross test
5.6 sweep of gaussian width of the wavepacket

Appendices
A. Full absorbing potential reflectivity curves
