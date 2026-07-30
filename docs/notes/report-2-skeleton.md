In this markdown, I am going to think and write down content as I am thinking. I am going to do this in different sections so that I can segregate the points I am thinking about and verifying and testing things before I write something. Here, I am almost going to spell out everything I want to write in terms of my ideas for the first draft. 

# Skeleton

## Outline
1. Introduction
    - Why Jellium (Penn paper) and experimental relevance established (need to research about this). 
    - Connection to previous report
    - Resolution failure of the loss function route (need to illustrate using a cost argument, also tradeoff between the wavepacket spreading and sampling time.)
    - periodic boundary condition problem (can mention this normally I think), 
    - objectives of the report (need to refer to the original project description document)
2. Methods
    (Here I need to describe the new simulation systems carefully, again mentioning the reason for everything here. )
    - complex absorbing potentials ()
    - jellium slab ()
    - classical projectile
    - mass of electrons
    - Stopping power definitions (In this section, we go over all t)
3. Results and Discussion 
    (Depending on the results I want to include, I need to modify methods)
    - Bulk Jellium studies
        (this can be thought of as KS orbital dependent definitions of stopping power). (yao-schielfe paper)
        - Comparison of the two stopping power definitions and bechmark against classical particles
        (Here I might need to explain the differences between the analytical, classical and the bulk wavepaket results. I already have a few simulations here. Need to re run wavepacket)

    - Jellium Slab studies 
        - Wavepacket through the jellium slab
            - The slab and its ground state
            - Anatomy of a single run
        - Classical versu wavepacket stopping power 
        - Approach to classical limit
            - Structural signature in induced density
        
    - Thin films (using suggested stopping power definition)


4. Conclusion
    (I don't want this section to be just a relay of information, I want this to add knowledge. So, this is where I bring all that I've learnt together)

5. Future Directions


Appendices

## Changs in report 1
1. Change the introduction of the report using the argument in the presentations
2. Add a figure in the appendix showcasing the classical projecile induced density patterns. Will use this as an argument for the introduction


References
1. Stopping Cross Sections for Protons Across Different Phases of Water - F. Matias
Phys. Rev. Lett. 135, 148003 – Published 1 October, 2025
2. TDDFT-Penn approach: F. Matias, T. F. Silva, N. E. Koval, J. J. N. Pereira, P. C. G. Antunes, P. T. D. Siqueira, M. H. Tabacniks, H. Yoriyaz, J. M. B. Shorto, and P. L. Grande, Efficient computational modeling of electronic stopping power of organic polymers for proton therapy optimization, Sci. Rep. 14, 9868 (2024).

References to consider
1. Radiation track, DNA damage and response—a review
    H Nikjoo et al 2016 Rep. Prog. Phys. 79 116601
    - This paper is a review about the radiation, and how they cause DNA damage. 
    - Key point: the electornic stopping power of these materials sets the initial conditions for molecular mechnaism responsible for radiation damage. 

### Introduction

1. Introduction
    - Why Jellium (Penn paper) and experimental relevance established (need to research about this). 
        - 1. Jellium is good because, as Penn paper suggests, we can readily approximate the stopping power of metals readily using this approach (as done in this study). 
        - 2. Moreover, given the experimental loss function of a certain material, we can use the loss function to then construct a tddft-penn approach to find its stopping power. 
        - 3. As we are concerned with the quantum stopping power, we will be able to find out the quantum stopping power for any material concerned. 
        - 4. Moreover, relevant to bilogical systems, for which water is usually taken as the target, such tddft-penn has been conducted. 
        - 5. Combined with the fact that, previous studies showed that 10 eV is the regime where there is a lot of DNA damage, where quantum effect of electrons become highly prominent, we can motivate the purpose of the project. 
        - 4. Furthermore, EELS connection? LEED connection? etc?
    - Connection to previous report
        - 1. Improving upon the previous suggestions of using the KS orbtial energy as an indicator of stopping power. However, there are some quetions left to answer, what is the right gaussian width of the wavepacket? Secondly, some computational artefacts produced with using a coulombic potential are to be considered. Furthermore, the question about the position of the projectile using the wavepacket approach is yet to be looked at carefully. 
        - 2. The suggestion of using the loss function as to define stopping power is challenged and understood. Resolution failure of the loss function route (need to illustrate using a cost argument, also tradeoff between the wavepacket spreading and sampling time.)
        - 3. Also, there are periodic boundary effects that corrupts both the wavepacket and the classical simulations. So, periodic boundary condition problem (can mention this normally I think), 
    - objectives of the report (need to refer to the original project description document)
        - Analysis of stopig power definitions. 
        - Quantum to classical cross over and deviation quantification. (I think, I can make deviation quantification for the different definitions I have. However, I also have to )
        - need to link points from the presentation. 

Narrative: 
1. Results from the previous studies show qualitative trends, prompt questions, and more importantly point towards areas of improvement in the methodolgy.
2. Describe the improvement points
3. This body of work aims at improve the method in an attempt to obtain achieve the objectives of the project. 
4. Meaningful


Preliminary studies using the wavepacket driven rt-TDDFT framework reproduced analytical benchmarks accurately and produced qualitative trends consistent with the underlying physics. However, these early results also raise open questions about the target system to consider and the specific methodology needed to meet the objective of quantifying the quantum contribution to stopping power.

Jellium is considered to be the ideal target system for this purpose, given the interpretability of its results, its experimental relevance, and its ability to generalise to a range of materials and biological systems. Although jellium's relevance to metals is apparent, several studies show that TDDFT combined with the Penn model can approximate the electronic stopping power for a much broader range of materials. In this approach, a material's experimental energy loss function is used to decompose it into a weighted combination of jellium systems of different densities [TDDFT-Penn], each of whose stopping powers is computed using TDDFT. A weighted average of these jellium stopping powers then approximates the stopping power of the material itself. This method has proven especially relevant for biological systems, where water serves as the primary target of interest [Stopping Cross Sections]. Given that the quantum contribution to radiation damage, especially DNA damage in cases of radiation therapy, makes jellium the natural system on which to refine the methodology and arrive at a physically grounded description of electronic stopping power [Radiation track]. These results can be compared directly to experimental probes such as EELS for further validation [reference].

Previously, two main approaches for calculating the stopping power were considered, one based on the Kohn-Sham (KS) orbitals of the wavepacket and the other on the loss function of the system. The KS orbital based stopping power depends on the Gaussian width of the wavepacket. Figure 11 shows that the total stopping power decreases as the width of the wavepacket increases. Taking this trend to its logical extreme, an infinitely wide wavepacket amounts to a plane wave, and the interaction with the jellium target would then be expected to vanish. This is because the target itself is composed of electrons occupying plane wave states, so a projectile in a plane wave state would traverse it without any scattering events occurring. The trend is therefore physically consistent. However, this dependence raises an obvious question of which Gaussian width represents a classical point particle, and therefore which choice makes the results comparable to a classical treatment.

The loss function offers a complementary route. It carries rich information about the mechanisms through which transitions occur in the system, specifically electron-hole transitions and plasmon excitations. Combined with the overlap of the propagated KS orbitals with the initial ones, it identifies which transitions are active and how strongly each contributes. This makes it a powerful interpretive tool, but as a route to extracting the stopping power it carries a fundamental limitation. The loss function $L(q,\omega)$ is obtained by Fourier transforming the density response over the finite time window of the simulation, so its frequency resolution is set by the total propagation time $T$ alone,

$$\Delta\omega = \frac{2\pi}{T}, \qquad \Delta E = \frac{h}{T}$$.

For the system considered previously, a 20 eV wavepacket in a 50 bohr cubic box, a prominent energy transfer feature is the plasmon at 3.5 eV. The momentum of the projectile sets a traversal time of only 41 a.u., which resolves the spectrum no more finely than about 15 eV and places the entire plasmon peak inside a single frequency bin. Resolving it requires a propagation of roughly 2000 a.u., close to two hundred times the traversal time, which makes the route computationally infeasible as a general extraction method.

Furthermore, the signal itself degrades on exactly these timescales. Under periodic boundary conditions the projectile wraps around after traversing the cell, and the transverse broadening of the wavepacket wraps it in the $x$ and $y$ directions as well, as seen in Figure 4. The resulting self-interference corrupts $L(q,\omega)$ before the required resolution is reached. Short runs therefore carry no resolution and long runs carry a corrupted spectrum.

Methodologically, further observations come to the front. To model classical projectiles, a bare Coulomb potential was used. The $1/r$ singularity of this potential cannot be represented in a plane wave basis with a finite energy cutoff, so the potential is effectively smoothed at short range. This distorts precisely the close collisions that dominate large energy transfers, and the physics of the high momentum transfer regime is therefore not captured reliably. An improved computational treatment of the classical projectile would consequently be required to build on these results.

The methods developed in this work address each of these issues directly. This work provides the development of a robust methodology that establishes the conditions under which the calculations remain reliable and sets out how the approach extends to the specific systems and regimes of interest. The objective throughout remains the same, which is to quantify the quantum stopping power and to isolate the quantum contribution within it, so that the point at which it departs from a classical description of the projectile can be identified.


### Methods

2. Methods
    (Here I need to describe the new simulation systems carefully, again mentioning the reason for everything here. )
    - complex absorbing potentials ()
    - jellium slab ()
    - classical projectile
    - mass of electrons
    - Stopping power definitions (In this section, we go over all t)

Refined method
- So new methods section can look as such - 
    - First, miscellaneous 
    - Orbital Free definition and system
        - Bulk jellium
    - Orbtial dependent definition and system
        - jellium slab
        - absorbing potentials (correction of error in K.E orbtials)
    - Classical projectiles
        - perturbation and pseudopotential method as a gaussian.
    - Further work
        - mass of electrons changed. 

### Results and Discussion






# Thinking section

I now want to think carefully about what results I want to showcase. By the suggestion of professor, 

1. Bulk jellium results (might requrie re-runs? or can I calcualte the stopping power from here? I would need the wavepacket KS orbtial. Need to think carefully about how I can get the momentum from the information I have. Make some tests, compile the results and make the required plots.) I have done some checking already here. Perhaps, I can extend this to include two different density results that i can tabulate. So, we only would have one definition here. 
    - Professor suggested that we have two things we can make $<p^2>/2m$ and $<p>^2/2m$. So, we should compute these stopping powers from here too. For this, I need to think about what INQ does internally when it gives the KS orbtial energy. Then, I should also be able to compute the <p>^2/2m component, that avoids the sigma_p dependent term. This can give us some information relative to each other. I don't know if this warrants other runs to be made (fresh one, which I would prefer so that I have clean resutls. For this, I can do a bunch of things, energy decomosition, and others. This way, I will have a clean dataset I am working with). 
    -  Now, 

2. Jellium slab resutls as a way to showcase that we are going beyond the KS orbtia


## Points to address
1. Need to adress the point that a wavepacket with the same momentum as a classical particle has a loclaisation energy. this is a fundamental difference between the two runs. So, we need to apply this lens to the simulations. 
2. I need to address the point that K.E of the wavepacket can be decomposed into <p>^2/2m term, which is dependent on momentum and the other term that depends on the variance of the wavpacket. In our definitions, we need to be careful about, what each definition means, and what information we can extract out of each. 
3. Some feedback ive to me was that, I also ned to establish a clearer pathway between experimental relvane and why these particualr systems are important. 
4. In comparisons, benchmarking against classical projectile is probably a better approach. This is because, the analytical expection (from bulk jellium), assumes that all channels of plasmons etc are available. however, when they are not available, then, we might face problems of the stopping power being under estimated. So, this can be the motivation for using classical projectile.
5. Before re-drafting the structure, I need to refer to papers that I respect in the field, and use them as a guide to re-formatting the structure.  
6. Perhaps, I need to organise my methods section as such - each system I consider aids a certain definition of stopping power. Orbtial dependent and orbital free definitions. Then, other stuff, such as classical projectiles and chanigng the mass of electrons to control specific effects. Classical projectiles have been used to give a better benchmark. This is because, the systems I've considered might not allow for specific plasmons or other effects. so, to quantify quantum effects, we need to compare with classical projecitles and understand the results. 
    - So new methods section can look as such - 
        - First, miscelleaneous 
        - Orbital Free definition and system
            - Bulk jellium
        - Orbtial dependent definition and system
            - jellium slab
            - absorbing potentials (correction of error in K.E orbtials)
        - Classical projectiles
            - perturbation and pseudopotential method as a gaussian.
        - Further work
            - mass of electrons changed. 
    
    - Points to mention:
        - Also, I need to talk about the chosen parameters for each of the things here. For example, paramter for classical projectile
        - Need to mention periodicity of the axes in each system
        - Need to mention the width w of the gaussian smoothening for the jellium slab. Rounded to the closest even number (changed the number of electrons, sometimes width). (Need to cite paper here) 
        - 
    - Thinking points
        - For method section, I need to organise everything logically. Forget about being succinct. I need to be clear first, and ensure everything that has to be said is said clearly. 

7. I need to mention that codebase has been changed in this project. Code can be perused to look at the changes. 
8. General thought, I want to selective with plots here. Any plot that can be explained in analysis, such as time animated energy decomposition must be produced in report-2 folder. But, only the analysis must make it into the report.
9. I want to go through all the work I've done, and check if there is anything else that is relevant in here.
    - Some stuff I would like to include. For the orbital free and orbital dependent definitions of stopping power, it would be good to compare the effects in different densities of slab and bulk jellium. This is because, I want to understand the effect on materials that are relevant for us to understand. And also understand if there is a density dependence of this stopping power. 
    - Nazarov Gross check would be something crucial I can do. I think, hopefully with the fix in WP run, I can run this study and conclude something important. 

10. Here, I think about what is relevant to include in the appendix. I would have a better clarify after completing the first draft of report-2. This is because, I would understand, what is essential to mention. 
    - Reflectivity curves of CAP. This can justify the selection of the CAP parameters for specific runs. 
11. Checklist of all the objectives to meet in this dissertation project - 
    - rt-TDDFT code development: 
        - systematically investaigate the influence of different wavepacket parameters on the ensuing quantum dyunamics
            - including wavepacket spreading (done) with electron density  
            - splitting with electron density (not done)
            - and entanglement with target electron density (not within scope)
    - analysis of stopping power definitions
        I have enough content to pass this checkpoint
        - orbital dependent definitions
        - orbital free definitions
    - quantum to classical crossover and deviation quantification
        - establish a link between wavepacket and classical point particle limit (the meaning of this is important to understand)
        - wave packet coherence, exchange and diffraction
        - finding the results for different regimes
    
    What objectives that we set out to achieve were achieved in this project?
    1. coming up with sensible definitions of stopping power
    2. able to compare classical and wavepacket effects
    3. (need to establish correspondence limit, but I think, once I have both the classical and wavepacket jellium slab runs ready, I can then understand the impact of these differences. This might lead me to correspondence limit.)
    4. (On the point of systematic explanation of the stopping power with different wavepacket parameters, i think I have some information here. I have some intuition here, but I need to quantify it and perhaps add to my explanation)
