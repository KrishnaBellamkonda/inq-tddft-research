In this markdown file, I am going to mind dumo the ideas I have for different campaigns. We will deliberate on each of them individually for a bit. Then, with the preliminary ideas, write rough drafts for these campaigns. Each campiagn idea will be enclosed in <campaign_idea></campiagn_idea> tags. For each of the campaign ideas, a prose of text is given to summarise the background and the problems this is trying to solve. When we deal with each campaign idea individually, while running the /campaigns command, we need to effectively come up with a plan, after careful thinking and researching, that can be executed to help through light on the questions or tasked posed.  

## Fixing the classical projectile. 
<campaign_idea>
In this campaign, we are talking about the localised jellium system.

<brief>
In this campaign, the aim is to come up with a classical version of the localised jellium system run, such that stopping power can be calculated effectively. There needs to be some benchmarking that needs to be done of the classical campaign with known results. These known results have to be sought and compared against for this type of simulation. 
</brief>
    
<mind_dump>

In the previus campaign that I ran, called quantum stopping power, the clasical projectile runs had a few problems. The most obvious way of claculating the stopping power in the classical projectile scenario, where ehrenfest is at work, is to use the delta E over time, and then find its gradient to find the stopping power. The second method, is the method used in the ase of quantum wavepacket, where we find out the total energy gained by the localised jellium system and then divide this by the length of traversal. Doing this, we find the stopping power. Both these definitions have problems.  

The problems are - 
1. Unlike the wavepacket, the CAP does not work with the classical projectile. So, in a few runs, I had the classical projectile re-enter simulation box. This changed the total energy dynamics of the run. I belive this was fixed by parking the classical electron on one of the edges. However, this produced artefacts that are not real (due to the abrupt stopping). 
2. Using the delta E vs time graph, and finding the slope, method is hard, because, the delta E, looks like an oscillation. I don't quite understand why this happens. Another campaign, that work on the understanding the localised jellium system better is probably to be done first, and then look at this system. i think, combining the knowledge from that campaign, and any screening, wake effects, and long range interaction effects in the box due to the coulomb potential (perhpas, need to have a long range cutoff of the classical projectile radial potential too, to avoid having loop around the box effects, need to think and brainstorm).
3. Using the quantum stopping like method, we would want the system to come to an equilibrium, before we calculate the E_final. In the classical runs I've seen, I don't think this has happened, chiefly due to the problems mentioned above. 

</mind_dump>

</campiagn_idea>

## Wide wavepacket campaign at the jellium slab system
<campaign_idea>

<brief>
To come up with a simulation system that is capable of identifying the artefacts that arise due to pure quantum effects that are missing in classical projectiles. To ensure this, we select a wave packet and classical projectile with the same gaussian width. Then, we ensure that this width ensures that the wavepacket does not change its width appreciably. This makes it convenient to identify the purely quantum effects from here.  
</brief>
    
<mind_dump>
There already might be a draft plan made, so I need to you to check if this already exists. The core idea is that, in the quantum stopping power run that we've made, the gaussian width of the wave packet is 0.5 bohr. So, by the time it reaches the jellium slab, it would have expanded in size appreciably. Now, its comparison is again the classical projectile. In this clasical projectile, the sigma of the radial potential does not change. So, I believe, there is a change in the coulomb interaction between the two. Also, I think, the effective interaction radius of the wavepacket increases as compared to the radial particle. This might also cause some problems. So, I want to plan a run where, the gaussian wavepacket does not spread appreciably, and its gaussian width is almost the same as its original gaussian width. In free propagation wavepacket, the final sigma depends on the original sigma and time only. So, we need to choose these two variables correctly in order to avoid gaussian broadening. Now, to achieve this, we can use the gaussian broadening equation, set to it a maximum sigma_final, and fix the sigma. I expect a large sigma where this happens, and this would mean that we need to adjust the cell sizes, the energy (and hence the total time of the simulation)
of the wavepacket amongst other parameters. 

The core idea is to isolate the purely "quantum effects" from the wave packet interaction and make it directly comparable to the classical trajectory. Meaning, as we ensure that the interaction range of the two is the same (or on a similar scale), we are able to understand what the purely quantum effects are. 


We might have to think out of the box to make design this simulation, and we would need to conduct an extensive research on different methods that can applied to achieve this purpose. 

I think, the best way to work on this is performing one simulation experiment at a time. Meaning, we start by looking at a single system setup, run it. I look at the results, and give you feedback. Then, we improve the simulation design considering the points, and then retry. This iteration would go on until we have a system that we belive works - meaning quantum stopping power can be calculated clearly from the data. Once, we are sure this works, we can then run a sweep among the energy of the wavepacket, and get a S(E)  plot for the quantum case (with varying simga sizes). 

We need an extensive suit of observables. I think the the observable suite, both raw and processed observables from here can be used here directly. I also like the current cadence of data recording. However, while deliberating this task, we are going to go through the observables too. 

</mind_dump>

</campiagn_idea>



## System of a cylinder jellium with a projectile passing throug hit
<campaign_idea>

<brief>
Make a system of a cylindrical jellium, placed along the z axis, through which a positive ion travels. Record important observables, such as the stopping power amongst others. Start with a simple system, then slowly increase complexity.  
</brief>
    
<mind_dump>
TDDFT-PENN paper introduced an idea that essentially says that a bunch of jelliums (with different r_s), can be usef to mimic the behaviour of different materials in different channels. The background for this task is that, it has been found that, in some systems (I don't remember), where the layer of material and the water are in close contact, when there is a flow of water, a current is produced, and when a curent flows, flow of water is produced. This is due to the dipole nature of water. I want you to find me core literature that help illustrate this. Some TDDFT methods have been used to quantify this. However, quantum effects have still not been accounted for. 

Water restricted in in carbon nanotubes, in my understanding (which is not very strong), can be used to emulate such systems. We can model such systems using layers of jellium. The idea is not quite clear in my mind as to how the the layers of jellium can be constructed; I would have to think about this. However, we can consider a tube of jellium for starters, parral to the z axis, through which a positive ion travels. Here, we can measure the stopping power of the ion. We can do this at different densities of jellium to derive some understanding. We would also need to use all the other observables I have in my armory and deploy all that are relevant. 

</mind_dump>

</campiagn_idea>




## Muon Campaign


<campaign_idea>


I belive the basics of the muon campaign are already written. I need to revisit it carefully, and check if all the specifications that are required from this campaign are right. I want you to help me understand all that is written in this rough draft. Then, I want to build an effective plan to execute all that is given in the specification. 

</campiagn_idea>



## Thoroguh exploration of the ground state of jellium slab system. The idea is to make analytical mental models to understand this. 
<campaign_idea>

<brief>
The idea of this campaign is to start from small building blocks, and gain a strong intuition for what's happening in this localised jellium slab system. I would need to make mental models using analogous ideas, for example (does not have to be this) of charged plated, or capacitors etc. to effectively understand. Then, the dynamics of what happens when a classical projectile goes through is also important for me to understand, again building from small controllable pieces. 
</brief>
    
<mind_dump>

One potential way to look at the localised jellium system, is to consider the ground state charge distribution. Then, by considering the mean, we can assign positive and negative sign to plates. Then, we can say that the localised jellium system is an alternating positive and negative signed plates. to be more specific, just before the actual extent of the jellium slab, there is smoothening potential, that makes it a net positive plate, then comes the tall peak in charge distribution near the boundary which is a negative plate, then come the internal region of the localised jellium system which is a net positive plate, and then comes in a net negative plate (at the boundary), and then a net positive plate beyond the limits of the localised jellium. When we think about the interaction of the classical projectile with this ground state, think about interaction with these systems instead in mind. As time changes, this this strucutre might change, but this picture can be useful. 

Also, there is reason to suspec the self interaction error magnitude in the simulation system. Now, currently, this is calculated using the E_total of the localised jellium system + the waevpacket at a certain distance. We then, subtract the E_GS of just the localised jellium system, T_Wavepacket (which is its initial kinetic energy plus the localisation cost). We call this the self interaction error. However, this is not quite right. I would assume the classical repulsion between the projectile and the jellium system is not accounted for. 

Also something to consider deeply is that coulomb interaction between the localised jellium and the wavepacket at the initial instant. I should use classical electrostatic to find the magnitude of the coulomb repulsion effect on the system. Perhaps, for different distances of the wavepacket from the localised jellium slab, I need to measure the E_total difference between the wave packet and the classical runs. This would indicate me what what's happening. This would give us how this interaction effect would behave as a function of the radial distance r from the jellium system. Also I think, the contribtion that does not disappear as we elongate the box can be thought of as the Self interaction effect, if we say that E_total(t=0) - E_GS - T_WP = E_SIE + E_WP_localised_jellium_repulson. I need to create the right system to test this on. Perhpas, a very elognated z axis, with the same x and y axes dimensions. This is to ensure that we far away. 

Along that point also comes the question of the cutoff for long range interactions in the clasical projectile. This is an important consideration I have not thought very deeply about yet. I need to conduct a search over the internet, find literature that has grappled with this problem, and find effective solutions for my system. 

Another open question to consider is that when we compare the classical projectile run to the wavepacket run, does the localisation cost, which can also be called the zero point error, make a huge difference? I assumed no, but we should consider this question deeply.

I also need to examine the effect of different configuration options of the localised jellium system and take them into account. Meaning, how does changing the x, y and z lengths (one change at a time), changing the w  (boundary smoothening term where erfc is used) is important. For each configuration, it is also important to understand how the surface effects change. If there is any literature on localised jellium, it would be worth it to collect it, read and understand it (you read it and I read it individually). Then, I might start making sense of any other effects that are taking place. 

Furthermore, there might be other clues in the already run classical projectile simulations on localised jellium run that can inform us something. So, I would also want you to scout for them using the VTI file, and perhaps, some algorithm that helps identify behaviours. This could help us understanding what's happening in these simulations better. 

I need to start building from simple systems, understanding them well, and then adding complexity step by step to build strong mental models for what's happening. Meaning, perhpas, I can consider simpler versions of the system to understand specific effects better. Then, slowly, I increase the complexity of these systems by adding one piece at a time. This would help answer some of the questions presented above. Then, I would be able to get a strong intuition for what's happening here. Also, some analytical work would also be good. If you point me to some literature that is relevant to what I mentioned here, then I can start reading it. I need to make a worksheet to understand the basics of jellium slab, and what it might mean. So, point me towards some literature that I can start working with.  

Now, I also have to come up with an effective plan that is capable of answering these questions. This has to be thought through carefully. 
</mind_dump>

</campiagn_idea>

## Changes to the campaigns skill

A few things I want to add to the campaigns skill-
1. I want to be able to research and sanity check plans while I am making it. So, I think, we need to restructure the different phases of a campaign to accomodate this. This is because, in many a cases, a concrete idea would not have taken root in my mind. So, a deal of checking with literature, brainstorming with you is reqruired to concretise the idea. We need to add this as one of the phases, or allow for such discussion in different stages of the skill. 
2. When making a plan for a campaign, it is always good to build simple experiments first to conclude and validate building blocks. Then, add complexity step by step, until a system is reached, and is to my liking, which can then be autonomously be run. Smoke test, and other tests to test simple hypotheses can be used. 