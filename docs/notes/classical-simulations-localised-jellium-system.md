I want to run classical simulations for the localised jellium system. This, way, we can monitor the actual stopping power. We do not need CAPs here. This establishes a baseline which we can then compare with the wavepacket runs. We use perturbations instead. So, after we have master the setup, and understand what's happening, we can then use this to set our expectations for what benchmark the quantum stopping should adhere to. 

The core idea is as such. We let the perturbation pass representing a particle propagate in the cell. The perturbation projectile propagates in the simulation cell until it reaches the end. I have to think carefully about how I would think of the edges. Here is my current thought process. I think, because we have a perturbation (represented by a gaussian charge particle), we can think of the charge distribution it represents moving out of the simulation cell (without wrapping). This means that, at a certain tim, the bulk of the perturbation potential would have moved out of the cell. We need to think carefully about how we can implement this and if this is feasible. 

Now, the system's core principle is the same as that of the localised system. We have a localised system. Then, we monitor the rise in the total energy in this entire interaction. We would expect, as the projectile moves out of the simulation cell that the total energy would plateau at a maximum. The delta energy would represent the energy gained by the slab. Then, we can formulate the E_gained/L_slab_z to be a measure of the quantum stopping power. The other way to do this is, look at the delta E_total graph, and then use the stopping power skill to come up with a metric for the stopping power. 

I want to start with one simple setup for which this can be run. Consider 82 electrons in the jellium slab density. The same Lz as we've been using for the localised jellium. 



In this run, we consider a localised jellium system. We want to establish the expectation of what the stopping power must be using these runs. 

### Using a perturbation