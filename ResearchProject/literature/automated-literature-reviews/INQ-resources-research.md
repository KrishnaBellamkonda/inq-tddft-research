# Deep research roadmap for TDDFT wave‑packet–jellium simulations using INQ and Octopus

## INQ at a glance

INQ is positioned as an “engine” for electronic-structure calculations that can be used (i) as a standalone electronic-structure code, (ii) as a library to build higher-level methods, or (iii) as a proxy-app for performance studies on HPC platforms. citeturn21view0turn30view0 In practice, this “library-first” framing matters for your plan to *fork and add functionality*: you can often build your new capability *on top of* INQ (via its C++ and/or Python interfaces) without needing to deeply rework the core. citeturn23view0turn30view0

Key technical signals, from INQ’s public descriptions and documentation:

- INQ targets DFT and TDDFT (including real-time propagation), with explicit emphasis on GPUs and multi-level parallelism (MPI + threads + GPU backends). citeturn30view0turn29search15  
- The documentation set is clearly “in motion” (draft / actively developed), but there is a substantial amount of *practical* content: machine-specific build recipes, a command-line interface, a Python interface, and multi-step tutorials (molecules, optical spectra, and an electronic stopping workflow that is conceptually close to “projectile + electron gas”). citeturn21view0turn22view0turn24view2turn51view0  
- Core developers are publicly identified and provide direct contact emails in the repository README: entity["people","Xavier Andrade","inq developer, tddft"] (xavier@llnl.gov) and entity["people","Alfredo A. Correa","inq developer, computational physicist"] (correaa@llnl.gov). citeturn8view0

INQ development is associated with the entity["organization","Center for Non-Perturbative Studies of Functional Materials Under Non-Equilibrium Conditions","doe cms centre"] and entity["organization","Lawrence Livermore National Laboratory","livermore ca, us"] (LLNL). citeturn5search0turn8view0turn30view0


## Resource map for INQ

This section is organised as a “where to look / what you get” map. It is not literally *every* resource on the internet (no search can guarantee that), but it captures the main public surfaces: code + docs + tutorials + talks + papers + auxiliary repos.

### Canonical code and mirrors

The canonical project is hosted on entity["company","GitLab","devops platform"] under the NPNEQ group. citeturn16search0 A public mirror exists on entity["company","GitHub","code hosting"] (llnl/inq), explicitly labelled as a mirror pointing back to GitLab. citeturn30view0

Why you should care: sometimes the mirror loads more reliably (for browsing files), while GitLab is better for issues/MRs and the “source of truth” history. citeturn30view0turn16search0

### Documentation (two parallel doc sites)

There are two prominent doc builds:

- A entity["organization","Read the Docs","documentation hosting"] site (“draft documentation”) with a conventional Sphinx structure (Intro, Compilation/Installation, Units, Tutorials, Theory, Code), and a visible build revision hash. citeturn21view0  
- A more extensive “INQ 0.95 documentation” site (alphataubio.com) that contains a very large navigation tree: build instructions (including many supercomputers), interface reference (shell + Python), tutorials (shell/Python + C++), and a developer guide. citeturn22view0turn23view0

The alphataubio docs are particularly valuable for day-to-day productivity because they include practical recipes like: local builds, cluster builds, dependencies (MPI HDF5, Boost), and workflow-level tutorials that show what the code *expects you to do* as a user. citeturn22view0turn23view0

### Tutorials that are directly relevant to your “projectile + electron gas” thinking

INQ includes an **electronic stopping** tutorial (C++) that simulates a proton moving through aluminium and measures forces/energy loss during real-time propagation. It explicitly frames the phenomenon as electronic stopping (fast projectile interacting with electrons), uses TDDFT to allow electronic excitations, and demonstrates a “compute ground state → propagate” workflow. citeturn24view2turn24view1

This tutorial is conceptually close to your jellium goal because it is the same *class of nonequilibrium coupling* (a moving perturbation interacting with an electron sea), even though the host is atomistic Al rather than a uniform positive background. citeturn24view2

A few implementation-level details from the tutorial that are useful as “mental anchors”:

- Metals: it uses extra unoccupied states and an electronic temperature, and saves/reloads a ground-state restart directory to avoid repeating SCF. citeturn24view0  
- Real-time propagation: it uses an explicit propagation call (`real_time::propagate`) and an “impulsive dynamics” ion propagator so the ions move at fixed velocity (appropriate for fast electronic stopping where the host lattice motion is neglected). citeturn24view2turn24view1  
- It explicitly notes that INQ “doesn’t use defined units for velocity at the moment”, which is a small but important footgun when you’re trying to design clean, reproducible workflows. citeturn24view2

### Templates / “front-end” repo for building your own features

The docs recommend an auxiliary project, **inq_template**, which is a CMake project that downloads INQ + dependencies and compiles example programs, i.e. a “starter kit” for writing your own C++ driver code against the INQ library. citeturn23view0

If your “additional functionality” can be implemented in a driver layer (new observables, new external potentials, new time-dependent protocols, new analysis), **starting from inq_template is usually a lower-risk path than forking the core immediately**. That is exactly aligned with your fork plan: prototype outside, then upstream what truly belongs in core. citeturn23view0

### Talks, slides, and workshop material

INQ’s README points to at least two introductory videos (including a MolSSI workshop talk). citeturn8view0turn30view0 A set of public slides (e.g., “slides_inq_vesw2020”) also exists and includes compilation workflow hints and tutorial pointers. citeturn29search15

(Pragmatically: these talks are often where “missing documentation” is implicitly explained—data layouts, intended workflows, design philosophy—so they’re worth watching even if you never quote them.)

### Related auxiliary projects / proxy apps

An Exascale Computing Project proxy application **minq** is described as “a DFT and DFT-MD mini-app” built from INQ. citeturn16search14 This is useful if you want to understand performance-relevant kernels or a “minimal viable” usage pattern.

### Publicly provided contact points

The INQ README provides direct developer emails (xavier@llnl.gov, correaa@llnl.gov). citeturn8view0  
A public slide deck also displays xavier@llnl.gov in its header context. citeturn29search15

If you email: include (a) a 3–5 line description of your physical problem, (b) the smallest “reproducible” pseudo-input sketch you imagine, and (c) the specific capability gap you think requires a fork. That tends to get much more actionable replies than “how do I use the code?”. (This is process advice rather than a fact claim.)


## Codebase evolution and change tracking for INQ

You explicitly noted that “the GitLab has changed a bit”. The public release and commit metadata supports that: there are clear feature jumps and refactors.

### Releases and their headline changes

INQ publishes release entries with short bullet summaries:

- **v0.5.0 (Initial Release; 2021‑06‑10)**: “Basic DFT”, “Basic TDDFT”, “MPI and GPU Parallel”, installation scripts/tests, and installation notes for specific HPC sites (including LLNL Lassen and NERSC Cori). citeturn51view0  
- **v0.9 (2023‑01‑13)**: adds (or at least highlights) DFT mixing, TDDFT, hybrid functionals including exact exchange, ACE (“Lin Lin’s algorithm”) for fast hybrids, GPU‑GPU direct communication via NCCL, velocity‑gauge laser pump for TDDFT, and non‑orthogonal cells. citeturn51view0  
- **0.95 (2024‑09‑30)**: “Fixes in pseudopotential reader”, “Colinear spin”, “Python interface”, and “Command line interface”. citeturn51view0  

Two additional technical clues embedded in the 0.95 release metadata:

- The release commit is a merge of a branch named **force_reorganization** and explicitly mentions creating a **forces_stress** object. citeturn51view0  
- That suggests a structural refactor around forces/stress responsibilities—relevant if your extension involves custom forces, external potentials, or coupling to projectile degrees of freedom.

### Recent (2026) change signals

The current master branch has commits in **late March / early April 2026** that indicate non-trivial internal movement:

- A merge of **inq_paw_dev** described as “PAW related refactoring”. citeturn52view0  
- Changes around the Boost/CMake linkage: removing the Boost “system” component because it became “part of headers after Boost 1.69”, merged via a “fix-boost-cmake” branch. citeturn52view0  

This matters for you because it implies:

- Build instructions and dependency assumptions can drift; prefer doc pages that are clearly tied to a version, and consider pinning a tag (e.g., 0.95) for reproducibility while you prototype. citeturn22view0turn51view0turn52view0  
- If your planned feature touches wavefunction representations or pseudopotential/PAW machinery, you should expect active churn and thus higher merge-conflict risk. citeturn52view0

### A practical “history reconstruction” workflow

Given the above, the most robust way to “document changes the codebase has undergone” (for yourself, for a thesis appendix, or for lab onboarding) is:

- Use tagged releases (v0.5.0 → v0.9 → 0.95) as your primary milestones because they carry explicit bullet summaries. citeturn51view0  
- Treat large merge commits with descriptive branch names (e.g., *force_reorganization*, *inq_paw_dev*) as “architecture events” and read the associated merge request threads for rationale. citeturn51view0turn52view0  
- Cross-check any doc set you rely on against its stated version (alphataubio “git info: 0.95”; ReadTheDocs “inq 0.9” style labelling and revision hash). citeturn22view0turn21view0

### Licensing note for your future fork

INQ is distributed under the Mozilla Public License 2.0 (MPL‑2.0) in its public materials. citeturn8view0turn23view0turn30view0  
MPL has reciprocal requirements at the file level (i.e., modifications to MPL-covered files must remain available under MPL when distributed). The INQ tutorial documentation explicitly flags that INQ itself is MPL and contrasts it with the permissive licence of inq_template. citeturn23view0

If you anticipate private development, eventual publication, or collaboration with industry partners, it is worth making an early “licence boundary” decision: write as much new logic as possible in new files / external driver layers unless core changes are strictly required. citeturn23view0


## Methodology for wave‑packet and jellium coupling in RT‑TDDFT

You described the physical aim as a *quantum wave packet interacting with an electronic jellium system* and asked “how one would go about making this simulation”, including whether to compute ground states separately and then evolve together.

The key methodological decision is what “wave packet” actually is:

- If it is a **charged projectile (ion/proton) wave packet**: you are in the domain of coupled electron–nuclear quantum dynamics (beyond standard electronic TDDFT with classical nuclei).  
- If it is an **electron wave packet**: you are trying to inject/propagate an additional electron excitation interacting with a metallic electron gas, which raises open-boundary/continuum issues and Pauli/exchange consistency.

Below are the most defensible routes, anchored to existing workflows that are documented and/or widely used.

### Route A: Classical (or prescribed) projectile + RT‑TDDFT electrons

This is the “electronic stopping power” paradigm. The core workflow is:

1. Build the *target* electronic ground state (jellium slab, bulk metal, etc.), usually with smearing and extra states when metallic. citeturn24view0turn33view0turn33view2  
2. Place the projectile at an initial position (often in a big supercell to avoid self-interaction with periodic images). citeturn24view2  
3. Start real-time propagation of the Kohn–Sham orbitals with the projectile moving along a prescribed trajectory (constant velocity is a common first step; then more sophisticated Ehrenfest schemes are possible). citeturn24view1turn49search3  
4. Measure energy transfer or force on the projectile, and average appropriately.

INQ’s electronic stopping tutorial is exactly an instance of this design: it sets up an Al supercell, inserts a proton, computes a ground state (saved/reloaded), then propagates in real time, outputting forces and total energy and interpreting the drag force as stopping power. citeturn24view2turn24view1

For the *jellium* variant of Route A, Octopus provides explicit machinery (see Route C), but conceptually it is still: “ground state of electron gas + moving external perturbation → RT propagation”. citeturn33view0turn49search3

This route answers your “do we compute GS separately?” question as follows:

- You compute the target GS first (optionally with the projectile present but stationary at the starting location). citeturn24view0turn24view2  
- You then time-propagate the electronic system under a time-dependent external potential created by the projectile trajectory. citeturn24view2turn49search3  
- You do **not** “exchange KS orbitals with a wave packet”; instead, the projectile influences the electrons through its potential, and the electrons back-react via forces (if you enable it). citeturn24view2turn49search3

A particularly good conceptual reference for this whole class of simulations is entity["people","Alfredo A. Correa","electronic stopping review author"]’s review article on first-principles electronic stopping, which is written explicitly to help graduate students understand modern RT‑TDDFT stopping workflows and modelling choices. citeturn49search3

### Route B: Quantum nuclear wave packet + TDDFT electrons

If your projectile is genuinely a **quantum wave packet** (nuclear), a purely classical “impulsive dynamics” treatment is no longer faithful. You then need a formalism such as:

- “Exact factorization” style electron–nuclear equations of motion, or  
- “Nuclear–electronic orbital” (NEO) frameworks (more common for quantum protons) combined with time dependence.

INQ’s theory documentation tree explicitly includes an “Exact factorization” section with nuclear and electronic equations of motion listed as subtopics. citeturn22view0 This is a strong hint that the INQ ecosystem is at least *thinking in that direction*, even if you will still need to confirm (in the code) what is implemented versus what is theoretical background.

A practical strategy if you want to pursue Route B without getting buried:

- Prototype the coupling in the smallest dimensionality and basis that can still show the effect (often 1D/2D model potentials, or very small 3D systems).  
- Only when the physics protocol is stable do you lift it into large-scale, production geometries.

### Route C: Explicit jellium in Octopus as a “method sandbox” (highly recommended)

Even if you ultimately want INQ, Octopus is unusually good as a *method sandbox* for jellium because the documentation includes explicit **jellium and jellium slab** tutorials.

Octopus’ jellium tutorial shows how to compute the ground state of a uniform electron gas and a jellium slab, including how to parameterise the electron density by the Wigner–Seitz radius \(r_s\), how to set up periodic dimensions, and how to define jellium through the Species block (including `species_jellium_slab`). citeturn33view0

It also explicitly notes that jellium slab calculations show Friedel oscillations in the electron density near the surface, which is an immediate validation target for your jellium setup. citeturn33view0

Octopus’ manual further emphasises that “species” can be nuclei (with pseudopotentials), a jellium sphere, or a user-defined potential. citeturn33view3 This flexibility is exactly what you need for a “wave packet interacts with electron gas” prototype: you can represent the projectile as a custom time-dependent external potential (or a moving Gaussian charge distribution) and test your protocol before committing to a deeper fork.

### A concrete “minimal simulation design” you can implement quickly

This is a method-first recipe (software-agnostic) that maps cleanly onto both INQ (as seen in the stopping tutorial) and Octopus (via jellium slabs):

- Target: jellium bulk (3D periodic) or jellium slab (2D periodic + vacuum). Use \(r_s\) as the primary density knob. citeturn33view0  
- Ground state: include smearing/finite temperature and sufficient empty states because you are modelling a metallic system. citeturn33view0turn24view0  
- Projectile model: start with a classical prescribed trajectory and represent the projectile as an external potential (Coulomb-like with smoothing / pseudisation to avoid grid singularities). This mirrors standard electronic stopping workflows. citeturn24view2turn49search3  
- Real-time propagation: short time steps (attosecond-scale is common in real-time electron dynamics) and careful convergence tests. (INQ’s BPVE paper using INQ states a 0.5 as time step in one setup, which gives a sense of magnitude.) citeturn44view0  
- Observables: induced density, energy drift, force on projectile, and—if you care about “stopping”—the average drag force along the trajectory. citeturn24view2turn24view1turn49search3  
- Validation: reproduce known jellium slab density profiles (Friedel oscillations) and basic qualitative stopping curves (linear at low velocity, peak near a characteristic velocity scale, then decay), before attempting any sophisticated wave-packet quantum treatments. citeturn33view0turn24view1

### Evidence that INQ is used for large-scale nonequilibrium RT‑TDDFT workflows

A recent example where authors explicitly state they used INQ for real-time TDDFT is a paper on ballistic photocurrents in monolayer GeS. In the arXiv HTML version, they write that rt‑TDDFT is performed (PBE, norm-conserving pseudopotentials, plane-wave basis, specified cutoff and timestep) “as implemented in the INQ code,” within a fixed-ion supercell workflow. citeturn44view0turn39view0

This matters less for your jellium project directly and more as a confirmation that “ground state → real-time propagation” workflows are actively used in the INQ ecosystem beyond toy inputs.


## Octopus assessment for your use case

Octopus is a mature, general-purpose real-space electronic-structure code designed for DFT and TDDFT “virtual experimentation”. Its documentation explicitly states: electrons are treated with DFT/TDDFT, nuclei are treated classically as point particles, and electron–nucleus interactions are typically through pseudopotentials. citeturn33view1

### Documentation quality and support

Octopus has a large, structured online manual and tutorials, including a “Getting started” guide that explains the input format, default pseudopotentials, and how to use the variable reference / `oct-help`. citeturn33view4turn33view1

For community support, Octopus explicitly advertises three mailing lists (announce, users, devel) with subscription links, including a user list intended for “getting help regarding the use of Octopus.” citeturn33view1

### HPC and performance posture

Octopus is described as parallelised with MPI and OpenMP, scaling to “tens of thousands” of processors, and also supporting GPUs via OpenCL and CUDA. citeturn33view1

The HPC tutorial on parallelisation is unusually explicit about strategies: distributing work over k-points, Kohn–Sham states, and domain decomposition (with ghost-point communication for finite-difference stencils). citeturn33view2

### Suitability for jellium and wave‑packet coupling

Octopus is strongly suited to your *method development* phase because:

- It includes dedicated tutorials for “Jellium and jellium slabs”, with explicit input templates and physical parameterisation. citeturn33view0  
- Its “species” concept explicitly supports jellium and user-defined potentials, which is the natural handle for implementing a moving perturbation or a model projectile potential. citeturn33view3turn33view0  
- The tutorial navigation explicitly mentions adjacent model-system topics such as “Absorbing boundaries” and “e-H scattering,” which are commonly needed when you try to represent wave packets and continuum-like behaviour. citeturn33view0

Putting that together: Octopus is an excellent platform for **prototyping the physics protocol** of “wave packet ↔ jellium” even if you later move the final simulation into INQ for GPU scaling or architectural reasons. citeturn33view0turn33view3turn33view1


## Roadmap and side quests that will materially de-risk your main project

These are “high-leverage tasks” that (a) teach you the right mental model, (b) generate reusable code/scripts, and (c) reduce the probability that you spend months debugging a conceptual mistake.

### Side quest: Reproduce a documented INQ end-to-end workflow

Run one of the molecule tutorials (e.g., N₂ or benzene spectrum) and then the electronic stopping tutorial. The learning outcome is not the chemistry; it is internalising INQ’s intended control flow: define system → ground state initialisation → SCF → restart handling → real-time propagation → observable extraction. citeturn23view0turn24view2turn22view0

### Side quest: Treat the INQ stopping tutorial as your “project skeleton”

Even if you ultimately do jellium (not atomistic Al), the stopping tutorial already answers the hardest software-engineering questions you will face:

- how they structure “moving projectile” logic (velocity + ion propagator choice), citeturn24view2turn24view1  
- how they compute/record forces and energy, citeturn24view2  
- how they handle metallic occupations (extra states + temperature). citeturn24view0  

Your future wave-packet/jellium driver will likely look structurally similar even if the underlying “potential” and “host” definitions change.

### Side quest: Build a jellium slab in Octopus and validate it

Do the Octopus “jellium slab” tutorial and validate that you can reproduce the density profile including Friedel oscillations. citeturn33view0  
This gives you (i) a verified jellium setup, and (ii) a baseline grid spacing / smearing / empty-state handling pattern.

### Side quest: Implement a moving external potential in Octopus first

Because Octopus supports user-defined potentials and already has explicit jellium species, it is a fast route to test:

- whether your projectile potential needs smoothing,  
- whether your boundary conditions cause spurious reflections,  
- whether your observable (force, energy loss, induced density) behaves sensibly with timestep and grid spacing.

This is the fastest way to uncover “I thought I needed a quantum wave packet, but actually a prescribed trajectory answers my physics question” (or the opposite). citeturn33view3turn33view0

### Side quest: Decide early whether you truly need a quantum projectile

A surprisingly common outcome in electronic-stopping and related nonequilibrium simulations is that a classical (or semiclassical) projectile already gives the electronic physics you care about (screening, wake formation, energy loss, induced currents), and the quantum nature of the projectile is a second-order effect for many regimes. The review literature on first-principles electronic stopping is a good way to calibrate when a full quantum treatment is necessary. citeturn49search3

If you do need a quantum projectile, treat that as a separate research problem (and likely a publishable methods contribution). INQ’s inclusion of “exact factorization” in its theory docs suggests a plausible conceptual direction, but you will want to confirm implementation status in the code before committing. citeturn22view0turn52view0

### Side quest: Postpone the “fork” by using inq_template as your extension layer

Because INQ encourages a library usage model and provides inq_template specifically to compile your own C++ drivers, you can often implement “additional functionality” (new perturbations, new observables, even new propagation protocols) as an external layer first. citeturn23view0turn21view0

Forking becomes safer when you can point to a working external prototype and say: “this single hook is missing in the core; here is the minimal change.”

### Side quest: Maintain your own mini‑changelog tied to INQ tags

Given the concrete release milestones (v0.5.0 → v0.9 → 0.95) and the active 2026 refactors (PAW-related changes, build-system updates), keeping a personal changelog that maps “my workflow depends on X” to “X changed in commit/tag Y” will save you substantial time. citeturn51view0turn52view0