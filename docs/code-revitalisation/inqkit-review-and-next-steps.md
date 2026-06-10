# Inqstack and review

Here are the observations I've made after completing a review of the inq-stack's 
inqkit library.  I am yet to complete a review of the Python library - inqview. 



## Where to find my reviews, comments & criticisms
I have left all of the reviews, comments and criticisms starting with TODO: in these files. In some cases, there are multi line comments formatted as /* TODO: */, and
in the other cases // TODO: . However, in most cases, these comments are multi-lined. I want you to find all of these comments from the inqkit part of the project. We are going to use these comments I've recoreded to move on to the next part, by grilling (interviewing) me, we are going to come up with specific tests to make for each case, and a proposed new architecture. We are going to implement these changes incrementally. Also, before we start changing the code, we need to be able to compare results from before and after the re-factor to ensure that the performance has not changed. We will have these tests too. For such tests, it is required to store data before moving to the next phase. So it must be kept in mind. Again, as instructed in the plan, we are going to start writing these tests, and running them, and modifying only if the tests fail (without an explanation, and confirmed no issues in the tests)


## Some overall comments
1. I think, there are minor issues in naming variables. It is fine to use x, y, z, and if field or vectors involed, vx, vy, vz etc. However, when a better name can be used, that is better. For example, norm could have been used for the integration of the density field. 

2. Compartmentalise the code better. For example, in the case of the wavepacket.hpp file, the Gram-Schmidt function could have been extracted out. This would have helped test it better (and writing unit tests for it). 

3. Utility functions have been named using the "function_", where the "_" represents that this is a utility function. We need to come up with some guidelines for coding (which should be turned into a skill, and rules, and hooks if required). We need to maintain a high standard of code up to scientific standards, and ensure it does the job well. 

4. Some code could be moved around to new files, to ensure better structure and modularisation. For example, the fft_shift and are defined in one grid_layout file, but have been used in multiple files. We can streamline this, and have a commonly re-used code in a utils file of some kind. 

5. We also need to standardise, for the user using the inqkit library, if we use the fft index pattern for the density, orbtial and other arrays we have or the usual (-L/2 is 0 index) convention. After fixing on the convention, we should then go ahead and standardise it throughout the codebase. 

6. Also, in this aim to standardise the outputs, we should have well thoughout out defaults that are outputted by the tddft runs. So, the default observable set, the post processing done etc. This is currently present in a skill, but, we can change the defaults to ensure that the is done automatically, and deterministically. 

7. Utility classes such as a vector that can have three components can be used to store vector quantities. Now, we have been storing these quantities in current_x, current_y, current_z etc. This is not clean. This can be cleaned up where appropriate. Need to carefully think such data structure, that could be use widely enough to make the code much more readable and easy to follow. 


8. Some other ideas, are brainstorming the potential analysis that can be done using the L2 signal using densities, and delta densities. Also, extending these to include total, system and wavepacket l2 for density and delta density. Recording this data and using it can be useful. 

9. Also, the coarse grained delta density is at a default value. I have not tried tuning it and playing around with it. This can be a powerful tool. So, for a given run, try using different values, comparing them, and come up with a default value. (we brainstorm and decide the run) 

10. Orthogonalisation and its impact has to be studied carefully. Here, it is necesary to compare the momentum space and real space implementations of the orthogonalisation algorithm. 

11. The GPU code, and the MPI code has to be ensured to be understood well, and implemented well in the codebase. This can be done using an independent task, that can be verified. 

12.  Mapping of coordintes of the output array to the correct directions is to be written down rigorously, perhaps grid_layout.hpp is the culprit or perhaps the pytho postporcessing. For the potential test for grid layout, run a simple simulation with an ion at  certain position (not center, at an arbitrary point). then, use the grid_layout.hpp module to get the 3D vector of the ion positions. Then take xz, xy and yz slices, and chec the result if it agress with the initially specified position. This is a test to see i the parsing of the 3D vector is right (in that the FFT shift that is taking place map to the right positions). I believe VTI is the format the post processing sripts use to make these density in 2D slices visualisation. Use that and check if the convention o reading these coordinates is tight in Python. Check the VTI output too.

In all, we need to come up with one-off tests, that answer my questions and test my hypotheses. The other tests would be unit and integration tests which are going to be implemented in the CI / CD pipeline. Both these ideas can be found in the comments in the files. 

## Some general principle for writing tests

### Unit tests
1.  One logical assertion per test. Each test should have a single, clear reason to fail. Multiple assertions are fine if they all verify the same concept, but split tests that check unrelated things. If multiple assertions are to be made for a single file or function, then they can be clubbed together in the same file. 
2.  Use descriptive names should_return_empty_list_when_user_has_no_orders() beats test_user_3(). A good name is a specification — it explains the scenario and expected outcome.
3. Test behavior, not implementation Tests should verify what code does, not how it does it. Avoid asserting on internal state or private methods; if you refactor and tests break without changing behavior, the tests are too brittle.
4. Follow the AAA Patter. Arrange (set up), Act (execute), Assert (verify) — keep each phase clearly separated and each test focused on one behavior.

### Integration tests
1. Be selective and targeted. Integration tests are slower and more expensive. Focus on critical paths and high-risk boundaries (auth flows, payment processing, DB transactions) rather than trying to cover every edge case.
2. Manage test data carefully. Use a dedicated test database or schema. Seed data before each test and clean up after, so tests don't bleed into each other. Tools like database transactions rolled back after each test work well here.
3. Keep them deterministic. Flaky integration tests are worse than no tests — they erode trust. Avoid depending on timing, external live services, or test-order assumptions.


## Materials required for the Python inqkit code-review
The python part of the library mainly has post processing scripts. However, 
there are important formulae that were used to calulate for example the loss function, the overlap between the occupation analysing of the evolved KS orbitals using ground state KS orbitals etc. I need to have a markdown table made, clearly documenting a metric / formula implemented in a certain file. I am going to review
these files and specifically these implementations to be sure that the code is doing
what I want it to do. For the review itself, i will be using understand everything plugin to create a website that I can browse easily and logically. So, these instructions have to be fed into understand-anything (dashboard). 