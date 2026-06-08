# Project: Write tests, validate and restructure the codebase and claude ecosystem
<identity>

</identity>

<description>
The aim of this task is to rejuvinate the codebase such that it is verified to scientific standards and is ready to be used and extended for further tasks. Broadly speaking, in this task, we are going to write unit tests for each component of the inq-stack library, make suitable chages based on the unit tests, restructure the library if necessary and run integration tests to check if all the changes made work well in tandem. 
</description>


<workspace>
The codebase in the folder tddft/ is where the changes will be made, Specifically, a new folder containing the tests will be added to the inq-stack library, and an automatic CI / CD pipeline is written to ensure that all future changes can be easily tested. 
</workspace>


## Tasks


<task>

Structurally, in a step by step manner, commit all of the active changes in git to github. Brainstorm which files can be captured by which commits. Ensure the commits are short, succinct, and are instantly recallable. Before writing commits, using interviw style questions to understand my commit writing style (using best practices). Once the voice and tone locked in, these commits can be made and pushed. 

</task>



<task>
This is the code unit testing task. This task has many subtasks

1. The first task is to map the entire inq-stack library so that the I can look at it and understand it. I intend to use a claude plugin named understand-anything which is at /plugin marketplace add Lum1104/Understand-Anything. To pass through this subtask, I need to go through the codebase in a first pass. While doing so, I would come up ideas for code restructing or reformatting. The ideas should be stored in a specified file. These will be used in a later subtask. 

2. Using the context gained in the previous task, the individual components of the library that need unit tests are mapped. Plans for each test are given to be in brief (using any mathematical expressions or equations) required to test. A simple method of how the test is written is also made. I review all of the plan, and lock in the plans for implementation. 

3. The unit tests locked in are written. A validation agent (with a completely new context) is used to verify the writing of the unit tests. The outcome of the unit tests are reviewed. In a loop, it is first checked if the unit tests are performing as expected. If any errors produced in the unit tests cannot be attributed to the unit tests, then the codebase is at fault. The specific errors in the codebase are found and documented. 

4. I review the set of errors in the codebase. Also, I consider the stored ideas that I have listed donw from subtask 1. Then, these have to be brainstormed critically thinking about the implications. In this brainstorming session, I want to go to and fro using interview style questions to which I give my answers. The outcome of this step is a locked in plan for changes in the codebase that is to be implemented. 

5. The implementation of the plan is done. After each component is changed, the user code reviews the outputs (using the understand plugin). The outcome of this stage is a restructured and validated codebase is modular and ready to be extended for future work. 

</task>

<task>
This task is about revitalising the claude ecosystem that exists around the codebase. I belive, I have a set of rules, docs, skills and subagents that are scattered and not structured to be readily usable. I want to modularise all of this content into well made rules, hooks, skills and workflows so that these can be readily implemented by me. Here are the subtasks in this task

1. Firstly, an evaluation set is created to test the skills, rules, hooks and any other claude tools on. This evaluation for hooks can be simple, and be written by you to validate them. Similarly, for skills, the evaluation set and tasks entailed in the evaluation set are to be interviewed with the user and then implemented. 

2. Form the first subtask, a general idea of the state of the entire claude ecosystem would be gained by me. Then, I will brainstorm alongside you to using interview style questions, a plan to modularise these skills, hooks, rules and subagents. After review by the user, each decision is locked in. 

3. The locked in decision are then implemented. Then, these are rigorously tested using the evaluation set. 

4. Then, finetuning of these skills, hooks, subagents and rules is performed using iteration. 

</task>

<task>
</task>

## Principles
1. The user controls the flow of all of the tasks. In a given task, after each sub-task, the user is to understand the work that is done completely. Here,  
2. Incremental changes to the codebase is carried out. All the changes are accepted by the user. 
3. Each feature change is accomponied by a test that is run to ensure that the code is correct, the output is satisfactory, and the simulation runs. 
4. The mathematics of formulae and their implementations are also checked rigorously. In these checks, the formula used in the codebase is given to the user who manually checks if it is right. After the formula is verified, the code implementation of it also has to be verified. 
5. Changes are marked as under-review and locked. Only locked changes are implemented.

## Rules (need to be brainstormed)
### ALWAYS
1. 
### NEVER
1. Never edit a file that is not required to be changed according to the plan
2. Do not edit or make changes to the results in runs/ in ResearchProject/, Tutorial/ and QuantumKickExtension/. If required to make re-runs after finding bugs, these are to be made as new runs of their own right. 
3.  