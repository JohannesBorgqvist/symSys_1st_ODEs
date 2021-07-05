# TASK LIST FOR THE SYM_SYS PROJECT
- [x] Implement a linear algebra approach in the solution algorithm
- [x] Add Loetka Volterra,
- [x] Add all models listed below under the section "Biological models"
- [x] Run symmetry calculations on all models under the section "Biological models"
- [x] Make sure code work for the general linear case,
- [x] Make the code work for the DBH model,
- [x] Add branches github,
- [ ] Fix sources for the Introduction,
- [ ] Fix sources for the Discussion,
- [ ] Shorten the code for the algorithm,
- [ ] Create a notebook for Hydon's model in the Github,
- [ ] Create notebooks for linear + DBH models
- [ ] Calculate the most general ODE system from the Lie algebra of the DBH model,
- [ ] Make sure that the code can run on Windows,
- [ ] Make sure that the code can run on Mac,
- [x] Find more system of ODEs with known symmetries from mathematical biology.
- [x] Fix error of "from sympy import *" where we need to import each part of sympy individually,
- [x] Fix so that the generator is saved with align instead of as equation,
- [x] Decide where to publish (see the second section below),
- [x] This is how you check a box! Just and an x in the box.
- [x] Calculate symmetries for LV-model
- [ ] Add charateristic to theory section
- [ ] Split defs in theory section
- [ ] Name parts of determining equations (currently "pseudo")
- [ ] Name type of anzats (currently "projective")
- [ ] Change name of algebraic equations from solution of ODEs to "constraints"
- [ ] Add block matrix derivation before algorithm
- [ ] Re-write algorithm section; polynomials, notation etc + structure
- [ ] Simpler model worked as detailed example
- [ ] Move full calc of Hydon's model to supplementary
- [ ] Include partial example Hydon's model in example section
- [ ] Derive DBH model from symmetries
- [ ] Complete derivation of Hydon's to standard form
- [ ] Change to tau=0 everywhere in the DBH model
- [x] Compute symmetries for oscillator models (Brusselator + Oregonator)
- [x] Compute SIR AIDS model (Murray p. 338)
- [ ] Check calculations for Hydon's model derivation
- [ ] Derive linear model from symmetries
- [x] Understand method of deriving ODEs from symmetries and invariant
- [x] Parallelise code using MPIs (each model takes up one core each, so we can compute the symmetries using sympy in parallel)
- [ ] Access the computers at the office using ssh, and set up the github project there
- [ ] Launch the simulations on the work computers 
- [ ] Transferring code to Julia: Read in models using Julia



## Branches:
REALLY GOOD LINK FOR WORKING WITH BRANCHES IN GIT:
https://thenewstack.io/dont-mess-with-the-master-working-with-branches-in-git-and-github/

One branch for each participant in the project:
1. Johannes's branch,
2. Fredrik's branch,
3. Ruth's branch.
  

## Potential journals where we will publish
Some potential candidate journals:
1. https://collections.plos.org/collection/compbiol-methods/ 
2. https://journals.plos.org/ploscompbiol/s/other-article-types 
3. https://www.springer.com/journal/11538/submission-guidelines#Instructions%20to%20Authors_Article%20Types 
4. https://royalsocietypublishing.org/action/doSearch?SeriesKey=rsif&sortBy=cited 


## Biological models
The following models have been added to the input folder:

1. AIDS_epidemic.xlsx
2. BZ_model.xlsx
3. Goodwin_1.xlsx
4. hydons_model.xlsx
5. Lotka_Volterra.xlsx
6. Lotka_Volterra_realistic.xlsx
7. DBH_model.xlsx
8. Brusselator.xlsx 
9. Goodwin_2.xlsx
10. Lactose_operon.xlsx
11. Lorenz_model.xlsx
12. SIR.xlsx 

A large run of symmetry calculation was started Thursday the 21st of June, and we will see how long it takes for these calculations to be finished. 


## List of computers at Oxford

1. https://www.maths.ox.ac.uk/members/it/remote-access 
2. https://www.maths.ox.ac.uk/members/it/machines 
3. https://www.maths.ox.ac.uk/members/it/faqs/connection


## Transferring code to Julia
As I (Johannes) launched all simulations on Midsummers eve (25/6-2021) it turned out that the simulations got stuck on the BZ model and it was not finished after three days (on Monday that is). So I had heard that symbolic calculations are slow before, but I had really not appreciated it until now. So what I realised is that we need faster symbolic solvers, and luckily that are people who have developed a faster symbolic solver in Julia called [Symbolic.jl](https://symbolics.juliasymbolics.org/dev/tutorials/symbolic_functions/ ). In fact, I read a discussion on a forum, and it turned out that many Python programmers left Python a few years back and started using Julia for performance. In particular, they discussed the difference between python and the Julia package, and sympy was written for solving small symbolic problems
conveniently but it is very slow as it is written purely in Python. However, Symbolic.jl is written with the sole purpose of performance, and it converts the code to C meaning that it can use powerful packages such as BLAS and LAPACK in order to solve matrix equations efficiently. So, after a lot of thinking, I think we need to convert our code to Julia.

In particular, there is one step of the algorithm which I think constitutes the bottleneck. This step involves computing the matrix exponential of the matrix B, and here there is actually *another* Julia implementation called [ExponentialUtilities.jl](https://github.com/SciML/ExponentialUtilities.jl ) and it is specialised at calculating (among other things) the matrix exponential. I am not sure that this latter package can be combined with the *Symbolic.jl* framework but what gives me hope on this front is that one of the developers of *ExponentialUtilities.jl* was also involved in developing *Symbolic.jl* so hopefully they can be combined. 

I am confident that we can implement our current Python code in Julia entirely. This is on account of the fact that Sympy exists in Julia as well (although it is still written in Python and therefore slow). However, we would like to avoid to use sympy as much as possible as we would like to ramp up the performance as much as possible. 
