Hypotheses
=======

When Schrödinger's equation is discretised over Hilbert space using coherent states, and variationally over time:

1. The discretisations of Hilbert space and time do not commute.  Discretising time first is stable, discretising to coherent states first is not.

2. The unstable method will be much more sensitive to how the initial state is discretised than the stable method.

3. This is due to stiffness of the unstable discretisation.

Tasks
====
5. Julia: lazily initialized matrices.  latrix.jl.  Trap on index errors.

1. Compute the direction where σ_R increases most rapidy with respect to the change in |ψ(z)〉.  From a random initial sampling, drift in that direction at a constant rate of change of |ψ(z)〉.  See what happens.

0. Look at the conditioning of G and |Dψ〉 sampled on random grids.  Compare to conditioning on an even grid, as that grid becomes sparser, and more random points are added.

1. Look at the difference between large singular directions of the expansion and variational problems.  (Say, operator distances between low-rank approximations.)

1. Try to solve the H=0 problem, while improving the condition of Dψ

1. Investigate the conditioning of the variational problem as a phase space grid shrinks smaller than the Fock space.

2. Prove hypothesised SVE of G.  Note that polynomials in z and z* span functions of a complex variable: could show that the non-constant z* ones annihilate G.  (Check Bargmann first.)

1. Find a way to automatically set the weight of the condition number

2. Why does truncating the Fock space make things better, while ignoring the distance components along those axes makes things worse?

1. Find out how people usually handle stiffness in the method of lines.

3. Check that Methods 2 and 3 actually differ.

4. Think about the assymetry between state and time discretisations in Method 2.

2. Derive the brackets 〈N|H|N〉 and 〈N|H|ψ〉 for the Morse oscillator.

5. Derive momentum-kicked ground state of Morse oscillator

Topics to write
=====

1. Details of comparing stability by varying timestep and measuring error.

2. What can be learnt by varying energy of initial state.

3. What can be learnt by varying number of components.

4. Stability of different time discretisations.

Chores
=====

4. Install AMPL—use the [web version](http://www.ampl.com/TRYAMPL/startup.html)

