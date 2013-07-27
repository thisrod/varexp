Hypotheses
=======

When Schrödinger's equation is discretised over Hilbert space using coherent states, and variationally over time:

1. The discretisations of Hilbert space and time do not commute.  Discretising time first is stable, discretising to coherent states first is not.

2. The unstable method will be much more sensitive to how the initial state is discretised than the stable method.

3. This is due to stiffness of the unstable discretisation.

Tasks
====

1. Regularise the expansion problem using the norm of <A|A>c.

3. In the initial state problem, monitor the deviation from a straight-line trajectory, perhaps as an angle.  Also monitor the angle between the direction of motion and the straight line, and how this changes with regularisation.

0. Watch the change from the variational problem to the expansion one as the Tychonov cost of |dα|² increases.

1. Think about the continuous expansion operator on a truncated phase space.

5. Julia: lazily initialized matrices.  latrix.jl.  Trap on index errors.

0. Look at the conditioning of G and |Dψ〉 sampled on random grids.  Compare to conditioning on an even grid, as that grid becomes sparser, and more random points are added.

1. Look at the difference between large singular directions of the expansion and variational problems.  (Say, operator distances between low-rank approximations.)

2. Prove hypothesised SVE of G.  Note that polynomials in z and z* span functions of a complex variable: could show that the non-constant z* ones annihilate G.  (Check Bargmann first.)

1. Find a way to automatically set the regularisation parameters.

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

