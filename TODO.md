Hypotheses
=======

When Schrödinger's equation is discretised over Hilbert space using coherent states, and variationally over time:

1. The discretisations of Hilbert space and time do not commute.  Discretising time first is stable, discretising to coherent states first is not.

2. The unstable method will be much more sensitive to how the initial state is discretised than the stable method.

3. This is due to stiffness of the unstable discretisation.

Tasks
====
1. Find a way to automatically set the weight of the condition number

2. Test large superpositions propagated step to step, not randomly drawn

2. Test how the condition number of the jacobian varies for states near those fitted to the exact oscillator states.

2. Test the Picard condition for H|ψ(z)〉=|Dψ(z)〉dz, for states fitted to the exact oscillator states.

1. Average particle number of singular vectors vs singular value

2. Why does truncating the Fock space make things better, while ignoring the distance components along those axes makes things worse?

0. Track the states and condition number gradients in follow.

0. Play with optimisation options to speed it up, or give up earlier.

2. Vary the singular-value truncation threshold

2. Test if chasing the exact solution converges when more components are added.

3. Tidy up dry.py

0. Derive details of backward Euler and semi-implicit method 2.

1. Find out how people usually handle stiffness in the method of lines.

6. Experiment with speed of sparse diagonal matrices

6. Experiment with speed of generators vs. matrix multiplication

3. Check that Methods 2 and 3 actually differ.

4. Think about the assymetry between state and time discretisations in Method 2.

2. Derive the brackets 〈N|H|N〉 and 〈N|H|ψ〉 for the Morse oscillator.

3. Choose time discretisations to compare.

4. Find sensible discretisations of the initial state for the quartic oscillator using AMPL.

5. Derive momentum-kicked ground state of Mores oscillator

4. Find sensible discretisations for the Morse oscillator using AMPL.

Topics to write
=====

1. Details of comparing stability by varying timestep and measuring error.

2. What can be learnt by varying energy of initial state.

3. What can be learnt by varying number of components.

4. Stability of different time discretisations.

Chores
=====

4. Install AMPL—use the [web version](http://www.ampl.com/TRYAMPL/startup.html)

