# Analysis of the effect of the element size

- **maxsz** is the set max size for the whole phantom
- **max_el_sz** is the max set size for the electrodes

Tested parameters

- maxsz [0.0025, 0.0050, 0.0075]
- max_el_size [0.001 - 0.005] for each value of maxsz
- presence of a ball in the middle of the phantom

The effect was evaluated on the value of B at

- a random point **B** (B_position = [0.05, models(1).eit.phantom.radius + 0.05, models(1).eit.phantom.elec_vert_position+0.05]) which is away from any symmetry plane so that every B component is nonzero.
- a point **B0** between the current carrying electrodes (B_position = [0, models(1).eit.phantom.radius + 0.01, models(1).eit.phantom.elec_vert_position]) which is where it also measured on a real simulation. This was tested only with **maxsz = 0.005**.

## Instruction

Run **elem_size_analysis.m** changing maxsz and the presence of the ball (you need to comment in and out some text) and B_position (commenting).

For each variation create the figures **maxsz_y_bb_ee.fig**.

where:

- **y** is the maxsz value (0.0025 -> 0_0025)
- **bb** is either *ball* or *nothing*
- **ee** is either *var* which gives the values of the B components and their relative difference or *elem* which gives the number of elements.

For the symmetry point **B0** the results are in figure **maxsz_0_0050_ball_var_sym_point.fig**

## Results

- The largest effect is visible on the **x component of B** with a max variation of 0.2 - 0.3 %. The other components have a variation of less than 0.1 %.
- The largest effect is played by **maxsz**.
- Changing **max_el_size** once **maxsz** is set has a smaller effect than changing **maxsz**.
- Best option is to use **maxsz = 0.0050** and **max_el_sz = 0.001**. It gives a reasonable error, but still keeps the number of element to a reasonable value (~1.2e5).
- The minimal reliable value of B from the simulation is **~1e-11 T**. In fact this value appears when one or more components of **B should be zero**.

