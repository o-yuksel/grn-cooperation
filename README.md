# The plasticity that helps gene networks escape the Prisoner's Dilemma also makes them vulnerable to mutation

This repository holds the MATLAB code behind the paper of this title by Özgür Yüksel: the simulations of small gene regulatory networks that play the Prisoner's Dilemma, and the analyses that turn them into the paper's figures, tables and numbers.

Plasticity, a trait's responsiveness to its inputs, is tied to its sensitivity to mutation.
In these model gene networks each gene sums its inputs, so a mutation in the wiring acts on a gene as an incoming signal does.
In a one-round Prisoner's Dilemma, networks that receive each other's output as a signal can evolve to cooperate conditionally: they defect back against partners that defect on them.
Populations that reach cooperation can lose it to mutations that turn a network's response into defection.
In a non-social task, networks tracked a target.
With the target held still, so that responding gained nothing, we lowered responsiveness by a constant input to the output genes and a change in their self-activation, set once for every individual.
These populations lost less fitness to mutation than untouched ones.

The paper is not yet published; a link to it will be added here on publication.
All figure, table, equation and section numbers below refer to its main text and its SI Appendix.

## Before you start

Every figure and table script reads stored simulation results from a separate data archive, so nothing in the paper can be drawn or printed without it.
It is archived at Zenodo: https://doi.org/10.5281/zenodo.24504968 (about 5.9 GB).

**Run times and hardware.**
Once the archive is available, every figure script and most table scripts run in seconds or minutes; the tables below give each one's time.
From scratch, without the archive, one chain is short enough to run: the three scripts behind Table S13, about 3 h on five workers (`REPRODUCING.md`, Table S13). They write into the data folder, so create `data/publication` first.
Every other result goes back to simulations of several hours each.

## The model in brief

Each individual is a network of four genes: two outputs, cooperation (C) and defection (D), and two sensors (S_C and S_D) that receive the partner's output.
Gene j regulates gene i with weight W_ij; these 16 weights are the wiring that evolves.
Each gene's expression x_i follows (Methods, Eq. 3)

```text
dx_i/dt = f(sum_j W_ij x_j + S_i) - mu x_i,    f(u) = max(u, 0) / (K_M + max(u, 0))
```

where S_i is the input from the partner, zero for the two output genes.
When two individuals meet, their networks develop together: at every step each network's sensors receive the other's current output with gain g (Methods, Eq. 4),

```text
separate sensing:  S_C = g * out_C,  S_D = g * out_D
contrast sensing:  S_C = -S_D = g * (out_C - out_D)
```

where out_C and out_D are what the partner passes on of its x_C and x_D: as they are (graded output, `hillCoefficient` Inf in the code) or as shares of their sum (saturating output, `hillCoefficient` 1).
The two choices cross into four **signaling architectures**; `sensingMode` is `'separate'` or `'contrast'` (data files store the older `'absolute'` and `'differential'`), and saturating contrast is the standard setting.
In the evolutionary runs, a pair starts at a random low level (0 to 0.01 per gene) and develops until expression stops changing or for at most 300 steps; each individual then cooperates if its C is higher than its D.
The paper calls this **the runs' encounter rule**, because a second procedure, the membership test below, plays the same pairs from a fixed start and for longer.

**Game and evolution.** Each encounter is one round of the Prisoner's Dilemma (T = 5, R = 3, P = 1, S = 0). Every generation each individual plays one random partner; parents are drawn in proportion to payoff, and offspring may carry mutations of their wiring. A population **succeeds** if its mean cooperation at generation 500 reaches 0.95.

**Members.** A wiring is a **member** of the **conditional genotype network**, which is a set of wirings rather than a gene network, if it cooperates with its clone (an identical copy) and no single-weight mutant (a copy that differs in one weight) earns more against it than its clone does; it meets the second condition by defecting back against each mutant that defects on it.
The **membership test** checks both on one wiring, trying 10 changes to each of the 16 weights (160 mutants), each pair developed from every gene at 0.005 for 600 steps, so the test and the runs' encounter rule can disagree on some wirings (SI S4.5).
The **share of members** is the fraction of a population's wirings that pass the test.

**Responsiveness and mutation-sensitivity.** The readout is the expression of C minus that of D.
Responsiveness (alpha) is the square root of the summed squares of the settled readout's sensitivities to an input at each of the four genes; mutation-sensitivity (beta) is the same over the 16 weights.
Because a weight change enters a gene's input sum as a signal does, beta equals alpha times ||x*||, the size of the settled expression state (the congruence identity; Methods, and SI S1). `model/alpha_and_beta.m` computes both.

**The tracking task.** A second model evolves the same network to follow a target, moving or held still (Methods, "Non-social tracking task"). Its four designs differ in the output genes' self-activation, fixed and never mutated; designs 1 and 2 of SI S5.1 are the strong (0.70) and weak (0.20) ones.

### Parameters

| Parameter | Value |
| --- | --- |
| Population size N | 1,000 |
| Generations per run | 500 |
| Replicate populations per setting | 100 |
| Mutations per offspring (Poisson mean, lambda) | 0.005 |
| Share of mutations that toggle a connection (switch it off if present, on if absent, its weight kept) | 0.2 (the rest change one weight by N(0, 0.2^2)) |
| Selection power zeta | 1.0 (parents drawn with weight max(payoff, 0)^zeta, at least 1e-10) |
| Payoffs T, R, P, S | 5, 3, 1, 0 |
| Degradation rate mu | 0.10 |
| Michaelis constant K_M | 1.0 |
| Sensor gain g (Eq. 4) | 1.0 |
| Euler step; development steps | 1.0; up to 300 in the runs, 600 in the membership test |
| Starting wiring | fully connected, weights N(0, 0.1^2) |

The game runs take these values from `model/publication_config.m`, which loads `model/setting_config.m` (which in turn loads `model/base_config.m`) and then sets the paper's values; `publication_params()` returns the result as one structure.
The mutation rate is written in all three files, and the line `config.mutationRate` in `model/publication_config.m` overrides the other two, so that is the line to change; name-value pairs passed to `publication_config` override it in turn.
The three runners of the repair test (`experiments/run_repair_*.m`; the test is defined under "Terms used in the tables") and three analyses (`analysis/s3_one_mutation/share_ending_membership.m`, `analysis/s3_one_mutation/share_ending_membership_spaced_snapshots.m`, `analysis/s4_losing_cooperation/measure_defector_encounter_frequency.m`) carry their own copy of the value.
The tracking task sets its own population size, generations and mutation size in `model/tracking_config.m`.

## Getting started

### Requirements

- MATLAB R2025a.
- Statistics and Machine Learning Toolbox.
- Parallel Computing Toolbox, for the scripts that run on several workers (most simulations).

### The data

Unpack the data archive into a folder named `data` beside `startup.m`.
If you unpack it elsewhere, set the environment variable `GRN_DATA` to that folder before starting MATLAB.
Every script, except the few named under "Your own reruns" below, finds the data through `grn_data_file.m`; if a stored file is missing, it stops with a message that names the file.

### One figure

```matlab
cd path/to/the/downloaded/folder   % the folder that holds startup.m
startup                            % puts the code folders on the MATLAB path
fig5_refusal_curve                 % Fig. 5, under a minute
```

The figure is saved as `output/figures/Fig5.pdf`, `.png` and `.tif`.

### All eight figures

Each takes a few minutes or less from the stored data.
They are saved in `output/figures/` as `Fig1` to `Fig6`, `FigS1` and `FigS2`.

```matlab
startup
fig1_encounter_pair; fig2_signaling_architecture; fig3_success_vs_share_of_members;
fig4_tipping_input; fig5_refusal_curve; fig6_tracking_load;
figS1_repair_test; figS2_moving_target;
```

## Finding the script behind a result

### Folders

| Folder | Contents |
| --- | --- |
| `model/` | The model: an evolving population (`evolve_population.m`), one encounter (`simulate_pair.m`, `encounter_batch.m`), the membership test (`membership_test.m`), alpha and beta (`alpha_and_beta.m`), the refusal of mutations (`mutate_with_refusal.m`), the tracking network (`tracking.m`), the parameter files, and `stored_value.m`, which reads the older spellings of option values stored in the data files |
| `experiments/` | The simulations; every name starts with `run_` |
| `shared/` | Small functions that more than one script calls (one file each, e.g. `steady_batch.m`, `member_batch.m`) |
| `analysis/methods/` | Checks in the SI's Extended Materials and Methods |
| `analysis/s1_congruence_identity/` | SI S1, the congruence identity |
| `analysis/s2_share_of_members/` | SI S2, and the share of members in the main text |
| `analysis/s3_one_mutation/` | SI S3 and Fig. 4 |
| `analysis/s4_losing_cooperation/` | SI S4 |
| `analysis/s5_tracking_task/` | SI S5 and Fig. 6 |
| `figures/` | One script per figure (`figN_*`, `figSN_*`) and the helpers they call |

Every script begins with a header that lists what it produces (the paper item), its inputs and outputs, the call that runs it with its run time where measured, and its random seeds.
The tables below name, for each paper item, the script that prints or draws it from the stored data; `REPRODUCING.md` lists, item by item, the scripts to run before it to recompute those data from scratch.
Arguments the paper's numbers need follow the path. Times: **seconds** (under a minute), **minutes** (under an hour), **hours**, or a measured time from the header; **not timed**: the header gives none.

### Terms used in the tables

These follow the paper.

- **The 11 settings** (Fig. 3, Table S2): the four signaling architectures; populations of 50, well mixed or on rings with four neighborhood sizes; and populations of 1,000 in 20 groups of 50, where a group splits to replace another either by payoff (group selection) or at random.
- **Tipping input**: the smallest input at a gene that ends a member's cooperation, or its retaliation against its worst-case mutant (Fig. 4A).
- **Worst-case mutant**: of a member's single-weight mutants that defect on it, the one against which the member's own readout lies nearest zero.
- **Selecting members for safety**: a search that changes a member one mutation at a time, keeping a change only if the wiring stays a member and no fewer of its own random mutations leave it a member that cooperates reliably with its clone (Methods, "Mutation tests in the game").
- **Refusing** a mutation: in runs continued from the generation at which cooperation first exceeds 0.90, a mutation that would end a network's cooperation with its own clone is rejected with a set probability and drawn again (Methods, "Failure route"). A run **collapses** when its cooperation falls below 0.5.
- **Repair test**: populations about to fail are restarted with every cooperator that does not retaliate replaced by its nearest member; the **sham** moves the same individuals about as far, to other cooperators that do not retaliate (Methods, "Failure route").
- **Arm**: one treatment of an experiment, such as mutation off, or the sham of the repair test.
- **Mutation load**: the fitness a population loses to mutation, measured in the tracking task with the target held still.
- Script names use `cooperate_only` for cooperators that do not retaliate, and `response_dial` for scaling a network's response to its partner (SI S3.5).

### Main-text figures

| Item | What it shows | Script |
| --- | --- | --- |
| Fig. 1 | (A) Two networks coupled through their sensors; (B) the readout of an evolved conditional network when the other network is a copy of it and when it is an unconditional defector; (C) what the two networks of one encounter gain by defecting on each other, by kind of pair | `figures/fig1_encounter_pair.m` (not timed) |
| Fig. 2 | (A) Mean cooperation over 500 generations in the four signaling architectures, with the direction of selection on cooperation above it; (B) cooperation at generation 500 in each population, and success | `figures/fig2_signaling_architecture.m` (minutes) |
| Fig. 3 | Success against the share of members in the 11 settings | `figures/fig3_success_vs_share_of_members.m` (minutes) |
| Fig. 4 | (A) Which random mutations end a member's cooperation or retaliation, by the input each adds over the member's tipping input; (B) what selecting members for safety leaves; (C) losses as the response to the partner is scaled | `figures/fig4_tipping_input.m` (minutes) |
| Fig. 5 | Collapses as a growing share of the mutations that end cooperation is refused. Each of the 100 saturating-contrast populations was continued once; populations 1 to 50 were then continued a second time with new seeds, so each arm counts 150 runs | `figures/fig5_refusal_curve.m` (seconds) |
| Fig. 6 | (A) Measured against predicted per-mutation load in the tracking task; (B) the cut in mutation load when responsiveness is lowered, against mutation size | `figures/fig6_tracking_load.m` (minutes) |

### SI Appendix, in the order of its sections

| Item | What it shows | Script |
| --- | --- | --- |
| Table S1 | Standard parameters | none: the values are in the four files of `model/` named under "Parameters" |
| Extended Materials and Methods, "Developmental integration and convergence" | Success against the number of development steps | `analysis/methods/development_steps_summary.m` (not timed) |
| Extended Materials and Methods, "Spatial reproduction" | Assortment on the rings | `analysis/methods/ring_assortment.m` (not timed) |
| S1.3 | The identity checked on members for finite mutations | `analysis/s1_congruence_identity/test_finite_step_identity.m`, with no arguments and with ([], [], 7700) (not timed) |
| Table S2 (S2.3, S2.7) | Share of members, success and the correlation in the 11 settings; the shares when the test also tries toggling connections | `analysis/s2_share_of_members/share_of_members_vs_success.m` (about 1 min); the toggles: `analysis/s2_share_of_members/scope_checks_summary.m` (about 2 min) |
| Table S16 (S2.3) | Whether populations reach cooperation from the starting wirings alone (mutation off) or from mutation alone (no starting response: the weights from the sensors to the outputs start at zero), against the standard runs | `experiments/run_member_source.m` ('summary') (not timed); it writes the table to `publication/member_source/summary_member_source.txt` in the data folder |
| S2.3, other numbers | The two conditions apart: `analysis/s2_share_of_members/cooperation_condition_shares.m` (about 10 min). The cooperation and retaliation margins (how far the readout lies from zero against the clone, and against the worst-case mutant): `analysis/s2_share_of_members/two_margin_plane.m` (fullfile(grn_data_file('publication'), 'retaliation_margin_vs_rho.mat'), struct()) (about 1 min). Members whose clone state is silent (both outputs below 0.001) or still changing when the test stops: `analysis/s2_share_of_members/silent_and_moving_members.m` (about 4 min on 3 workers). The interval on the share of members: `analysis/s2_share_of_members/share_interval_effective_n.m` (about 1 min) | as named |
| Table S3 (S2.3, S2.5, S3.2) | The four signaling architectures. Failed populations and the generation-1 share of members: `analysis/s2_share_of_members/tableS3_failed_and_gen1_rows.m` (about 30 s). Member pairs joined: `analysis/s2_share_of_members/member_connectivity.m` (about 5 min on 3 workers). Mutations ending membership: `analysis/s3_one_mutation/share_ending_membership.m` (struct('nReps', 100)) (not timed). The slope of the share of members after generation 300: `analysis/s3_one_mutation/share_of_members_drift.m` (about 15 s) | as named |
| Table S4 (S2.4) | When selection turns in cooperation's favor, and against it again as populations lose cooperation | `analysis/s2_share_of_members/selection_inherited.m` (not timed); counts per population: `analysis/s2_share_of_members/selection_around_collapse.m` (not timed) |
| S2.4, other numbers | How bimodal final cooperation is; when the ranking of the architectures settles | `analysis/s2_share_of_members/bimodality_by_architecture.m`, `analysis/s2_share_of_members/ranking_settle_generation.m` (not timed) |
| Table S5 (S2.6) | Success under one-parameter sweeps | `analysis/s2_share_of_members/tableS5_parameter_sweeps.m` (not timed); it writes the table to `publication/table_s3_sweeps.txt` in the data folder |
| S2.7 | Group selection against random splitting and undivided populations | `analysis/s2_share_of_members/group_selection_arms.m` (not timed) |
| Table S6 (S3.2) | Share of mutations that end membership, cooperation or retaliation; the main text's "more than a third end membership" | `analysis/s3_one_mutation/tableS6_per_mutation_shares.m` (about 20 s) |
| S3.2, other numbers | Per architecture, the share of mutations that end membership (0.591 under graded separate to 0.381 under saturating contrast): `analysis/s3_one_mutation/share_ending_membership.m` (struct('nReps', 100)) (not timed); re-measured on snapshots at least ten generations apart (0.340 under saturating contrast): `analysis/s3_one_mutation/share_ending_membership_spaced_snapshots.m` (about 1 h) | as named |
| S3.3 | The worked example; tipping inputs; mutations that add almost no input; cooperators that do not retaliate | `analysis/s3_one_mutation/worked_example.m` (seconds); `analysis/s3_one_mutation/tipping_input_flag.m` ('analyze', 7000) (not timed); `analysis/s3_one_mutation/tipping_flag_vs_input_size.m` (seconds); `analysis/s3_one_mutation/measure_silent_source.m` (about 1 min); `analysis/s3_one_mutation/cooperate_only_losses.m` ('full') (seconds) |
| Table S7 (S3.4) | Selecting members for safety | `analysis/s3_one_mutation/safety_search_cooperate_only.m` (seconds); `analysis/s3_one_mutation/safety_search_losses.m` ('long') (4.4 min) |
| S3.4, other numbers | Whether the searched wirings retaliate when started at their resting state with their clone; the search without its two rules | `analysis/s3_one_mutation/fig4b_retaliation_from_rest.m` (seconds); `experiments/run_safety_search_without_rules.m` ('report') (calls in its header) |
| Table S8 (S3.5) | Losses as the response to the partner is halved or scaled up | `analysis/s3_one_mutation/dial_loss_shares.m` (seconds): the first three rows; `analysis/s3_one_mutation/response_dial_recount.m` (seconds); `analysis/s3_one_mutation/response_dial_from_rest.m` (seconds): the last row |
| S3.5, other numbers | The scaled response treated as an input | `analysis/s3_one_mutation/response_dial_as_input.m` (seconds) |
| Table S9 (S4.1, S4.2) | How often members meet a defector, selection on retaliation, and the census before failure | `analysis/s4_losing_cooperation/measure_defector_encounter_frequency.m` (about 20 min); `analysis/s4_losing_cooperation/membership_exits_split.m` (10 to 15 min); `analysis/s4_losing_cooperation/census_before_failure.m` (struct('reportOnly', true)) (about 1 min) |
| Table S10 (S4.3), Fig. S1 | Collapse in the repair test, by arm | `analysis/s4_losing_cooperation/repair_test_summary.m` ('main') (about 1 min); `analysis/s4_losing_cooperation/repair_test_without_defector.m` (about 1 min); Fig. S1: `figures/figS1_repair_test.m` (about 1 min) |
| S4.3, other numbers | Failing populations with no member at the restart point; checks of the repair against the sham | `analysis/s4_losing_cooperation/repair_test_observed_arm.m` ('main') (about 1 min); `analysis/s4_losing_cooperation/repair_vs_sham_checks.m` (about 2.3 min) |
| Table S11 (S4.4) | What a member's defecting mutants earn | `analysis/s4_losing_cooperation/tableS11_mutant_encounters.m` ('full') (not timed) |
| S4.4, other numbers | How members answer their own mutants | `analysis/s4_losing_cooperation/defect_back_rates.m` (about 1 min); `analysis/s4_losing_cooperation/measure_joint_loss.m` (about 2 min; a supporting analysis, no number of the paper comes from it); `analysis/s4_losing_cooperation/measure_mutation_payoff_spectrum.m` (not timed) |
| S4.5 | The membership test against the runs' encounter rule | `analysis/s4_losing_cooperation/membership_test_vs_encounter_rule.m` (struct('stage', 'report')) (about 8 min) |
| Table S12 (S4.6) | Collapses by share of mutations refused, for the first 100 runs, the further 50 and the pooled 150 | `experiments/run_partial_refusal.m` ('report_batch2') (not timed) |
| Table S13 (S5) | Tracking task: shift and load per mutation, target moving and held still | `experiments/run_tracking_task.m` (struct('nPop', 1000)): the shift; `analysis/s5_tracking_task/tracking_responsiveness.m` (struct('nPop', 1000)): the size of x* and the share of responsiveness carried by the input genes; `analysis/s5_tracking_task/tracking_load_vs_identity.m` (struct('nPop', 1000)): measured and predicted load (about 3 h together) |
| S5, other numbers | The excess of the measured load over the prediction | `analysis/s5_tracking_task/tracking_load_excess.m` (struct('nDraws', 6000)) (about 5 min) |
| Table S14 (S5.1) | Tracking task: population tests of the load | `analysis/s5_tracking_task/tracking_cut_after_150.m` (seconds): the cut with the target held still; the other rows are printed by the runs that make them (`REPRODUCING.md`) |
| Table S15 (S5.1) | Tracking task: the cut in mutation load against mutation size | `experiments/run_tracking_cut_vs_mutation_size.m` ('final') (seconds); `experiments/run_tracking_cut_vs_shift.m` ('final') (seconds) |
| Fig. S2 (S5.1) | Tracking task with the target moving: populations whose responsiveness was lowered partly restore it | `figures/figS2_moving_target.m` (minutes) |
| S5.1, other numbers | Checks of the load with the target held still | `analysis/s5_tracking_task/tracking_held_still_checks.m`, `analysis/s5_tracking_task/tracking_load_exponents.m`, `analysis/s5_tracking_task/tracking_load_rank_correlation.m` (seconds each) |

S1.1, S1.2, S1.4, S1.5, S2.1, S2.2, S2.8, S2.9, S3.1 and S3.6 hold derivations, settings and background, not simulated results.

## Running the slow parts

- **Workers.** Most simulations run on a pool of parallel workers. Most scripts open a pool only when none is open, so to limit the number of workers, open one first, for example `parpool('Processes', 2)`; the `response_dial*.m` scripts and the 'np_*' phases of `conditional_cooperation_measures.m` close an open pool of another size and open their own. Some scripts take the worker count as an argument (see the header).
- **Stored reports.** Some scripts stop rather than overwrite a report already in the data folder; the header says which. To rerun one on the archive, rename the stored report first.
- **Your own reruns.** To keep them apart from the archive, the headers of a few scripts name the environment variables `GRN_RESULTS` and `GRN_PUBDIR`; leave both unset to read the archive.

## Reproducibility

- Seeds are set in the code, and each header names them.
- In the evolution runs, each replicate has its own seed: replicate r runs under `rng(base + r, 'twister')`. The base is 12345 in `run_signaling_architectures.m` and is stated in the header of every other run. A replicate's result therefore does not depend on how many replicates or workers a run uses.
- Scripts that take a structure of settings can run fewer replicates: `run_tracking_task(struct('nPop', 1000, 'nReps', 2))` runs replicates 1 and 2 of each design and target, as in the full run. It writes to the file the full run writes, so use a data folder of your own.
- `membership_test.m` seeds each wiring itself and restores the caller's random stream, so its verdict depends only on the wiring.
- A few scripts seed MATLAB's default generator without naming it (their headers say so). A pool worker's default generator differs from the desktop session's, so run these exactly as the header's Run line does.

## License

MIT License; see `LICENSE`.

## How to cite

Please cite the paper:

Yüksel, Ö. The plasticity that helps gene networks escape the Prisoner's Dilemma also makes them vulnerable to mutation. Journal details will be added on publication.

Code: this repository, version 1.0.0, archived at https://doi.org/10.5281/zenodo.24505011.
Data: the data archive, version 1.0.0, https://doi.org/10.5281/zenodo.24504968.
`CITATION.cff` holds the same information in machine-readable form.

## Contact

Özgür Yüksel, <ozguryuksel@live.com>, ORCID [0009-0005-7431-8075](https://orcid.org/0009-0005-7431-8075)
