# Electrode Balancing & OCP Reconstruction

This project reconstructs physically consistent **half-cell OCP curves**
for graphite (anode) and the cathode, while simultaneously **inferring
the stoichiometric limits** of both electrodes in a Li-ion full cell.

The workflow uses the experimental **full-cell OCV curve** as ground
truth, extends the cathode OCP to the full stoichiometric range, infers
the usable stoichiometric windows of both electrodes, and rebuilds a
consistent full-cell OCV model.

The final outputs are half-cell OCP curves defined over **SoL ∈
\[0,1\]**, suitable for ECM, DFN, or pack-level simulations.

------------------------------------------------------------------------

## Required Data Structure

Place the following file in:

Data/Electrode_Balancing_Data.xlsx

The workbook must contain:

### Sheet: SOC_Fullcell

-   SOC --- Full-cell state of charge\
-   OCV --- Experimental full-cell open-circuit voltage

### Sheet: Cathode_Relative

-   Relative_SoL --- Cathode OCP measured on a relative (0→1) SoL axis\
-   OCP --- Cathode open-circuit potential

### Sheet: Graphite_Literature

-   SoL --- Absolute graphite stoichiometry\
-   OCP --- Literature graphite OCP

Results are automatically written to:

Results/

All paths are built relative to `pwd` for portability.

------------------------------------------------------------------------

## Method Overview

### 1. Load Data

-   Experimental full-cell OCV vs SOC\
-   Cathode OCP vs relative SoL\
-   Literature graphite OCP

The full-cell OCV curve is treated as the reference to match.

------------------------------------------------------------------------

### 2. Convert Cathode Relative SoL → Absolute SoL

Initial bootstrap values:

y_min = 0.13\
y_max = 1.00

Mapping:

SoL_abs = y_min + (y_max − y_min) \* SoL_rel

These limits are later optimized.

------------------------------------------------------------------------

### 3. Extend Cathode OCP to SoL ∈ \[0,1\]

Because measured data usually covers only part of the window:

-   Low SoL edge → MSMR-style logarithmic fit\
-   High SoL edge → Linear fit\
-   Middle region → Interpolation

Edge model:

OCP(s) = a + b · log((1 − s)/s)

This stabilizes the extension without overfitting.

------------------------------------------------------------------------

### 4. Infer Stoichiometric Limits (Electrode Balancing)

Optimization variables:

\[x_min, x_max, y_min, y_max\]

Stoichiometry mapping:

x(SOC) = x_min + (x_max − x_min) \* SOC\
y(SOC) = y_min + (y_max − y_min) \* (1 − SOC)

Full-cell model:

OCV_model = OCP_cathode(y) − OCP_graphite(x)

A constrained least-squares solver adjusts the limits to minimize:

\|\| OCV_exp − OCV_model \|\|

------------------------------------------------------------------------

### 5. Reconstruct Graphite OCP

After balancing:

OCP_graphite_inferred = OCP_cathode(y_SOC) − OCV_exp

Then: - MSMR fits at both graphite edges\
- Interpolation inside \[x_min, x_max\]\
- Enforced strictly decreasing behavior

This ensures physical consistency.

------------------------------------------------------------------------

### 6. Validation

The script: - Recomputes full-cell OCV\
- Plots experimental vs reconstructed curves\
- Prints RMSE

Low RMSE indicates a successful balancing.

------------------------------------------------------------------------

## How to Run

1.  Ensure Excel file and sheet names match exactly.
2.  Confirm Data/ and Results/ folders exist.
3.  Run the MATLAB script.

Expected output: - Inferred limits: x_min, x_max, y_min, y_max\
- Two validation figures\
- RMSE printed in console\
- Two Excel files in Results/

------------------------------------------------------------------------

## Outputs

Results/OCP_anode_halfcell.xlsx\
- Columns: SoL, OCP_anode\
- Monotonic\
- Edge-regularized\
- SoL ∈ \[0,1\]

Results/OCP_cathode_halfcell.xlsx\
- Columns: SoL, OCP_cathode\
- Extended to full window\
- Smooth and physically consistent

------------------------------------------------------------------------

## Assumptions

-   Initial cathode limits are bootstrap values only.\
-   Graphite literature data provides shape guidance.\
-   Edge MSMR fits stabilize extrapolation.\
-   Monotonic enforcement avoids non-physical oscillations.

------------------------------------------------------------------------

## Rationale

This workflow keeps the physics simple and the numerics robust:

-   Full-cell OCV is the ground truth\
-   Edge behavior is lightly regularized\
-   Stoichiometric limits are inferred, not assumed\
-   Monotonicity ensures physical realism

The result is a stable and practical electrode balancing tool suitable
for research and industrial battery modeling.
