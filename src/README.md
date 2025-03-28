1. Main Script (main.m)
   |
   ├─ Calculate Global ISPC (calculateGlobalISPC)
   |  └─ Compute mean ISPC across entire dataset for normalization reference
   |
   ├─ Process each stimulation protocol
   |  ├─ Calculate ISPC (calculateISPC)
   |  ├─ Normalize ISPC values (normalizeISPC)
   |  └─ Generate visualizations
   |
   └─ Save subject-level results
      
2. Interpolation (interpolateNormalizedSubjects)
   |
   ├─ Load subject data
   |
   ├─ Interpolate to standard electrode layout
   |  └─ Use EEGLAB's pop_interp function
   |
   └─ Save interpolated data

3. Group Analysis (runGroupAnalysisWithInterpolatedData)
   |
   ├─ Load interpolated data for all subjects
   |
   ├─ Calculate group averages
   |
   └─ Generate group-level visualizations

4. Protocol-Specific Analysis (runProtocolGroupAnalysisWithInterpolation)
   |
   ├─ Process each stimulation protocol separately
   |
   ├─ Load interpolated data for all subjects
   |
   ├─ Calculate protocol-specific group averages
   |
   └─ Generate protocol-specific visualizations


Data Collection: The script processes EEG data for multiple subjects exposed to a 1Hz electrical stimulus, with subjects divided into ACTIVE (real stimulation) and SHAM (control) conditions.
ISPC Calculation: The pipeline calculates Inter-Site Phase Coherence (ISPC) between each EEG channel and the stimulus signal, measuring how synchronized brain activity is with the stimulus frequency.
Global Normalization: A global ISPC value is calculated for each subject across their entire dataset, serving as a baseline reference. All ISPC values are then normalized relative to this global value.
Z-score Standardization: Normalized ISPC values are converted to Z-scores using each subject's global standard deviation, allowing for statistical interpretation (how many standard deviations from baseline).
Segment Analysis: The recording is divided into segments (pre-stim, early-stim, late-stim, post-stim) to track changes in entrainment over time.
Comparison Metrics:

Normalized ISPC: How much entrainment occurs relative to baseline
Percent Change: Relative increase/decrease between segments
Z-scores: Statistical significance of entrainment changes


Group Analysis: Subject data is interpolated to a standard electrode layout and averaged by condition (ACTIVE vs SHAM).
Contrast Analysis: Direct subtraction of SHAM from ACTIVE data reveals treatment effects, showing where and when real stimulation produces significantly different entrainment than placebo.


1. SUBJECT-LEVEL PROCESSING
   ┌───────────────────────────────────────┐
   │            main.m                     │
   └───────────────────────┬───────────────┘
                           ▼
   ┌────────────────────────────────────────────────┐
   │ Calculate Global ISPC (calculateGlobalISPC)    │
   │ ├─ Compute mean ISPC across entire dataset     │
   │ └─ Store global reference & standard deviation │
   └──────────────────────┬─────────────────────────┘
                          ▼
   ┌────────────────────────────────────────────────┐
   │ Process Stimulation Protocols                  │
   │ ├─ Calculate ISPC (calculateISPC)              │
   │ │  └─ Measure EEG-stimulus phase coherence     │
   │ │                                              │
   │ ├─ Normalize ISPC (normalizeISPC)              │
   │ │  └─ Generate Z-scores using global std       │
   │ │                                              │
   │ └─ Visualize Results                           │
   │    ├─ Time course plots                        │
   │    └─ Topographic maps                         │
   └──────────────────────┬─────────────────────────┘
                          ▼
   ┌────────────────────────────────────────────────┐
   │ Save Subject-Level Results                     │
   └──────────────────────┬─────────────────────────┘
                          │
                          │
2. CROSS-SUBJECT STANDARDIZATION
                          │
                          ▼
   ┌────────────────────────────────────────────────┐
   │ Interpolation (interpolateNormalizedSubjects)  │
   │ ├─ Load subject data                           │
   │ ├─ Interpolate to standard electrode layout    │
   │ └─ Save interpolated data                      │
   └──────────────────────┬─────────────────────────┘
                          │
                          │
3. GROUP-LEVEL ANALYSIS
                          │
                          ▼
   ┌────────────────────────────────────────────────┐
   │ Group Analysis                                 │
   │ ├─ Load all subjects' interpolated data        │
   │ ├─ Calculate ACTIVE and SHAM group averages    │
   │ └─ Generate group-level visualizations         │
   └──────────────────────┬─────────────────────────┘
                          │
                          ▼
   ┌────────────────────────────────────────────────┐
   │ Protocol-Specific Analysis                     │
   │ ├─ Process each stimulation separately         │
   │ ├─ Calculate protocol-specific averages        │
   │ └─ Generate protocol visualizations            │
   └──────────────────────┬─────────────────────────┘
                          │
                          ▼
   ┌────────────────────────────────────────────────┐
   │ ACTIVE vs SHAM Comparison                      │
   │ ├─ Calculate condition differences             │
   │ │  ├─ Normalized (raw effect size)             │
   │ │  └─ Z-scored (statistical significance)      │
   │ │                                              │
   │ └─ Generate contrast visualizations            │
   └────────────────────────────────────────────────┘
