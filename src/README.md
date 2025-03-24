# EEG Entrainment Analysis Pipeline

This pipeline analyzes entrainment of EEG data to exogenous electrical stimulation across multiple subjects and conditions (ACTIVE vs SHAM).

## Pipeline Overview

The analysis pipeline consists of the following steps:

1. **Subject-Level Analysis** (`main.m`):
   - Process individual subject's EEG data
   - Calculate ISPC (Inter-Site Phase Coherence) for different segments
   - Create visualizations for each subject and stimulation protocol 
   - Save subject-level average data

2. **Channel Interpolation** (`interpolateSubjectData.m`):
   - Interpolate missing/bad channels in subject data
   - Create consistent electrode structures across subjects
   - Exclude bad channels based on exclude.json

3. **Group-Level Analysis** (`runGroupAnalysisWithInterp.m`):
   - Create global averages across all ACTIVE and SHAM subjects
   - Use consistent color scaling for visualizations
   - Generate topoplots with electrode labels
   - Save group average data and visualizations

4. **Protocol-Specific Analysis** (`runProtocolGroupAnalysis.m`):
   - Create averages for each stimulation protocol (Stim_1, Stim_2, etc.)
   - Use same color scaling as global analysis for consistency
   - Generate protocol-specific visualizations

## Directory Structure

- `src/` - Main source code directory
  - `assets/` - Configuration files
    - `exclude.json` - List of channels to exclude
    - `regions.json` - Region definitions
    - `subject_condition.json` - Subject condition assignments (ACTIVE/SHAM)
    - `GSN-HydroCel-256.sfp` - Electrode positions file
  - `eeglab/` - Custom EEGLAB functions
    - `eeg_topoplot.m` - Custom topoplot function
    - `pop_interp.m` - Interpolation function
    - `pop_topoplot.m` - Topoplot wrapper
  - `utils/` - Utility functions
    - `egi256_GSN_HydroCel.sfp` - Alternate electrode positions file
    - Other visualization utilities

## Configuration Files

### exclude.json
Contains a list of channels to exclude from analysis (e.g., reference electrodes, bad channels).

### subject_condition.json
Specifies which subjects received active stimulation vs. sham:
```json
{
  "subjects": [
    {"id": "101", "condition": "ACTIVE"},
    {"id": "107", "condition": "SHAM"},
    ...
  ]
}
```

## Running the Pipeline

### Prerequisites
- MATLAB with EEGLAB installed
- EEGLAB path set to `/Users/idohaber/Documents/MATLAB/eeglab2024.0/` (modify as needed)

### Option 1: Run the Complete Pipeline Automatically

The easiest way to run the full pipeline is to use the `runFullPipeline` function:

```matlab
runFullPipeline('/path/to/experiment', 'N1', 5);
```

Parameters:
- First parameter: Path to the experiment data directory
- Second parameter: Session ID (e.g., 'N1')
- Third parameter: Maximum number of stimulation protocols to process (default: 5)

This will:
1. Validate the pipeline setup
2. Run channel interpolation for all subjects
3. Perform global group analysis
4. Perform protocol-specific group analysis
5. Generate a detailed log file

### Option 2: Run Each Step Manually

If you prefer to run each step individually, follow these steps:

#### Step 1: Validate the Pipeline
Run the validation script to ensure all required dependencies are available:
```matlab
validateAnalysisPipeline();
```

#### Step 2: Run Channel Interpolation
This creates interpolated versions of all subject data files:
```matlab
interpolateSubjectData('/path/to/experiment', 'N1');
```

#### Step 3: Run Group Analysis with Interpolated Data
This creates global group averages for ACTIVE and SHAM conditions:
```matlab
runGroupAnalysisWithInterp('/path/to/experiment', 'N1');
```

#### Step 4: Run Protocol-Specific Analysis
This creates protocol-specific averages for ACTIVE and SHAM conditions:
```matlab
runProtocolGroupAnalysis('/path/to/experiment', 'N1', 5);
```
The last parameter (5) specifies the maximum number of stimulation protocols to process.

## Output Files

### Subject-Level Output
- `subject_<ID>_average.mat` - Subject average data (raw)
- `subject_<ID>_interp.mat` - Subject interpolated data

### Group-Level Output
- `group_active_average_interp.mat` - ACTIVE group average data
- `group_sham_average_interp.mat` - SHAM group average data
- `group_active_topoplots_interp.png` - ACTIVE group topoplot visualizations
- `group_sham_topoplots_interp.png` - SHAM group topoplot visualizations
- `color_scales.mat` - Color scale limits for consistent visualization

### Protocol-Specific Output
Located in subdirectories for each protocol (stim_1, stim_2, etc.):
- `group_active_average.mat` - ACTIVE group average for this protocol
- `group_sham_average.mat` - SHAM group average for this protocol
- `group_active_topoplots.png` - ACTIVE visualizations for this protocol
- `group_sham_topoplots.png` - SHAM visualizations for this protocol

## Troubleshooting

### Electrode position issues:
The scripts will look for electrode position files in the following order:
1. First in `src/assets/GSN-HydroCel-256.sfp`
2. Then in `src/utils/egi256_GSN_HydroCel.sfp` as a fallback

Make sure at least one of these files exists. The scripts print which file they're using for transparency.

### Missing subject data:
If subject data files cannot be found, check:
1. The experiment path is correct (`/Volumes/Ido/analyze` by default)
2. Subject directories follow the expected structure: `<experiment_path>/<subject_id>/<session_id>/`
3. Files named `Strength_<subject>_<session>_forSW.set` exist

### Channel exclusion:
The excluded channels are defined in `exclude.json` and are consistently applied across:
1. Subject-level interpolation
2. Global group analysis
3. Protocol-specific analysis

### EEGLAB dependencies:
Make sure EEGLAB is installed and properly initialized. The scripts need the following EEGLAB functions:
- `readlocs`: For reading electrode positions
- `topoplot`: For creating topographical plots
- `eeglab`: Main EEGLAB function (called with the 'nogui' option)

## Notes on Implementation

### Path handling:
- All file paths use `fullfile()` for cross-platform compatibility
- Most paths are derived from a base code path and experiment path

### Performance considerations:
- Calculation-intensive operations use efficient MATLAB vectorization
- Large files use MATLAB's v7.3 format for proper handling of large arrays
