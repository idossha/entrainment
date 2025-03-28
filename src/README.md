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
