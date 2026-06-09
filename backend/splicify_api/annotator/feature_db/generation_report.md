# Feature DB Generation Report

- Generated: 2026-03-20T12:41:22.483557+00:00
- Score: **-28**
- Grade: **failed**
- Release readiness: **not_ready**
- QC overall: **FAIL**

## Coverage

- Total Canonical Features: 19697
- Included: 3092
- Included Categories: 21
- Seed Features Detected: 3
- Snapgene Kb Records: 160
- Snapgene Kb Mapped: 156
- Snapgene Kb Curated: 17

## Clustering

- Total Clusters: 38390
- High Cv Clusters: 0
- Multimodal Clusters: 0
- Boundary Drift Features: 126

## Output Integrity

- Anchor Fasta Records: 405
- Anchor Meta Rows: 405
- Motif Fasta Records: 68
- Motif Meta Rows: 68
- Anchor Fasta Meta Match: True
- Motif Fasta Meta Match: True

## QA

- Check status counts: {'PASS': 3, 'SKIP': 1, 'WARN': 2, 'FAIL': 1}

## Biosecurity

- Enabled: False
- Source Plasmids Screened: 0
- Source Plasmids Blocked: 0
- Total Flagged Excluded: 0
- Errors: 0
- Uncertain: 0

## Recommendations

- Address FAIL checks in qc_report.json before release.
- Review boundary drift features and manually curate problematic boundaries.
- Enable biosecurity screening for production releases.
