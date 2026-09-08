# Steps that need a browser (cannot be done from the command line)

1. **Zenodo integration** — **DONE 2026-09-08.** The repository was made
   public, the Zenodo GitHub switch for `dmai287/csa-audit` was flipped, and
   release `v0.1.0` was cut, which minted the DOIs below. The concept DOI is
   recorded in `CITATION.cff` and in the manuscript's data availability
   statement.

   - Concept DOI (always the latest version): `10.5281/zenodo.22648224`
   - Version DOI for `v0.1.0`: `10.5281/zenodo.22648225`

   Every later release is archived automatically; cut a fresh tagged release
   before submission so the archived code matches the reported analysis.

2. **fastMRI data-use agreement.** Request access at
   https://fastmri.med.nyu.edu/ (brain multi-coil). The download links arrive
   by e-mail and expire; store the data outside this repository (`data/` is
   git-ignored). Status 2026-09-08: **done for this paper's needs.** val 0-2,
   test_full 0-2, train batch 0 and test 0-2 (514 GB, 10 files) downloaded and
   verified on `/Volumes/T9/fastMRI_brain/archives/` (see its README.txt).
   DICOM and train batches 1-9 were skipped (T9 has no more free space and
   Paper 2's protocol does not need them: no model is retrained here, and
   DICOM is vendor-reconstructed, not raw k-space). To fetch them anyway,
   attach a larger drive, copy the folder there, and re-run
   `download_fastmri.sh` before 2026-12-07 (link expiry).
   fastMRI+ annotations are public at
   https://github.com/microsoft/fastmri-plus (CSV files).

3. **OSF pre-registration.** Create the registration at https://osf.io from
   the analysis plan (paper Section 4.8) and `configs/thresholds.yaml` once the
   thresholds are chosen; save the PDF and its SHA-256 under `prereg/`.

4. **Ethics** for the optional reader spot-check, through the RMIT HREC portal.
