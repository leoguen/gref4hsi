# Otter / Stavøya RGB mosaic: assumptions, findings and open points

Status of the COLMAP RGB orthomosaic pilot on the Otter USV transect
`uhi_20241024_140601` (24 Oct 2024, Stavøya). Data, config, logs and results live in
`<mission>/gref4hsi_test/colmap/` on NESP_1; the scripts are in this folder. This
file records what was assumed along the way so that the results are not mistaken for
more than they are.

## Dataset

| Item | Value |
|---|---|
| RGB frames | 217, 648 x 486 px, 3 fps, 72 s, from three H5 chunks |
| Track length | ~33 m, vehicle speed ~0.45 m/s |
| Altitude above seabed | 1.17 to 1.29 m (altimeter) |
| Navigation | ECEF positions at 25 Hz from gref4hsi `Intermediate/pose.csv`, vehicle fixed at ellipsoidal height 0 |
| Seabed surface | gref4hsi altimeter DEM, 5 cm, ellipsoidal height about -1.21 m |
| Vehicle motion | heading oscillates between ~65 and ~130 deg with a period of ~3.3 s (fishtailing), pitch/roll within ~2 deg |

## Assumptions made

1. **Camera intrinsics are fixed during SfM.** COLMAP's self-calibration diverges on
   the flat seabed under pure forward motion (focal lengths from 43 to 4300 px across
   fragments). The intrinsics used are the ones COLMAP refined on its largest
   free-intrinsics fragment (f ~487 px), with the principal point at the image centre.
   The RGB camera has not been calibrated for this deployment. The in-water
   calibration in Løvås et al. (2022, Table II) is for the UHI-4 unit of that paper
   and gave a worse fit here.

2. **Projection focal length is 605 px, not 487.** Chosen because the overlap
   consistency between neighbouring frames is best there, and because an independent
   estimate from optical flow against navigation gives ~605 px. Mapping with 605 px
   fixed inside COLMAP was worse (11 frames dropped), so 605 is applied only when
   projecting. The two numbers disagreeing is a symptom of the missing calibration.

3. **The SfM model is used only locally.** Over the whole chain the model drifts
   17-fold in scale and the optical axis rotates by more than 90 deg (no loop closure,
   forward motion, planar scene). Every image is therefore georeferenced with its own
   similarity transform fitted on its +-6 neighbours. Relative orientation between
   neighbours comes from the images; absolute position comes from navigation.

4. **The camera looks straight down.** The navigation has no usable height and the
   track is nearly straight, so camera centres alone leave the roll about the track
   undetermined. A near-nadir constraint resolves it. Logged pitch/roll are within
   2 deg, and the fitted tilts come out at 3 to 4 deg median, so the assumption is
   consistent with the data but it does replace a measured boresight.

5. **The camera is offset from the navigation point by a fixed lever arm and a
   clock offset.** Estimated from the data: camera = navigation(t + 1.5 s) shifted
   0.69 m backwards along the image-up direction and 0.02 m to the right. The
   estimate halves the misfit between image-derived and navigation-derived camera
   paths (12.8 to 6.4 cm RMS) and improves overlap consistency (4.49 to 4.05 DN).
   Caveats:
   * The lever arm and the clock offset trade off against each other, and the yaw
     is periodic, so near-equivalent solutions exist every ~3.3 s of clock offset
     (e.g. -1.5 s with the camera 0.46 m ahead). +1.5 s is the deepest minimum, but
     the choice is not unique.
   * The fit over-straightens. The images alone see ~10 cm of sideways wobble per
     13-frame window, the antenna track has ~16 cm, the corrected track ~4 cm. The
     truth is between the corrected and antenna tracks.
   * gref4hsi uses a zero lever arm for the HSI (`HSI_2b.xml`: tx = ty = tz = 0), so
     the hyperspectral products are placed on the same navigation point. If the lever
     arm is real, they inherit the same offset. The RGB camera is ~5 cm ahead of the
     HSI along image-up (UHI housing geometry, cf. Løvås et al. 2022), which is
     negligible against the 0.7 m.

6. **The seabed is the altimeter DEM.** Relief within a frame (rocks, kelp) is not
   modelled; a 10 cm height error at the frame edge displaces that pixel by about
   5 cm on the ground.

7. **Image preprocessing.** COLMAP sees grey-world normalised frames (per-pixel
   statistics over the transect) with the housing rim masked; the mosaic uses the raw
   frames. On this transect the normalisation quadrupled the keypoints and more than
   doubled the inliers but changed neither the poses nor the estimated offsets (forward
   -0.688 vs -0.694 m, clock 1.501 vs 1.505 s), a useful reproducibility check.
   **One global flat-field and one global colour stretch.** The flat field is the
   heavily smoothed mean of all frames; the stretch is common to all channels so
   colour ratios are preserved. No water-column correction.

8. **Nadir blending.** Per cell the frame whose centre is nearest wins. Seams remain
   visible where frames disagree by a few centimetres.

## What the results are good for

* Visual interpretation of the seabed at ~2 mm (camera pixel ~2 mm on the ground).
* Checking the navigation: the track overlays in `output/tracks/` show antenna,
  corrected camera and image-derived camera positions.
* Absolute position is correct to the lever-arm uncertainty (decimetres along the
  track), relative position within the strip to a few centimetres.

## Possible improvements, roughly in order of payoff

1. **Measure the lever arm from the GNSS antenna to the UHI on the Otter** and fix it
   in the alignment (only the clock offset is then estimated). This removes the
   trade-off and settles the absolute position along the track. Apply the same
   lever arm to the HSI in gref4hsi.
2. **Calibrate the RGB camera in water** (checkerboard, as in Løvås et al. 2022). A
   known focal length and distortion removes the 487 vs 605 px ambiguity and most of
   the remaining seam offsets. Alternatively, run the calibration on a transect with
   loop closure so COLMAP can self-calibrate reliably.
3. **Find the clock offset independently.** Compare the heading from the images
   with the heading from the navigation at a sharp turn and read off the lag, or log
   a synchronisation event. 1.5 s is large for a logger offset and may be latency
   inside the navigation solution.
4. **Use pose priors in bundle adjustment** (COLMAP >= 3.10 or pycolmap) instead of
   the sliding-window fit. With navigation as a prior, drift and scale are handled
   inside the optimisation and the per-image fit becomes unnecessary.
5. **Use the SfM surface instead of the altimeter DEM** once the model is metric: a
   Poisson or Delaunay mesh from the sparse/dense points captures relief.
6. **Reduce matching cost** for longer transects: sequential overlap 8 instead of
   15 (overlap vanishes after 5 to 6 frames) cuts matching time by about two thirds.
7. **Better blending** (multi-band or graph-cut seams) once registration is at the
   centimetre level; before that, blending only hides errors.
8. **Vehicle behaviour.** The 30 deg yaw oscillation drives most of the trouble.
   Slower speed or better heading control on the Otter would help every step above.

## Closing the loop to the hyperspectral data

The RGB poses have been written into a copy of the mission
(`gref4hsi_test_sfmpose`) as `raw/nav_rgbsfm` and gref4hsi was re-run on them with
the HSI 5 cm behind the camera. Assumptions specific to this step:

* The RGB image-up direction is the vehicle's forward axis. Check: image-up heading
  minus compass heading is -1.1 deg median (IQR -4 to +2 deg) over 217 frames.
* The HSI boresight from `HSI_2b.xml` (rz = -90 deg, rx = ry = 0) still applies when
  the body frame is defined from the RGB camera. Residual boresight errors of a few
  degrees are expected and would need the in situ luminance calibration.
* Poses at 3 Hz are interpolated (slerp) to the ~25 Hz HSI lines. With yaw rates up to
  30 deg/s this smooths out sub-frame motion; the navigation's attitude rates could be
  blended in for the high-frequency part.
* The DEM is still the altimeter surface.

## Open issue: HSI and RGB products do not coincide (not worked on yet)

After the loop closure, the hyperspectral composite and the RGB mosaic follow the
same strip but are offset against each other: the HSI strip sits to one side of the
RGB strip and a sea urchin visible in both appears about half a swath apart. Both
products use the same camera poses, so the pose is not the cause. Candidates, most
likely first:

1. **Mirrored HSI slit.** The boresight in `HSI_2b.xml` (rz = -90 deg) came from the
   vehicle-frame preprocessing; with the body frame defined from the RGB camera the
   sign may be wrong, which flips left and right across the track.
2. **Roll/pitch tilt of the HSI relative to the camera** (rx, ry are 0 in the
   calibration file). ~2 cm of swath shift per degree at 1.2 m altitude.
3. **Clock offset between HSI lines and RGB frames.** Both are on the UHI clock, but
   the timestamps may mark different moments of the exposure/readout; 1 s = 45 cm
   along track.
4. Not likely: the 5 cm lever arm, the shared DEM, a 180 deg mounting error (the
   image-up heading matches the compass heading to ~1 deg).

Planned test: cross-correlate the two rasters along the strip and split the shift
into along-track and across-track parts using the heading. A sign change across the
strip means a mirror, a constant across-track part a tilt, an along-track part a
clock offset. Then re-run gref4hsi (~4 min) with the corrected boresight sign/tilt or
clock offset to confirm. Still open: in the raw RGB frames, does the seabed move
down the image when the vehicle moves forward (i.e. is image-up really forward)?

## Altitude check: SfM versus altimeter (measured, not acted on)

Per-frame camera height above the SfM seabed points (horizontal scale from the
navigation, f = 487 px) against the 1 Hz altimeter, over the whole transect:

| Quantity | Value |
|---|---|
| Altimeter | 1.17 to 1.29 m, std 3 cm |
| SfM height above points | median 1.35 m, std 11 cm |
| Ratio SfM / altimeter | median 1.08, 20-frame block means 1.02 to 1.14 |
| Difference SfM - altimeter | median +10 cm |
| Correlation with window scale | 0.01 (not an alignment artefact) |

Reading: a constant ~8 % (10 cm) offset, i.e. either the focal length is ~8 % too
high (~450 px) or the altimeter sits ~10 cm below the camera. The SfM scatter (11 cm)
is four times the altimeter variation, so terrain or kelp-canopy effects cannot be
separated on this transect. The effect scales the swaths about their centrelines by
a few centimetres at the edges and does not explain the HSI/RGB offset above. No
correction applied; revisit together with the focal-length calibration.

## References

* Løvås, H. S., Mogstad, A. A., Sørensen, A. J., Johnsen, G. (2022). A Methodology
  for Consistent Georegistration in Underwater Hyperspectral Imaging. IEEE J.
  Oceanic Eng. 47(2), doi:10.1109/JOE.2021.3108229.
* Mission README with the full log of runs: `<mission>/gref4hsi_test/colmap/README.md`.
