# Directory structure

- Exposition/
  - material presented for the exhibition
- formExploration/
  - Openmole/electra: contain .oms files to run the model in openmole and perform different analyses/explorations : directSampling, NSGA (is it possible to produce a given trajectory with electra?), and PSE (explore the diversity in trajectories)  
  - R electra: We have also implemented the model in R.
    - R model: R code to plot trajectories (we indicate some parameters that produce simple figures) and to save them (to create a movie)
    - postProc_OM: R code to analyse results from openmole, sorted by explicit name directory (NSGA,PSE,direcSampling,1run). The final figure is obtained in the PSE folder.
  - scala: contain the scala code of the model. In particular there is a function to produce the stationary trajectory (also a function for the non stationary bu not used yet),  It also contains functions that compute measures on the trajectories (loop points, Moran index, density index, curvature and singular points.
- shadow/
  - interactive interface for running and visualising the model. See `shadow/README.md`.
