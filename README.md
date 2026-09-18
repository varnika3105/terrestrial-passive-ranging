# Terrestrial Passive Ranging Using Angle-Only Measurements

A MATLAB-based simulation and estimation pipeline for terrestrial passive target localization using only azimuth and elevation measurements from multiple sensors.

## Overview

The project combines classical estimation with deep learning:

**Angle Measurements → WLS Triangulation → Kalman Filtering → BiLSTM Residual Correction**

The simulation models moving targets, moving sensors, sensor-position uncertainty, angular measurement noise/bias, and geometry quality. The BiLSTM learns the residual error remaining after Kalman filtering.

## Methods

- Angle-only 3D line-of-sight measurement generation
- Weighted least-squares / geometric triangulation
- Kalman filtering for trajectory smoothing
- BiLSTM residual learning for position correction
- Quantitative evaluation using position RMSE and error statistics

## Repository Structure

```text
matlab/
├── passive_ranging_baseline.m
├── moving_sensor_triangulation.m
├── triangulation_with_kalman.m
├── passive_ranging_dataset_generator.m
└── passive_ranging_bilstm.m
results/
└── training process.png
```

## Requirements

- MATLAB
- Statistics/Deep Learning functionality required by the BiLSTM script (e.g. `sequenceInputLayer`, `bilstmLayer`, `trainNetwork`)

## Running the project

1. Open MATLAB and set the repository root as the working directory.
2. Run `passive_ranging_baseline.m` for the basic triangulation experiment.
3. Run `moving_sensor_triangulation.m` or `triangulation_with_kalman.m` for moving-sensor estimation and filtering.
4. Run `passive_ranging_dataset_generator.m` to generate the multi-trajectory dataset used for learning.
5. Run `passive_ranging_bilstm.m` to train and evaluate the BiLSTM residual model.

> **Note:** The learning script expects the generated CSV dataset with the feature columns used in the code. Large generated datasets and trained model binaries are intentionally not committed to keep the repository lightweight. Generate them locally using the provided scripts.

## Experimental setup

The research simulations use multiple sensors observing a moving 3D target. Measurement uncertainty is introduced through sensor-position errors and azimuth/elevation noise and bias. The learning experiment uses separate trajectory runs for training and testing.

## Results

The project evaluates localization quality using RMSE and position-error comparisons. A training-progress figure is included under `results/`. Reported experimental results should be interpreted as simulation results rather than field measurements.

## Reproducibility

The MATLAB scripts use fixed/randomized seeds in different experiments. For exact reproduction of a particular experiment, use the corresponding script and preserve its parameter settings and random seed.

## Important note

This repository is a cleaned portfolio version of a research project. It intentionally excludes reports, internal documents, proprietary material, and large generated datasets. Only non-sensitive implementation and simulation material should be published.
