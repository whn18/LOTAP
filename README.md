# LOTAP
A tensor singular value decomposition (t-SVD) based third-order time-series forecasting algorithm, which tactically incorporates the unique advantages of t-SVD and AR into a unified framework.  
More details (including parameter settings) refer to [the original paper](https://arxiv.org/abs/2403.02835).

### Paper
- [Low-rank Tensor Autoregressive Predictor for Third-Order Time-Series Forecasting](https://arxiv.org/abs/2002.12135)

### Datasets
  
Synthetic (SYN) dataset. The SYN dataset is a low-rank, third-order tensor time series that generated .More details refer to [the original paper](https://arxiv.org/abs/2403.02835). We here select the first **80** time points to forecast the following **20** time points.

### Getting Started

#### Prerequisites  

- Matlab
  -   <a href="https://www.tensorlab.net/" class="textlink">Tensolab Toolbox</a>

#### Run

``` MATLAB
main.m
```
