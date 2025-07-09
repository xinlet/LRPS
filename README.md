# Large Multi-response Linear Regression Estimation Based on Low-rank Pre-smoothing

## Abstract

Pre-smoothing is a technique aimed at increasing the signal-to-noise ratio in data to improve subsequent estimation and model selection in regression problems. Motivated by the many scientific applications in which multi-response regression problems arise, particularly when the number of responses is large, we propose here to extend pre-smoothing methods to the multiple outcomne setting. Specifically, we introduce and study a simple technique for pre-smoothing based on low-rank approximation. We establish theoretical results on the performance of the proposed methodology, which show that in the large-response setting, the proposed technique outperforms ordinary least squares estimation with the mean squared error criterion, whilst being computationally more efficient than alternative approaches such as reduced rank regression. We quantify our estimator's benefit empirically in a number of simulated experiments. We also demonstrate our proposed low-rank pre-smoothing technique on real data arising from the environmental and biological sciences.

## Authors

[Xinle Tian](https://xinlet.github.io/), [Alex Gibberd](https://sites.google.com/view/gibberd/) , [Matthew Nunes](https://people.bath.ac.uk/man54/homepage.html), [Sandipan Roy](https://roysandipan.github.io/)

## ArXiv link
Preprint paper link can be found at [https://arxiv.org/abs/2411.18334]<br />

## Datasets
UK air pollution dataset can be found at [https://uk-air.defra.gov.uk/]<br />
Beijing air pollution dataset can be found at [http://archive.ics.uci.edu/ml/]<br />
US air pollution dataset can be found at [https://www.epa.gov/outdoor-air-quality-data)]<br />
Gene association dataset can be found at [https://genomebiology.biomedcentral.com/articles/10.1186/gb-2004-5-11-r92]<br />

## Code and output
[UK dataset Output](R/US-output.md)<br />
[Beijing dataset Output](R/Beijing-output.md)<br />
[US dataset Output](R/US-output.md)<br />
[Geno dataset Output](R/Geno-output.md)<br />
