We imagine that the expected log chance species $i$ is detected at site $j$ is a
function of the number of individuals sampled at site $j$ $n_j$, and the
(unknown) per-individual detection probabilities of species $i$ at site $j$,
$p_{i,j}$. We will define $q_{i,j} = \log(-\log(1-p_{i,j}))$, explained below.
Our model would look something like this: $$\log(d_{i,j}) \sim c_{j} +
\log(n_{j}) + q_{i,j}$$ We imagine that the per-individual detection
probabilities of species at each site sum to 1 (though not all species will
necessarily be detected in a given sample, so the sum of the detection
probabilities of the observed species may be <1), and that these probabilities
are normal on the log-hazard scale, i.e.,  $$q_{i,j} \sim N(0, \sigma_{j})$$,
and $$\sum_{i}{p_{i,j} =1}$$ at each site $j$.
