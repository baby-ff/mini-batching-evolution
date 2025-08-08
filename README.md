# mini-batching-evolution

Wright Fisher dynamics with a stochastic genotype-to-fitness map. 

At each generation, each individual picks a random fitness value from a genotype-dependent distribution (*individual mini-batching*). This generates an effective noise that can alter the staionary distribution of population states. 

**Question**: can stochasticity deriving from mini-batching help generalization? The goal of evolution is to learn about the statistics of environments/prepare for future environment perturbations.

To answer this question, we compare results with the reference case of a fixed environment given by the average genotype-to-fitness map.
