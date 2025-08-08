# mini-batching-evolution

Wright Fisher dynamics with a stochastic genotype-to-fitness map. 

At each generation, each individual picks a random fitness value from a genotype-dependent distribution (*individual mini-batching*). This generates an effective noise that can alter the stationary distribution of population states. We compare the results with a reference process where evolution takes place in a fixed environment given by the average genotype-to-fitness map.

The simulations can be run from model.ipynb notebook. Modify simulation and computing parameters as desired.

The analysis reported in the manuscript's figures can be run from the model_interface/analysis folder. 

Simualtions adn analysis related to the 3-species case can be run frm the 3species folder.