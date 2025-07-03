# Numerical_Simulations_Confined_Langevin_Equation

This is my PhD code to simulate the Brownian motion confined between two rigid walls.
It is make with `Python`, `Cython` and I use `Jupyter Notebook` to analyse datas.

## The problem :
![Schematic of the system. A particle of radius *a* diffuses in two dimensions, between two walls separated by a distance *2H*<sub>p </sub>.](https://media.springernature.com/lw685/springer-static/image/art%3A10.1140%2Fepje%2Fs10189-023-00281-y/MediaObjects/10189_2023_281_Figa_HTML.png?as=webp) \
This is the schematic of the system. \
A particle of radius *a* diffuses in two dimensions *(x<sub>t</sub>, z<sub>t</sub>)*, immersed in a fluid of viscosity *η*<sub>0</sub>, undergoing acceleration gravity *g* and confined between two walls separated by a distance *2H*<sub>p</sub>.

The diffusion is affected by the walls presence and depend of the high *z*<sub>t</sub>. \
[🔗Here](Double_Walls_Overdamped_Langevin_Python/Diffusion.pdf) is the plot the parrallèle diffusion *D*<sub>∥</sub>*(z*<sub>t</sub>*)* $${\color{#009de0}(in \space blue)}$$ and perpendicular diffusion *D*<sub>⊥</sub>*(z*<sub>t</sub>*)* $${\color{#008000}(in \space green)}$$ in fonction of the high *z*<sub>t</sub>.

[🔗Here](Optimisations/Figures/Traj_Peq.pdf) there is `(a)` an example of numerical tajectory *(x<sub>t</sub>, z<sub>t</sub>)* in fonction of time *t*, with `(b)` the Gibbs-Boltzmann equilibrium distribution of it altitude *z<sub>t</sub>* (the probability *P*<sub>eq</sub> to find the particle at a certain high *z*<sub>t</sub>).

### For more informations about the project, 
- you can find my PhD manuscript about this project (french) [🔗here]([Optimisations/Figures/Traj_Peq.pdf](https://theses.hal.science/tel-04583730)).

- and my article research (english) [🔗here](https://link.springer.com/article/10.1140/epje/s10189-023-00281-y).

## ©️License and Citation
This code is released under the Creative Commons Attribution 4.0 International License (CC BY 4.0).
You are free to use, modify, and share this code for any purpose, provided that proper credit is given.

🔗 License details: [Creative commons licenses 4.0/](https://creativecommons.org/licenses/by/4.0/)

If you use this code or any part of it in your work, please cite:

> **@ElodieMillan**, "[Numerical_Simulations_Confined_Langevin_Equation]" (2023), Github, https://github.com/ElodieMillan/Numerical_simulations_of_confined_brownian_motion,

## 📰References: 
[1] Arthur Alexandre, Maxime Lavaud, Nicolas Fares, **Elodie Millan** et al. _« Non-Gaussian Diffusion
Near Surfaces »_. In : Phys. Rev. Lett. 130 (7 (2023)), p. 077101, 🔗https://doi.org/10.1103/PhysRevLett.130.077101 ,

[2] **Elodie Millan**, Maxime Lavaud, Yacine Amarouchene et Thomas Salez. _« Numerical simulations of
confined Brownian-yet-non-Gaussian motion. »_ In : Eur. Phys. J. E 46 (4 (2023)), p. 24, 🔗https://doi.org/10.1140/epje/s10189-023-00281-y.

Thank you for your interest and collaboration!

![Enjoy 'Python'!](https://myoctocat.com/assets/images/base-octocat.svg)
