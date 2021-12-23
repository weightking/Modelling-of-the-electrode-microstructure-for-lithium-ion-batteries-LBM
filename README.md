# Modelling-of-the-electrode-microstructure-for-lithium-ion-batteries-LBM

A three-dimensional pore-scale lattice Boltzmann model (LBM) was developed to simulate LIBs electrodes, capturing the transport of ions, electrons, and liquid electrolyte species, coupled with the electrochemical reaction at the interface between the active material and the liquid electrolyte with complex boundary conditions. The simulated discharge curves of a realistic nickel-manganese-cobalt (NMC) battery electrode are in good agreement with the experimental measurement, demonstrating the validity of the model prediction. 

A 3D pore-scale LBM model was previously developed at the University of Surrey to simulate electrodes in proton exchange membrane fuel cells with gas-liquid two-phase flow in the porous gas diffusion layer, coupled with electrochemical reactions at the catalyst layer [1]. This model was later developed to simulate redox flow batteries [2-3] with liquid flow coupled with electrochemical reactions at the electrode-electrolyte interface within the 3D porous electrodes to study electrode wetting and effect of different porous structure on the electrolyte flow and electrochemical performance. In the current study, this model is further developed to simulate LIB electrodes with transport of ions, and electrolyte species within the 3D porous space, coupled with electrochemical charge transfer reactions at the active material particle-electrolyte interface, and Li ion and electron transport through the active material particles.  

The 3D half-cell model developed based on the LBM model is verified by comparing the model prediction with the experimentally measured electrochemical performance [4] for the NMC electrode. To implement the NMC electrode microstructure provided by Ebner et al. [5] in our LBM simulation, an additional coarsening step is used to reduce the computational cost by reducing the voxel resolution to 1.48 μm from 0.37 μm in the raw tomography structure. A resolution sensitivity study is carried out by comparing the model prediction based on the structure of 1.48 μm resolution with that of 0.74 μm. It is shown that the resolution of 1.48 μm gives the same electrochemical performance as that of the higher resolution and thus the accuracy of the LBM simulation on the structure of 1.48 μm resolution can be guaranteed. In addition, the resolution choice of 1.48 μm is the same as the previous numerical study [6]. 
Table 1 summarize the parameters that are adopted for this verification study. 

Table 1. Input parameters for the LBM simulations

Electrolyte
Initial concentration ce0	1000 mol m-3 [75]
transference number t+	0.4 [76]
Diffusivity De	2.72×10-10 m2 s-1[76]
Electrode
Initial concentration cs0	7976.16 mol m-3
Electrical conductivity σs	1.0 S m-1 [52]
Diffusivity Ds	2.0×10-15 m2 s-1[77]
Reaction rate k0	6.325×10-8 m2.5 mol-0.5 s-1 [a]
Maximum capacity csmax	36224 mol m-3 [52]
[a] Determined by fitting to the electrochemical data.
 
1. D. Zhang, A. Forner-Cuenca, O. O. Taiwo, V. Yufit, F. R. Brushett, N. P. Brandon, S. Gu, Q. Cai, Understanding the role of the porous electrode microstructure in redox flow battery performance using an experimentally validated 3D pore-scale lattice Boltzmann model, J. Power. Sources. 447(2020) 227249. https://doi.org/10.1016/j.jpowsour.2019.227249.

2. D. Zhang, Q. Cai, O. O. Taiwo, V. Yufit, N. P. Brandon, S. Gu, The effect of wetting area in carbon paper electrode on the performance of vanadium redox flow batteries: A three-dimensional lattice Boltzmann study, Electrochimica Acta. 283(2018) 1806-1819. https://doi.org/10.1016/j.electacta.2018.07.027.

3. D. Zhang, Q. Cai, S. Gu, Three-dimensional lattice-Boltzmann model for liquid water transport and oxygen diffusion in cathode of polymer electrolyte membrane fuel cell with electrochemical reaction, Electrochimica Acta. 262(2018) 282-296. https://doi.org/10.1016/j.electacta.2017.12.189.

4. M. Singh, J. Kaiser, H. Hahn, Thick electrodes for high energy Lithium ion Batteries, J. Electrochem. Soc. 162 (2015) A1196-A1201. doi:10.1149/2.0401507jes.

5. M. Ebner, F. Geldmacher, F. Marone, M. Stampanoni, V. Wood, X‐ray tomography of porous, transition metal oxide based lithium ion battery electrodes, Adv. Energy Mater. 3 (2013) 845–850. https://doi.org/10.1002/aenm.201200932.

6. T. Danner, M. Singh, S. Hein, J. Kaiser, H. Hahn, Thick electrodes for Li-ion batteries: A model based analysis, J. Power Sources. 334 (2016) 191-201. http://dx.doi.org/10.1016/j.jpowsour.2016.09.143.