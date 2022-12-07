# =========================== Interdependent Superconducting Nnetworks =========================== #

The files in this folder refer to the thermally interdependent resistively shunted Josephson junction (RSJJ) networks studied in the article entitled "Interdependent Superconducting Networks" and freely accessible at https://doi.org/10.48550/arXiv.2207.01669. 

The folder "Video files" contains the link to the Supplementary Videos S1 and S2 (referred in the article, see Fig.2 and details therein) where the emergence of hotspots and the heterogeneous redistribution of the currents in the two arrays during the electro-thermal cascading processes can be visually observed. 

The folder "MATLAB code" contains instead the codes needed to execute the thermally coupled Kirchhoff equations solving the superconductor-normal metal transitions and the electro-thermal propagation of overheating cascades within and across the layers. The code presented in this folder adopts a global thermal coupling obtained by increasing adaptively the cryostat temperatures of the two layers via the Joule heating dissipation measured via the macroscopic resistive behavior of the two arrays. 

The zip folder "RSJJ_local.zip" provides a similar set of algorithms as in the above with the major difference that now the local thermal coupling is set locally, i.e. by means of the Joule heating produced by single junctions. Qualitatively, global and local thermal couplings result into the same phase diagram for the mutual SN-phase transitions of the model (for details about this point, see the Methods section of the article, Secs. M5 and M6). 
