time ./map3r << EOF
yes -----------------change abundances to asplund et al abundances 
file                 ***************************
abund/1sol.abn
no ------------------change abund : abundance file now asplund so no further changes needed.
no ------------------disable react : ignore this - this is always no 
yes -----------------include dust : dust can be included or expluded depending on models
yes -----------------change depletions
file
abund/1sol.dpl
no ------------------change depletions : these depletions are based on the solar neighborhood, may be excessive for ionized regions.
Shattering grain profile ---- note: change DUSTDAT to full Here you can choose MRN for simplicity or a power-law form
-3.3 ----------------power law
40 ------------------graphite grain min. radius
1600 ----------------graphite grain max. radius
40 ------------------silicate grain min. radius
1600 ----------------silicate grain max. radius
1.8 -----------------graphite grain density
3.5 -----------------silicate grain density
yes -----------------Include PAH molecules
0.05 -----------------fraction of Carbon Dust Depletion in PAHs
Q -------------------Q-value QHDH < Value PAH photodissociation parameter
1000 ----------------PAH switch on Value : PAHs appear more sensitive to ionizing photons so we have a switch for them
no -----------------Destroy all C grains with PAHs
no -----------------Evaluate dust temperatures and IR flux This calculates the IR emission - it will make the computation ~2-4x longer
P5 model	: this is the important split - for most people P5 (photoionization) or S4 (shocks) is sufficient
Default preion. balance : -this sets up the initial step. the closer to reality ie 100% ionized hydrogen for HII regions the quicker and better the model 
Blackbody :here you enter the ionizing source, which cna be a model (eg SB99) or functional forms as listed
20000 : temperature in Kelvin
0 -------------------turn-on energy (in Rydberg) (depends on model
4.0 ---------------cut-off energy (in Rydberg)
1.0 ---------------scaling factor : can add fractional components
E  ------------------Bremsstrahlung component (c*exp(-nu/k2e7))
3e6 -----------------Temperature
3e18 ----------------Flux
13.6 -------------------turn-on energy (in eV)
15000 ---------------cut-off energy 
1.0 -----------------scaling factor
Xit with current ionizing spectrum
No ------------------cosmic ray heating : causes continuous low leverl ionization and heating
s -------------------Spherical geometry (can have spherical or plane-parallel geometry)
L -------------------Luminosity of ionizing Source
T -------------------Total Luminosity
40.146                (in log (erg/s))
C -------------------constant density : can also have constant pressure, which includes radiation pressure
100.0  --------------density of H (cm^-3)
1.0 -----------------Filling Factor
d -------------------Distance sets initial radius from ionizing source of HII region (either by distance or Flux)
1.97E+20             (cm)************************
yes -----------------Integrate over whole sphere (in spherical systems)
Equilibrium ---------ionisation balance
0.04 ----------------photon absorption fraction (this should always be from 0.02-0.05)
A  :   Radiation bounded, HII < 1.00%
G			: Output format [G=Everything]
2 6 8 10		: Atomic numbers of four elements to monitor
HII_5
HII_5
HII Z=1sol, P/k=1e4, t=5Myr
X
EOF

