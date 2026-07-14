## Why we still doing CMB
- $H_0$ tension models
- $\Delta N_{\rm eff}$ 
- DM interactions
- B modes
- Secondaries: reconstructions and cross correlations (classical ones: lensing, tSZ, kSZ)
	- Neutrino mass
	- PNL
	- Baryonic feedback

LCDM agree wells, no matter Planck, ACT or SPT
### Typical data: e.g. ACT DR6.02
- maps, nilc, beams, passbands, chains/likelihood, depth1 (time domain data), simulations, pspipe (measure PS of the maps)

### Detector noise:
NET (Noise Equivalent Temperature): temperature fluctuation RMS gives SNR=1 in 1 second of integration. Pixel noise: $\sigma_p=NET_{\rm arr}/\sqrt{t_p}$ 
	- Ground: around $300 \mu K s^{1/2}$ for individual detectors, space is 10x lower
	- Can convert to survey property: $t_p = t_{\rm eff} \Omega_{\rm pix}/A$ , but this leads to noise depends on size of pixel
	- Map depth ("sensitivity") independent of pixerl size: Define  $\Delta_T = \sigma_p \sqrt{\Omega_{\rm pix}}$      
	- $\Delta_T = NET_{\rm arr}\sqrt{A/t_{\rm eff}}$  
	- Note all have units: ($\mu K$- arcmin) as sqrt(angular power), interpret as standard deviation
	- All applied to white noise, correlated atm noise does not avg down in time and detector counts
	- atm noise is correlated, but not polarized (except ice crystals), red noise (1/f)
### Noise properties
- Correlated noise
- Inhomogeneous (hit variation)
- Anisotropic (2D Fourier space)
	- Each scan has a perpendicular direction, plane wave -> DC mode
	- Cross-linking alleviates this
- Inhomogeneity of anisotropy
	- Cross-linking varies across the map

### Telescope scans
Constant elevation scans: Constant elevation, scan in azimuth
	- But projects in RA and Dec, tracks cross-link, sky drifts (look rising and setting part of the sky)

### Pixell and rectangular pixelization
- Healpix is analyzed by healpy, Healpix maps are defined by $N_{\rm side}$ 
- Rectangular pixel maps can be pixell, unified by FITS world coordinate system (WCS)
- Rectangular pixel does NOT mean "Flat-sky", that only kicks in in Fourier transform 

### Secondaries:
- Moving lens: CMB photons lensed by different part of potential as it enters and exits moving lens halo
	- Dipole-like signal probing transverse velocity field
	- Lensing + kSZ
- Patchy screening: reionization, CMB photons scattered into or out of LOS by Thomson scattering off ionized electrons, proportional to optical depth, suppressed by primary anisotropy.
- Polarized SZ: tSZ is also polarized, CMB scatter off ionized gas, quadrupole leads to polarization just like at recombination
- Baxter effect: rotational kSZ