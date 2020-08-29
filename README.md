# ATgeomPSFsim
Simulate geometrical PSF for Auxtel, for Ronchi and Hologram

- author : Sylvie Dagoret-Campagne

- affiliation : IJClab/IN2P3/CNRS (ex LAL)
- creation date August 27th 2020
- update August 29th 2020


## Generate Beamfour Ray files

     
        >python generateB4rayfile.py --config default.ini
        >python generateB4rayfile.py --config config/conf_CTIO_2020_08_29.ini
        >python generateB4rayfile.py --config config/conf_AuxTel_2020_08_27.ini

- the two last config files are used for production

- the output file may be found in **beamfour\_data/** to work with beam_four. 

- **beamfour\_data/** contains *.OPT and *.RAY files used by Beamfour


## tools : Bac à sable