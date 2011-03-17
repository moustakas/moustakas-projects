pro sirtfz_initialize, omega_0=omega_0, omega_lambda=omega_lambda, h_100=h_100, $
          zmin=zmin, zmax=zmax, dz=dz, alpha=alpha, beta=beta, templates=templates, $
          filters=filters
;+
; NAME:
;	SIRTFZ_INITIALIZE
;
; PURPOSE:
;	Initialize the defaults parameters and the common blocks for
;	the photometric redshift program SIRTFZ.
;
; INPUTS:
;
;
; OPTIONAL INPUTS:
;	omega_0      - matter density
;	omega_lambda - vacuum energy density
;	h_100        - Hubble constant normalized by 100 km/s/Mpc
;	zmin         - minimum redshift (when simulating photometric
;                      catalogs) 
;	zmax         - maximum redshift
;	dz           - redshift interval
;	alpha        - luminosity evolution power
;	beta         - density evolution power
;	templates    - SED templates to use (Devriendt, CWW)
;	filters      - default SIRTFz filters
;	
; COMMON BLOCKS:
;	sirtf_simulations
;	cosmology
;
; COMMENTS:
;
; PROCEDURES USED:
;	READ_SEDCUBE()
;	READ_BANDCUBE()
;	READ_LF()
;	
; MODIFICATION HISTORY:
;	J. Moustakas, 2001 August, U of A
;-

    common sirtf_simulations, sirtf
    common cosmology, cosmo_params

    if not keyword_set(omega_0) then omega_0 = 0.3D
    if not keyword_set(omega_lambda) then omega_lambda = 0.7D
    if not keyword_set(h_100) then h_100 = 0.65D
    if not keyword_set(zmin) then zmin = 0.03D
    if not keyword_set(zmax) then zmax = 6.0D
    if not keyword_set(dz) then dz = 0.03D
    if not keyword_set(alpha) then alpha = 0.0D
    if not keyword_set(beta) then beta = 0.0D
;   if not keyword_set(templates) then templates = 'CWW'
    if not keyword_set(templates) then templates = 'DEVRIENDT'
    if not keyword_set(filters) then filters = ['IRAC 3.6', 'IRAC 4.5','IRAC 5.8','IRAC 8',$
                                                'MIPS 24','MIPS 70','MIPS 160']
    
; clean memory from a previous call

    sirtfz_shutdown, /all

; read
    
    sedcube = read_sedcube(templates=templates) ; SED templates
    bandcube = read_bandcube()                  ; filter structure
    lf = read_lf()                              ; IRAS 60mu luminosity function

    cosmo_params = {name: 'Cosmology Parameters', $ ; cosmology structure
                    omega_0: omega_0, $
                    omega_lambda: omega_lambda, $
                    h_100: h_100}

    zarray = (findgen((zmax-zmin+dz)/dz+1.0))*dz
;   zarray = (findgen((zmax-zmin+dz)/dz))*dz+zmin
    nz = n_elements(zarray)
    
    redshift = {name: 'Redshift Grid', $            ; redshift structure
                zmin: zmin, $
                zmax: zmax, $
                dz: dz, $
                zarray: ptr_new(zarray), $
                nz: nz}

    evolution = {name: 'Evolution Parameters', $    ; luminosity/density evolution
                 alpha: alpha, beta: beta}
    
; structure containing all the above information    
    
    sirtf = {name      : 'SIRTFz',   $
             evolution : evolution,  $
             sedcube   : sedcube,    $
             templates : templates,  $
             bandcube  : bandcube,   $
             lf        : lf,         $
             filters   : ptr_new(filters), $
             redshift  : redshift,   $
             catalog   : '',         $ ; photometric catalog to analyze (none by default)
             baseid    : lonarr(1)}    ; SIRTFz widget base ID number

return
end
