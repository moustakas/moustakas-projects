function read_sedcube, templates=templates
;+
; NAME:
;	READ_SEDCUBE()
;
; PURPOSE:
;	Read in all the model galaxy SEDs.
;
; INPUTS:
;
; OPTIONAL INPUTS:
; 	models - select the SED models to read
; 
; OUTPUTS:
;	sedcube - structure containing the fields described below
;    
; COMMENTS:
;	Need to add Bruzual-Charlot templates.
;    
; PROCEDURES USED:
;	READFAST
;
; MODIFICATION HISTORY:
;	John Moustakas, 2000 October 25, U of A
;	jm01may8uofa, added templates keyword
;-

    if not keyword_set(templates) then templates = 'devriendt'

    path = '/home/ioannis/synthesis/devriendt/'
;   path = filepath('',root_dir=getenv('SIRTFZ_DIR'),subdirectory='seds/'+strlowcase(templates))
    readcol, path+'sed.dat', nr, sednames, names, gtypes, format='I,A,A,A', /silent, skip=2
    nseds = n_elements(sednames)

; parse the true galaxy name

    truenames = strarr(nseds)
    for i = 0L, nseds-1L do truenames[i] = strjoin(strsplit(names[i],'_',/extract),' ')
    
    template = {galaxy :          '', $	; model galaxy name
                gtype  :          '', $	; galaxy type
                lbolir :         0.0, $ ; bolometric IR luminosity [L_sun]
                lambda :   ptr_new(), $ ; wavelength (microns)
                mlum_nu:   ptr_new(), $ ; monochromatic luminosity [W/Hz]                
                mlum_mu:   ptr_new()}   ; monochromatic luminosity [W/mu]
    sedcube = replicate(template,nseds)

    lsun = 3.826D26    ; [W]
    light = 2.99793D14 ; [mu/s]
    
    for k = 0L, nseds-1L do begin
       
       readfast, path+sednames[k], seddata ; read the SED

       lambda = seddata[0,*]  ; [mu]
       mlum_nu = seddata[1,*] ; [W/Hz]
       
       mlum_mu = light * mlum_nu / lambda / lambda ; [W/mu]

; bolometric IR luminosity [1-1000 mu]
       
       ir = where((lambda gt 1.0) and (lambda lt 1000.0),nir) ; IR SED
       if nir ne 0L then lbolir = int_tabulated(lambda[ir],mlum_mu[ir],/double)/lsun else lbolir = 0.0

; fill the structure

       sedcube[k].galaxy = truenames[k]
       sedcube[k].gtype = gtypes[k]
       sedcube[k].lbolir = lbolir
       sedcube[k].lambda = ptr_new(lambda)
       sedcube[k].mlum_nu = ptr_new(mlum_nu)
       sedcube[k].mlum_mu = ptr_new(mlum_mu)

    endfor              

return, sedcube
end

