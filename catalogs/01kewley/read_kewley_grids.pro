;+
; NAME:
;       READ_KEWLEY_GRIDS()
;
; PURPOSE:
;       Read in the Kewley et al. photoionization model grids. 
;
; INPUTS:
;
; OPTIONAL INPUTS:
;       model - model number (see below)
;
; KEYWORD PARAMETERS:
;       oversamp - oversample the metallicity and ionization parameter
;                  spacing 
;       nolog    - return the linear line ratios (default is to return
;                  the logarithmic line ratios)
;
; OUTPUTS:
;       grids - output data structure (see below)
;
; OPTIONAL OUTPUTS:
;       Z - metallicity vector
;       U - ionization parameter vector
;
; COMMENTS:
;       This routine assumes that in the published models
;       [OIII] and [NII] corresponds to the 5007 and 6584 lines,
;       respectively. 
;
; EXAMPLES:
;
; MODIFICATION HISTORY:
;       J. Moustakas, 2003 Aug 06, U of A - written
;       jm04jun27uofa - added the NOLOG and OVERSAMP keywords
;       jm05mar21uofa - documented; added additional models
;       jm08may13nyu - documentation tweaks
;
; Copyright (C) 2003-2005, 2008, John Moustakas, U of A
; 
; This program is free software; you can redistribute it and/or modify 
; it under the terms of the GNU General Public License as published by 
; the Free Software Foundation; either version 2 of the License, or
; (at your option) any later version. 
; 
; This program is distributed in the hope that it will be useful, but 
; WITHOUT ANY WARRANTY; without even the implied warranty of
; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
; General Public License for more details. 
;-

function read_kewley_grids, model=model, Z=Z, U=U, $
  oversamp=oversamp, nolog=nolog
    
    light = 2.99792458D10        ; speed of light [cm/s]

; correction factors    

    bb = im_branch_ratios()
    nratio = bb.n_ii  ; 3.054
    oratio = bb.o_iii ; 2.984
    
    ncor = 1.0+1.0/nratio
    ocor = 1.0+1.0/oratio

    logncor = alog10(ncor)
    logocor = alog10(ocor)

    path = filepath('',root_dir=getenv('CATALOGS_DIR'),subdirectory='01kewley')

    if n_elements(model) eq 0L then model = 1
    
    case model of

       1L: file = 'starburst_pegase2_cont.dat'
       2L: file = 'starburst_pegase2_inst.dat'

       3L: file = 'starburst_SB99_cont.dat'
       4L: file = 'starburst_SB99_inst.dat'

       5L: file = 'HII_pegase2_cont.dat'
       6L: file = 'HII_pegase2_inst.dat'

       7L: file = 'HII_SB99_cont.dat'
       8L: file = 'HII_SB99_inst.dat'

       else: begin
          splog, 'Model '+strn(model)+' not supported yet.'
          return, -1L
       endelse

    endcase

; read the data

    if file_test(path+file) then begin
       readfast, path+file, data, skipline=10, ncols=ncols, nlines=nlines
    endif else begin
       splog, 'File '+path+file+' not found.'
       return, -1L
    endelse

    Z = reform(data[0,uniq(data[0,*])]) ; unique metallicities [Z_sun]
    q = [5E6,1E7,2E7,4E7,8E7,1.5E8,3E8] ; *fixed* ionization parameter [dimensionless]
    U = alog10(q/light) 

    nZ = n_elements(Z)
    nU = n_elements(U)
    nratios = ncols-1L

    Zindx = lindgen(nZ)
    Uindx = lindgen(nU)

    bigZ = Z # (fltarr(nU)+1)
    bigU = U # (fltarr(nZ)+1)
    
; initialize the 3D data structure: 0th = metallicity; 1st =
; ionization parameter; for example, grids[0,*].oiii_oii is the ratio
; [O III]/[O II] over the range of ionization parameter U at a fixed
; metallicity (Z=0.01 for MODEL=1)

    grids = {$
      nii_6584_h_alpha:                  0.0, $ ; [N II] 6584 / Ha
      oiii_5007_h_beta:                  0.0, $ ; [O III] 5007 / Hb
      sii_h_alpha:                       0.0, $ ; [S II] 6716,6731 / Ha
      oi_6300_h_alpha:                   0.0, $ ; [O I] / Ha
      nii_6584_oii:                      0.0, $ ; [N II] 6584 / [O II] 3727
      oiii_5007_nii_6584:                0.0, $ ; [O III] 5007 / [N II] 6584
      oiii_5007_oii:                     0.0, $ ; [O III] 5007 / [O II] 3727
      oii_h_alpha:                       0.0, $ ; additional line ratios begin here
      oii_h_beta:                        0.0, $ ; ---------------------------------
      nii_6584_sii:                      0.0, $ ; [N II] 6584 / [S II] 6716,6731
      nii_6584_oiii_5007:                0.0, $ ; [N II] 6584 / [O III] 5007
      nii_6584_h_beta:                   0.0, $
      oii_sii:                           0.0, $ ; [O II] / [S II] 6716,6731
      sii_h_beta:                        0.0, $
      oii_h_beta_oii_sii:                0.0, $ ; ([O II] / Hb) ([O II] / [S II] 6716,6731)
      oiii_5007_h_beta_nii_6584_h_alpha: 0.0, $ ; Pettini & Pagel 2004
      P:                                 0.0, $ ; {[O III] 4959,5007 / ([O III] 4959,5007 + [O II] 3726,3729)}
      O32:                               0.0, $ ; ([O III] 4959,5007 / [O II] 3726,3729)
      R23:                               0.0}   ; ([O II] 3726,3729 + [O III] 4959,5007) / Hb
    gridstemplate = grids
    grids = replicate(gridstemplate,nZ,nU)

    Uindx = lindgen(nU)
    for i = 0L, nZ-1L do begin
       indx = i*nU + Uindx
       for k = 1L, nratios do grids[i,Uindx].(k-1L) = data[k,indx] ; ignore the first (metallicity) column
    endfor

; correct for doublets

;   grids.nii_oii = grids.nii_oii - logncor
;   grids.oiii_nii = grids.oiii_nii + logncor - logocor
    
; compute new, useful ratios
    
    grids.oii_h_alpha                       = -grids.nii_6584_oii + grids.nii_6584_h_alpha
    grids.oii_h_beta                        =  grids.oii_h_alpha + alog10(2.859)
    grids.nii_6584_sii                      =  grids.nii_6584_h_alpha - grids.sii_h_alpha
    grids.nii_6584_oiii_5007                = -grids.oiii_5007_nii_6584
    grids.nii_6584_h_beta                   =  grids.nii_6584_h_alpha + alog10(2.859)
    grids.oii_sii                           =  grids.oii_h_alpha - grids.sii_h_alpha
    grids.sii_h_beta                        =  grids.sii_h_alpha + alog10(2.859)
    grids.oii_h_beta_oii_sii                =  grids.oii_h_beta + grids.oii_sii
    grids.oiii_5007_h_beta_nii_6584_h_alpha =  grids.oiii_5007_h_beta - grids.nii_6584_h_alpha
    grids.P                                 =  10^(grids.oiii_5007_h_beta+logocor) / (10^(grids.oiii_5007_h_beta+logocor) + 10^grids.oii_h_beta)
    grids.O32                               =  (grids.oiii_5007_h_beta+logocor) - grids.oii_h_beta
    grids.R23                               =  alog10(10^grids.oii_h_beta + (10^(grids.oiii_5007_h_beta+logocor)))

; "undo" the logarithm of the line ratios

    if (n_elements(nolog) ne 0L) then begin

       for itag = 0L, n_tags(grids[0])-1L do grids.(itag) = 10^grids.(itag)
    
    endif
    
    if n_elements(oversamp) ne 0L then begin ; over-sample if requested

; code to use INTERPOLATE(), below
       
;      dZ = 0.025 & nnewZ = fix((2.0-0.1)/dZ+1) & newZ = findgen(nnewZ)*dZ+0.1
;      dU = 0.051 & nnewU = fix((max(U)-min(U))/dU+1) & newU = findgen(nnewU)*dU+min(U)
;      bigZ = rebin(newZ,nnewZ,nnewU) & bigU = rebin(newU,nnewU,nnewZ)
;
;      bigZindx = rebin(findex(Z,newZ),nnewZ,nnewU)
;      bigUindx = transpose(rebin(findex(U,newU),nnewU,nnewZ))
;      newgrids = replicate(gridstemplate,nnewZ,nnewU)

; code to use TRIGRID() or KRIG2D(), below

       bigZ = rebin(Z,nZ,nU)
       bigU = transpose(rebin(U,nU,nZ))
       ZUspacing = [0.05,0.05]
       ZUlimits = [0.1,min(U),2.0,max(U)]
;      ZUlimits = [min(Z),min(U),max(Z),max(U)]
          
       tags = tag_names(gridstemplate)
       ntags = n_tags(gridstemplate)
       
       for itag = 6L, ntags-1L do begin
;      for itag = 0L, ntags-1L do begin

          oldfluxes = grids.(itag)

; TRIGRID() or KRIG2D() code
          
          triangulate, bigZ, bigU, tr, b
          newfluxes = trigrid(bigZ, bigU, oldfluxes, tr, ZUspacing, ZUlimits, $
            xgrid=newZ, ygrid=newU, extra=b, /quintic)
          nnewZ = n_elements(newZ) & nnewU = n_elements(newU)

          if (n_elements(newgrids) eq 0L) then newgrids = replicate(gridstemplate,nnewZ,nnewU)

; INTERPOLATE() or BILINEAR() code
          
;          newfluxes = bilinear(oldfluxes,bigZindx,bigUindx)
;;         newfluxes = interpolate(oldfluxes,bigZindx,bigUindx);,cubic=-0.1)

          newgrids.(itag) = newfluxes

          get_element, newU, U, indx
          for iU = 0L, nU-1L do begin
             plot, Z, alog10(oldfluxes[*,iU]), xsty=3, ysty=3, charsize=2, $
               title=tags[itag]+', U = '+string(U[iU],format='(G0)'), $
               xtitle='Z', ytitle='log Flux Ratio'
             oplot, newZ, alog10(newfluxes[*,indx[iU]]), ps=4, syms=2
             cc = get_kbrd(1)
          endfor

       endfor

       Z = newZ & nZ = n_elements(newZ)
       U = newU & nU = n_elements(newU)
       grids = newgrids
       
    endif

return, grids
end    
