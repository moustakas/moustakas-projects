pro nfgs_specfit, integrated=integrated
; jm04jan6uofa

; some notes:

; CR in the blue, but does not affect the fitting, just the plotting:
;
;       A01187-0048, NGC825

; low S/N:
;
;       A02493-0122, NGC4034, NGC5356
;

; possibly broader than the instrumental resolution: 
;       IC195, NGC193
;

;   otherbad = [$
;     'NGC3213',     $          ; the red continuum is very wrong
;     'NGC5267',     $          ; funny red emission lines
;     'A08567+5242', $          ; funny red emission lines
;     'A10389+3859', $          ; funny red emission lines
;     'A13281+3153', $          ; funny red emission lines
;     'NGC695'       $          ; the continuum is weird
;     ]
      
; blatant AGN

;   agn = [$
;     'MRK421',$          ; BL Lac [A11017+3828W]
;     'A00510+1225', $    ; Sy1
;     'A12195+7535', $    ; Sy1
;     'A15016+1037'  $    ; Sy1
;     ]
    
    datapath = nfgs_path(/spec1d)
    specfitpath = nfgs_path(/specfit)
    nfgs = mrdfits(nfgs_path(/analysis)+'nfgs_info.fits.gz',1,/silent)

    if keyword_set(integrated) then begin
       
       drift = where(nfgs.drift and (nfgs.drift_agnflag eq 0L))

       suffix = 'nfgs_int'
       speclist = strtrim(nfgs[drift].drift_file,2)

       specdata = ispeclinefit(speclist,specres=6.0,snrcut=0.0,dustmodel=0,$ ; dustmodel=3,$
         datapath=datapath,linepath=specfitpath,suffix=suffix,Zmulti=Zmulti,$
         /charlot,/zcrosscor,/postscript,/write,vmaxshift=vmaxshift,/nologfile,$
         starvdisp=100.0,nback=0L)

       nfgs_parse_specfit, datapath=specfitpath

    endif
    
return
end    
