function read_mockages_sample, evolve=evolve, zbin1=zbin1, $
  zbin2=zbin2, zbin3=zbin3, zbin4=zbin4, zbin5=zbin5, zbin6=zbin6
; jm09mar19nyu - read the output from BUILD_MOCKAGES_SAMPLE

    common mz_mockages, zbin1_evolve, zbin2_evolve, zbin3_evolve, $
      zbin4_evolve, zbin5_evolve, zbin6_evolve, zbin1_noevolve, $
      zbin2_noevolve, zbin3_noevolve, zbin4_noevolve, zbin5_noevolve, $
      zbin6_noevolve

; note: the evolution parameters must match BUILD_MOCKAGES
    datapath = ages_path(/projects)+'mz/mockages/'

    if keyword_set(evolve) then begin
       suffix = 'q1.50-a0.00' 
; ZBIN1
       if keyword_set(zbin1) then begin
          thisfile = datapath+'mockages_zbin1_'+suffix+'.fits.gz'
          if (size(zbin1_evolve,/type) ne 8) then begin
             if (not keyword_set(silent)) then splog, 'Reading '+thisfile
             zbin1_evolve = mrdfits(thisfile,1,silent=0)
          endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
          return, zbin1_evolve
       endif 
; ZBIN2
       if keyword_set(zbin2) then begin
          thisfile = datapath+'mockages_zbin2_'+suffix+'.fits.gz'
          if (size(zbin2_evolve,/type) ne 8) then begin
             if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
             zbin2_evolve = mrdfits(thisfile,1,silent=0)
          endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
          return, zbin2_evolve
       endif 
; ZBIN3
       if keyword_set(zbin3) then begin
          thisfile = datapath+'mockages_zbin3_'+suffix+'.fits.gz'
          if (size(zbin3_evolve,/type) ne 8) then begin
             if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
             zbin3_evolve = mrdfits(thisfile,1,silent=0)
          endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
          return, zbin3_evolve
       endif 
; ZBIN4
       if keyword_set(zbin4) then begin
          thisfile = datapath+'mockages_zbin4_'+suffix+'.fits.gz'
          if (size(zbin4_evolve,/type) ne 8) then begin
             if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
             zbin4_evolve = mrdfits(thisfile,1,silent=0)
          endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
          return, zbin4_evolve
       endif 
; ZBIN5
       if keyword_set(zbin5) then begin
          thisfile = datapath+'mockages_zbin5_'+suffix+'.fits.gz'
          if (size(zbin5_evolve,/type) ne 8) then begin
             if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
             zbin5_evolve = mrdfits(thisfile,1,silent=0)
          endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
          return, zbin5_evolve
       endif 
; ZBIN6
       if keyword_set(zbin6) then begin
          thisfile = datapath+'mockages_zbin6_'+suffix+'.fits.gz'
          if (size(zbin6_evolve,/type) ne 8) then begin
             if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
             zbin6_evolve = mrdfits(thisfile,1,silent=0)
          endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
          return, zbin6_evolve
       endif 
    endif else begin
       suffix = 'noevol'
; ZBIN1
       if keyword_set(zbin1) then begin
          thisfile = datapath+'mockages_zbin1_'+suffix+'.fits.gz'
          if (size(zbin1_noevolve,/type) ne 8) then begin
             if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
             zbin1_noevolve = mrdfits(thisfile,1,silent=0)
          endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
          return, zbin1_noevolve
       endif 
; ZBIN2
       if keyword_set(zbin2) then begin
          thisfile = datapath+'mockages_zbin2_'+suffix+'.fits.gz'
          if (size(zbin2_noevolve,/type) ne 8) then begin
             if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
             zbin2_noevolve = mrdfits(thisfile,1,silent=0)
          endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
          return, zbin2_noevolve
       endif 
; ZBIN3
       if keyword_set(zbin3) then begin
          thisfile = datapath+'mockages_zbin3_'+suffix+'.fits.gz'
          if (size(zbin3_noevolve,/type) ne 8) then begin
             if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
             zbin3_noevolve = mrdfits(thisfile,1,silent=0)
          endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
          return, zbin3_noevolve
       endif 
; ZBIN4
       if keyword_set(zbin4) then begin
          thisfile = datapath+'mockages_zbin4_'+suffix+'.fits.gz'
          if (size(zbin4_noevolve,/type) ne 8) then begin
             if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
             zbin4_noevolve = mrdfits(thisfile,1,silent=0)
          endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
          return, zbin4_noevolve
       endif 
; ZBIN5
       if keyword_set(zbin5) then begin
          thisfile = datapath+'mockages_zbin5_'+suffix+'.fits.gz'
          if (size(zbin5_noevolve,/type) ne 8) then begin
             if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
             zbin5_noevolve = mrdfits(thisfile,1,silent=0)
          endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
          return, zbin5_noevolve
       endif 
; ZBIN6
       if keyword_set(zbin6) then begin
          thisfile = datapath+'mockages_zbin6_'+suffix+'.fits.gz'
          if (size(zbin6_noevolve,/type) ne 8) then begin
             if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
             zbin6_noevolve = mrdfits(thisfile,1,silent=0)
          endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
          return, zbin6_noevolve
       endif 
    endelse
    
end
