;+
; NAME:
;   READ_MZ_SAMPLE()
;
; PURPOSE:
;   Read the various AGES/MZ samples.
;
; INPUTS: 
;
; KEYWORD PARAMETERS: 
;   sdss - 
;   silent - 
;   parent -
;   mz_ancillary - 
;   mz_ispec - 
;   mzhii_ancillary - 
;   mzhii_ispec - 
;   mzhii_log12oh - 
;   nodust - 
;
; OUTPUTS: 
;
; COMMENTS:
;
; MODIFICATION HISTORY:
;   J. Moustakas, 2009 Mar 19, NYU
;   jm10may08ucsd - generalized and documented
;
; Copyright (C) 2010, John Moustakas
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

function read_mz_sample, sdss=sdss, silent=silent, nodust=nodust, parent=parent, $
  mass=mass, sfrs=sfrs, mz_ancillary=mz_ancillary, mz_mass=mz_mass, mz_ispec=mz_ispec, $
  mzhii_ancillary=mzhii_ancillary, mzhii_mass=mzhii_mass, mzhii_ispec=mzhii_ispec, $
  mzhii_log12oh=mzhii_log12oh, mzagn_ancillary=mzagn_ancillary, mzagn_mass=mzagn_mass, $
  mzagn_ispec=mzagn_ispec, mzagn_log12oh=mzagn_log12oh

    common mz_parent, $
      ages_parent, ages_mass, ages_sfrs, ages_mz_ancillary, ages_mz_mass, ages_mz_ispec, $
      ages_mzhii_ancillary, ages_mzhii_mass, ages_mzhii_ispec, ages_mzhii_log12oh, $
      ages_mzhii_log12oh_nodust, ages_mzagn_ancillary, ages_mzagn_mass, $
      ages_mzagn_ispec, ages_mzagn_log12oh, ages_mzagn_log12oh_nodust, $

      sdss_parent, sdss_mass, sdss_sfrs, sdss_mz_ancillary, sdss_mz_mass, sdss_mz_ispec, $
      sdss_mzhii_ancillary, sdss_mzhii_mass, sdss_mzhii_ispec, sdss_mzhii_log12oh, $
      sdss_mzhii_log12oh_nodust, sdss_mzagn_ancillary, sdss_mzagn_mass, $
      sdss_mzagn_ispec, sdss_mzagn_log12oh, sdss_mzagn_log12oh_nodust

    mzpath = mz_path()
    
; ##################################################
; SDSS
    if keyword_set(sdss) then begin
; see BUILD_MZ_PARENT_SAMPLE
       if keyword_set(parent) then begin
          thisfile = mzpath+'sdss_parent.fits.gz'
          if (size(sdss_parent,/type) ne 8) then begin
             if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
             sdss_parent = mrdfits(thisfile,1,silent=0)
          endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
          return, sdss_parent
       endif 
; see SDSS_ISEDFIT
       if keyword_set(mass) then begin
          thisfile = mzpath+'sdss_parent_mass.fits.gz'
          if (size(sdss_mass,/type) ne 8) then begin
             if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
             sdss_mass = mrdfits(thisfile,1,silent=0)
          endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
          return, sdss_mass
       endif 
; see BUILD_MZ_SFRS
       if keyword_set(sfrs) then begin
          thisfile = mzpath+'sdss_parent_sfrs.fits.gz'
          if (size(sdss_sfrs,/type) ne 8) then begin
             if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
             sdss_sfrs = mrdfits(thisfile,1,silent=0)
          endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
          return, sdss_sfrs
       endif 
; MZ/ancillary - see BUILD_MZ_EMLINE_SAMPLE
       if keyword_set(mz_ancillary) then begin
          thisfile = mzpath+'sdss_mz_ancillary.fits.gz'
          if (size(sdss_mz_ancillary,/type) ne 8) then begin
             if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
             sdss_mz_ancillary = mrdfits(thisfile,1,silent=0)
          endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
          return, sdss_mz_ancillary
       endif
; MZ/mass - see BUILD_MZ_EMLINE_SAMPLE
       if keyword_set(mz_mass) then begin
          thisfile = mzpath+'sdss_mz_mass.fits.gz'
          if (size(sdss_mz_mass,/type) ne 8) then begin
             if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
             sdss_mz_mass = mrdfits(thisfile,1,silent=0)
          endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
          return, sdss_mz_mass
       endif
; MZ/ispec - see BUILD_MZ_EMLINE_SAMPLE
       if keyword_set(mz_ispec) then begin
          thisfile = mzpath+'sdss_mz_ispec.fits.gz'
          if (size(sdss_mz_ispec,/type) ne 8) then begin
             if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
             sdss_mz_ispec = mrdfits(thisfile,1,silent=0)
          endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
          return, sdss_mz_ispec
       endif
; MZ/HII/ancillary - see BUILD_MZ_EMLINE_SAMPLE
       if keyword_set(mzhii_ancillary) then begin
          thisfile = mzpath+'sdss_mz_hii_ancillary.fits.gz'
          if (size(sdss_mzhii_ancillary,/type) ne 8) then begin
             if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
             sdss_mzhii_ancillary = mrdfits(thisfile,1,silent=0)
          endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
          return, sdss_mzhii_ancillary
       endif
; MZ/HII/mass - see BUILD_MZ_EMLINE_SAMPLE
       if keyword_set(mzhii_mass) then begin
          thisfile = mzpath+'sdss_mz_hii_mass.fits.gz'
          if (size(sdss_mzhii_mass,/type) ne 8) then begin
             if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
             sdss_mzhii_mass = mrdfits(thisfile,1,silent=0)
          endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
          return, sdss_mzhii_mass
       endif
; MZ/HII/ispec - see BUILD_MZ_EMLINE_SAMPLE
       if keyword_set(mzhii_ispec) then begin
          thisfile = mzpath+'sdss_mz_hii_ispec.fits.gz'
          if (size(sdss_mzhii_ispec,/type) ne 8) then begin
             if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
             sdss_mzhii_ispec = mrdfits(thisfile,1,silent=0)
          endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
          return, sdss_mzhii_ispec
       endif
; MZ/HII/log12oh - see BUILD_MZ_LOG12OH_SAMPLE
       if keyword_set(mzhii_log12oh) then begin
          if keyword_set(nodust) then begin
             thisfile = mzpath+'sdss_mz_hii_log12oh_nodust.fits.gz'
             if (size(sdss_mzhii_log12oh_nodust,/type) ne 8) then begin
                if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
                sdss_mzhii_log12oh_nodust = mrdfits(thisfile,1,silent=0)
             endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
             return, sdss_mzhii_log12oh_nodust
          endif else begin
             thisfile = mzpath+'sdss_mz_hii_log12oh.fits.gz'
             if (size(sdss_mzhii_log12oh,/type) ne 8) then begin
                if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
                sdss_mzhii_log12oh = mrdfits(thisfile,1,silent=0)
             endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
             return, sdss_mzhii_log12oh
          endelse
       endif 
; MZ/AGN/ancillary - see BUILD_MZ_EMLINE_SAMPLE
       if keyword_set(mzagn_ancillary) then begin
          thisfile = mzpath+'sdss_mz_agn_ancillary.fits.gz'
          if (size(sdss_mzagn_ancillary,/type) ne 8) then begin
             if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
             sdss_mzagn_ancillary = mrdfits(thisfile,1,silent=0)
          endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
          return, sdss_mzagn_ancillary
       endif
; MZ/AGN/mass - see BUILD_MZ_EMLINE_SAMPLE
       if keyword_set(mzagn_mass) then begin
          thisfile = mzpath+'sdss_mz_agn_mass.fits.gz'
          if (size(sdss_mzagn_mass,/type) ne 8) then begin
             if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
             sdss_mzagn_mass = mrdfits(thisfile,1,silent=0)
          endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
          return, sdss_mzagn_mass
       endif
; MZ/AGN/ispec - see BUILD_MZ_EMLINE_SAMPLE
       if keyword_set(mzagn_ispec) then begin
          thisfile = mzpath+'sdss_mz_agn_ispec.fits.gz'
          if (size(sdss_mzagn_ispec,/type) ne 8) then begin
             if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
             sdss_mzagn_ispec = mrdfits(thisfile,1,silent=0)
          endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
          return, sdss_mzagn_ispec
       endif
; MZ/AGN/log12oh - see BUILD_MZ_LOG12OH_SAMPLE
       if keyword_set(mzagn_log12oh) then begin
          if keyword_set(nodust) then begin
             thisfile = mzpath+'sdss_mz_agn_log12oh_nodust.fits.gz'
             if (size(sdss_mzagn_log12oh_nodust,/type) ne 8) then begin
                if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
                sdss_mzagn_log12oh_nodust = mrdfits(thisfile,1,silent=0)
             endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
             return, sdss_mzagn_log12oh_nodust
          endif else begin
             thisfile = mzpath+'sdss_mz_agn_log12oh.fits.gz'
             if (size(sdss_mzagn_log12oh,/type) ne 8) then begin
                if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
                sdss_mzagn_log12oh = mrdfits(thisfile,1,silent=0)
             endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
             return, sdss_mzagn_log12oh
          endelse
       endif 
; ##################################################
; AGES
    endif else begin
; see BUILD_MZ_PARENT_SAMPLE
       if keyword_set(parent) then begin
          thisfile = mzpath+'ages_parent.fits.gz'
          if (size(ages_parent,/type) ne 8) then begin
             if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
             ages_parent = mrdfits(thisfile,1,silent=0)
          endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
          return, ages_parent
       endif 
; see MZ_ISEDFIT
       if keyword_set(mass) then begin
          thisfile = mzpath+'ages_parent_mass.fits.gz'
;         thisfile = isedfitpath+'UBwRIzJHKs_bc03_chab_calzetti_sfhgrid02.fits.gz'
;         thisfile = isedfitpath+'BwRIJHKsirac_bc03_chab_calzetti_sfhgrid02.fits.gz'
          if (size(ages_mass,/type) ne 8) then begin
             if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
             ages_mass = mrdfits(thisfile,1,silent=0)
          endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
          return, ages_mass
       endif 
; see BUILD_MZ_SFRS
       if keyword_set(sfrs) then begin
          thisfile = mzpath+'ages_parent_sfrs.fits.gz'
          if (size(ages_sfrs,/type) ne 8) then begin
             if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
             ages_sfrs = mrdfits(thisfile,1,silent=0)
          endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
          return, ages_sfrs
       endif 
; MZ/ancillary
       if keyword_set(mz_ancillary) then begin
          thisfile = mzpath+'ages_mz_ancillary.fits.gz'
          if (size(ages_mz_ancillary,/type) ne 8) then begin
             if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
             ages_mz_ancillary = mrdfits(thisfile,1,silent=0)
          endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
          return, ages_mz_ancillary
       endif
; MZ/mass
       if keyword_set(mz_mass) then begin
          thisfile = mzpath+'ages_mz_mass.fits.gz'
          if (size(ages_mz_mass,/type) ne 8) then begin
             if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
             ages_mz_mass = mrdfits(thisfile,1,silent=0)
          endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
          return, ages_mz_mass
       endif
; MZ/ispec
       if keyword_set(mz_ispec) then begin
          thisfile = mzpath+'ages_mz_ispec.fits.gz'
          if (size(ages_mz_ispec,/type) ne 8) then begin
             if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
             ages_mz_ispec = mrdfits(thisfile,1,silent=0)
          endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
          return, ages_mz_ispec
       endif
; MZ/HII/ancillary
       if keyword_set(mzhii_ancillary) then begin
          thisfile = mzpath+'ages_mz_hii_ancillary.fits.gz'
          if (size(ages_mzhii_ancillary,/type) ne 8) then begin
             if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
             ages_mzhii_ancillary = mrdfits(thisfile,1,silent=0)
          endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
          return, ages_mzhii_ancillary
       endif
; MZ/HII/mass
       if keyword_set(mzhii_mass) then begin
          thisfile = mzpath+'ages_mz_hii_mass.fits.gz'
          if (size(ages_mzhii_mass,/type) ne 8) then begin
             if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
             ages_mzhii_mass = mrdfits(thisfile,1,silent=0)
          endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
          return, ages_mzhii_mass
       endif
; MZ/HII/ispec
       if keyword_set(mzhii_ispec) then begin
          thisfile = mzpath+'ages_mz_hii_ispec.fits.gz'
          if (size(ages_mzhii_ispec,/type) ne 8) then begin
             if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
             ages_mzhii_ispec = mrdfits(thisfile,1,silent=0)
          endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
          return, ages_mzhii_ispec
       endif
; MZ/HII/log12oh - see BUILD_MZ_LOG12OH_SAMPLE
       if keyword_set(mzhii_log12oh) then begin
          if keyword_set(nodust) then begin
             thisfile = mzpath+'ages_mz_hii_log12oh_nodust.fits.gz'
             if (size(ages_mzhii_log12oh_nodust,/type) ne 8) then begin
                if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
                ages_mzhii_log12oh_nodust = mrdfits(thisfile,1,silent=0)
             endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
             return, ages_mzhii_log12oh_nodust
          endif else begin
             thisfile = mzpath+'ages_mz_hii_log12oh.fits.gz'
             if (size(ages_mzhii_log12oh,/type) ne 8) then begin
                if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
                ages_mzhii_log12oh = mrdfits(thisfile,1,silent=0)
             endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
             return, ages_mzhii_log12oh
          endelse
       endif 
; MZ/AGN/ancillary
       if keyword_set(mzagn_ancillary) then begin
          thisfile = mzpath+'ages_mz_agn_ancillary.fits.gz'
          if (size(ages_mzagn_ancillary,/type) ne 8) then begin
             if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
             ages_mzagn_ancillary = mrdfits(thisfile,1,silent=0)
          endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
          return, ages_mzagn_ancillary
       endif
; MZ/AGN/mass
       if keyword_set(mzagn_mass) then begin
          thisfile = mzpath+'ages_mz_agn_mass.fits.gz'
          if (size(ages_mzagn_mass,/type) ne 8) then begin
             if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
             ages_mzagn_mass = mrdfits(thisfile,1,silent=0)
          endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
          return, ages_mzagn_mass
       endif
; MZ/AGN/ispec
       if keyword_set(mzagn_ispec) then begin
          thisfile = mzpath+'ages_mz_agn_ispec.fits.gz'
          if (size(ages_mzagn_ispec,/type) ne 8) then begin
             if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
             ages_mzagn_ispec = mrdfits(thisfile,1,silent=0)
          endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
          return, ages_mzagn_ispec
       endif
; MZ/AGN/log12oh - see BUILD_MZ_LOG12OH_SAMPLE
       if keyword_set(mzagn_log12oh) then begin
          if keyword_set(nodust) then begin
             thisfile = mzpath+'ages_mz_agn_log12oh_nodust.fits.gz'
             if (size(ages_mzagn_log12oh_nodust,/type) ne 8) then begin
                if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
                ages_mzagn_log12oh_nodust = mrdfits(thisfile,1,silent=0)
             endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
             return, ages_mzagn_log12oh_nodust
          endif else begin
             thisfile = mzpath+'ages_mz_agn_log12oh.fits.gz'
             if (size(ages_mzagn_log12oh,/type) ne 8) then begin
                if (keyword_set(silent) eq 0) then splog, 'Reading '+thisfile
                ages_mzagn_log12oh = mrdfits(thisfile,1,silent=0)
             endif else if (keyword_set(silent) eq 0) then splog, 'Restoring '+file_basename(thisfile)
             return, ages_mzagn_log12oh
          endelse
       endif 
    endelse
end
