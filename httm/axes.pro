pro axes,dttms,rmin,rwidth,firstfile,ctr,ang,outsz,tms,alts, $
               cadence,rotate=rotate,interp_r=interp_r

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                            ;;
;; Purpose: Makes height-time axes for r_strp images.         ;;
;;                                                            ;;
;; Inputs: dttms - array holding date/time of each image      ;;
;;         rmin - minimum radius used, in pixels              ;;
;;         rwidth - max radius used, in pixels                ;;
;;         firstfile - filename of first image in sequence    ;;
;;         ctr - center point from which to calculate radius  ;;
;;         ang - angular position of radial strip             ;;
;;         outsz - size of images, in pixels                  ;;
;;                                                            ;;
;; Outputs: tms - times of images in seconds                  ;;
;;          alts - altitudes of radial points, in solar radii ;;
;;                                                            ;;
;; Optional Outputs: cadence - cadence of image sequence,     ;;
;;                             according to firstfile header  ;;
;;                                                            ;;
;; Keywords: rotate - if set, keywords rotate_on and          ;;
;;                    rotinterp_on were used in prepping      ;;
;;                    image files                             ;;
;;           interp_r - number of radial points interpolated  ;;
;;                      by r_strp                             ;;
;;                                                            ;;
;; Dependencies: secchi_prep                                  ;;
;;                                                            ;;
;; Created: 10/15/07                                          ;;
;;                                                            ;;
;; Created by: Shaela Jones                                   ;;
;;                                                            ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; error check
if N_ELEMENTS(outsz) eq 0 then begin
  print,'Problem with variable outsz.'
  return
endif




;initializations
if KEYWORD_SET(interp_r) then proxy=interp_r else proxy=rwidth


; set time axis
tms=ANYTIM2TAI(dttms)



; get position info from image file header
if N_ELEMENTS(wcshead) eq 0 then begin
  if keyword_set(rotate) then begin
    secchi_prep,firstfile,hdr,im,/calfac_off,/calimg_off,/silent, $
               outsize=outsz,/smask_on,/rotate_on,/rotinterp_on
  endif else begin
    secchi_prep,firstfile,hdr,im,/calfac_off,/calimg_off,/silent, $
               outsize=outsz,/smask_on
  endelse
  struct=FITSHEAD2WCS(hdr)
  convert=hdr.dsun_obs/206265.0/1000.0/696000.  
endif else begin
  struct=wcshead
  convert=struct.position.dsun_obs/206265.0/1000.0/696000.
endelse



; get image cadence
cadence=hdr.cadence
;if cadence eq 0 then 



;make image axes
radius=rmin+FINDGEN(proxy)*rwidth/proxy
;x=FLTARR(2,proxy)
;x[0,*]=ctr[0]+radius*COS(ang*!DTOR)
;x[1,*]=ctr[1]+radius*SIN(ang*!DTOR)
polar=TRANSPOSE([[REPLICATE(ang,proxy)],[radius]])
x=CV_COORD(from_polar=polar,/to_rect,/degrees)
x[0,*]=x[0,*]+ctr[0]
x[1,*]=x[1,*]+ctr[1]
pos=WCS_GET_COORD(struct,x)*convert
pol=CV_COORD(from_rect=pos,/degrees,/to_polar)
alts=pol[1,*]



end



; obsolete

;date=strmid(dttms[0],0,4)+strmid(dttms[0],5,2) $
;      +strmid(dttms[0],8,2)
;time=strmid(dttms[0],11,2)+strmid(dttms[0],14,2) $
;      +strmid(dttms[0],17,2)
