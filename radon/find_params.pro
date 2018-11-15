pro find_params,httm,xs,ys,param1,param2,param3

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                            ;;
;; Purpose: Calculates an optimal array of theta values for   ;;
;;          the Radon transform of a given height-time image. ;;
;;                                                            ;;
;; Inputs: httm - the original height-time image to be trans- ;;
;;           formed                                           ;;
;;         xs,ys - x,y coordinates of httm image axes         ;;
;;                                                            ;;
;; Outputs: param1 - recommended thetas for given image, axes ;;
;;          param2 - recommended rs for given image, axes     ;;
;;                                                            ;;
;; Optional Outputs: param3 - recommended as for given image, ;;
;;                     axes                                   ;;
;;                                                            ;;
;; Created: 04/06/08                                          ;;
;;                                                            ;;
;; Created by: Shaela Jones                                   ;;
;;                                                            ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; parameters
nas=50
maxa=5./(696000.*2.)



; input checks
szhttm=SIZE(httm)
if szhttm[0] ne 2 then begin
  print,'Variable httm must be a 2D image array.'
  param1=-1
  param2=-1
  return
endif
szxs=SIZE(xs)
szys=SIZE(ys)
if szxs[0] ne 1 or szys[0] ne 1 then begin
  print,'Variables xs and ys must be 1-D vectors.'
  param1=-1
  param2=-1
  return
endif
if szys[1] ne szhttm[2] or szxs[1] ne szhttm[1] then begin
  print,'Variables xs and ys must have same length as ', $
            'height-time image horizontal and vertical lengths.'
  param1=-1
  param2=-1
  return
endif



; initializations
yrange=szhttm[2]
xrange=szhttm[1]
maxx=MAX(xs,min=minx)
maxy=MAX(ys,min=miny)
;rexs=(xs-minx)/(maxx-minx)-0.5
;reys=(ys-miny)/(maxy-miny)-0.5
thetas=FLTARR(xrange)
deltx=maxx-minx
delty=maxy-miny



; main
if yrange gt xrange then deltr=1./yrange else deltr=1./xrange
for i=1,xrange do thetas[i-1]=ATAN(yrange,i)
thetas=thetas+!dpi/2
param1=REVERSE(thetas)
nrs=CEIL(1/deltr/SQRT(2.))
if nrs MOD 2 ne 0 then nrs++
temp=(INDGEN(nrs)+1)*deltr
rs=[REVERSE(-temp),0.,temp]
param2=rs
if N_PARAMS() eq 6 then begin
;  nas=FLOOR(maxamax/ares)
  maxap=maxa*(yrange/delty)*(deltx^2/xrange^2)
  ares=maxap/nas
  temp=(INDGEN(nas)+1)*ares
  param3=[REVERSE(-temp),temp]
endif


end



; obsolete

;for i=1,xrange-1 do thetas[i-1]=ATAN(yrange,i)
;for i=0,yrange-1 do thetas[i+xrange-1]=ATAN(yrange-i,xrange)
;for i=0,xrange-2 do thetas[i]=ATAN(yrange,i+1)

;thetas=FLTARR(yrange+xrange-1)
