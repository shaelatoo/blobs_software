pro adj_pix,intensities,times,degree,adjust,fit

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                            ;;
;; Purpose: Fits the provided intensities and times to a poly-;;
;;           nomial for "background" subtraction.             ;;
;;                                                            ;;
;; Inputs: intensities - 1D array of pixel intensities        ;;
;;         times - corresponding TAI-formatted date-time stamp;;
;;         degree - degree of polynomial data is fit to       ;;
;;                                                            ;;
;; Outputs: adjust - fit values at times given in times array ;;
;;          fit - array of four fit parameters                ;;
;;                                                            ;;
;; Keywords: none                                             ;;
;;                                                            ;;
;; Dependencies: none                                         ;;
;;                                                            ;;
;; Created: 07/17/07                                          ;;
;;                                                            ;;
;; Created by: Shaela Jones                                   ;;
;;                                                            ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; parameters
n_sig=3.
resfactor=10


; initializations
num=N_ELEMENTS(intensities)
x=times-times[0]
mag=MEAN(intensities)
if mag lt 10.^(-8) then lowmag=1 else lowmag=0
if lowmag eq 1 then y=reform(intensities)*10.^9 $
   else y=REFORM(intensities)



; main
list=WHERE(FINITE(y))
if list[0] ne -1 then begin
  xp=x[list]
  yp=y[list]
  fit=poly_fit(xp,yp,degree,status=stat,yfit=adjust)
  if stat ne 0 then stop
  res=ABS(yp-adjust)
  list=WHERE(res gt n_sig*MEAN(res),nlist,complement=comp)
  if nlist gt 1 and nlist lt num/resfactor then begin
    fit=POLY_FIT(xp[comp],yp[comp],degree,status=stat,adjust)
    if stat ne 0 then stop
  endif
  if lowmag eq 1 then fit=fit*10.^(-9)
  adjust=fit[0]
  for i=1,degree do adjust=adjust+fit[i]*x^FLOAT(i)
;  adjust=fit[0]+fit[1]*x+fit[2]*x^2+fit[3]*x^3
endif



end


; obsolete

;  adjust=fit[0]+fit[1]*x+fit[2]*x^2+fit[3]*x^3
