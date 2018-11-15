function centroid,array 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                            ;;
;; Purpose: Calculates the center of mass (COM) of an         ;;
;;          array.                                            ;;
;;                                                            ;;
;; Created: ?                                                 ;;
;;                                                            ;;
;; Created by: David Foster                                   ;;
;;                                                            ;;
;; Modified by: Shaela Jones                                  ;;
;;              04/29/08 - generalized to n dimensions        ;;
;;                                                            ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

s=SIZE(array,/dimensions)
ndim=N_ELEMENTS(s)
inds=FLTARR(ndim)
totalmass=TOTAL(array)
for i=0,ndim-1 do begin
  temp=array
  for j=0,ndim-1 do begin
    if j ne i then temp=TOTAL(temp,j)
  endfor
  inds[i]=TOTAL(temp*INDGEN(s[i]))/totalmass
endfor
RETURN,inds
END


; obsolete:

; xcm=TOTAL(TOTAL(array,2)*INDGEN(s[0]))/totalmass
; ycm=TOTAL(TOTAL(array,1)*INDGEN(s[1]))/totalmass
; RETURN,[xcm,ycm]


