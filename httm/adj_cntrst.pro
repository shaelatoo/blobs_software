function adj_cntrst,image,zeromean=zeromean

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                            ;;
;; Purpose: Adjusts pixel intensities in image such that each ;;
;;          horizontal row has the same mean value.           ;;
;;                                                            ;;
;; Inputs: image - 2D image array to be adjusted              ;;
;;                                                            ;;
;; Keywords: zeromean - if set, mean of each row is scaled to ;;
;;             zero; if not, mean is scaled to image mean     ;;
;;                                                            ;;
;; Created: 09/18/07                                          ;;
;;                                                            ;;
;; Created by: Shaela Jones                                   ;;
;;                                                            ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


; input checks
szimage=SIZE(image)
if szimage[0] ne 2 then begin
  print,'Variable image must be a two-dimensional array.'
  return,-1
endif



; adjust mean row-by-row
sig=FLTARR(szimage[2])
mn=AVERAGE(image,1)
locs=WHERE(mn ne 0,complement=locs0)
temp=image
szimage=SIZE(temp)
avgmn=MEAN(mn)
for i=0,szimage[2]-1 do temp[*,i]=temp[*,i]*avgmn/mn[i]
if KEYWORD_SET(zeromean) then begin
  mn=MEAN(temp)
  temp=temp-mn
endif


return, temp
end
     

; obsolete

;if locs0[0] ne -1 then begin
;  if N_ELEMENTS(locs0) gt 1 OR locs0[0] ne $
;             szimage[2]-1 then begin 
;    print,'Input image has more than one blank row.'
;    return,-1
;  endif else begin
;    mn=mn[0:szimage[2]-2]
;    temp=image[*,0:szimage[2]-2]
;  endelse
;endif else temp=image
