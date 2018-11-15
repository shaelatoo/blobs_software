
function enlarge_image2,image,factor1,factor2

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                      ;;
;; Purpose: Enlarges the input image by factor1 in the  ;;
;;           horizontal and factor2 in the vertical     ;;
;;           directions.                                ;;
;;                                                      ;;
;; Inputs: image - 2D image to be enlarged              ;;
;;         factor1 - horizontal enlargement factor      ;;
;;         factor2 - vertical enlargement factor        ;;
;;                                                      ;;
;; Notes: factor1, factor2 must be ge 1                 ;;
;;                                                      ;;
;; Created by: Shaela Jones                             ;;
;;                                                      ;;
;; Created: ?                                           ;;
;;                                                      ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

szimage=size(image)
if szimage[0] ne 2 then begin
  print,'Input image is not a 2D array.'
  return,-1
endif



sz=szimage
new_image=fltarr(sz[1]*factor1,sz[2]*factor2)

for i=0,sz[1]-1 do begin
  for j=0,sz[2]-1 do begin
    new_image[i*factor1:(i+1)*factor1-1,j*factor2:$
                    (j+1)*factor2-1]=image[i,j]
  endfor
endfor


return,new_image

end
