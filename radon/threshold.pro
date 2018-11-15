pro threshold,trans,nmean,szstruc,fun2=fun2,fun3=fun3,fun4=fun4, $
   fun5=fun5,fun6=fun6,fun7=fun7

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                            ;;
;; Purpose: Thresholds Radon transform several different ways ;;
;;                                                            ;;
;; Inputs: trans - transform (2D or 3D)                       ;;
;;         nmean - number of means at which to threshold the  ;;
;;           transform                                        ;;
;;         szstruc - size of morphological closing structure  ;;
;;                                                            ;;
;; Keywords: fun2 - transform thresholded at median           ;;
;;           fun3 - fun2 thresholded at nmean*median          ;;
;;           fun4 - fun3 closed                               ;;
;;           fun5 - smoothed fun4                             ;;
;;           fun6 - fun3 closed greyscale                     ;;
;;           fun7 - smoothed fun6                             ;;
;;                                                            ;;
;; Created: 04/29/08                                          ;;
;;                                                            ;;
;; Created by: Shaela Jones                                   ;;
;;                                                            ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

; make closing structure
vals=DIST(szstruc)
struc=vals
struc[where(struc ge 1.)]=1.
struc[where(struc lt 1.)]=0.




fun2=trans
fun2[where(fun2 le MEDIAN(fun2))]=0.
fun3=fun2
fun3[where(fun2 le nmean*MEAN(fun2))]=0.
fun4=MORPH_CLOSE(fun3,struc)
fun5=SMOOTH(fun4,2)
fun6=MORPH_CLOSE(fun3,struc,/gray,values=vals)
fun7=SMOOTH(fun6,2)


end
