function rep2d,mask,trans,noaves=noaves 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;                                                            ;;
;; Purpose: Selects a representative point for each numbered  ;;
;;          point cluster in mask.                            ;;
;;                                                            ;;
;; Keywords: noaves - don't use average cluster values as rep-;;
;;             resentative values in output array             ;;
;;                                                            ;;
;; Created: 01/29/08                                          ;;
;;                                                            ;;
;; Created by: Shaela Jones                                   ;;
;;                                                            ;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


; initializations
sztrans=SIZE(trans)
nclusts=MAX(mask)
edit1=FLTARR(sztrans[1],sztrans[2])



; main
for i=1,nclusts do begin
  locs=WHERE(mask eq i)
  if locs[0] eq -1 then begin
    print,'Problem with cluster_mask resuilt'
    return,-1
  endif
  if KEYWORD_SET(noaves) then ave=1. else ave=AVERAGE(trans[locs])
  ind=ARRAY_INDICES(trans,locs)
  minthet=MIN(ind[0,*],max=maxthet)
  minrho=MIN((ind[1,*]),max=maxrhp)
  subarray=trans[minthet:maxthet,minrho:maxrho]
  submask=mask[minthet:maxthet,minrho:maxrho]
  nthet=maxthet-minthet+1
  locs=WHERE(submask ne i)
  if locs[0] ne -1 then subarray[locs]=0.
  szsub=SIZE(subarray)
;stop
;  if szsub[0] gt 1 then cent=CENTROID(subarray) else cent=[0,0]
  if szsub[0] gt 1 then maxpos=ARRAY_INDICES(subarray, $
      WHERE(subarray eq MAX(subarray))) else maxpos=[0,0]
  edit1[maxpos[0]+minthet,maxpos[1]+minrho]=ave
endfor


return,edit1
end
