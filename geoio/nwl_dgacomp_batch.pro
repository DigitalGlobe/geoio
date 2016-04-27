;==========================================================
;Script to call idl DGAComp from the command line - did 
;this so that I could put it in Python subprocess.
;
;This script is called like:
;
;idl -e ".run nwl_batch.pro" -args '/path/to/DGAComp'
; '/path/to/in_file' '/path/to/out_file'
; '/path/to/ms_aod_map'
;
;
;The args are all positional without call names. 
;If 'ms_aod_map' is passed, it is applied to 'in_file'.
;This is useful for runnin DGAComp on pan or swir.
;============================================================

compile_opt idl2

;pro nwl_batch

; Grab the command line arguments
args_in=COMMAND_LINE_ARGS(count=nargs)

;;;;;;;;;;;;; Start ENVI ;;;;;;;;;;;;;; 
;;; These are the commands for a headless setup of envi 5+
;e = ENVI(/HEADLESS)
;;; These are the commands for a headless setup of envi classic
; For some reason the headless version causes acomp to have calculation problems
ENVI ;Run this for normal windowed envi
;ENVI, /RESTORE_BASE_SAVE_FILES
;ENVI_BATCH_INIT

print, '###################################################'
print, '### Ignore all the errors that are about to pop ###'
print, '###################################################'

; Set error catch to hack a double compile
CATCH, Error_status
IF Error_status NE 0 THEN BEGIN
  RESOLVE_ROUTINE, 'dgacomp_v023', /COMPILE_FULL_FILE, /EITHER
  CATCH, /CANCEL
ENDIF

; Set the path to compile from and the routine
; then do a double compile with the catch defined above.
resolve_path = args_in[0]
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;The DGAComp pro file MUST be in lower case
; i.e. DGAComp_v023.pro -> dgacomp_v023.pro
; -- just manually rename it if you have a
; version with upper case, everything else
; should work.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
CD, resolve_path
RESOLVE_ROUTINE, 'dgacomp_v023', /COMPILE_FULL_FILE, /EITHER

print, '###################################################'
print, '### Ignore all the errors that just popped ###'
print, '###################################################'

;inFile=swirFile
;outFile=swirOutFile
;aodFile=msAodFile
if nargs EQ 3 then begin
  inFile=args_in[1]
  ;print, inFile
  outFile=args_in[2]
  ;print, outFile
  dgacomp_v023, IN_FILE=inFile, OUT_FILE=outFile
endif else if nargs EQ 4 then begin
  inFile=args_in[1]
  ;print, inFile
  outFile=args_in[2]
  ;print, outFile
  aodFile=args_in[3]
  ;print, aodFile
  dgacomp_v023, IN_FILE=inFile, OUT_FILE=outFile, MS_AOD_MAP=aodFile
endif else stop

; Exit ENVI
;ENVI
ENVI_BATCH_EXIT, /EXIT_IDL, /NO_CONFIRM ;Force in case preferences are set oddly

end
