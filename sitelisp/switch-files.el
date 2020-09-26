;; swich-files.el: a method of switching between matched pairs of
;; files and for following include directives.
;;
;; Copyright (C) 2002-2004  Wes Hardaker <elisp@hardakers.net>
;;
;; Modified by Tom 7, 19 Dec 2009.
;;
;; This program is free software; you can redistribute it and/or modify
;; it under the terms of the GNU General Public License as published by
;; the Free Software Foundation; either version 2, or (at your option)
;; any later version.
;;
;; This program is distributed in the hope that it will be useful,
;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;; GNU General Public License for more details.
;;
;; A copy of the GNU General Public License can be obtained from this
;; program's author (send electronic mail to psmith@BayNetworks.com) or
;; from the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA
;; 02139, USA.
;;
;; $Revision: 1.8 $

(defvar switch-files-paths '("." "/usr/include" "/usr/local/include")
  "the list of paths to look through for matching files.")

;; XXX there need to be multiple possible destinations for
;; certain suffixes, for example .h (.c, .cc, .cpp, etc.).
;; just return a list below.
(defvar switch-files-list '(
			    ("\\.cc" ".h")
			    ("\\.h" ".cc")

                            ("\\.c" ".h")
			    ("\\.h" ".c")

			    ("\\.C" ".H")
			    ("\\.H" ".C")

			    ("\\.cpp" ".h")
			    ("\\.h" ".cpp")

			    ("-sig\\.sml" ".sml")
			    ("\\.sml" "-sig.sml")
			    )
  "A list of file regexp matches and replacements switch between.")

(defvar switch-file-backto nil
  "A buffer local variable containing a buffer to jump back to for switch-file.")

(require 'cl)
(require 'seq)

(defun switch-files ()
  (interactive)
  (let* ((curbuffer (current-buffer))
	 (buffername (buffer-file-name))
	 (pathlist switch-files-paths)
	 ;; single file based on history or what we're looking at
	 (already
	  ;; are we staring at an include directive?
	  ;; XXX do I need this? I don't think I ever use C-o for this. -tom7
	  (or
	   (and
	    (looking-at
	     "#include [<\\\"]\\([^/>\\\"]*\\|.*/[^/>\\\"]*\\)[>\\\"]")
	    (buffer-substring (match-beginning 1) (match-end 1)))
           ;; are we in a buffer that we know where to go back to already?
	   ;; XXX if so, just prioritizing it as the first element might make
	   ;; more sense thatn returning a singleton list...
	   switch-file-backto))

	   ;; can we guess at a file name from the current buffer name?

	 
	 ; list of files to try switching to, in priority order
	 (startfile
	  (if already
	      ;; if we already have a single file, return that
	      (list already)
	    ;; otherwise, get all associations that match the buffer
	    ;; (XXX more convoluted than it needs to be. we want mapPartial, but
	    ;; instead we first map to the filename or nil, then filter out nils.
	    ;; note that due to the history of this code, we need to do the
	    ;; transform right after the match (uses regex match state))
	    (let* ((matches-with-nils
		    (seq-map
		     (function (lambda (row) 
				 ;; (message "regex %s" (car key))
				 (if (string-match
				      (concat "\\(.*\\)" (car row))
				      buffername)
				     ;; matched
				     (file-name-nondirectory
				      (concat (substring buffername (match-beginning 1)
							 (match-end 1))
					      (cadr row)))
				   ;; not matched
				   nil)))
		     switch-files-list))
		   (targets (seq-filter
			     (function (lambda (elt) (not (null elt))))
			     matches-with-nils)))

	      (message "%s" targets)
	      ;; XXX no return the list!
	      (if (null targets) nil (car targets)))))
	 )

    (if (not startfile)
	(message "Can't switch. Buffer '%s' name '%s' startfile '%s'" 
		 curbuffer buffername startfile)
      (progn
	(if (bufferp startfile)
	    (switch-to-buffer startfile)
	  (setq done nil)
	  (while (and (not done) pathlist)
	    (setq thefile (expand-file-name startfile (car pathlist)))
	    (if (file-exists-p thefile)
		(setq done t)
	      (setq pathlist (cdr pathlist))))

	  (if done
	      (progn
		(find-file thefile)
		(make-local-variable 'switch-file-backto)
		(setq switch-file-backto curbuffer))
	    (message "Can't find file: %s" startfile)
	    ))))
    )
  )

(global-set-key "\C-x\M-f" 'switch-files)

; (switch-files)

(provide 'switch-files)
