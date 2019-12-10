(in-package :membrane-packer)

(defun build-num-lipids (solute num-lipids &key (input-bounding-box (chem:make-bounding-box '(60.0 60.0 60.0)))
                                             (lipid-selector (list (cons 1.0 *popc*)))
                                             (z-height *z-height*)
                                             (lipid-radius *lipid-radius*))
  (let (top-lipids bottom-lipids)
    (let* ((lipid-selector (make-optimized-lipid-selector lipid-selector))
           (nonbond-db (chem:compute-merged-nonbond-force-field-for-aggregate solute))
           (x-width (chem:get-x-width input-bounding-box))
           (y-width (chem:get-y-width input-bounding-box))
           (z-width (chem:get-z-width input-bounding-box))
           (z-offset *z-space*)
           (half-x-width (/ x-width 2.0))
           (half-y-width (/ y-width 2.0))
           (lipid-id 0))
      (let* ((hex-angle (* 0.0174533 60.0)) ; 60 degrees in radians
             (xstep (* 2.0 lipid-radius))
             (xdir (geom:v* (geom:vec 1.0 0.0 0.0) xstep))
             (ypdir (geom:v* (geom:vec (cos hex-angle) (sin hex-angle) 0.0) xstep))
             (ymdir (geom:v* (geom:vec (- (cos hex-angle)) (sin hex-angle) 0.0) xstep))
             (ystep (* (sin hex-angle) xstep))
             (xnum (floor (/ x-width xstep)))
             (ynum (let ((simple-ynum (floor (/ y-width ystep))))
                     (if (evenp simple-ynum) simple-ynum (1+ simple-ynum)))) ; need even number of cells in y dir
             (x-width (* xnum xstep))
             (y-width (* ynum ystep))
             (xstart (- (* x-width 0.5)))
             (ystart (- (* y-width 0.5)))
             (lipid-id 0)
             (bounding-box (chem:make-bounding-box (list x-width y-width z-width)
                                                   :angles-degrees (chem:get-bounding-box-angles-degrees input-bounding-box)
                                                   :center (chem:get-bounding-box-center input-bounding-box))))
        (let (top-solute-colliders
              bottom-solute-colliders)
          (loop for yindex from (/ ynum 2) to (/ ynum 2)
                do (loop for xindex from (- (/ xnum 2) (/ num-lipids 2)) below (+ (/ xnum 2) (/ num-lipids 2))
                         for wrapped-xypos = (wrapped-xypos xindex yindex xstart ystart xdir ypdir ymdir xnum)
                         for top-ga-lipid = (make-random-ga-lipid wrapped-xypos z-offset lipid-selector (incf lipid-id) :top)
                         do (push top-ga-lipid top-lipids)))
          (format t "There are ~a top-lipids~%" (length top-lipids))
          (let* ((lipids (make-array (+ (length top-lipids) (length bottom-lipids))
                                     :initial-contents (append top-lipids bottom-lipids)))
                 (membrane (make-instance 'ga-membrane
                                          :bounding-box bounding-box
                                          :z-height z-height
                                          :lipids lipids
                                          :lipid-selector lipid-selector)))
            #+(or)(optimize-lipid-placement membrane)
            membrane
            ))))))


(defun rigid-body-coordinates-from-lipid-transforms (membrane)
  (let ((rigid-body-coords (make-array (* 7 (length (lipids membrane))) :element-type 'double-float)))
    (loop for index from 0 by 7
          for lipid across (lipids membrane)
          do (multiple-value-bind (qw qx qy qz tx ty tz)
                 (geom:matrix-to-quaternion-translation (transform lipid))
               (setf (elt rigid-body-coords (+ 0 index)) qw
                     (elt rigid-body-coords (+ 1 index)) qx
                     (elt rigid-body-coords (+ 2 index)) qy
                     (elt rigid-body-coords (+ 3 index)) qz
                     (elt rigid-body-coords (+ 4 index)) tx
                     (elt rigid-body-coords (+ 5 index)) ty
                     (elt rigid-body-coords (+ 6 index)) tz)))
    rigid-body-coords))


(defun rigid-body-coordinates-to-lipid-transforms (membrane rigid-body-coords)
  (loop for index from 0 below (length rigid-body-coords) by 7
        for lipid across (lipids membrane)
        for qw = (elt rigid-body-coords (+ 0 index))
        for qx = (elt rigid-body-coords (+ 1 index))
        for qy = (elt rigid-body-coords (+ 2 index))
        for qz = (elt rigid-body-coords (+ 3 index))
        for tx = (elt rigid-body-coords (+ 4 index))
        for ty = (elt rigid-body-coords (+ 5 index))
        for tz = (elt rigid-body-coords (+ 6 index))
        for matrix = (geom:make-matrix nil)
        do (geom:set-from-quaternion matrix qw qx qy qz tx ty tz)
        do (setf (transform lipid) matrix)))


(defun reset-lipid-transforms (membrane) 
  (loop for lipid across (lipids membrane)
        for matrix = (geom:make-matrix t)
        do (setf (transform lipid) matrix)))


#+(or)
(defun transform-lipids-using-rigid-body-coordinates (membrane rigid-body-coordinates)
  "Build an aggregate for the membrane and return it as well as a vector of vectors of atoms for
testing the scoring."
  (let* ((aggregate (chem:make-aggregate)) ; make an empty aggregate
         (lipids (lipids membrane))
         (index -1)
         (bounding-box (bounding-box membrane))
         (lipid-molecules (make-array 1024 :adjustable t :fill-pointer 0))
         (atom-vectors (make-instance 'atom-vectors
                                      :thing-array (make-array (length lipids) :element-type t)
                                      :bounding-box bounding-box)))
    (chem:set-bounding-box aggregate (bounding-box membrane))
    (let ((molecule-atoms (make-array 256 :fill-pointer 0 :adjustable t)))
      (loop for index below (length lipids)
            for ga-thing = (elt lipids index)
            unless (typep ga-thing 'ga-solute)
              do (let* ((ga-lipid ga-thing)
                        (transform (transform ga-lipid))
                        (lipid-info (lipid-info ga-lipid))
                        (lipid-index (lipid-index ga-lipid)))
                   (multiple-value-bind (lipid-mol lipid-atom-vector)
                       (make-molecule-from-ga-lipid ga-lipid :debug debug)
                     (chem:add-molecule aggregate lipid-mol)
                     (let ((lipid-mol-index (vector-push-extend lipid-mol lipid-molecules)))
                       (setf (elt (thing-array atom-vectors) index) lipid-mol-index))
                     (vector-push-extend lipid-atom-vector molecule-atoms))))
      ;; add in the solute atoms at the end
      (when (slot-boundp membrane 'ga-solute)
        (let ((solute (chem:matter-copy (solute (ga-solute membrane)))))
          (chem:apply-transform-to-atoms solute (transform (ga-solute membrane)))
          (let ((solute-atoms (make-array (chem:number-of-atoms solute)))
                (si -1))
            (cando:do-molecules (mol solute)
              (let ((mol-list (list mol)))
                (chem:add-molecule aggregate mol)
                (cando:do-residues (res mol)
                  (let ((res-list (list* res mol-list)))
                    (cando:do-atoms (atm res)
                      (let ((atm-list (list* atm res-list)))
                        (setf (elt solute-atoms (incf si)) (make-instance 'atm-res-mol :atm atm :res res :mol mol))))))))
            (setf (solute-atom-vector-index atom-vectors) (length molecule-atoms))
            (vector-push-extend solute-atoms molecule-atoms))))
      (setf (atom-vectors atom-vectors) molecule-atoms)
      ;; (check-aggregate-atom-vector-atom-order aggregate atom-vectors)
      (values aggregate atom-vectors))))
