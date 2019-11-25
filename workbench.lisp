(eval-when (:compile-toplevel :execute :load-toplevel)
  (ql:quickload "membrane-packer"))

(setupAmberPaths)

(source "leaprc.protein.ff14SB")

(source "leaprc.lipid14")

(source "leaprc.water.tip3p")

(source "leaprc.gaff")

(defparameter prot (loadMoe "/Users/meister/Development/cando-dev/extensions/cando/src/lisp/membrane-packer/data/merck_suvorexant_derivatives.moe"))

(defparameter ligands (loop for index from 3 to 20 collect (removeMatter prot (molid prot index)) ))

(progn
  (leap.commands:leap-set (resid (molid prot 1) 1) :name :npro)
  (leap.commands:leap-set (resid (molid prot 1) 205) :name :ngln)
  (leap.commands:leap-set (resid (molid prot 2) 206) :name :nlys)
  (leap.commands:leap-set (resid (molid prot 2) 293) :name :ccys))

(assignAtomTypes prot)

(move-geometric-center-to-origin prot)

(.= scored-membranes (membrane-packer:pack prot :number-of-generations 50))


(.= ef (chem:make-energy-function :bounding-box (chem:make-bounding-box '(20.0 20.0 20.0))))
(.= nb (chem:get-nonbond-component ef))
(.= at (chem:atom-table ef))

(.= tm (membrane-packer:build-ga-membrane prot))
(.= agg (membrane-packer::build-aggregate-from-ga-membrane tm))


#|
(.= memb  (membrane-packer::membrane (first  |scored-membranes|))

membagg = (membrane-packer::build-aggregate-from-ga-membrane |memb|)

(defparameter *rbmin* (membrane-packer::build-rigid-body-minimizer |memb| :debug t))

*rbmin*

(chem:minimize |*rbmin*|)
|#


