quickload "membrane-packer"

(membrane-packer::setup-system-membrane-force-field)

(list-force-fields)

setupAmberPaths

source "leaprc.protein.ff14SB"

source "leaprc.lipid14"

source "leaprc.water.tip3p"

source "leaprc.gaff"

prot = loadMoe /Users/meister/Development/cando-dev/extensions/cando/src/lisp/membrane-packer/data/merck_suvorexant_derivatives.moe

ligands = (loop for index from 3 to 20 collect (removeMatter prot (molid prot index)))

desc ligands


desc prot

set prot/1.1 name :NPRO
set prot/1.205 name :CGLN
set prot/2.206 name :NLYS
set prot/2.293 name :CCYS

assignAtomTypes prot

(move-geometric-center-to-origin prot)


memb = (membrane-packer:build-ga-membrane prot)

aggtop = (membrane-packer::build-simple-membrane-aggregate memb 10 :top)


ef-frozen = (multiple-value-list (membrane-packer::build-simple-membrane-energy-function memb 10 :top))

ef = (first ef-frozen)

(chem:enable-debug ef)

(chem:calculate-energy-and-force ef)


