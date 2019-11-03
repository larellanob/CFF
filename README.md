# CFF

Compact file format for CLAS data. Takes large data/simulation files
in UTFSM cluster and applies filters, generating a more manageable
file. Particle identification is made by Analyser.

## Installation

### Requirements
* ROOT
* ClasTool - https://github.com/orsosa/ClasTool
* Analyser (**install development branch at least**, or MW_branch, but not master as it is too outdated) - https://github.com/orsosa/Analyser/tree/development
* clas_pack - https://github.com/orsosa/clas_pack
* clas_lib - https://github.com/orsosa/clas_lib
* makedepend (for compiling some of these packages)

Compile with make when necessary


### Environment variables
* Add Analyser and ClasTool "include" folders to ROOT_INCLUDE_PATH
* Add ClasTool/slib/$OS_NAME and Analyser/slib to LD_LIBRARY_PATH

* see example .bashrc file - https://pastebin.com/bebvjcNa

## Running
* Compile running Makefile
* Make a directory named "local" in the main program directory.
* Generate Trees using short_tree and Ntuples using short_ntuple (development on short_tree is privileged)

### Simulation files

* You need a file named simulFiles.txt in the same directory as
  short_tree which contains the path to the full simulation files you
  want to filter

  e.g.
  ```
  [arellano@ui02 CFF]: head simulFiles.txt 
  /data/tsunami/user/o/orsosa/HSim/D2_pb/ROOT/simul_D2_Pb101.root
  /data/tsunami/user/o/orsosa/HSim/D2_pb/ROOT/simul_D2_Pb102.root
  /data/tsunami/user/o/orsosa/HSim/D2_pb/ROOT/simul_D2_Pb103.root
  /data/tsunami/user/o/orsosa/HSim/D2_pb/ROOT/simul_D2_Pb104.root
  /data/tsunami/user/o/orsosa/HSim/D2_pb/ROOT/simul_D2_Pb105.root
  /data/tsunami/user/o/orsosa/HSim/D2_pb/ROOT/simul_D2_Pb106.root
  /data/tsunami/user/o/orsosa/HSim/D2_pb/ROOT/simul_D2_Pb107.root
  /data/tsunami/user/o/orsosa/HSim/D2_pb/ROOT/simul_D2_Pb108.root
  /data/tsunami/user/o/orsosa/HSim/D2_pb/ROOT/simul_D2_Pb109.root
  /data/tsunami/user/o/orsosa/HSim/D2_pb/ROOT/simul_D2_Pb10.root  
  ```
* Run using
  ```
  ./short_tree mandatorytext
  ```

  Note that the argument in this case is mandatory. The argument will
  be added to the resulting filename.
  
* The resulting tree will be in ./local/CFFTree_mandatorytext.root

### Data files

* Put data files to filter in a file named dataFiles.txt

  e.g.
  ```
  [arellano@ui02 CFF]: head dataFiles.txt 
  /data/jlab/mss/clas/eg2a/production/Pass2/Clas/clas_42071_00.pass2.root
  /data/jlab/mss/clas/eg2a/production/Pass2/Clas/clas_42071_01.pass2.root
  /data/jlab/mss/clas/eg2a/production/Pass2/Clas/clas_42071_02.pass2.root
  /data/jlab/mss/clas/eg2a/production/Pass2/Clas/clas_42071_03.pass2.root
  /data/jlab/mss/clas/eg2a/production/Pass2/Clas/clas_42071_04.pass2.root
  /data/jlab/mss/clas/eg2a/production/Pass2/Clas/clas_42071_05.pass2.root
  /data/jlab/mss/clas/eg2a/production/Pass2/Clas/clas_42071_06.pass2.root
  /data/jlab/mss/clas/eg2a/production/Pass2/Clas/clas_42071_07.pass2.root
  /data/jlab/mss/clas/eg2a/production/Pass2/Clas/clas_42071_08.pass2.root
  /data/jlab/mss/clas/eg2a/production/Pass2/Clas/clas_42071_09.pass2.root  
  ```
* Run without arguments
  ```
  ./short_tree
  ```
  
* The resulting file will be in ./local/CFFTree_data.root

  * Remember to rename this output file if running multiple times, otherwise it
    will just be overwritten.


### TO DO

* Take target selection as command line input rather than it being hardcoded in short_tree.cxx

* Make particle selection selectable by user in a user-friendly way (e.g. config file, cli, etc.)

* Fix huge memory consumption, there's probably a leak

* Keep short_ntuple.cxx up to date or otherwise eliminate it
  * Maybe merge it with short_tree and make ntuples an option within it