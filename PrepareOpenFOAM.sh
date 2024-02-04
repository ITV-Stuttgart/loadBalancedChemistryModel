#!/bin/bash
cd ${0%/*} || exit 1    # run from this directory

cat << EOF
The TDACChemistryModel of OpenFOAM uses a private scope for the variables
and member functions. This however prevents other classes, such as the 
LoadBalancedTDACChemistryModel to derive from this class without replicating
all functionality. Further, the tabulation methods require a TDACChemistryModel.

Therefore, this script will add the protected keyword, making the private 
member variables protected in the TDACChemistryModel
EOF


sed -i '88i protected:' ${WM_PROJECT_DIR}/src/thermophysicalModels/chemistryModel/chemistryModel/TDACChemistryModel/TDACChemistryModel.H

# ----------------------------------------------------------------- end-of-file
