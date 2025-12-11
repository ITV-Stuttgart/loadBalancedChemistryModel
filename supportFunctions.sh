#======================================================================
# Support Functions
#======================================================================

__die()
{
    RED='\033[0;31m'
    echo -e "${RED}$*"
    exit 1
}


__banner()
{
    echo "==============================================================================="
    echo $*
    echo "==============================================================================="
}


__warningBanner()
{
    BYellow='\033[1;33m'
    Color_Off='\033[0m'
    >&2 echo -e "${BYellow}==============================================================================="
    >&2 echo -e "Warning: $*"
    >&2 echo -e "=============================================================================== ${Color_Off}"
}


__abort()
{
        cat <<EOF
***************
*** ABORTED ***
***************
An error occurred. Exiting...
EOF
        exit 1
}

__checkOpenFOAMEnvironment()
{
    if [[ -z ${WM_PROJECT_DIR} ]]; then
        __banner Checking OpenFOAM environment
        __die OpenFOAM environment is not sourced
    fi
}


compileTDACLibrary() {
    __checkOpenFOAMEnvironment

    FILE="${WM_PROJECT_DIR}/src/thermophysicalModels/chemistryModel/chemistryModel/TDACChemistryModel/TDACChemistryModel.H"
    if sed -n "88p" "$FILE" | grep -q "protected:"; then
        echo "OpenFOAM is ready to compile the load balanced TDAC model"
        true
    else
        if [[ -r $FILE && -w $FILE ]]; then 
            cat << EOF
    ================================================================================
    The TDACChemistryModel of OpenFOAM uses a private scope for the variables
    and member functions. This however prevents other classes, such as the 
    LoadBalancedTDACChemistryModel to derive from this class without replicating
    all functionality. Further, the tabulation methods require a TDACChemistryModel.

    Therefore, this script will add the protected keyword, making the private 
    member variables protected in the TDACChemistryModel
    ================================================================================
EOF
            sed -i '88i protected:' ${WM_PROJECT_DIR}/src/thermophysicalModels/chemistryModel/chemistryModel/TDACChemistryModel/TDACChemistryModel.H
            true
            return
        else
            __warningBanner "TDAC model is not compiled"
            echo "No write permissions for $FILE"
            cat << EOF
    ================================================================================
    The TDACChemistryModel of OpenFOAM uses a private scope for the variables
    and member functions. This however prevents other classes, such as the 
    LoadBalancedTDACChemistryModel to derive from this class without replicating
    all functionality. Further, the tabulation methods require a TDACChemistryModel.

    Without write permissions the load-balanced TDAC model cannot be compiled.
    ================================================================================
EOF
            false
            return
        fi
    fi
}



# ----------------------------------------------------------------- end-of-file
