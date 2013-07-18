export TERM=xterm-color

module load IDL
source ~/astronomy/lib/export_idl.sh

export PLATEDESIGN_DIR=~/astronomy/proposals/BOSSANOVA/platedesign/
export BOSSTILE_DIR=~/astronomy/proposals/BOSSANOVA/bosstile

export IDL_PATH=+$PLATEDESIGN_DIR/pro:+$BOSSTILE_DIR/pro:$IDL_PATH

