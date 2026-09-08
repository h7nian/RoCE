#!/bin/bash
#SBATCH --job-name=install_em
#SBATCH --time=00:40:00
#SBATCH --mem=10g
#SBATCH --cpus-per-task=4
#SBATCH --output=diagnosis/face_probe/validation/install_em_%j.out
module load R/4.2.2-gcc-8.2.0-vp7tyde
export R_LIBS_USER=$HOME/Rlibs
mkdir -p $HOME/Rlibs_em
cd "${ROCE_PROJECT_ROOT:-$(git rev-parse --show-toplevel)}"
echo "Installing EM-enabled RoCE to ~/Rlibs_em (separate library) ..."
R CMD INSTALL --library=$HOME/Rlibs_em --no-multiarch --no-test-load . 2>&1 | tail -15
echo "INSTALL exit: ${PIPESTATUS[0]}"
ls -la $HOME/Rlibs_em/RoCE/libs/*.so 2>/dev/null && echo "[OK] .so present"
