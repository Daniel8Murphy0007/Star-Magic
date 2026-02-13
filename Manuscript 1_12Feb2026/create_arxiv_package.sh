#!/bin/bash
# Create arXiv submission package for UQFF Production Manuscript
# Usage: bash create_arxiv_package.sh

set -e  # Exit on error

echo "========================================"
echo "arXiv Submission Package Creator"
echo "UQFF Production Manuscript"
echo "========================================"
echo ""

# Configuration
PACKAGE_NAME="uqff_arxiv_submission"
MANUSCRIPT="uqff_production_arxiv.tex"
FIGURES_DIR="figures"
ANC_DIR="anc"

# Clean previous build
echo "[1/6] Cleaning previous builds..."
rm -rf ${PACKAGE_NAME}
rm -f ${PACKAGE_NAME}.tar.gz
rm -f uqff_production_arxiv.aux uqff_production_arxiv.log uqff_production_arxiv.out
echo "✓ Cleaned"

# Create package directory structure
echo "[2/6] Creating package structure..."
mkdir -p ${PACKAGE_NAME}
mkdir -p ${PACKAGE_NAME}/${FIGURES_DIR}
mkdir -p ${PACKAGE_NAME}/${ANC_DIR}
echo "✓ Directories created"

# Generate figures
echo "[3/6] Generating figures..."
cd manuscript
python generate_figures.py
if [ $? -ne 0 ]; then
    echo "✗ Figure generation failed!"
    exit 1
fi
cd ..
echo "✓ Figures generated"

# Compile manuscript (twice for references)
echo "[4/6] Compiling manuscript..."
cd manuscript
pdflatex -interaction=nonstopmode ${MANUSCRIPT} > /dev/null 2>&1
pdflatex -interaction=nonstopmode ${MANUSCRIPT} > /dev/null 2>&1
if [ ! -f uqff_production_arxiv.pdf ]; then
    echo "✗ Compilation failed!"
    exit 1
fi
cd ..
echo "✓ Manuscript compiled (31 pages)"

# Copy files to package
echo "[5/6] Copying files to package..."
cp manuscript/${MANUSCRIPT} ${PACKAGE_NAME}/
cp manuscript/${FIGURES_DIR}/*.pdf ${PACKAGE_NAME}/${FIGURES_DIR}/
cp manuscript/generate_figures.py ${PACKAGE_NAME}/${ANC_DIR}/
cp manuscript/README_ARXIV.md ${PACKAGE_NAME}/${ANC_DIR}/
cp manuscript/SUBMISSION_CHECKLIST.md ${PACKAGE_NAME}/${ANC_DIR}/
echo "✓ Files copied"

# Create tarball
echo "[6/6] Creating submission archive..."
tar -czf ${PACKAGE_NAME}.tar.gz ${PACKAGE_NAME}/
SIZE=$(du -h ${PACKAGE_NAME}.tar.gz | cut -f1)
echo "✓ Archive created: ${PACKAGE_NAME}.tar.gz (${SIZE})"

# Summary
echo ""
echo "========================================"
echo "✓ PACKAGE READY FOR SUBMISSION"
echo "========================================"
echo ""
echo "Contents:"
echo "  - Main manuscript: ${MANUSCRIPT}"
echo "  - Figures: 5 PDF files"
echo "  - Ancillary: 3 supplementary files"
echo ""
echo "Archive: ${PACKAGE_NAME}.tar.gz (${SIZE})"
echo ""
echo "Next steps:"
echo "  1. Review: manuscript/uqff_production_arxiv.pdf"
echo "  2. Upload: ${PACKAGE_NAME}.tar.gz to arxiv.org/submit"
echo "  3. Categories: gr-qc (primary), astro-ph.CO, hep-ph, physics.comp-ph"
echo ""
echo "Submission URL: https://arxiv.org/submit"
echo ""

# Optional: Open PDF for review
if command -v xdg-open &> /dev/null; then
    read -p "Open PDF for review? [y/N] " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        xdg-open manuscript/uqff_production_arxiv.pdf
    fi
elif command -v open &> /dev/null; then
    read -p "Open PDF for review? [y/N] " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        open manuscript/uqff_production_arxiv.pdf
    fi
fi

echo "Done!"
