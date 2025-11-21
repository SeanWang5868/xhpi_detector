# XH-$\pi$ Interaction Detector

**xhpi** is designed for detection of XH-$\pi$ interactions in protein structures (PDB/mmCIF). It features automated hydrogen addition and geometric validation based on Hudson/Plevin criteria


## Installation

Requires Python $\ge$ 3.9.

```bash
# 1. Clone the repository
git clone [https://github.com/SeanWang5868/xhpi_detector](https://github.com/SeanWang5868/xhpi_detector)
cd xhpi
pip install .
```

## Configuration

The detection of xh $/pi$ interactions depends on the position of H atoms in the protein structure. In order to add H to the structure to be detected, the path to the monomer library (such as the CCP4 monomer library) needs to be specified first.

```bash
xhpi --set-mon-lib /path/to/ccp4/monomers
```

## Usage

### Basic Analysis
Recursively scan a directory for `.cif` files and output JSON results to the same directory.

```bash
xhpi ./pdb_data/
```

## Detection Logic

The tool classifies interactions based on two distinct geometric systems:

### 1. Plevin System
* **Distance ($d$)**: $< 4.3 \text{\AA}$
* **Angle ($XH\dots\pi$)**: $> 120^\circ$
* **Angle ($XPCN$)**: $< 25^\circ$

### 2. Hudson System
* **Distance ($d$)**: $\le 4.5 \text{\AA}$
* **Angle ($\theta_{norm}$)**: $\le 40^\circ$
* **Distance ($d_{proj}$)**: 
    * $\le 1.6 \text{\AA}$ for **HIS** and **TrpA**.
    * $\le 2.0 \text{\AA}$ for **TRP**, **TYR**, **PHE**.

## Contact
sean.wang@york.ac.uk

York Structural Biology Laboratory (YSBL)
Department of Chemistry, University of York
Heslington, York, YO10 5DD, UK
