# Philips-to-MRD Converter

[![License](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Python](https://img.shields.io/badge/Python-3.7%2B-blue.svg)](https://www.python.org/downloads/)

## Overview

This tool converts Philips raw MRI data into the [ISMRMRD](https://ismrmrd.github.io/) (International Society for Magnetic Resonance in Medicine Raw Data) format. ISMRMRD is an open standard for MRI data that enables compatibility with a wide range of open-source reconstruction and processing tools.

### Key Features
- Converts Philips raw data from various scanner versions and sequence types
- Supports both Cartesian and non-Cartesian (radial, spiral) trajectories
- Handles multiple data formats: `.data/.list` and `.raw/.lab/.sin`
- Preserves important metadata including encoding parameters, sequence timing, and system information
- Includes specialized converter for Xenon gas exchange imaging
- Supports multiple trajectory ordering schemes: standard interleaved Archimedean spiral, golden mean ordering, and Halton randomized Archimedean spiral

### Why Convert to ISMRMRD?
ISMRMRD provides a standardized format that enables:
- Interoperability with popular MRI reconstruction frameworks (e.g., [Gadgetron](https://gadgetron.github.io/))
- Use of open-source reconstruction tools and algorithms
- Simplified data sharing and collaboration
- Reproducible MRI data processing pipelines

## Installation

### Prerequisites
- Python 3.7 or higher (tested with Python 3.12)
- pip package manager

### Steps

1. Clone the repository:
   ```bash
   git clone https://github.com/Xe-MRI-CTC/philips-to-mrd-converter.git
   cd philips-to-mrd-converter
   ```

2. Set up a virtual environment (recommended):
   ```bash
   python -m venv .venv
   source .venv/bin/activate  # On Windows: .venv\Scripts\activate
   ```

3. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

## Usage

### Basic Conversion

The main conversion functionality is provided by the `Ph2Mrd` class:

```python
from pathlib import Path
import philips2mrd as p2m

# Set file paths
dl_name = Path("path/to/your/file.data")  # .data/.list files
rls_name = Path("path/to/your/file.sin")  # .raw/.lab/.sin files
out_dir = Path("path/to/output/directory")

# Create converter instance
converter = p2m.Ph2Mrd(dl_name, rls_name)

# Run conversion
mrd_file, rls_data, dl_data = converter.convert(out_dir)
print(f"Data converted to: {mrd_file}")
```

### Advanced Usage

#### Trajectory Ordering

The `trajorder` parameter specifies the trajectory ordering pattern for non-Cartesian reconstructions:

```python
converter = p2m.Ph2Mrd(dl_name, rls_name)
converter.trajorder = 2  # 0: Standard interleaved Archimedean spiral, 1: Golden mean, 2: Halton randomized
````

#### Manual Gradient Delay

For non-Cartesian reconstructions, you may need to specify a manual gradient delay:

```python
converter.delay = -1.25  # microseconds
```

#### Conversion with Different File Combinations

The converter supports different combinations of Philips data files:

```python
# Only .data/.list files
converter = p2m.Ph2Mrd(dl_name, None)

# Only .raw/.lab/.sin files
converter = p2m.Ph2Mrd(None, rls_name)

# Both file types (recommended for complete metadata)
converter = p2m.Ph2Mrd(dl_name, rls_name)
```

### Xenon Gas Exchange Imaging

This repository also includes specialized scripts for converting Xenon gas exchange data to the format required by the [Xenon Clinical Trials Consortium](https://www.129xectc.org/). See the [`Scripts/XeGasExchange2XeCTCMRD.py`](Scripts/XeGasExchange2XeCTCMRD.py) script for details.

### Example Script

The [`Example.py`](Example.py) file demonstrates how to:
- Run the converter
- Load the resulting ISMRMRD file
- Extract metadata and data
- Visualize trajectories and acquisition data

Run it with:
```bash
python Example.py
```

## Data Formats

### Supported Input Formats

| Format | Extensions | Description |
|--------|------------|-------------|
| Data/Lists | `.data`, `.list` | Raw data and acquisition list files |
| Raw/Lab/Sin | `.raw`, `.lab`, `.sin` | Philips raw data files (.sin files are *not* sinogram files, but Philips-specific acquisition data) |

### Output Format

The converter outputs data in ISMRMRD format (HDF5 file with `.h5` extension) containing:
- **Header**: Full ISMRMRD header with experimental conditions, system information, encoding parameters, and acquisition metadata
- **Acquisitions**: All MRI acquisitions with k-space data and trajectory information

## Testing

To run the test suite:

```bash
python test_philips2mrd.py
```

The test suite includes conversion tests for various data types:
- 2D Spiral
- 3D Radial (CTC gas exchange)
- 3D Radial with bonus spectroscopy
- 3D Radial with 2 echoes
- 3D FLORET with gas exchange

Test data is included in the `testdata/` directory.

## API Reference

### Ph2Mrd Class

The main converter class.

#### Constructor
```python
Ph2Mrd(dlName=None, rlsName=None)
```
- `dlName`: Path to .data file (can also be .list file)
- `rlsName`: Path to .sin file (can also be .raw or .lab file)

#### Attributes
- `trajorder`: Trajectory ordering pattern (0: Standard interleaved Archimedean spiral, 1: Golden mean, 2: Halton randomized)
- `delay`: Manual gradient delay in microseconds

#### Methods
```python
convert(outDir)
```
- `outDir`: Output directory for the ISMRMRD file
- Returns: Tuple of (MRD filename, RLS data, DL data)

## Limitations

- The converter currently focuses on standard MRI acquisition types
- Some sequence-specific parameters may require manual adjustment
- Non-Cartesian trajectory support is functional but may require calibration for specific sequences

## Contributing

Contributions are welcome! Please follow these steps:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Contact

- Repository: [https://github.com/Xe-MRI-CTC/philips-to-mrd-converter](https://github.com/Xe-MRI-CTC/philips-to-mrd-converter)
- Issues: [Report bugs or request features](https://github.com/Xe-MRI-CTC/philips-to-mrd-converter/issues)
- Xenon Clinical Trials Consortium: [https://www.129xectc.org/](https://www.129xectc.org/)

---
