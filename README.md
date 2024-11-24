
# DAC 2025 Repository for ROAR (Robust Optimal Analog Reuse) Software, Scripts, Results, and Images
## Author: Alec S. Adair 
For questions or comments e-mail alecadair1@gmail.com
ROAR (Robust Optimal and Analog Reuse) - Turku, Finland 2024

This repository contains the implementation of the Robust and Efficient Analog IC Design Automation framework described in the paper under review, **"Analytical Optimization for Robust and Efficient Analog IC Design Automation"**. The repository includes the tools to run the code, generate figures, and reproduce results presented in the paper.

## Prerequisites
Code Intended to run with Centos/Rocky/Redhat and Ubuntu operating systems
tcsh or csh based shell should be used to run setups and installation

To run this software, ensure the following prerequisites are met:

1. **Python 3.x** is installed with the following libraries:
    - `matplotlib`
    - `pandas`
    - `numpy`
    - Standard Python libraries: `os`, `sys`, `math`

2. A **Make** utility is installed on your system.

3. You must use either the `csh` or `tcsh` shell to execute the commands.

## Installation Instructions

1. Clone or download the repository:
    ```bash
    git clone https://github.com/alecadair/DAC2025_ROAR.git
    cd DAC2025_ROAR
    ```

2. Run `Make` at the top level of the repository:
    ```csh
    make
    ```

   This will generate a file called `roar_env.csh`.

3. Source the `roar_env.csh` file to set up all necessary environment variables:
    ```csh
    source roar_env.csh
    ```

## Running the Example

The repository includes an example design of a **N-Input Current Mirror OTA**, as described in the paper. To run this example:

1. Navigate to the `design` directory:
    ```bash
    cd design
    ```

2. Run the script `cm_ota.py` using Python:
    ```bash
    python3 cm_ota.py
    ```

## Additional Notes

- The scripts will automatically generate the figures and results presented in the paper.
- Ensure all environmental variables are correctly set by sourcing `roar_env.csh` before running any scripts.

## Contributing

Feel free to submit issues or pull requests to enhance the functionality or documentation of this repository.

## License

This project is licensed under the terms specified in the repository's LICENSE file.
